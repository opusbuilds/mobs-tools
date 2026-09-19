#!/usr/bin/env python3
"""
Pre-flight an EXOTIC inits file before spending an hour reducing with it.

Written 2026-07-27 after finding that two of five inits files here still carried
the shipped HAT-P-32 b sample's ephemeris (P=2.1500082 d) because the sample was
used as a template and the planetary block never edited. One of them also aimed
its target aperture 4.5 arcmin off the star. Both nights were written off with
confident, wrong diagnoses. Every check below would have failed loudly.

    ../venv/bin/python check_inits.py ../inits_corot2_resolved.json

Checks, in order of how badly they bite:
  1. planetary parameters against the NASA Exoplanet Archive
  2. the predicted mid-transit actually falls inside the observing window
  3. the target pixel actually lands on the target, via a local plate solve
     of three frames, because these fields drift ~100 px across a night
  4. comparison stars are inside the frame and not saturated
"""
import csv, io, json, os, subprocess, sys, tempfile, urllib.parse, urllib.request
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def _urlopen_retry(url, timeout, tries=2):
    """The NASA Exoplanet Archive occasionally stalls for a full timeout and then
    answers the next request in two seconds (first seen 2026-09-13 from the home
    box). One retry covers it; a second failure is a real outage."""
    import time
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=timeout)
        except (TimeoutError, OSError) as e:
            if i == tries - 1:
                raise
            time.sleep(3)


warnings.filterwarnings('ignore')
TAP = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'
FAILED = []


def check(ok, label, detail=''):
    mark = 'ok  ' if ok else 'FAIL'
    print(f'  [{mark}] {label}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAILED.append(label)


def warn(ok, label, detail=''):
    """Report without failing: things that need a judgment call, not a rule."""
    print(f"  [{'ok  ' if ok else 'look'}] {label}" + (f'  {detail}' if detail else ''))


def archive(planet):
    cols = ('pl_orbper,pl_tranmid,pl_ratror,pl_ratdor,pl_orbincl,ra,dec,sy_vmag,pl_trandur,default_flag')
    q = f"select {cols} from ps where pl_name='{planet}'"
    url = TAP + urllib.parse.urlencode({'query': q, 'format': 'csv'})
    rows = list(csv.DictReader(io.StringIO(_urlopen_retry(url, 60).read().decode())))
    if not rows:
        return None
    merged = {}
    for row in sorted(rows, key=lambda r: r['default_flag'] != '1'):
        for k, v in row.items():
            if v not in ('', 'null'):
                merged.setdefault(k, v)
    return merged


def window(fits_dir):
    files = sorted(f for f in os.listdir(fits_dir) if f.lower().endswith(('.fits', '.fit', '.fts')))
    if not files:
        raise SystemExit(f'no FITS files in {fits_dir}')
    picks = [os.path.join(fits_dir, files[i]) for i in (0, len(files) // 2, -1)]
    jd = [fits.getheader(f)['MJD-OBS'] + 2400000.5 for f in (picks[0], picks[-1])]
    return picks, jd[0], jd[1]


def solve(path, workdir):
    hdr = fits.getheader(path)
    scale = float(hdr['IM_SCALE']) if 'IM_SCALE' in hdr else 206.265 * float(hdr['XPIXSZ']) / float(hdr['FOCALLEN'])
    cmd = ['solve-field', '--no-plots', '--overwrite', '--dir', workdir, '--downsample', '2', '--cpulimit', '120',
           '--scale-units', 'arcsecperpix', '--scale-low', f'{scale * 0.85:.3f}', '--scale-high', f'{scale * 1.15:.3f}',
           '--ra', str(float(hdr['RA'])), '--dec', str(float(hdr['DEC'])), '--radius', '2', path]
    subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    stem = os.path.join(workdir, os.path.splitext(os.path.basename(path))[0])
    return WCS(fits.getheader(stem + '.wcs')) if os.path.exists(stem + '.solved') else None


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print('usage: check_inits.py <inits.json>'); sys.exit(2)
    path = sys.argv[1]
    d = json.load(open(path))
    u, p = d['user_info'], d['planetary_parameters']
    planet = p['Planet Name']
    print(f'{os.path.basename(path)}  ->  {planet}')

    a = archive(planet)
    if a is None:
        check(False, 'planet is in the NASA Exoplanet Archive', planet)
        return report()

    print('\narchive agreement')
    pairs = [('Orbital Period (days)', 'pl_orbper', 1e-4),
             ('Ratio of Planet to Stellar Radius (Rp/Rs)', 'pl_ratror', 0.15),
             ('Ratio of Distance to Stellar Radius (a/Rs)', 'pl_ratdor', 0.15),
             ('Orbital Inclination (deg)', 'pl_orbincl', 0.02)]
    for key, col, tol in pairs:
        if a.get(col) in (None, '', 'null'):
            # pscomppars does not always carry the column (TOI-2570 b, 2026-09-10: no pl_ratror);
            # nothing to compare against is a skip, not a failure and not a crash.
            print(f'  [skip] {key.split(" (")[0]}  archive has no {col}; inits {p.get(key)}')
            continue
        mine, theirs = p.get(key), float(a[col])
        rel = abs(mine - theirs) / theirs
        check(rel <= tol, key.split(' (')[0], f'inits {mine}  archive {theirs}  ({rel * 100:.2f}% off)')

    # Internal consistency, no network needed: a/Rs and the inclination fix the
    # impact parameter, and b >= 1 + Rp/Rs is a planet that misses its star. A
    # WBoM run on WASP-53 b (2026-09-15) reported inc 83.40 and a/Rs 11.78, b =
    # 1.35, and timed a mid-transit on a flat line; the archive rows were fine.
    rprs, ars, inc = (p.get(k) for k in ('Ratio of Planet to Stellar Radius (Rp/Rs)',
                                          'Ratio of Distance to Stellar Radius (a/Rs)',
                                          'Orbital Inclination (deg)'))
    if None not in (rprs, ars, inc):
        import math
        b = abs(ars * math.cos(math.radians(inc)))
        detail = f'b = a/Rs cos(i) = {b:.3f}  (1 + Rp/Rs = {1 + rprs:.3f})'
        check(b < 1 + rprs, 'given geometry transits', detail)
        if b < 1 + rprs:
            warn(b <= 1 - rprs, 'geometry is not grazing', detail)

    # A null in the planetary parameters is a FAIL, not a skip: EXOTIC's own
    # pre-flight rejects the file (TOI-5300 b, 2026-09-19, a/Rs None) and its
    # reduction would crash later. Cheaper to catch here than in the archive
    # comparison, which has nothing to compare a null against.
    nulls = [k for k, v in p.items() if v is None and 'Uncertainty' not in k]
    check(not nulls, 'no null planetary parameters', ', '.join(nulls) if nulls else 'all present')

    print('\ntiming')
    fits_dir = u['Directory with FITS files']
    picks, jd0, jd1 = window(fits_dir)
    first = picks[0]
    T0, P = p['Published Mid-Transit Time (BJD-UTC)'], p['Orbital Period (days)']
    n = round((jd0 - T0) / P)
    tmid = T0 + n * P
    dur = float(a['pl_trandur']) / 24 if a.get('pl_trandur') else 0.0
    print(f'  window {jd0:.5f} -> {jd1:.5f}  ({(jd1 - jd0) * 24:.2f} h)')
    check(jd0 <= tmid <= jd1, 'predicted mid-transit inside the window',
          f'Tmid {tmid:.5f} ({(tmid - jd0) * 24:+.2f} h from start)')
    if dur:
        check(tmid - dur / 2 >= jd0 and tmid + dur / 2 <= jd1, 'full transit inside the window',
              f'ingress {tmid - dur / 2:.5f}  egress {tmid + dur / 2:.5f}')

    # Solve three frames, not one. MicroObservatory does not guide well: CoRoT-2
    # moved 97 px across one night and KELT-20 about 100, so a seed taken from a
    # mid-night frame looks badly wrong against the first frame alone. Checking
    # only the first frame produced exactly that false failure on 2026-07-27.
    print('\npointing (local plate solve of first, middle and last frame)')
    ra, dec = float(a['ra']), float(a['dec'])
    got = np.array(u['Target Star X & Y Pixel'], dtype=float)
    offsets = {}
    with tempfile.TemporaryDirectory(prefix='preflight-') as work:
        for label, path in zip(('first', 'middle', 'last'), picks):
            wcs = solve(path, work)
            if wcs is None:
                print(f'  [    ] {label} frame did not solve  {os.path.basename(path)}')
                continue
            want = np.array(wcs.all_world2pix(ra, dec, 0), dtype=float)
            offsets[label] = (want, float(np.hypot(*(want - got))))
            print(f'  [    ] {label:6} frame: target at ({want[0]:6.1f}, {want[1]:6.1f}), '
                  f'{offsets[label][1]:5.1f} px from the inits pixel')
    if not offsets:
        check(False, 'at least one frame plate-solves')
    else:
        best = min(offsets.items(), key=lambda kv: kv[1][1])
        check(best[1][1] <= 8, 'inits pixel is on the target somewhere in the night',
              f'closest is the {best[0]} frame at {best[1][1]:.1f} px')
        if 'first' in offsets:
            warn(offsets['first'][1] <= 8, 'inits pixel matches the FIRST frame',
                 f'{offsets["first"][1]:.1f} px; EXOTIC seeds from the first frame, '
                 f'so a mid-night seed can still cost the run')

    print('\ncomparison stars')
    data = fits.getdata(first).astype(float)
    ny, nx = data.shape
    sat = float(fits.getheader(first).get('DATAMAX', data.max()))
    bg = float(np.median(data))
    noise = float(np.std(data[data < np.percentile(data, 95)]))

    def peak_over_bg(cx, cy):
        if not (0 <= cx < nx and 0 <= cy < ny):
            return None
        return float(data[max(0, int(cy) - 5):int(cy) + 6, max(0, int(cx) - 5):int(cx) + 6].max()) - bg

    tgt = peak_over_bg(*u['Target Star X & Y Pixel'])
    print(f'  background {bg:.0f}, noise {noise:.0f}, target sits {tgt:.0f} ADU above background')
    for cx, cy in u['Comparison Star(s) X & Y Pixel']:
        pk = peak_over_bg(cx, cy)
        if pk is None:
            check(False, f'comp ({cx}, {cy})', 'outside the frame')
            continue
        # Too faint is as disqualifying as too bright: the 2026-06-30 run used four
        # "comparison stars" that were within noise of blank sky, and the 2026-07-11
        # run used comparisons 49x brighter than the target.
        check(pk + bg < sat * 0.95, f'comp ({cx}, {cy}) unsaturated', f'peak {pk + bg:.0f} / DATAMAX {sat:.0f}')
        # Brightness is a judgment call, not a rule, and one frame is a small sample:
        # on a clouded night the first frame can be the worst one. Reported, not enforced.
        ratio = pk / tgt if tgt else float('inf')
        warn(pk > 5 * noise and 0.1 <= ratio <= 20, f'comp ({cx}, {cy}) brightness',
             f'{pk:+.0f} ADU over background, {ratio:.1f}x the target')

    report()


def report():
    print()
    if FAILED:
        print(f'{len(FAILED)} check(s) failed. Do not reduce with this file yet.')
        sys.exit(1)
    print('All checks passed.')


if __name__ == '__main__':
    main()
