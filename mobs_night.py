#!/usr/bin/env python3
"""
One command from a MicroObservatory target name and date to a pre-flighted
EXOTIC inits file, with the reduction as an option.

    python3 tools/mobs_night.py Qatar-1 260907 --planet "Qatar-1 b"
    python3 tools/mobs_night.py Qatar-1 260907 --planet "Qatar-1 b" --run

Steps, each of which is an existing tool in this directory:
  1. list the night's frames on the MObs Image Directory and download the
     missing ones (mo-www's certificate expired 2026-09-06; verification falls
     back to off for that host, and says so)
  2. archive parameters for the planet (NASA Exoplanet Archive, pscomppars)
  3. observing window from the frame headers; epoch, predicted Tmid, ingress
     and egress; baseline minutes either side; stop if no transit is inside
  4. locate_target.py on the first frame: target pixel and comparison candidates
  5. night_triage.py over the night: cloud, drift, comp box, phase split
  6. choose comparisons inside the comp box within a brightness factor of the
     target; write the inits file from the archive values
  7. check_inits.py on the result
  9. a pre-registration scaffold (PROCEED only) with the bar DERIVED from the V-magnitude
     scatter calibration (0.84% at V 11.57 on 2026-09-05, photon scaling), for
     editing and committing BEFORE any fit
  8. a triage verdict: proceed or reject, with the numbers. With --run and a
     proceed verdict, launch EXOTIC detached (setsid) with post_run_check.py
     appended, exactly as the hand-written run scripts did.

The verdict is advisory. The inits file is written either way; --run refuses on
a reject verdict unless --force. Nothing here submits anything anywhere.
"""
import argparse, csv, io, json, math, os, re, ssl, subprocess, sys, urllib.parse, urllib.request
import numpy as np
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from night_triage import frame_time, detect, voted_shift  # noqa: E402

def _urlopen_retry(url, timeout, tries=2):
    """The NASA Exoplanet Archive occasionally stalls for a full timeout and then
    answers the next request in two seconds (first seen 2026-09-13 from the home
    box). One retry covers it; a second failure is a real outage."""
    import time
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=timeout)
        except urllib.error.HTTPError as e:
            # 2026-09-17: the archive put a Cloudflare browser challenge in front
            # of TAP/sync (403, cf-mitigated: challenge). A real browser passes it;
            # fetch through headless Chromium and hand back a file-like object.
            if e.code == 403 and e.headers.get('cf-mitigated') == 'challenge' and 'exoplanetarchive' in url:
                import subprocess
                # Optional: a browser-backed fetcher (see browser/tap-fetch.mjs in
                # opus-infra; set MOBS_TAP_FETCH to its path). Without one, the 403
                # is raised as before.
                fetcher = os.environ.get('MOBS_TAP_FETCH', '/opt/opus-infra/browser/tap-fetch.mjs')
                r = subprocess.run(['node', fetcher, url], capture_output=True, timeout=150) if os.path.exists(fetcher) else None
                if r is not None and r.returncode == 0 and r.stdout.strip():
                    print('  archive: API behind a browser challenge; fetched through headless Chromium instead')
                    return io.BytesIO(r.stdout)
            if i == tries - 1:
                raise
            time.sleep(3)
        except (TimeoutError, OSError) as e:
            if i == tries - 1:
                raise
            time.sleep(3)


LISTING = 'https://waps.cfa.harvard.edu/microobservatory/MOImageDirectory/ImageDirectory.php?SortBy=Filename&SortPos=DESC'
FITS_URL = 'https://mo-www.cfa.harvard.edu/ImageDirectory/{name}.FITS'
TAP = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'
SITE = {'lat': '+31.68', 'lon': '-110.88', 'elev': 1268}   # MObs, Whipple Observatory, Arizona
# Calibration from 2026-09-05: WASP-11 b (archive V 11.57) gave 0.84% residual scatter at 60 s;
# scatter/depth 0.42-0.43 gave Tmid bars of 6.6 and 9.8 min (KELT-23A, WASP-11) at ~3 min cadence.
CAL_SCATTER, CAL_V, CAL_RATIO, CAL_BAR_MIN = 0.84, 11.57, 0.425, 8.2
# Scatter law, revised 2026-09-18. The first law scaled the 09-05 point by
# photon noise, 10**(0.2*(V-11.57)), and predicted 1.31% at V 12.54; the night
# delivered 0.90%. Three clean nights now: 0.84% at V 11.57 (09-05), 0.85% at
# 11.55 (TOI-3693 b, 09-13), 0.90% at 12.54 (WASP-67 b, 09-18). That is a
# systematics floor with almost no magnitude slope on this instrument, so the
# law is now linear and shallow. The bar model was never the problem: fed the
# observed scatter it predicted 9.1 min for WASP-67 b against 8.6 delivered.
CAL_SLOPE = 0.05   # percent scatter per magnitude, from the three points above
CAL_V_KNEE = 12.6  # beyond this the target approaches the 200 ADU floor and photon
                   # noise takes over: HAT-P-54 b (V 13.40, ~60 ADU) gave another
                   # observer 2.6% scatter and a 24 min bar, not the flat law's 0.93%.
# Seed floor, peak pixel above background in the seed frame. Set at 300 from the nights that failed
# (WASP-80 34 ADU, Qatar-1 39, TrES-5); lowered to 200 on 2026-09-13 when WASP-52 b (197 ADU in a
# 76%-transmission frame 1, ~260 clear, V 12.19) ran under --force to a genuine QC PASS, KTMF 4.20,
# 8 min Tmid bar. One night; the number moves again on the next one that bears on it.
SEED_MIN_ADU = 200
FORCE_MAX_RATIO = 0.5              # predicted scatter/depth above which a forced no-seed night can't time the transit (TOI-2570 b 09-25)
# The CEILING. Set 2026-09-22 from HD 189733 b (V 7.67, 8 s exposures): PROCEED on
# every other gate with its core at DATAMAX in 35/75 alignable frames (45%), 17 of the
# 37 in-transit frames, and above 90% of full well in 47/75. 90% because a CCD leaves
# its linear regime before full well, so clipping is the last symptom rather than the
# first; EXOTIC's own reject_overexposed_stars uses the same 0.9 fraction, which is a
# check on the number and not a coincidence I should take credit for. This gate exists
# to say NO BEFORE the reduction, since EXOTIC would discard those frames itself and
# leave a thin, unevenly sampled transit -- the WASP-52 b row-56 mechanism. One night;
# like the floor, the number moves again on the next night that bears on it.
SAT_NONLINEAR_FRAC = 0.90   # of DATAMAX: above this a frame's core is not trustworthy
SAT_FRAME_FRAC_MAX = 0.10   # reject when more than this fraction of frames are there
# KAF-1402ME as listed for MicroObservatory on science.nasa.gov/citizen-science/exoplanet-watch/how-to-contribute/how-to-submit-your-data/
MOBS_NOISE = {'gain': 53.6, 'read_noise': 15.0, 'dark': 15.0}


def say(s=''):
    print(s, flush=True)


# ---------------------------------------------------------------- 1. frames
def listing_names(target, yymmdd):
    html = urllib.request.urlopen(LISTING, timeout=90).read().decode(errors='ignore')
    pat = re.compile(r'ImageDirectory/(' + re.escape(target) + yymmdd + r'\d{6})\.FITS')
    return sorted(set(pat.findall(html)))


def download(names, ddir):
    os.makedirs(ddir, exist_ok=True)
    ctx_ok = ssl.create_default_context()
    ctx_no = ssl._create_unverified_context()
    ctx, warned, got = ctx_ok, False, 0
    for n in names:
        dest = os.path.join(ddir, n + '.FITS')
        aside = os.path.join(ddir, 'excluded', n + '.FITS')   # set aside on an earlier run; still present
        if any(os.path.exists(d) and os.path.getsize(d) > 100000 for d in (dest, aside)):
            continue
        url = FITS_URL.format(name=n)
        for attempt in range(2):
            try:
                data = urllib.request.urlopen(url, timeout=120, context=ctx).read()
                break
            except urllib.error.URLError as e:
                if 'CERTIFICATE' in str(e).upper() and ctx is ctx_ok:
                    ctx = ctx_no
                    if not warned:
                        say('  mo-www certificate did not verify (expired 2026-09-06); continuing without verification for this host')
                        warned = True
                    continue
                raise
        open(dest, 'wb').write(data)
        got += 1
    return got


# --------------------------------------------------------------- 2. archive
def archive(planet):
    cols = ('pl_name,hostname,pl_orbper,pl_orbpererr1,pl_tranmid,pl_tranmiderr1,pl_trandur,pl_ratror,pl_ratrorerr1,'
            'pl_ratdor,pl_ratdorerr1,pl_orbincl,pl_orbinclerr1,pl_orbeccen,pl_orblper,st_teff,st_tefferr1,st_tefferr2,'
            'st_met,st_meterr1,st_meterr2,st_logg,st_loggerr1,st_loggerr2,sy_dist,sy_pmra,sy_pmdec,sy_vmag,ra,dec,'
            'pl_radj,st_rad,st_mass,pl_trandep')
    q = f"select {cols} from pscomppars where pl_name='{planet}'"
    url = TAP + urllib.parse.urlencode({'query': q, 'format': 'csv'})
    rows = list(csv.DictReader(io.StringIO(_urlopen_retry(url, 60).read().decode())))
    if not rows:
        raise SystemExit(f'archive: no pscomppars row for {planet!r} (try the alias: HAT-P-10 b is WASP-11 b)')
    r = rows[0]
    f = lambda k, d=None: float(r[k]) if r.get(k) not in (None, '', 'null') else d
    # pscomppars does not always carry pl_ratror (TOI-2570 b, 2026-09-10, had radii and a
    # depth but no ratio); derive it rather than fail, and say where it came from.
    rprs, rprs_from = f('pl_ratror'), 'pl_ratror'
    if rprs is None and f('pl_radj') and f('st_rad'):
        rprs, rprs_from = f('pl_radj') * 0.10045 / f('st_rad'), 'pl_radj / st_rad'
    if rprs is None and f('pl_trandep'):
        rprs, rprs_from = math.sqrt(f('pl_trandep') / 100.0), 'sqrt(pl_trandep)'
    if rprs is None:
        raise SystemExit(f'archive: no Rp/Rs, radii or depth for {planet!r} in pscomppars')
    # Same gap on a/Rs (TOI-5300 b, 2026-09-19: no pl_ratdor, so my inits carried a
    # null and EXOTIC's own pre-flight FAILED on it while check_inits skipped it).
    # Kepler's third law from the stellar mass and radius is what EXOTIC does
    # internally; it reproduces EXOTIC's 9.763 for TOI-5300 b as 9.748.
    ars, ars_from = f('pl_ratdor'), 'pl_ratdor'
    if ars is None and f('st_mass') and f('st_rad') and f('pl_orbper'):
        a_au = (f('st_mass') * (f('pl_orbper') / 365.25) ** 2) ** (1.0 / 3.0)
        ars, ars_from = a_au / (f('st_rad') * 0.00465047), 'Kepler III from st_mass, st_rad'
    if ars is None:
        raise SystemExit(f'archive: no a/Rs and no stellar mass+radius for {planet!r} in pscomppars')
    return {
        'planet': r['pl_name'], 'host': r['hostname'],
        'P': f('pl_orbper'), 'Perr': f('pl_orbpererr1', 1e-6), 'T0': f('pl_tranmid'), 'T0err': f('pl_tranmiderr1', 1e-3),
        'T14h': f('pl_trandur'), 'rprs': rprs, 'rprs_from': rprs_from, 'rprserr': f('pl_ratrorerr1', 0.005),
        'ars': ars, 'ars_from': ars_from, 'arserr': f('pl_ratdorerr1', 0.2), 'inc': f('pl_orbincl'), 'incerr': f('pl_orbinclerr1', 0.5),
        'ecc': f('pl_orbeccen', 0.0), 'omega': f('pl_orblper', 90.0),
        'teff': f('st_teff'), 'teffp': f('st_tefferr1', 100.0), 'teffm': f('st_tefferr2', -100.0),
        'met': f('st_met', 0.0), 'metp': f('st_meterr1', 0.1), 'metm': f('st_meterr2', -0.1),
        'logg': f('st_logg', 4.5), 'loggp': f('st_loggerr1', 0.1), 'loggm': f('st_loggerr2', -0.1),
        'dist': f('sy_dist'), 'pmra': f('sy_pmra', 0.0), 'pmdec': f('sy_pmdec', 0.0), 'V': f('sy_vmag'),
        'ra': f('ra'), 'dec': f('dec'),
    }


def sexa(ra, dec):
    h = ra / 15.0
    hh = int(h); mm = int((h - hh) * 60); ss = ((h - hh) * 60 - mm) * 60
    sgn = '+' if dec >= 0 else '-'; d = abs(dec)
    dd = int(d); dm = int((d - dd) * 60); ds = ((d - dd) * 60 - dm) * 60
    return f'{hh:02d}:{mm:02d}:{ss:05.2f}', f'{sgn}{dd:02d}:{dm:02d}:{ds:04.1f}'


def ut(jd):
    return (datetime(2000, 1, 1, 12, tzinfo=timezone.utc)
            .__class__.fromtimestamp((jd - 2440587.5) * 86400, tz=timezone.utc)).strftime('%H:%M')


# --------------------------------------------------------------- 3. timing
TWILIGHT_SKY = 1000    # median sky above this is dusk or dawn, not a night frame. Calibrated on
                       # WASP-80 2026-09-08: dusk frames read 4095, 3197, 1692; the first usable
                       # frame read 718 with 85 stars; the night settled at 430-500
MIN_SEED_STARS = 10    # a seed frame needs a star field to plate-solve and to vote with
TWILIGHT_SUN_ALT = -12.0  # deg. Twilight is a statement about the SUN, not the sky: a bright sky with
                          # the Sun well below this is moonlight or lit cloud, not dusk. WASP-177 b
                          # 2026-09-25: all 74 frames read above 1000 ADU under a ~98% moon at local
                          # midnight, the sky-only rule set every one aside as 'twilight', and the
                          # tool then crashed on an empty frame list instead of giving a verdict.


def set_aside_twilight(ddir):
    """Move LEADING and TRAILING twilight frames (median sky above TWILIGHT_SKY) to <ddir>/excluded,
    stopping at the first night frame. These are not observations of anything and must
    not seed a reduction or be counted in the triage. WASP-80 2026-09-08: frames 1-3
    were dusk; the real first frame was 16 minutes in.

    A leading frame with a night-dark sky and few or no stars is NOT set aside here:
    that is cloud or an empty low field, and it belongs in the triage as a lost frame.
    WASP-50 2026-09-09: the first 41 frames (07:33-09:33 UT, two hours) had 0-9 stars
    on a 402-414 ADU sky and were set aside as "twilight" by the old star-count rule,
    which hid the whole clouded pre-ingress stretch from the triage counts.

    TRAILING (dawn) frames are set aside on the Sun alone: every frame at the end with the Sun
    at or above -12 deg, whatever its sky. TOI-2570 b 2026-09-25: the last 8 frames had the Sun
    above -12 (sky 453 -> 4095 ADU); kept, they carried a noise dip over the QC line (+29 min,
    depth 2.3% against 1.25%); dropped, the fit failed QC as the floor had predicted. The two
    dimmest of them alone were enough to restore it (prereg_toi2570_20260925.md, tests D-F).

    Returns (leading, trailing) lists of moved filenames."""
    import shutil
    from astropy.io import fits
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz, get_sun
    import astropy.units as u
    loc = EarthLocation(lat=float(SITE['lat']) * u.deg, lon=float(SITE['lon']) * u.deg, height=SITE['elev'] * u.m)
    files = sorted(f for f in os.listdir(ddir) if f.upper().endswith('.FITS'))
    lead = twilight_run(ddir, files, loc)
    trail = twilight_run(ddir, [f for f in reversed(files) if f not in lead], loc, need_bright=False)
    for f in lead + trail:
        os.makedirs(os.path.join(ddir, 'excluded'), exist_ok=True)
        shutil.move(os.path.join(ddir, f), os.path.join(ddir, 'excluded', f))
    return lead, sorted(trail)


def twilight_run(ddir, ordered, loc, need_bright=True):
    """Frames from the start of `ordered` with the Sun at or above TWILIGHT_SUN_ALT (and, when
    need_bright, a sky above TWILIGHT_SKY), stopping at the first frame that is not. Moves nothing.

    The trailing end does not need a bright sky (TOI-2570 b 2026-09-25, tests E and F): two dawn
    frames at 453 and 607 ADU with ~600 stars, Sun at -12.0 and -10.7 deg, carried the whole
    forced fit, while dropping two night frames in their place changed nothing. A dim dawn
    frame is still a rising, changing sky."""
    from astropy.io import fits
    from astropy.time import Time
    from astropy.coordinates import AltAz, get_sun
    moved = []
    for f in ordered:
        data, med, xy, _ = detect(os.path.join(ddir, f), 4.0, 8.0)
        if need_bright and med <= TWILIGHT_SKY:
            break
        # Bright, but is it twilight? Only if the Sun is near the horizon at this frame.
        try:
            t = Time(frame_time(fits.getheader(os.path.join(ddir, f))), format='jd', scale='utc')
            sun_alt = get_sun(t).transform_to(AltAz(obstime=t, location=loc)).alt.deg
        except Exception:
            sun_alt = None
        if sun_alt is None or sun_alt < TWILIGHT_SUN_ALT:
            break                       # night-time bright sky (moon, lit cloud): the triage judges it
        moved.append(f)
    return moved


def first_seedable_frame(ddir, files):
    """(index of the first frame with at least MIN_SEED_STARS detections, per-frame star counts).

    The index is None when no frame qualifies: a night with 0-2 stars in every
    frame (TrES-5 2026-09-10) is a verdict, not a plate-solving problem.
    """
    counts, first = [], None
    for i, f in enumerate(files):
        _, med, xy, _ = detect(os.path.join(ddir, f), 4.0, 8.0)
        counts.append((len(xy), med))
        if first is None and med <= TWILIGHT_SKY and len(xy) >= MIN_SEED_STARS:
            first = i
    return first, counts


def window(ddir):
    from astropy.io import fits
    files = sorted(f for f in os.listdir(ddir) if f.upper().endswith('.FITS'))
    if not files:
        return [], None, None, None     # the caller turns this into a verdict, not a traceback
    h0, h1 = fits.getheader(os.path.join(ddir, files[0])), fits.getheader(os.path.join(ddir, files[-1]))
    return files, frame_time(h0), frame_time(h1) + float(h1.get('EXPTIME', 60)) / 86400.0, h0


# ---------------------------------------------------------------- helpers
def run_tool(script, args, capture=True):
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def parse_triage(text):
    g = {}
    m = re.search(r'comp box.*?x in \[(-?\d+), (-?\d+)\]\s+y in \[(-?\d+), (-?\d+)\]', text)
    if m:
        g['box'] = tuple(int(v) for v in m.groups())
    m = re.search(r'clear reference (\d+) ADU.*?clear (\d+), partial (\d+), lost (\d+)', text)
    if m:
        g['clear_ref'], g['clear'], g['partial'], g['lost'] = (int(v) for v in m.groups())
    m = re.search(r'shift range dx \[(-?\d+), (-?\d+)\]\s+dy \[(-?\d+), (-?\d+)\]', text)
    if m:
        g['shift'] = tuple(int(v) for v in m.groups())
    for ph in ('pre', 'in', 'post'):
        m = re.search(ph + r'\s+(\d+) frames: clear (\d+), partial (\d+), lost (\d+)', text)
        if m:
            g[ph] = tuple(int(v) for v in m.groups())
    g['steps'] = len(re.findall(r'^\s+\S+: \([+-]\d+, [+-]\d+\)$', text, re.M))
    m = re.search(r'target core: peak median (\d+) of DATAMAX (\d+) \((\d+)%\); '
                  r'clipped (\d+)/(\d+), over 90% of full well (\d+)/(\d+)', text)
    if m:
        g['core'] = {'peak_median': int(m.group(1)), 'datamax': int(m.group(2)), 'pct': int(m.group(3)),
                     'clipped': int(m.group(4)), 'n': int(m.group(5)), 'nonlinear': int(m.group(6))}
    return g


def seed_above_bg(path, x, y):
    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats
    d = fits.getdata(path).astype(float)
    _, med, _ = sigma_clipped_stats(d, sigma=3)
    xi, yi = int(round(x)), int(round(y))
    return float(d[max(0, yi - 2):yi + 3, max(0, xi - 2):xi + 3].max() - med), float(med)


def choose_comps(cands, target_flux, box, lo, hi, n):
    x0, x1, y0, y1 = box
    inbox = [c for c in cands if x0 <= c['x'] <= x1 and y0 <= c['y'] <= y1]
    ratio = lambda c: c['flux'] / target_flux if target_flux > 0 else float('inf')
    good = [c for c in inbox if lo <= ratio(c) <= hi]
    if good:
        return good[:n], True
    # nothing in range: nearest in log-brightness, flagged
    inbox.sort(key=lambda c: abs(math.log10(max(ratio(c), 1e-6))))
    return inbox[:n], False


# ------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('target', help='MObs target name as in the listing, e.g. Qatar-1, WASP-2, HAT-P-10')
    p.add_argument('yymmdd', help='UT night as in the filenames, e.g. 260907')
    p.add_argument('--planet', required=True, help='archive planet name, e.g. "Qatar-1 b"')
    p.add_argument('--data-root', default=os.path.join(ROOT, 'data'))
    p.add_argument('--out-root', default=os.path.join(ROOT, 'output'))
    p.add_argument('--ncomps', type=int, default=6)
    p.add_argument('--ratio', default='0.15,8', help='comparison/target brightness range accepted')
    p.add_argument('--run', action='store_true', help='launch EXOTIC detached if the verdict is proceed')
    p.add_argument('--no-pf', action='store_true', help="skip EXOTIC's own pre-flight (exotic -pf), the slow plate-solving step")
    p.add_argument('--force', action='store_true', help='run even on a reject verdict, and carry a night with no transit in its window past step 3')
    p.add_argument('--no-download', action='store_true')
    p.add_argument('--data-dir', help='use this existing frame directory instead of data/<target>_<date>')
    a = p.parse_args()
    lo, hi = (float(v) for v in a.ratio.split(','))

    date = '20' + a.yymmdd
    iso = f'{date[:4]}-{date[4:6]}-{date[6:]}'
    slug = re.sub(r'[^a-z0-9]', '', a.planet.lower().replace(' b', '').replace(' ', ''))
    ddir = a.data_dir or os.path.join(a.data_root, f'{a.target}_{date}')
    odir = os.path.join(a.out_root, f'{slug}_{date}')
    inits_path = os.path.join(ROOT, f'inits_{slug}_{date}.json')
    prereg_path = os.path.join(ROOT, f'prereg_{slug}_{date}.md')
    log_path = os.path.join(ROOT, f'{slug}_{date[4:]}_run.log')

    # 1. frames
    say(f'== {a.planet}: MObs {a.target} night {iso}')
    if not a.no_download and not a.data_dir:
        names = listing_names(a.target, a.yymmdd)
        say(f'  listing: {len(names)} frames for {a.target}{a.yymmdd}')
        if not names and not os.path.isdir(ddir):
            raise SystemExit('  nothing to do')
        got = download(names, ddir)
        say(f'  downloaded {got} new frame(s) to {ddir}')
    lead, trail = set_aside_twilight(ddir)
    for moved, which in ((lead, 'leading'), (trail, 'trailing')):
        if moved:
            say(f'  set aside {len(moved)} {which} twilight frame(s) (' + (f'sky above {TWILIGHT_SKY} ADU, ' if which == 'leading' else '') + f'Sun above {TWILIGHT_SUN_ALT} deg) to excluded/: {moved[0]}' + (f' .. {moved[-1]}' if len(moved) > 1 else ''))
    files, jd0, jd1, h0 = window(ddir)
    if not files:
        raise SystemExit('  no usable frames')
    say(f'  {len(files)} frames on disk, {ut(jd0)}-{ut(jd1)} UT, {h0.get("EXPTIME")} s {h0.get("FILTER")}')

    # 2. archive
    ar = archive(a.planet)
    depth = 100 * ar['rprs'] ** 2
    # Cache the archive ephemeris for publish.ts (its vsArchiveMin column) if this planet is not cached yet.
    # Never overwrite: existing entries may be hand-checked. Added 2026-09-26: 6 of 38 targets had no entry
    # because nothing wrote the cache automatically (Gaia-2 b was the one that showed it).
    try:
        ep = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ephemerides.json')
        cache = json.load(open(ep)) if os.path.exists(ep) else {}
        if ar['planet'] not in cache:
            cache[ar['planet']] = {'P': ar['P'], 'Perr': ar['Perr'], 'T0': ar['T0'], 'T0err': ar['T0err'], 'T14h': ar['T14h'],
                                   'depthPct': round(depth, 3), 'fetched': datetime.now(timezone.utc).strftime('%Y-%m-%d'),
                                   'rprs': ar['rprs'], 'source': 'NASA Exoplanet Archive pscomppars', 'vmag': ar['V']}
            with open(ep + '.new', 'w') as fh:
                json.dump(cache, fh, indent=1, sort_keys=True)
            os.replace(ep + '.new', ep)
    except Exception as e:
        say(f'  (ephemeris cache not updated: {e})')
    say(f'  archive: P {ar["P"]:.8f} d, T0 {ar["T0"]:.6f}, T14 {ar["T14h"]:.3f} h, Rp/Rs {ar["rprs"]:.4f} (depth {depth:.2f}%), V {ar["V"]}' + ('' if ar['rprs_from'] == 'pl_ratror' else f' (Rp/Rs derived from {ar["rprs_from"]}; no pl_ratror in pscomppars)')
        + ('' if ar['ars_from'] == 'pl_ratdor' else f' (a/Rs {ar["ars"]:.3f} derived: {ar["ars_from"]}; no pl_ratdor in pscomppars)'))

    # 3. timing
    n = round(((jd0 + jd1) / 2 - ar['T0']) / ar['P'])
    tmid = ar['T0'] + n * ar['P']
    half = ar['T14h'] / 48.0
    ing, egr = tmid - half, tmid + half
    pre_min, post_min = (ing - jd0) * 1440, (jd1 - egr) * 1440
    eph_err_min = (ar['T0err'] + abs(n) * ar['Perr']) * 1440
    say(f'  epoch {n}: Tmid {tmid:.5f} = {ut(tmid)} UT, ingress {ut(ing)}, egress {ut(egr)}; '
        f'baseline pre {pre_min:.0f} min, post {post_min:.0f} min; ephemeris bar {eph_err_min:.1f} min')
    if egr < jd0 or ing > jd1:
        if not a.force:
            raise SystemExit('  no transit inside the window; nothing to reduce')
        say('  NO TRANSIT inside the window by the archive ephemeris; continuing under --force '
            '(a fit here is an extrapolation and the verdict will say so)')
        coverage = 'none'
    else:
        coverage = 'full' if ing > jd0 and egr < jd1 else ('ingress-only' if ing > jd0 else 'egress-only')
    geometry_note = 'NO transit in the window. ' if coverage == 'none' else f'{coverage} transit. '

    # 4. locate: frame 1, or the first of the next few that solves, with the pixel
    #    carried back to frame 1 by the star-pair vote (EXOTIC seeds from frame 1)
    first = os.path.join(ddir, files[0])
    start, counts = first_seedable_frame(ddir, files)
    if start is None:
        stars = [c[0] for c in counts]; skies = [c[1] for c in counts]
        say(f'  no frame has {MIN_SEED_STARS} stars on a night-dark sky: {len(files)} frames, '
            f'{min(stars)}-{max(stars)} stars (median {int(np.median(stars))}), sky {min(skies):.0f}-{max(skies):.0f} ADU/px; '
            f'nothing to seed from, so nothing to solve')
        # Name the actual reason. WASP-177 b 2026-09-25 had a median of 14 stars (above the
        # threshold) on a 1978-2430 ADU moonlit sky, and the old clause called it 'star-poor'.
        dark = [c for c in counts if c[1] <= TWILIGHT_SKY]
        why = (f'every frame star-poor ({min(stars)}-{max(stars)} stars)' if dark else
               f'no night-dark frame to seed from: sky {min(skies):.0f}-{max(skies):.0f} ADU/px all night with the Sun down '
               f'(moonlight or lit cloud), {min(stars)}-{max(stars)} stars')
        say(f'== verdict: REJECT: {why}; '
            f'transit window {ut(ing)}-{ut(egr)} UT inside a {ut(jd0)}-{ut(jd1)} night')
        raise SystemExit(1)
    if start:
        say(f'  first {start} frame(s) have fewer than {MIN_SEED_STARS} stars on a night-dark sky (cloud or an empty field); '
            f'seeding from frame {start + 1} ({files[start]}) and grading them in the triage')
    loc, solved_on = None, None
    for i in range(start, min(start + 6, len(files))):
        rc, out, err = run_tool('locate_target.py', [os.path.join(ddir, files[i]), '--ra', str(ar['ra']), '--dec', str(ar['dec']), '--json', '--comps', '14'])
        if rc == 0:
            loc, solved_on = json.loads(out), i
            break
    if loc is None:
        raise SystemExit(f'  no frame among {start + 1}-{min(start + 6, len(files))} would plate-solve; last error: {err.strip()[-200:]}')
    tx, ty = loc['target']['x'], loc['target']['y']
    ref_frame = 0
    if solved_on:
        _, _, xy0, fl0 = detect(first, 4.0, 8.0)
        _, _, xyi, fli = detect(os.path.join(ddir, files[solved_on]), 4.0, 8.0)
        ref = xy0[np.argsort(fl0)[::-1][:25]] if len(xy0) else xy0
        cur = xyi[np.argsort(fli)[::-1][:25]] if len(xyi) else xyi
        sh, votes = voted_shift(cur, ref)
        if sh is None:
            # Frame 1 cannot be aligned (too few stars). Keep the pixel on the solved
            # frame, triage relative to it, and set the unseedable leading frames aside
            # AFTER the triage so they are counted as lost rather than hidden.
            ref_frame = solved_on
            say(f'  frame 1 could not be aligned to frame {solved_on + 1} (votes {votes}); the pixel refers to frame {solved_on + 1}')
        else:
            tx, ty = tx - float(sh[0]), ty - float(sh[1])
            for c in loc['comparisons']:
                c['x'], c['y'] = c['x'] - float(sh[0]), c['y'] - float(sh[1])
            say(f'  frame 1 did not solve; solved frame {solved_on + 1} ({files[solved_on]}) and carried the pixel back by ({-sh[0]:+.1f}, {-sh[1]:+.1f}) on {votes} votes')
    seed_file = os.path.join(ddir, files[ref_frame])
    above, bg = seed_above_bg(seed_file, tx, ty)
    say(f'  target in frame {ref_frame + 1}: ({tx:.1f}, {ty:.1f}), {above:.0f} ADU above a {bg:.0f} background'
        + ('  SATURATED' if loc['target']['saturated'] else ''))

    # 5. triage
    rc, out, err = run_tool('night_triage.py', [ddir, '--x', str(tx), '--y', str(ty), '--tmid', f'{tmid:.5f}', '--t14', f'{ar["T14h"]:.4f}',
                                                '--ref-frame', str(ref_frame)])
    if rc != 0:
        raise SystemExit(f'  night_triage failed: {err.strip()[-300:]}')
    open(os.path.join(ddir, 'triage.txt'), 'w').write(out)
    tg = parse_triage(out)
    say('  triage: ' + '; '.join(l for l in out.strip().splitlines() if l.startswith(('frames', 'shift', 'target flux'))))
    # The seed floor is measured on the seed frame; say what it would be on the clearest frames, so a
    # verdict on a target AT the floor (WASP-52 09-13: 197 ADU in a 76%-transmission frame 1, ~260 clear)
    # is read as one, and a target UNDER it only by cloud is not mistaken for a faint one.
    seed_clear = None
    mref = re.search(r'clear reference (\d+) ADU', out)
    mseed = re.search(re.escape(files[ref_frame]) + r'\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+(\d+)', out)
    if mref and mseed and int(mseed.group(1)) > 0:
        frac = int(mseed.group(1)) / int(mref.group(1))
        seed_clear = above / frac if frac > 0 else None
        if seed_clear is not None and abs(frac - 1) > 0.05:
            say(f'  seed frame holds {frac:.0%} of the clear reference flux; the seed target would be ~{seed_clear:.0f} ADU above background on the clearest frames')
    for ph in ('pre', 'in', 'post'):
        if ph in tg:
            say(f'    {ph:4s} {tg[ph][0]:3d} frames: clear {tg[ph][1]}, partial {tg[ph][2]}, lost {tg[ph][3]}')
    if tg['steps']:
        say(f'    pointing steps > 25 px: {tg["steps"]}')
    if ref_frame:
        # Now that they have been counted, move the leading frames EXOTIC could not
        # seed from out of the way, so a reduction (if any) starts at the seed frame.
        import shutil
        os.makedirs(os.path.join(ddir, 'excluded'), exist_ok=True)
        for f in files[:ref_frame]:
            shutil.move(os.path.join(ddir, f), os.path.join(ddir, 'excluded', f))
        say(f'  set aside the {ref_frame} leading star-poor frame(s) to excluded/ for seeding only; they are counted above as lost')
        files = files[ref_frame:]

    # Leading clouded frames (2026-09-14, Qatar-1 09-14): frame 1 had 33 stars, so it
    # counted as seedable, yet the target held 5% of its clear flux (23 ADU against
    # ~446 clear): thin cloud that leaves stars detectable and a faint target gone.
    # EXOTIC seeds from the first file, so a target invisible there fails the
    # reduction regardless of the night. When the seed is under the floor only by
    # cloud, set aside the LEADING run of frames the triage grades lost (they carry
    # no photometry anyway, and they stay counted above), re-reference the seed
    # pixel to the first surviving frame by the triage's own shift, and measure the
    # floor there. Before this the tool only SAID "the seed frame is the problem".
    # Per-frame rows from the triage table: (file, dx, dy, tflux, grade), shifts relative to the seed frame.
    rows = []
    for line in out.splitlines():
        m = re.match(r'(\S+\.FITS)\s+\S+\s+\d+\s+(\S+)\s+(\S+)\s+\d+\s+(\S+)\s+\d+\s+(clear|partial|lost)$', line.strip())
        if m:
            num = lambda v: float(v) if v != '?' else None
            rows.append((m.group(1), num(m.group(2)), num(m.group(3)), num(m.group(4)), m.group(5)))
    if above < SEED_MIN_ADU and seed_clear is not None and seed_clear >= SEED_MIN_ADU:
        lead = 0
        while lead < len(rows) - 1 and rows[lead][4] == 'lost':
            lead += 1
        if lead and rows[lead][1] is not None:
            import shutil
            os.makedirs(os.path.join(ddir, 'excluded'), exist_ok=True)
            for name, *_ in rows[:lead]:
                if os.path.exists(os.path.join(ddir, name)):
                    shutil.move(os.path.join(ddir, name), os.path.join(ddir, 'excluded', name))
            dx, dy = rows[lead][1], rows[lead][2]
            tx, ty = tx + dx, ty + dy
            for c in loc['comparisons']:
                c['x'], c['y'] = c['x'] + dx, c['y'] + dy
            gone = {r[0] for r in rows[:lead]}
            files = [f for f in files if f not in gone]
            rows = [(r[0], None if r[1] is None else r[1] - dx, None if r[2] is None else r[2] - dy, r[3], r[4]) for r in rows[lead:]]
            seed_file = os.path.join(ddir, files[0]); first = seed_file; ref_frame = 0
            above, bg = seed_above_bg(seed_file, tx, ty)
            say(f'  set aside the {lead} leading lost frame(s) (stars present, target extinguished) to excluded/; '
                f'seeding from {files[0]} at ({tx:.1f}, {ty:.1f}), {above:.0f} ADU above a {bg:.0f} background')

    # 6. comps + inits
    from astropy.io import fits
    from photutils.aperture import CircularAperture, aperture_photometry
    # Ratios on the CLEAREST frame (2026-09-14): measured on a clouded seed frame the
    # target reads dim and a 25-ADU nothing passes as a "0.3x comparison". The triage
    # table names the frame with the most target flux and its shift from the seed.
    from astropy.stats import sigma_clipped_stats
    best = max((r for r in rows if r[3] is not None and r[1] is not None), key=lambda r: r[3], default=None)
    cx_off, cy_off, ratio_frame = (best[1], best[2], os.path.join(ddir, best[0])) if best else (0.0, 0.0, seed_file)
    d0 = fits.getdata(ratio_frame).astype(float)
    _, rbg, _ = sigma_clipped_stats(d0, sigma=3)
    tflux = float(aperture_photometry(d0 - rbg, CircularAperture((tx + cx_off, ty + cy_off), r=5))['aperture_sum'][0])
    for c in loc['comparisons']:
        c['flux'] = float(aperture_photometry(d0 - rbg, CircularAperture((c['x'] + cx_off, c['y'] + cy_off), r=5))['aperture_sum'][0])
    if best and ratio_frame != seed_file:
        # The seed-floor estimate for the clearest frames was a fraction of the seed
        # reading; measure it directly here instead (the fraction swung 446 -> 166 on
        # Qatar-1 09-14 when the seed frame changed).
        seed_clear, _ = seed_above_bg(ratio_frame, tx + cx_off, ty + cy_off)
        say(f'  comparison ratios measured on the clearest frame, {best[0]} (target {tflux:.0f} ADU in 5 px, {seed_clear:.0f} ADU peak above background)')
    box = tg.get('box', (16, 634, 16, 484))
    comps, in_range = choose_comps(loc['comparisons'], tflux, box, lo, hi, a.ncomps)
    say(f'  comps: {len(comps)} chosen in box x[{box[0]},{box[1]}] y[{box[2]},{box[3]}]'
        + ('' if in_range else f'  NONE within {lo}-{hi}x of the target; nearest taken, flagged'))
    for c in comps:
        say(f'    ({c["x"]:6.1f}, {c["y"]:6.1f})  {c["flux"] / tflux if tflux > 0 else float("nan"):6.2f}x target')

    ra_s, dec_s = sexa(ar['ra'], ar['dec'])
    inits = {
        'user_info': {
            'Directory with FITS files': ddir, 'Directory to Save Plots': odir,
            'Directory of Flats': None, 'Directory of Darks': None, 'Directory of Biases': None,
            'AAVSO Observer Code (blank if none)': '', 'Secondary Observer Codes (blank if none)': 'MOBS',
            'Observation date': iso, 'Obs. Latitude': SITE['lat'], 'Obs. Longitude': SITE['lon'], 'Obs. Elevation (meters)': SITE['elev'],
            'Camera Type (CCD or DSLR)': 'CCD', 'Pixel Binning': '2x2', 'Filter Name (aavso.org/filters)': 'CV',
            'Observing Notes': (f'MicroObservatory public Image Directory. {len(files)} frames {ut(jd0)}-{ut(jd1)} UT. '
                                f'Triage (tools/night_triage.py): clear {tg.get("clear")}, partial {tg.get("partial")}, lost {tg.get("lost")}; '
                                f'shift range dx {tg.get("shift", ("?",) * 4)[:2]} dy {tg.get("shift", ("?",) * 4)[2:]} px. '
                                f'Seed from a local plate solve of frame {solved_on + 1 if solved_on else 1}' + (f' ({ref_frame} leading star-poor frames set aside after triage)' if ref_frame else '') + '; comps chosen to stay on the chip through the full shift range. '
                                f'Archive ephemeris Tmid {tmid:.5f} ({ut(tmid)} UT), {geometry_note}'
                                f'Pre-registered (prereg_{slug}_{date}.md). Reduced by Opus (AI); not submitted.'),
            'Plate Solution? (y/n)': 'n', 'Add Comparison Stars from AAVSO? (y/n)': 'n',
            'Target Star X & Y Pixel': [int(round(tx)), int(round(ty))],
            'Comparison Star(s) X & Y Pixel': [[int(round(c['x'])), int(round(c['y']))] for c in comps],
            'Demosaic Format': None, 'Demosaic Output': None,
        },
        'planetary_parameters': {
            'Target Star RA': ra_s, 'Target Star Dec': dec_s, 'Planet Name': ar['planet'], 'Host Star Name': ar['host'],
            'Orbital Period (days)': ar['P'], 'Orbital Period Uncertainty': ar['Perr'],
            'Published Mid-Transit Time (BJD-UTC)': ar['T0'], 'Mid-Transit Time Uncertainty': ar['T0err'],
            'Ratio of Planet to Stellar Radius (Rp/Rs)': ar['rprs'], 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty': ar['rprserr'],
            'Ratio of Distance to Stellar Radius (a/Rs)': ar['ars'], 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty': ar['arserr'],
            'Orbital Inclination (deg)': ar['inc'], 'Orbital Inclination (deg) Uncertainty': ar['incerr'],
            'Orbital Eccentricity (0 if null)': ar['ecc'], 'Argument of Periastron (deg)': ar['omega'],
            'Star Effective Temperature (K)': ar['teff'], 'Star Effective Temperature (+) Uncertainty': ar['teffp'], 'Star Effective Temperature (-) Uncertainty': ar['teffm'],
            'Star Metallicity ([FE/H])': ar['met'], 'Star Metallicity (+) Uncertainty': ar['metp'], 'Star Metallicity (-) Uncertainty': ar['metm'],
            'Star Surface Gravity (log(g))': ar['logg'], 'Star Surface Gravity (+) Uncertainty': ar['loggp'], 'Star Surface Gravity (-) Uncertainty': ar['loggm'],
            'Star Distance (pc)': ar['dist'], 'Star Proper Motion RA (mas/yr)': ar['pmra'], 'Star Proper Motion DEC (mas/yr)': ar['pmdec'],
        },
        'optional_info': {
            'Pre-reduced File:': '', 'Pre-reduced File Time Format (BJD_TDB, JD_UTC, MJD_UTC)': 'BJD_TDB',
            'Pre-reduced File Units of Flux (flux, magnitude, millimagnitude)': 'flux',
            'Filter Minimum Wavelength (nm)': 350, 'Filter Maximum Wavelength (nm)': 850,
            'Image Scale (Ex: 5.21 arcsecs/pixel)': None, 'Exposure Time (s)': float(h0.get('EXPTIME', 60)),
            # MicroObservatory noise budget (KAF-1402ME), from NASA's Exoplanet Watch 'How to Submit Your Data'
            # page. MObs headers carry no gain key, and EXOTIC's default is 1 e-/ADU, which puts every ADU in
            # as one photon and inflates the per-point error bars (2026-09-13: WASP-52 b per-point 2.65% at
            # gain 1 against a 1.81% residual scatter).
            'gain_electrons_per_adu': MOBS_NOISE['gain'], 'read_noise_electrons': MOBS_NOISE['read_noise'],
            'dark_current_electrons_per_second_per_pixel': MOBS_NOISE['dark'],
        },
    }
    if os.path.exists(inits_path):
        inits_path = inits_path.replace('.json', '.auto.json')
        say('  an inits file for this night already exists; writing the generated one beside it')
    json.dump(inits, open(inits_path, 'w'), indent=2)
    say(f'  wrote {os.path.relpath(inits_path, ROOT)}')

    # 7. pre-flight: this directory's check, then EXOTIC's own (exotic -pf, upstream
    #    since 2026-09-10: archive agreement, transit in window, seed pixel against a
    #    plate solution, comps on-frame and unsaturated). The second one plate-solves
    #    a frame, so it is the slow part; skip it with --no-pf.
    retained = sum(1 for f in os.listdir(ddir) if f.upper().endswith('.FITS'))
    aside = len(os.listdir(os.path.join(ddir, 'excluded'))) if os.path.isdir(os.path.join(ddir, 'excluded')) else 0
    if aside:
        say(f'  pre-flight runs on the {retained} retained frame(s); the {aside} in excluded/ are outside its window, '
            'so a transit-in-window FAIL below can be the set-aside, not the sky')
    rc, out, err = run_tool('check_inits.py', [inits_path])
    say('  check_inits: ' + ('PASS' if rc == 0 else f'FAIL (exit {rc})'))
    for line in out.strip().splitlines():
        if '[FAIL]' in line or '[look]' in line:
            say('    ' + line.strip())
    pf_rc = 0
    if not a.no_pf:
        exotic_bin = os.path.join(ROOT, 'venv', 'bin', 'exotic')
        r = subprocess.run([exotic_bin, '-pf', inits_path], capture_output=True, text=True, cwd=ROOT)
        pf_rc = r.returncode
        say('  exotic -pf: ' + ('PASS' if pf_rc == 0 else f'FAIL (exit {pf_rc})'))
        for line in (r.stdout + r.stderr).splitlines():
            if re.search(r'\[(FAIL|look|skip)', line):
                say('    ' + line.strip())

    if ar['V']:
        scatter = CAL_SCATTER + CAL_SLOPE * (min(ar['V'], CAL_V_KNEE) - CAL_V)
        if ar['V'] > CAL_V_KNEE:
            scatter *= 10 ** (0.2 * (ar['V'] - CAL_V_KNEE))
    else:
        scatter = None
    bar = (scatter / depth) / CAL_RATIO * CAL_BAR_MIN if scatter else None

    # 8. verdict
    reasons = []
    if above < SEED_MIN_ADU:
        r = f'seed target {above:.0f} ADU above background (< {SEED_MIN_ADU})'
        if seed_clear is not None and seed_clear >= SEED_MIN_ADU:
            r += f'; ~{seed_clear:.0f} on the clearest frames, so the seed frame is the problem, not the target'
            # TOI-2570 b 2026-09-25: this line said '--force is reasonable' unconditionally. Forced,
            # a 1.25% transit at scatter/depth ~0.9 came back +29 min late and twice too deep with a
            # QC MARGINAL: a noise dip, confidently fitted. Forcing is only worth it when the night
            # could time the transit if it were clear (scatter/depth <= FORCE_MAX_RATIO predicted).
            if scatter and scatter / depth <= FORCE_MAX_RATIO:
                r += f': --force is reasonable (predicted scatter/depth {scatter / depth:.2f})'
            else:
                r += (f': --force would be a test, not a measurement (predicted scatter/depth '
                      + (f'{scatter / depth:.2f}' if scatter else 'unknown, no V') + f' > {FORCE_MAX_RATIO})')
        elif seed_clear is not None:
            r += f'; ~{seed_clear:.0f} on the clearest frames'
        reasons.append(r)
    if 'in' in tg and tg['in'][0] and tg['in'][3] / tg['in'][0] > 0.5:
        reasons.append(f'in-transit frames lost {tg["in"][3]}/{tg["in"][0]}')
    if 'in' in tg and (tg['in'][1] + tg['in'][2]) < 10:
        reasons.append(f'only {tg["in"][1] + tg["in"][2]} usable in-transit frames')
    if not in_range:
        reasons.append(f'no comparison within {lo}-{hi}x of the target')
    core = tg.get('core')
    if core and core['n'] and core['nonlinear'] / core['n'] > SAT_FRAME_FRAC_MAX:
        reasons.append(f"target core saturated: {core['clipped']}/{core['n']} frames clipped at "
                       f"DATAMAX {core['datamax']}, {core['nonlinear']}/{core['n']} above "
                       f"{SAT_NONLINEAR_FRAC:.0%} of full well (peak median {core['pct']}%)")
    if rc != 0:
        reasons.append('check_inits failed')
    if pf_rc != 0:
        reasons.append('exotic -pf failed')
    verdict = 'REJECT' if reasons else 'PROCEED'
    say(f'== verdict: {verdict}' + (': ' + '; '.join(reasons) if reasons else ''))
    # Night record for observatory.opusgarden.dev (2026-09-20): what the triage table
    # does not hold. publish.ts pairs it with triage.txt.
    json.dump({'verdict': verdict, 'clauses': reasons,
               'seed': {'file': os.path.basename(seed_file), 'x': round(tx, 1), 'y': round(ty, 1), 'above': round(above), 'bg': round(bg),
                        'clear': round(seed_clear) if seed_clear is not None else None},
               'floor': SEED_MIN_ADU, 'ratio_range': [lo, hi], 'in_range': bool(in_range),
               'comps': [{'x': round(c['x'], 1), 'y': round(c['y'], 1), 'ratio': round(c.get('ratio', 0), 2) if c.get('ratio') is not None else None} for c in comps],
               'window': {'ingress_ut': ut(ing), 'mid_ut': ut(tmid), 'egress_ut': ut(egr), 'coverage': coverage},
               'frames': {k: list(v) for k, v in tg.items() if k in ('pre', 'in', 'post')} if tg else {},
               'expected': {'scatter_pct': round(scatter, 2) if scatter else None, 'bar_min': round(bar, 1) if bar else None},
               'written': datetime.now(timezone.utc).isoformat(timespec='minutes')},
              open(os.path.join(ddir, 'night.json'), 'w'), indent=1)
    if bar:
        say(f'   derived expectation if clear: scatter ~{scatter:.2f}%, Tmid bar ~{bar:.0f} min')

    # 9. prereg scaffold (PROCEED only; never overwrite)
    if verdict == 'REJECT':
        say('  no prereg scaffold written for a rejected night')
    elif not os.path.exists(prereg_path):
        open(prereg_path, 'w').write(f"""# Pre-registration: {ar['planet']}, MicroObservatory night of {iso} UT

SCAFFOLD from tools/mobs_night.py. Edit the judgment lines, then commit BEFORE any fit.

Frames: {len(files)}, {ut(jd0)} to {ut(jd1)} UT, {h0.get('EXPTIME')} s {h0.get('FILTER')}.
Archive (NEA pscomppars, fetched {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC): P {ar['P']:.8f} +/- {ar['Perr']:.1e} d,
T0 {ar['T0']:.6f} +/- {ar['T0err']:.1e}, T14 {ar['T14h']:.3f} h, Rp/Rs {ar['rprs']:.4f} (depth {depth:.2f}%), a/Rs {ar['ars']},
inc {ar['inc']}, e {ar['ecc']}, Teff {ar['teff']}, V {ar['V']}. Epoch {n}: predicted Tmid {tmid:.5f} ({ut(tmid)} UT),
ingress {ut(ing)}, egress {ut(egr)}; {pre_min:.0f} min pre-ingress and {post_min:.0f} min post-egress in the window
({coverage}). Propagated ephemeris bar {eph_err_min:.1f} min.

Predictions:
- Sky judged by stars per frame and the target's aperture flux (night_triage), not the background.
  If the counts hold through {ut(ing)}-{ut(egr)} and the seed target is >{SEED_MIN_ADU} ADU above background:
  QC PASS or MARGINAL (not predicting which).   [EDIT: state which if there is a reason]
- Scatter, from the V-magnitude calibration ({CAL_SCATTER}% at V {CAL_V} on 2026-09-05, +{CAL_SLOPE}%/mag from three clean nights, revised 2026-09-18):
  ~{scatter:.2f}%. On a {depth:.2f}% depth that is scatter/depth {scatter / depth:.2f}; with {CAL_RATIO} giving
  {CAL_BAR_MIN} min at ~3 min cadence, expected Tmid bar ~{bar:.0f} min (range {0.7 * bar:.0f}-{1.5 * bar:.0f}; two clean nights gave 9.6 and 8.6 against 6.7 and 9.0 predicted, so the bar model is two-point-crude).
  Under {0.6 * bar:.0f} min: calibration wrong the other way. Over {1.8 * bar:.0f} min: worse night than the counts show.
- Tmid within 1 sigma (that bar) of {tmid:.5f}.   [EDIT: what a miss would and would not mean]
- Depth within 1 sigma of {depth:.2f}%; Rp/Rs within 1 sigma of {ar['rprs']:.4f}.
- If stars-per-frame falls by more than half at any point between {ut(ing)} and {ut(egr)}, or the target is under
  ~{SEED_MIN_ADU} ADU above background in the seed frame: expect FAIL or a rejected night, and the record will say so
  with no depth or timing claim.
- Any printed bar under 0.001 d is a replaced bar (#1401) and gets refit before the ledger.
- First-frame coordinates from a local plate solve; drift by star-pair voting; comps in the triage box within a
  brightness factor of the target; #1409 confirms the rest.
- Not submitted regardless; checkpoint holds.
""")
        say(f'  wrote {os.path.relpath(prereg_path, ROOT)} (scaffold; EDIT and commit before any fit)')
    else:
        say(f'  prereg exists, left alone: {os.path.relpath(prereg_path, ROOT)}')

    # --run
    if a.run:
        if verdict == 'REJECT' and not a.force:
            say('   --run refused on a REJECT verdict (use --force to override)')
            return 1
        script = os.path.join(ROOT, f'{slug}_{date[4:]}_run.sh')
        rel = os.path.relpath(inits_path, ROOT)
        open(script, 'w').write(f"""#!/bin/bash
cd {ROOT}
echo "=== {ar['planet']} {date} run start $(date -u +%FT%TZ) ===" >> {os.path.basename(log_path)}
venv/bin/exotic -red {rel} -ov >> {os.path.basename(log_path)} 2>&1
echo "=== {ar['planet']} {date} run end $(date -u +%FT%TZ) exit=$? ===" >> {os.path.basename(log_path)}
venv/bin/python tools/post_run_check.py {rel} >> {os.path.basename(log_path)} 2>&1
echo "=== post_run_check exit=$? ===" >> {os.path.basename(log_path)}
""")
        os.chmod(script, 0o755)
        subprocess.Popen(['setsid', 'nohup', script], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
        say(f'   launched detached: {os.path.basename(script)} -> {os.path.basename(log_path)}; '
            f'wait for "post_run_check exit" in the log')
    return 0 if verdict == 'PROCEED' else 1


if __name__ == '__main__':
    sys.exit(main())
