#!/usr/bin/env python3
"""
Triage a night of MicroObservatory frames before writing an inits file.

Built 2026-09-06 after two nights in a row taught the same two lessons by hand:

1. The sky BACKGROUND level is not a cloud detector at night. WASP-2 on 09-06
   held 400-423 ADU/px all night while thin cloud cut the target's flux by up to
   99% from 04:09 UT; I had called the sky "stable" from the background alone in
   the pre-registration. What detects cloud is the number of stars found per
   frame and the target's own aperture flux.
2. Pointing on this telescope moves by tens to a hundred pixels in a night,
   sometimes in steps and sometimes with an intermittent nod, and both phase
   correlation and centroid-following gave wrong answers on WASP-11 (09-05).
   Star-list translation voting against the first frame is what worked.

Given the target's pixel position in the FIRST frame (from tools/locate_target.py;
EXOTIC seeds from the first frame), this prints per frame: UT time, stars found,
the voted shift from frame 1, the target's background-subtracted flux in a 4 px
aperture at the shifted position, and the sky. Then a summary: the shift range,
the first-frame box a comparison star must sit in to stay on the chip with an
aperture margin through the whole night, the clear/partial/lost frame counts
(relative to the median flux of the clearest quarter), and, if --tmid and --t14
are given, those counts split by transit phase.

    python3 tools/night_triage.py DIR --x 224 --y 229 [--tmid 2461289.7007 --t14 1.76]
"""
import argparse, glob, warnings
from datetime import datetime, timezone
warnings.filterwarnings('ignore')
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry

APERTURE_MARGIN_PX = 16   # widest aperture EXOTIC will try, plus slack


def detect(path, fwhm, nsigma):
    data = fits.getdata(path).astype(float)
    _, med, sd = sigma_clipped_stats(data, sigma=3)
    found = DAOStarFinder(fwhm=fwhm, threshold=nsigma * sd)(data - med)
    if found is None:
        return data, med, np.zeros((0, 2)), np.zeros(0)
    xy = np.array([found['xcentroid'], found['ycentroid']]).T
    return data, med, xy, np.array(found['flux'])


def voted_shift(xy, ref, bin_px=3.0):
    """Translation that the most star pairs agree on (mode of pairwise offsets)."""
    if len(xy) < 5 or len(ref) < 5:
        return None, 0
    diffs = (xy[:, None, :] - ref[None, :, :]).reshape(-1, 2)
    keys = np.round(diffs / bin_px).astype(int)
    uniq, counts = np.unique(keys, axis=0, return_counts=True)
    k = np.argmax(counts)
    sel = np.all(keys == uniq[k], axis=1)
    return diffs[sel].mean(0), int(counts[k])


def frame_time(hdr):
    """Frame start as a UT JD. MObs writes DATE-OBS in local time with an offset
    ('2026-09-05T19:49:31.882-0700'); a bare timestamp is taken as UT."""
    raw = hdr.get('DATE-OBS') or hdr.get('DATE_OBS')
    if not raw:
        return np.nan
    try:
        dt = datetime.fromisoformat(raw.replace('Z', '+00:00'))
    except ValueError:
        return np.nan
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
    return Time(dt.isoformat(), format='isot').jd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('directory')
    p.add_argument('--x', type=float, required=True, help='target x in the reference frame (the first, unless --ref-frame)')
    p.add_argument('--y', type=float, required=True, help='target y in the reference frame (the first, unless --ref-frame)')
    p.add_argument('--ref-frame', type=int, default=0,
                   help='0-based index of the frame the target pixel and the shifts refer to. Use it when the '
                        'first frame has too few stars to vote with (cloud or an empty low field), so the '
                        'star-poor leading frames are graded lost instead of aborting the triage')
    p.add_argument('--tmid', type=float, help='predicted mid-transit, JD')
    p.add_argument('--t14', type=float, help='transit duration, hours')
    p.add_argument('--fwhm', type=float, default=4.0)
    p.add_argument('--nsigma', type=float, default=8.0)
    p.add_argument('--ref-stars', type=int, default=25, help='brightest stars used for voting')
    p.add_argument('--quiet', action='store_true', help='summary only')
    a = p.parse_args()

    files = sorted(glob.glob(f'{a.directory}/*.FITS') + glob.glob(f'{a.directory}/*.fits') + glob.glob(f'{a.directory}/*.fit'))
    if not files:
        raise SystemExit('no frames found')
    if not 0 <= a.ref_frame < len(files):
        raise SystemExit(f'--ref-frame {a.ref_frame} is outside the {len(files)} frames found')
    data0, _, xy0, fl0 = detect(files[a.ref_frame], a.fwhm, a.nsigma)
    ref = xy0[np.argsort(fl0)[::-1][:a.ref_stars]]
    if a.ref_frame:
        print(f'reference frame {a.ref_frame + 1} ({files[a.ref_frame].split("/")[-1]}); shifts are relative to it')
    shape = data0.shape
    tgt = np.array([a.x, a.y])
    rows = []
    for f in files:
        data, med, xy, fl = detect(f, a.fwhm, a.nsigma)
        jd = frame_time(fits.getheader(f))
        top = xy[np.argsort(fl)[::-1][:a.ref_stars]] if len(xy) else xy
        sh, votes = voted_shift(top, ref)
        if sh is None:
            rows.append((f, jd, len(xy), None, 0, np.nan, med))
            continue
        pos = tgt + sh
        flux = aperture_photometry(data - med, CircularAperture(pos, r=4))['aperture_sum'][0]
        rows.append((f, jd, len(xy), sh, votes, float(flux), med))

    fluxes = np.array([r[5] for r in rows])
    good = fluxes[np.isfinite(fluxes)]
    if len(good) == 0:
        raise SystemExit('no frame could be aligned to the reference frame')
    clear_ref = np.median(np.sort(good)[-max(3, len(good) // 4):])

    def grade(fl):
        if not np.isfinite(fl):
            return 'lost'
        if fl >= 0.75 * clear_ref:
            return 'clear'
        if fl >= 0.45 * clear_ref:
            return 'partial'
        return 'lost'

    if not a.quiet:
        print(f'{"frame":32s} {"UT":8s} {"stars":>5s} {"dx":>7s} {"dy":>7s} {"votes":>5s} {"tflux":>7s} {"sky":>5s} grade')
        for f, jd, n, sh, votes, flux, med in rows:
            ut = Time(jd, format='jd').iso[11:19] if np.isfinite(jd) else '?'
            dx, dy = (f'{sh[0]:7.1f}', f'{sh[1]:7.1f}') if sh is not None else ('      ?', '      ?')
            fl = f'{flux:7.0f}' if np.isfinite(flux) else '      ?'
            print(f'{f.split("/")[-1]:32s} {ut:8s} {n:5d} {dx} {dy} {votes:5d} {fl} {med:5.0f} {grade(flux)}')
        print()

    shifts = np.array([r[3] for r in rows if r[3] is not None])
    dxmin, dxmax = shifts[:, 0].min(), shifts[:, 0].max()
    dymin, dymax = shifts[:, 1].min(), shifts[:, 1].max()
    print(f'frames {len(rows)}, aligned {len(shifts)}, sky {min(r[6] for r in rows):.0f}-{max(r[6] for r in rows):.0f} ADU/px')
    print(f'shift range dx [{dxmin:.0f}, {dxmax:.0f}]  dy [{dymin:.0f}, {dymax:.0f}] px from frame 1')
    m = APERTURE_MARGIN_PX
    print(f'comp box (first-frame pixels, stays on chip all night): '
          f'x in [{m - dxmin:.0f}, {shape[1] - m - dxmax:.0f}]  y in [{m - dymin:.0f}, {shape[0] - m - dymax:.0f}]')
    grades = [grade(r[5]) for r in rows]
    print(f'target flux: clear reference {clear_ref:.0f} ADU (4 px aperture); '
          f'clear {grades.count("clear")}, partial {grades.count("partial")}, lost {grades.count("lost")}')
    steps = [(rows[i][0].split('/')[-1], rows[i][3] - rows[i - 1][3]) for i in range(1, len(rows))
             if rows[i][3] is not None and rows[i - 1][3] is not None and np.hypot(*(rows[i][3] - rows[i - 1][3])) > 25]
    if steps:
        print('pointing steps > 25 px between consecutive frames:')
        for name, d in steps:
            print(f'  {name}: ({d[0]:+.0f}, {d[1]:+.0f})')

    if a.tmid and a.t14:
        half = a.t14 / 48.0
        phase = {'pre': [], 'in': [], 'post': []}
        for r, g in zip(rows, grades):
            jd = r[1]
            if not np.isfinite(jd):
                continue
            ph = 'pre' if jd < a.tmid - half else 'in' if jd <= a.tmid + half else 'post'
            phase[ph].append(g)
        print(f'transit window {Time(a.tmid - half, format="jd").iso[11:16]} to {Time(a.tmid + half, format="jd").iso[11:16]} UT:')
        for ph in ('pre', 'in', 'post'):
            g = phase[ph]
            print(f'  {ph:4s} {len(g):3d} frames: clear {g.count("clear")}, partial {g.count("partial")}, lost {g.count("lost")}')


if __name__ == '__main__':
    main()
