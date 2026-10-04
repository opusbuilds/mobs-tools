"""Inject a known transit into a constant field star in real MObs frames, then reduce the edited frames as an ordinary night.

Usage: inject_transit.py SRC_DIR OUT_DIR X Y TMID_BJD_TDB INITS [--ld a,b,c,d] [--radius R] [--no-verify]

Why image level: a timing bias can live in the reduction itself (sky subtraction, aperture choice, comparison
selection). A fake transit added to a finished photometry file skips exactly those steps and tests only the fitter.
Here the transit goes into the pixels, and the whole pipeline has to find it again.

What it does, per frame: locates the star (star-pair voting against frame 1, then a centre-of-mass centroid), measures
the local sky (median of a 15-25 px annulus), and multiplies only the star's light ABOVE that sky, within R px
(default 10), by the model flux at the frame's mid-exposure time in BJD_TDB. Nothing else in the frame changes. Output
frames keep the header and file name, and the dtype for plain integer and float data (frames with BSCALE/BZERO
scaling come out as float32, values kept), so EXOTIC sees an ordinary night. Planet parameters (period, a/Rs,
inclination, Rp/Rs) and the site come from the inits file; the orbit is circular. Limb darkening is the four-term
Claret law, defaulting to HAT-P-32's coefficients; pass your own with --ld. A symmetric transit shape cannot shift Tmid,
so a mismatch with the fit's limb darkening matters for depth, not timing.

Verification (on by default): after writing, the star's flux above sky in each edited frame is compared with the
original frame times the model; the script exits 1 if any frame deviates by more than 1% (plus integer-rounding noise),
if any pixel outside R changed, or if ANY frame was left uninjected (star not located, within R px of the chip edge,
or with too little light above sky to check). This proves the
edit is on disk as intended; it cannot prove the star, time or model were the ones you meant. X, Y are pixel positions
in the FIRST file by sorted name, and a star must be detected within 3 px of them there. The sky ring is fixed at
15-25 px whatever R is.

To run the test: point an inits file at OUT_DIR with the injected star as the target and the usual comparisons, set
the prior mid-transit time to TMID, reduce, and compare the fitted Tmid with TMID. Use several constant stars per night,
and reduce each star once WITHOUT injection too (same inits, SRC_DIR): a star with structure of its own near TMID will
drag the injected fit. Written 2026-09-28 for a timing that came back 12 minutes early (the moonlit sky, it turned out).
"""
import sys, os, json, argparse, warnings
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from photutils.centroids import centroid_com
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from night_triage import detect, voted_shift, frame_time  # noqa: E402
import pylightcurve as plc  # noqa: E402

HATP32_LD = '2.2830070901439385,-4.194444639860928,4.598187757478915,-1.6867502077619263'
ap = argparse.ArgumentParser(description='Inject a known transit into a field star in FITS frames.')
ap.add_argument('src'); ap.add_argument('out'); ap.add_argument('x', type=float); ap.add_argument('y', type=float)
ap.add_argument('tmid', type=float, help='injected mid-transit time, BJD_TDB')
ap.add_argument('inits', help='EXOTIC inits file: planet parameters, target RA/Dec (for BJD) and site')
ap.add_argument('--ld', default=HATP32_LD, help='four Claret coefficients, comma-separated (default: HAT-P-32)')
ap.add_argument('--radius', type=float, default=10.0, help='edit radius in px (default 10)')
ap.add_argument('--no-verify', action='store_true')
a = ap.parse_args()
src, out, X, Y, TMID, R = a.src, a.out, a.x, a.y, a.tmid, a.radius
LD = [float(v) for v in a.ld.split(',')]
if len(LD) != 4:
    sys.exit('--ld needs exactly four coefficients')
rs, ro = os.path.realpath(src), os.path.realpath(out)
if rs == ro or (os.path.exists(out) and os.path.samefile(src, out)) or ro.startswith(rs + os.sep) or rs.startswith(ro + os.sep):
    sys.exit('OUT_DIR must be a separate folder from SRC_DIR, not the same one, a link to it, or nested in it')
d = json.load(open(a.inits)); U, P = d['user_info'], d['planetary_parameters']
per = float(P['Orbital Period (days)']); ars = float(P['Ratio of Distance to Stellar Radius (a/Rs)'])
inc = float(P['Orbital Inclination (deg)']); rprs = float(P['Ratio of Planet to Stellar Radius (Rp/Rs)'])
coord = SkyCoord(P['Target Star RA'], P['Target Star Dec'], unit=(u.hourangle, u.deg))
loc = EarthLocation(lat=float(U['Obs. Latitude']) * u.deg, lon=float(U['Obs. Longitude']) * u.deg,
                    height=float(U['Obs. Elevation (meters)']) * u.m)
files = sorted(f for f in os.listdir(src) if f.upper().endswith(('.FITS', '.FIT', '.FTS')))
if not files:
    sys.exit(f'no FITS files in {src}')
# Times first, before anything is written: a header the time parser misreads (e.g. a date-only DATE-OBS) would put every
# frame at the same instant, and the injection would quietly be no transit at all.
times = []
for f in files:
    h = fits.getheader(os.path.join(src, f))
    times.append(frame_time(h) + float(h.get('EXPTIME', 60)) / 2 / 86400)
if not np.all(np.isfinite(times)) or len(set(np.round(times, 7))) < len(times):
    sys.exit('frame times are missing or repeated; check DATE-OBS/EXPTIME before injecting')
_, _, ref, _ = detect(os.path.join(src, files[0]), 4.0, 8.0)
ref = np.asarray(ref)
if len(ref) == 0 or np.min(np.hypot(ref[:, 0] - X, ref[:, 1] - Y)) > 3:
    sys.exit(f'no star detected within 3 px of ({X}, {Y}) in {files[0]} (the reference frame: X, Y are in the FIRST file by name)')
os.makedirs(out, exist_ok=True)
log = []
for f, tmid_frame in zip(files, times):
    p = os.path.join(src, f)
    with fits.open(p) as hl:
        h = hl[0].header.copy(); raw = hl[0].data; dtype = raw.dtype; data = raw.astype(float)
    _, med, xy, _ = detect(p, 4.0, 8.0)
    sh, votes = voted_shift(xy, ref)
    t = Time(tmid_frame, format='jd', scale='utc', location=loc)
    bjd = (t.tdb + t.light_travel_time(coord, 'barycentric')).jd
    fl = float(plc.transit(LD, rprs, per, ars, 0.0, inc, 0.0, TMID, np.array([bjd]), method='claret')[0])
    H, W = data.shape
    if sh is not None and votes >= 5 and R <= X + sh[0] < W - R and R <= Y + sh[1] < H - R:
        x0, y0 = X + sh[0], Y + sh[1]
        b = 5; xi, yi = int(round(x0)), int(round(y0))
        cut = data[yi - b:yi + b + 1, xi - b:xi + b + 1] - med
        if cut.shape == (2 * b + 1, 2 * b + 1):
            cx, cy = centroid_com(np.clip(cut, 0, None)); x0, y0 = xi - b + cx, yi - b + cy
        yy, xx = np.mgrid[0:data.shape[0], 0:data.shape[1]]
        r = np.hypot(xx - x0, yy - y0)
        sky = float(np.median(data[(r >= 15) & (r <= 25)]))
        m = r <= R
        data[m] = sky + (data[m] - sky) * fl
        log.append((f, bjd, fl, float(x0), float(y0), sky))
    else:
        log.append((f, bjd, fl, None, None, None))  # star not located, or too near the edge: frame copied unchanged
    lo, hi = (np.iinfo(dtype).min, np.iinfo(dtype).max) if np.issubdtype(dtype, np.integer) else (-np.inf, np.inf)
    newd = np.clip(np.round(data), lo, hi).astype(dtype) if np.issubdtype(dtype, np.integer) else data.astype(dtype)
    fits.PrimaryHDU(newd, h).writeto(os.path.join(out, f), overwrite=True)
json.dump(log, open(os.path.join(out, 'injection_log.json'), 'w'))
fmin = min(l[2] for l in log); miss = sum(l[3] is None for l in log)
print(f'{out}: {len(files)} frames, injected Tmid {TMID}, min model flux {fmin:.4f} (depth {100*(1-fmin):.2f}%), '
      f'star not located in {miss} frames')

if not a.no_verify:
    # What this proves: the files on disk carry the edit the model asked for (no rounding, clipping, dtype or write
    # problem) in EVERY frame. It cannot prove the star was the right one, or the time or model were the ones you meant.
    worst, outside, n, skipped = 0.0, 0.0, 0, []
    for f, bjd, fl, x0, y0, sky in log:
        if x0 is None:
            skipped.append(f); continue
        newraw = fits.getdata(os.path.join(out, f)); is_int = np.issubdtype(newraw.dtype, np.integer)
        orig = fits.getdata(os.path.join(src, f)).astype(float); new = newraw.astype(float)
        yy, xx = np.mgrid[0:orig.shape[0], 0:orig.shape[1]]; r = np.hypot(xx - x0, yy - y0)
        core = r <= min(6.0, R)
        base = (orig[core] - sky).sum()
        rnd = 0.29 * np.sqrt(core.sum())  # sd of the summed rounding error of integer pixels
        if base <= 20 * rnd:  # no real star light here (or a centroid on noise): counts as not injected
            skipped.append(f); continue
        # integer frames are rounded after the edit: allow 5 sigma of that rounding noise on top of 1%
        tol = 0.01 + (5 * rnd / base if is_int else 0)
        worst = max(worst, abs((new[core] - sky).sum() / base - fl) / tol * 0.01)
        outside = max(outside, float(np.abs(new - orig)[r > R].max()))
        n += 1
    print(f'verify: {n} frames checked; worst deviation from the model {worst:.5f} (scaled to a 0.01 limit); '
          f'largest change outside r={R:g}: {outside:.0f}; frames NOT injected: {len(skipped)}')
    if skipped:
        print('  not injected (star not located, within R px of the edge, or too little light above sky): ' + ', '.join(skipped[:10]) + (' ...' if len(skipped) > 10 else ''))
    if n == 0 or skipped or worst > 0.01 or outside > 0:
        print('verify: FAIL (do not trust fits of these frames)')
        sys.exit(1)
    print('verify: PASS')
