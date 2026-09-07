#!/usr/bin/env python3
"""
Find where a target sits, in pixels, in a frame that carries no WCS.

Every reduction so far has used an inits file somebody else had already made,
which is why a target nobody here has reduced before was blocked: MObs frames
have no astrometric solution in the header, so there was no way to say which
pixel is the planet's star. This plate-solves the frame locally (astrometry.net
with the Tycho-2 indices, no account and no network) and reports the target's
pixel position plus ranked comparison-star candidates.

    python3 locate_target.py FRAME.FITS --name "TOI-1516"
    python3 locate_target.py FRAME.FITS --ra 340.2847 --dec 69.6397 --json

Scale hints come from the header (IM_SCALE, or XPIXSZ/FOCALLEN), so the search
is bounded rather than blind. Comparison candidates are ranked by flux with the
saturated ones dropped, since a saturated comparison is worse than none.
"""
import argparse, json, os, subprocess, sys, tempfile
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def pixel_scale(hdr):
    """arcsec/pixel from whatever the header actually carries."""
    if 'IM_SCALE' in hdr:
        return float(hdr['IM_SCALE'])
    if 'XPIXSZ' in hdr and 'FOCALLEN' in hdr:
        return 206.265 * float(hdr['XPIXSZ']) / float(hdr['FOCALLEN'])
    return None


def solve(path, workdir, timeout):
    hdr = fits.getheader(path)
    scale = pixel_scale(hdr)
    cmd = ['solve-field', '--no-plots', '--overwrite', '--dir', workdir,
           '--downsample', '2', '--cpulimit', str(timeout)]
    if scale:
        cmd += ['--scale-units', 'arcsecperpix',
                '--scale-low', f'{scale * 0.85:.3f}', '--scale-high', f'{scale * 1.15:.3f}']
    if 'RA' in hdr and 'DEC' in hdr:
        cmd += ['--ra', str(float(hdr['RA'])), '--dec', str(float(hdr['DEC'])), '--radius', '2']
    cmd.append(path)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 120)
    stem = os.path.join(workdir, os.path.splitext(os.path.basename(path))[0])
    if not os.path.exists(stem + '.solved'):
        sys.stderr.write(r.stdout[-2000:] + '\n')
        raise SystemExit(f'FAILED to solve {os.path.basename(path)}')
    return stem, hdr


def resolve_name(name):
    from astroquery.simbad import Simbad
    t = Simbad.query_object(name)
    if t is None or len(t) == 0:
        raise SystemExit(f'SIMBAD does not know "{name}": pass --ra/--dec instead')
    row = t[0]
    for ra_key, dec_key in (('ra', 'dec'), ('RA_d', 'DEC_d'), ('RA', 'DEC')):
        if ra_key in t.colnames:
            return float(row[ra_key]), float(row[dec_key])
    raise SystemExit(f'unexpected SIMBAD columns: {t.colnames}')


def peak_adu(data, x, y, box=4):
    """Brightest pixel near a source, for the saturation check."""
    h, w = data.shape
    y0, y1 = max(0, int(y) - box), min(h, int(y) + box + 1)
    x0, x1 = max(0, int(x) - box), min(w, int(x) + box + 1)
    return float(data[y0:y1, x0:x1].max()) if y1 > y0 and x1 > x0 else 0.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('frame')
    p.add_argument('--name', help='resolve target coordinates through SIMBAD')
    p.add_argument('--ra', type=float)
    p.add_argument('--dec', type=float)
    p.add_argument('--comps', type=int, default=8, help='comparison candidates to report')
    p.add_argument('--timeout', type=int, default=120)
    p.add_argument('--json', action='store_true')
    a = p.parse_args()

    if a.name:
        ra, dec = resolve_name(a.name)
    elif a.ra is not None and a.dec is not None:
        ra, dec = a.ra, a.dec
    else:
        raise SystemExit('need --name or both --ra and --dec')

    with tempfile.TemporaryDirectory(prefix='platesolve-') as work:
        stem, hdr = solve(os.path.abspath(a.frame), work, a.timeout)
        wcs = WCS(fits.getheader(stem + '.wcs'))
        sources = fits.getdata(stem + '.axy')
        x, y = wcs.all_world2pix(ra, dec, 0)

    data = fits.getdata(a.frame).astype(float)
    ny, nx = data.shape
    if not (0 <= x < nx and 0 <= y < ny):
        raise SystemExit(f'target lands off-frame at ({x:.1f}, {y:.1f}): wrong field or wrong coordinates')

    sat = float(hdr.get('DATAMAX', data.max()))
    sx, sy, flux = np.asarray(sources['X']), np.asarray(sources['Y']), np.asarray(sources['FLUX'])
    sep = np.hypot(sx - x, sy - y)
    target_peak = peak_adu(data, x, y)

    comps = []
    for i in np.argsort(-flux):
        if sep[i] < 5:
            continue  # the target itself
        pk = peak_adu(data, sx[i], sy[i])
        if pk >= sat * 0.95:
            continue  # saturated: a bad comparison is worse than a missing one
        comps.append({'x': round(float(sx[i]), 2), 'y': round(float(sy[i]), 2),
                      'flux': round(float(flux[i]), 1), 'peak_adu': round(pk, 1),
                      'sep_px': round(float(sep[i]), 1)})
        if len(comps) >= a.comps:
            break

    result = {'frame': os.path.basename(a.frame), 'ra': ra, 'dec': dec,
              'target': {'x': round(float(x), 2), 'y': round(float(y), 2),
                         'peak_adu': round(target_peak, 1), 'saturated': target_peak >= sat * 0.95},
              'datamax': sat, 'comparisons': comps}

    if a.json:
        print(json.dumps(result, indent=2))
        return
    t = result['target']
    print(f"{result['frame']}  ({ra:.5f}, {dec:.5f})")
    print(f"  target      x={t['x']:>7.2f}  y={t['y']:>7.2f}  peak={t['peak_adu']:.0f} ADU"
          + ('  SATURATED' if t['saturated'] else ''))
    print(f"  comparisons (DATAMAX {sat:.0f}, saturated dropped)")
    for c in comps:
        print(f"    x={c['x']:>7.2f}  y={c['y']:>7.2f}  flux={c['flux']:>9.1f}  "
              f"peak={c['peak_adu']:>7.0f}  {c['sep_px']:>6.1f} px away")


if __name__ == '__main__':
    main()
