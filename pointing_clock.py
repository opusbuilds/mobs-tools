#!/usr/bin/env python3
"""Where did the telescope ACTUALLY point? A mount aims by hour angle from its sidereal clock, so a clock
N minutes wrong misses its target by about N minutes of RA. Plate-solves the first solvable frame of each
night given (data/<night>/) and reports the miss between the solved centre (converted to coordinates of
date) and the commanded header RA/Dec. The header TELALT/LST-OBS are computed from the same clock and
cannot test it; this can, PROVIDED the mount points open-loop and DATE-OBS shares its clock.
Usage: pointing_clock.py <frames-dir> [...]   (written 2026-09-24 for HAT-P-32 b row 57)"""
import glob, subprocess, os, sys, tempfile, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, FK5
import astropy.units as u

def solve(path):
    d = tempfile.mkdtemp()
    cmd = ['solve-field','--no-plots','--overwrite','--dir',d,'--scale-units','arcsecperpix',
           '--scale-low','4.25','--scale-high','5.75','--downsample','2','--cpulimit','60',
           '--new-fits','none','--corr','none','--rdls','none','--match','none','--index-xyls','none',
           '--solved','none','--temp-axy', path]
    subprocess.run(cmd, capture_output=True, timeout=120)
    w = glob.glob(os.path.join(d,'*.wcs'))
    return WCS(fits.getheader(w[0])) if w else None

out = {}
args = sys.argv[1:]
if not args:
    sys.exit('usage: pointing_clock.py <frames-dir> [<frames-dir> ...]   (a dir of one night\'s .FITS)')
for night in args:
    d = night if os.path.isdir(night) else f'data/{night}'
    fs = sorted(glob.glob(os.path.join(d, '*.FITS')) + glob.glob(os.path.join(d, '*.fits')))
    night = os.path.basename(os.path.normpath(d))
    for f in fs[:6]:
        h = fits.getheader(f)
        w = solve(f)
        if w is None: continue
        ny, nx = fits.getdata(f).shape
        c = w.pixel_to_world((nx-1)/2, (ny-1)/2)                         # ICRS centre, physical
        t = Time(h['UT-OBS'].replace('-0000',''), scale='utc')
        cd = c.transform_to(FK5(equinox=t))                             # to coordinates of date
        dra = ((cd.ra.deg - h['RA'] + 180) % 360 - 180)                 # deg of RA
        ddec = cd.dec.deg - h['DEC']
        dra_sky = dra*np.cos(np.radians(h['DEC']))*60                  # arcmin on sky
        dra_time = dra*4                                                # minutes of time (1 deg RA = 4 min)
        out[night] = dict(frame=os.path.basename(f), ut=h['UT-OBS'][11:19],
                          dRA_arcmin=round(dra_sky,2), dDec_arcmin=round(ddec*60,2), dRA_minutes_of_time=round(dra_time,2))
        print(f"{night:22} {os.path.basename(f)[-11:-5]}  pointing miss: RA {dra_sky:+6.2f}' on sky = {dra_time:+6.2f} min of time, Dec {ddec*60:+6.2f}'", flush=True)
        break
    else:
        print(f"{night:22} no frame among the first 6 solved", flush=True)
