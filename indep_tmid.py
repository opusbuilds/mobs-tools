#!/usr/bin/env python
"""Independent Tmid posterior check for an EXOTIC FinalLightCurve CSV.

Refits the detrended light curve EXOTIC wrote (working_artifacts/FinalLightCurve_*.csv)
with EXOTIC's own transit model (exotic.api.elca.transit), the same nonlinear
limb darkening (exotic.api.ld), the same fixed baseline (a0=1, a2=0), and a plain
emcee sampler over {tmid, rprs, ars, inc} with the same priors the inits file
carries. Purpose: give the pipeline's reported Tmid uncertainty an external
reference on byte-identical points (EXOTIC #1401).

usage: indep_tmid.py <inits.json> <FinalLightCurve.csv> [--free-baseline] [--nwalk 40] [--nstep 4000]
"""
import argparse, json, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('inits'); ap.add_argument('csv')
ap.add_argument('--free-baseline', action='store_true', help='also sample a0 (scale) and a2 (airmass slope)')
ap.add_argument('--nwalk', type=int, default=40); ap.add_argument('--nstep', type=int, default=4000)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--bounds', help='JSON dict of {param:[lo,hi]} to override bounds (e.g. EXOTIC scout bounds)')
a = ap.parse_args()

import emcee
from exotic.api.elca import transit
from exotic.api.ld import LimbDarkening
from exotic.exotic import nonlinear_ld

pp = json.load(open(a.inits))['planetary_parameters']
oi = json.load(open(a.inits)).get('optional_info', {})
stellar = dict(teff=pp['Star Effective Temperature (K)'], teffUncPos=pp['Star Effective Temperature (+) Uncertainty'],
               teffUncNeg=pp['Star Effective Temperature (-) Uncertainty'], met=pp['Star Metallicity ([FE/H])'],
               metUncPos=pp['Star Metallicity (+) Uncertainty'], metUncNeg=pp['Star Metallicity (-) Uncertainty'],
               logg=pp['Star Surface Gravity (log(g))'], loggUncPos=pp['Star Surface Gravity (+) Uncertainty'],
               loggUncNeg=pp['Star Surface Gravity (-) Uncertainty'])
ld = LimbDarkening(stellar)
ld.set_filter('CV', 'custom', oi.get('Filter Minimum Wavelength (nm)', 350), oi.get('Filter Maximum Wavelength (nm)', 850))
ld.calculate_ld()
u = [ld.ld0[0], ld.ld1[0], ld.ld2[0], ld.ld3[0]]

d = np.genfromtxt(a.csv, delimiter=',', comments='#', names=None)
t, f, e, am = d[:, 0], d[:, 2], d[:, 3], d[:, 5]
ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
t, f, e, am = t[ok], f[ok], e[ok], am[ok]

P = pp['Orbital Period (days)']; T0 = pp['Published Mid-Transit Time (BJD-UTC)']
n = np.round((np.mean(t) - T0) / P); tmid0 = T0 + n * P
pri = dict(rprs=(pp['Ratio of Planet to Stellar Radius (Rp/Rs)'], pp['Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty']),
           ars=(pp['Ratio of Distance to Stellar Radius (a/Rs)'], pp['Ratio of Distance to Stellar Radius (a/Rs) Uncertainty']),
           inc=(pp['Orbital Inclination (deg)'], pp['Orbital Inclination (deg) Uncertainty']))

def model(th):
    tmid, rprs, ars, inc = th[:4]
    a0 = th[4] if a.free_baseline else 1.0
    a2 = th[5] if a.free_baseline else 0.0
    v = dict(rprs=rprs, ars=ars, inc=inc, tmid=tmid, per=P, ecc=pp.get('Orbital Eccentricity (0 if null)', 0) or 0,
             omega=pp.get('Argument of Periastron (deg)', 0) or 0, u0=u[0], u1=u[1], u2=u[2], u3=u[3])
    return transit(t, v) * a0 * np.exp(a2 * am)

# EXOTIC-style bounds: tmid +/- 0.25*duration-ish window, rprs/ars/inc broad uniform (data-driven), no Gaussian on geometry
lo = dict(tmid=tmid0 - 0.05, rprs=0.02, ars=2.0, inc=70.0, a0=0.9, a2=-1.0)
hi = dict(tmid=tmid0 + 0.05, rprs=0.4, ars=20.0, inc=90.0, a0=1.1, a2=1.0)
if a.bounds:
    for k, (l, h) in json.loads(a.bounds).items(): lo[k], hi[k] = l, h
names = ['tmid', 'rprs', 'ars', 'inc'] + (['a0', 'a2'] if a.free_baseline else [])

def lnp(th):
    for k, x in zip(names, th):
        if not (lo[k] <= x <= hi[k]): return -np.inf
    r = (f - model(th)) / e
    return -0.5 * np.sum(r * r)

rng = np.random.default_rng(a.seed)
p0 = np.array([tmid0, pri['rprs'][0], pri['ars'][0], min(pri['inc'][0], 89.9)] + ([1.0, 0.0] if a.free_baseline else []))
sc = np.array([0.003, 0.005, 0.2, 0.5] + ([0.003, 0.01] if a.free_baseline else []))
walk = p0 + sc * rng.standard_normal((a.nwalk, len(p0)))
walk[:, 3] = np.minimum(walk[:, 3], 89.99)
s = emcee.EnsembleSampler(a.nwalk, len(p0), lnp)
s.run_mcmc(walk, a.nstep, progress=False)
ch = s.get_chain(discard=a.nstep // 2, flat=True)
try: tau = s.get_autocorr_time(discard=a.nstep // 2, quiet=True)
except Exception: tau = None
print(f'# points={len(t)} LD={np.round(u,4).tolist()} tmid0={tmid0:.5f} baseline={"free" if a.free_baseline else "fixed a0=1,a2=0"}')
print(f'# autocorr tau={np.round(tau,1).tolist() if tau is not None else "n/a"}  samples={len(ch)}')
for i, k in enumerate(names):
    q = np.percentile(ch[:, i], [16, 50, 84])
    print(f'{k:5} = {q[1]:.5f}  -{q[1]-q[0]:.5f} +{q[2]-q[1]:.5f}   (std {ch[:,i].std():.5f})')
best = ch[np.argmax(s.get_log_prob(discard=a.nstep // 2, flat=True))]
r = (f - model(best)) / e
# Reproduce EXOTIC's _loglike_neighborhood_uncertainty on these samples (elca.py ~L2270):
lp = s.get_log_prob(discard=a.nstep // 2, flat=True)
for dchi2 in (1.0, 4.0, 9.0):
    m = 2.0 * (lp.max() - lp) <= dchi2
    v = ch[m, 0]
    if v.size >= 8:
        lo_, hi_ = np.percentile(v, [15.8655, 84.1345])
        print(f'# EXOTIC-neighborhood estimator on tmid, dchi2<={dchi2:g}: n={v.size} std={v.std():.5f} halfwidth={(hi_-lo_)/2:.5f} halfRANGE={(v.max()-v.min())/2:.5f} -> ratio posterior/neighborhood = {ch[:,0].std()/max(v.std(),(hi_-lo_)/2):.2f}')
        break
print(f'# best-fit chi2={np.sum(r*r):.2f} dof={len(t)-len(names)}  rms={np.std(f-model(best))*100:.4f}%')
