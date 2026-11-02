#!/usr/bin/env python3
"""Decide whether a reduced night may be used as a timing point, WITHOUT ever reading or printing its timing.

Why: if you decide which nights count after seeing their mid-transit times, the choice can bend toward the answer
you want, and nobody, including you, can tell afterwards. So the decision reads only the quality fields the rules
name, one by one, and never a line containing a mid-transit time. Run it, commit its verdict, and only then compute
the O-C. (Written for a pre-registered bet on MObs timing; the defaults are that bet's rules.)

Rules (each prints PASS or FAIL with its evidence):
  - EXOTIC's transit-detection QC is PASS or MARGINAL (the newest run log in <plots dir>/Diagnostics/);
  - post_run_check.py printed no partial-coverage and no saturation warning for the comparison stars;
  - the fitted Rp/Rs is within --rprs-sigma of the published value in the inits;
  - the Moon is not MOONLIT: more than --moon-pct lit AND above the horizon at ingress, mid or egress
    (computed for the MicroObservatory site from the transit window in <frames dir>/night.json, written by
    mobs_night.py);
  - optionally (--clocks FILE: a JSON list of planet names), the target has a verified ephemeris you trust.
Anything it cannot read makes the verdict UNDECIDED rather than passing quietly.

Exit 0 = counts, 1 = not a point, 2 = could not decide.
usage: inclusion_check.py <inits.json> [--clocks verified.json] [--rprs-sigma 3] [--moon-pct 50]
Run post_run_check.py on the night first. test_inclusion_check.py proves every rule can fail.
"""
import argparse, glob, json, os, re, sys

ap = argparse.ArgumentParser(usage=__doc__.split('usage: ')[1].split('\n')[0])
ap.add_argument('inits')
ap.add_argument('--clocks', help='JSON list of planet names with a verified ephemeris (optional rule)')
ap.add_argument('--rprs-sigma', type=float, default=3.0)
ap.add_argument('--moon-pct', type=float, default=50.0)
a = ap.parse_args()

inits = json.load(open(a.inits))
U, P = inits['user_info'], inits['planetary_parameters']
out = U['Directory to Save Plots']
target = P.get('Planet Name', '?')
checks, undecided = [], []


def rule(name, ok, detail):
    checks.append((name, ok, detail))


# QC verdict, from the newest run log: only the "Transit detection QC:" line is read.
logs = sorted(glob.glob(os.path.join(out, 'Diagnostics', 'EXOTIC_RunLog_*.log')), key=os.path.getmtime)
qc = None
if logs:
    for line in open(logs[-1], errors='ignore'):
        m = re.search(r'Transit detection QC: (PASS|MARGINAL|FAIL)', line)
        if m:
            qc = m.group(1)
if qc is None:
    undecided.append('no QC line in the run log')
else:
    rule('transit QC is PASS or MARGINAL', qc in ('PASS', 'MARGINAL'), qc)

# post_run_check.py output: the comparison-coverage line (with its WARNs) and the noise line, matched by leading words.
chk = os.path.join(out, 'TMID_BAR_CHECK.txt')
cov_lines, noise_line = [], None
if os.path.exists(chk):
    for line in open(chk):
        if line.startswith('transit-fit comparison'):
            cov_lines.append(line)          # every comparison line: a WARN on any one of them must count
        elif line.startswith('noise:'):
            noise_line = line
if not cov_lines:
    undecided.append('no comparison-coverage line from post_run_check (run it first)')
else:
    # A line without "coverage N/M frames" is post_run_check saying it could NOT read the coverage. It holds
    # neither warning word, so treating it as clean would be a quiet pass (found by review, 2026-11-02).
    unread = [l for l in cov_lines if not re.search(r'coverage \d+/\d+ frames', l)]
    if unread:
        undecided.append('post_run_check could not read the comparison coverage: ' + unread[0].strip()[:120])
    else:
        cov = ', '.join(re.search(r'coverage (\d+/\d+) frames', l).group(1) + ' frames' for l in cov_lines)
        rule('no partial-coverage WARN', not any('partial-coverage' in l for l in cov_lines), cov)
        rule('no saturation WARN', not any('saturates' in l for l in cov_lines), 'comparison saturation checked')
scale = None
if noise_line is None:
    undecided.append('no noise line from post_run_check')
else:
    m = re.search(r'errors x([\d.]+) to match', noise_line)
    scale = float(m.group(1)) if m else None
    if scale is None:
        undecided.append('noise scale factor not parsed')

# Rp/Rs against the published value.
fps = sorted(glob.glob(os.path.join(out, 'FinalParams_*.json')), key=os.path.getmtime)
def num(v, default):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


pub = num(P.get('Ratio of Planet to Stellar Radius (Rp/Rs)'), float('nan'))      # unreadable -> nan -> UNDECIDED
pub_err = num(P.get('Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty'), 0.0)
if fps:
    F = json.load(open(fps[-1])).get('FINAL PLANETARY PARAMETERS', {})
    m = re.search(r'([\d.]+)\s*(?:\+/-\s*([\d.]+)|-([\d.]+)/\+([\d.]+))', str(F.get('Ratio of Planet to Stellar Radius (Rp/R*)', '')))
    if m and pub == pub:
        k = float(m.group(1))
        e = float(m.group(2)) if m.group(2) else (float(m.group(3)) + float(m.group(4))) / 2
        z = abs(k - pub) / max((e ** 2 + pub_err ** 2) ** 0.5, 1e-9)
        rule(f'Rp/Rs within {a.rprs_sigma:g} sigma of the published value', z <= a.rprs_sigma,
             f'{k:.3f} +/- {e:.3f} vs {pub:.4f} ({z:.1f} sigma)')
    else:
        undecided.append('Rp/Rs not parsed from FinalParams or inits')
else:
    undecided.append('no FinalParams file')

# The Moon, from night.json's transit window (window_jd: ingress, mid, egress as JDs).
try:
    night = json.load(open(os.path.join(U['Directory with FITS files'], 'night.json')))
    jds = [float(v) for v in night['window_jd']]
    if len(jds) != 3:
        raise ValueError('window_jd must hold ingress, mid, egress')
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mobs_night import moon_during
    pct, alts = moon_during(jds)
except Exception as e:  # anything missing or broken: the Moon is UNKNOWN, never a silent pass
    undecided.append(f'Moon not computed ({type(e).__name__}: {e})')
else:
    flagged = pct > a.moon_pct and max(alts) > 0
    rule(f'Moon not MOONLIT (>{a.moon_pct:g}% lit and up during the transit)', not flagged,
         f'{pct:.0f}% lit, altitude {alts[0]:+.0f}/{alts[1]:+.0f}/{alts[2]:+.0f} deg at ingress/mid/egress')

if a.clocks:
    try:
        verified = json.load(open(a.clocks))
        rule('verified ephemeris for the target', target in verified, target)
    except Exception as e:
        undecided.append(f'clocks file unreadable ({type(e).__name__}: {e})')

print(f'inclusion check for {target} ({os.path.basename(os.path.normpath(out))}); no timing is read or shown')
for name, ok, detail in checks:
    print(f'  [{"PASS" if ok else "FAIL"}] {name}: {detail}')
for u in undecided:
    print(f'  [????] {u}')
if scale is not None:
    print(f'  noise scale factor (multiply the fit bar by this): x{scale:.2f}')
if undecided:
    print('VERDICT: UNDECIDED (fix the missing inputs; do not look at the timing yet)')
    sys.exit(2)
if all(ok for _, ok, _ in checks):
    print('VERDICT: COUNTS. Commit this decision, then compute the O-C.')
    sys.exit(0)
print('VERDICT: NOT A POINT under these rules.')
sys.exit(1)
