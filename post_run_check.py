#!/usr/bin/env python
"""Post-reduction check: was the reported Tmid uncertainty the posterior, or a replaced bar?

EXOTIC (WBoM, elca.py `_finalize_ultranest_fit_results`) can replace UltraNest's
posterior stdev with a delta-chi2<=1 dead-point spread that is biased low ~3x
(EXOTIC #1401, diagnosed 2026-09-02). The run log records it in one line:
  "UltraNest uncertainty fallback note: replaced posterior summary error(s) for ... tmid ..."
This script finds the SELECTED final fit in the run log, reports whether that
line fired for tmid, and if so runs tools/indep_tmid.py on the run's own
FinalLightCurve CSV (scout bounds from FinalParams when present) to get the
honest bar. Writes output/<dir>/TMID_BAR_CHECK.txt so the ledger step sees it.

usage: post_run_check.py <inits.json> [--no-refit]
"""
import argparse, glob, json, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument('inits'); ap.add_argument('--no-refit', action='store_true')
a = ap.parse_args()

inits = json.load(open(a.inits))
out = inits['user_info']['Directory to Save Plots']
logs = sorted(glob.glob(os.path.join(out, 'Diagnostics', 'EXOTIC_RunLog_*.log')), key=os.path.getmtime)
fps = sorted(glob.glob(os.path.join(out, 'FinalParams_*.json')), key=os.path.getmtime)
csvs = sorted(glob.glob(os.path.join(out, 'working_artifacts', 'FinalLightCurve_*.csv')), key=os.path.getmtime)
if not logs or not fps:
    print(f'post_run_check: no run log / FinalParams under {out}'); sys.exit(2)

fp = json.load(open(fps[-1])).get('FINAL PLANETARY PARAMETERS', {})
reported = fp.get('Mid-Transit Time (Tmid)', '?')
lines = open(logs[-1], errors='ignore').read().splitlines()

# The selected final fit: the LAST block that reports "parameters: fit_method=" with the
# same Tmid central value FinalParams carries; scan the 14 lines after it for the fallback note.
m = re.search(r'([\d]+\.[\d]+) \+/- ([\d.]+)', str(reported))
central = m.group(1) if m else None
cands = [i for i, l in enumerate(lines) if 'parameters: fit_method=' in l and (central is None or f'Tmid={central}' in l)]
matched = 'by Tmid'
if not cands:
    # A single-comparison run prints its evaluation-stage Tmid here and a different final
    # one after detrending (WASP-52 b 09-20: 2461303.7409 vs .7417), so match by position.
    cands = [i for i, l in enumerate(lines) if 'parameters: fit_method=' in l]
    matched = 'by position (last block; its Tmid differs from FinalParams)'
fired = None; keys = '-'
if cands:
    i = cands[-1]
    fb = [l for l in lines[i:i + 14] if 'uncertainty fallback note' in l]
    if fb:
        km = re.search(r'error\(s\) for ([a-z0-9, ]+) using', fb[0])
        keys = km.group(1) if km else '?'
        fired = 'tmid' in keys.split(', ')
    else:
        fired = False

verdict = []
verdict.append(f'reported Tmid: {reported}')
verdict.append(f'selected-final-fit fallback fired for tmid: {fired}  (replaced keys: {keys}; block matched {matched})')

# Coverage of the comparison the transit fit actually used. WASP-52 b 09-20: EXOTIC chose a
# star with 11 PSF-quality rejects (67/78 frames) over full-coverage stars 0.08% worse on its
# score; the fit ran on 60 points, lost 7 of them in the first half of the transit, went grazing
# and 7 minutes late. A full-coverage comparison put the same night at -0.9 min.
sel = [l for l in lines if 'Transit Fit Comparison Star: #' in l and '[' in l]
if sel:
    sm = re.search(r'#(\d+) - \[([\d.]+), ([\d.]+)\]', sel[-1])
    if sm:
        n, cx, cy = sm.group(1), sm.group(2), sm.group(3)
        # Use the last line that carries the reject counts. Since dev16x EXOTIC prints a second, later line per
        # comparison WITHOUT them, and taking cov[-1] made this whole check vanish silently (TOI-2570 b 2026-10-04:
        # comparison #6 at 63/101 frames, 37 overexposure rejects, and no warning printed).
        cov = [l for l in lines if re.search(rf'Comp {n}\b.*\(x={cx}, y={cy}\).*coverage=', l)]
        # EXOTIC omits psf_quality_rejects / overexposure_rejects when they are zero, so both are optional; a line that
        # carries 'reason=' is a complete per-comparison line even without them.
        def parse(l):
            m = re.search(r'coverage=(\d+) valid frame\(s\) out of (\d+) total', l)
            if not m or ('psf_quality_rejects=' not in l and 'reason=' not in l):
                return None
            pr = re.search(r'psf_quality_rejects=(\d+)', l); ov = re.search(r'overexposure_rejects=(\d+)', l)
            return int(m.group(1)), int(m.group(2)), int(pr.group(1)) if pr else 0, ov.group(1) if ov else None, bool(pr)
        parsed = [r for r in (parse(l) for l in cov) if r]
        # Prefer the last line that carries reject counts: later stages repeat the comparison WITHOUT them, which reads
        # as 'zero rejects' and hid a 37-frame overexposure (TOI-2570 10-04) twice. Count-free lines only as a fallback.
        cms = [r for r in parsed if r[4]] or parsed
        if cms:
            a_, b_, rej, over, _ = cms[-1]
            line = f'transit-fit comparison #{n} ({cx},{cy}): coverage {a_}/{b_} frames, {rej} PSF rejects' + (f' ({over} of them overexposure)' if over else '')
            if a_ < 0.9 * b_:
                line += ' -- WARN partial-coverage comparison; the fit lost frames. Refit with a full-coverage comparison before quoting Tmid.'
            if over and int(over) > 0.1 * b_:
                line += ' -- WARN the comparison saturates in over a tenth of the frames.'
            verdict.append(line)
        else:
            verdict.append(f'transit-fit comparison #{n} ({cx},{cy}): COULD NOT READ its coverage from the log (format changed?); check it by hand before quoting Tmid.')
    else:
        verdict.append('transit-fit comparison: COULD NOT PARSE the selected-comparison line from the log; check coverage by hand.')
else:
    verdict.append('transit-fit comparison: no "Transit Fit Comparison Star" line in the log; coverage not checked.')

if fired and not a.no_refit and csvs:
    bounds = fp.get('Pre-UltraNest LM boundary scout final bounds')
    cmd = [sys.executable, os.path.join(os.path.dirname(__file__), 'indep_tmid.py'), a.inits, csvs[-1], '--nwalk', '40', '--nstep', '4000']
    if bounds: cmd += ['--bounds', bounds if isinstance(bounds, str) else json.dumps(bounds)]
    # If the pipeline sampled a0/a2 (not fixed), refit with the baseline free too.
    if 'fixed' not in str(fp.get('Baseline flux (a0)', '')): cmd.append('--free-baseline')
    r = subprocess.run(cmd, capture_output=True, text=True)
    keep = [l for l in r.stdout.splitlines() if l.startswith(('tmid', '# EXOTIC', '# best-fit', '# points'))]
    verdict.append('INDEPENDENT REFIT (indep_tmid.py):')
    verdict += ['  ' + l for l in keep] or ['  refit failed: ' + r.stderr[-400:]]
    verdict.append('QUOTE THE INDEPENDENT POSTERIOR BAR, NOT THE REPORTED ONE.')
elif fired:
    verdict.append('replaced bar detected; refit skipped (--no-refit or no CSV). Do NOT quote the reported bar.')
elif fired is False:
    verdict.append('reported bar is the posterior summary; safe to quote.')
else:
    verdict.append('could not locate the selected final fit in the log; inspect manually.')

text = '\n'.join(verdict)
open(os.path.join(out, 'TMID_BAR_CHECK.txt'), 'w').write(text + '\n')
print(text)
sys.exit(1 if fired else 0)
