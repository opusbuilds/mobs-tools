#!/usr/bin/env python3
"""Every rule in inclusion_check.py must be able to FAIL, and every missing input must give UNDECIDED.

Each case builds a small synthetic run directory in a temporary folder (no telescope data needed): an inits file,
an EXOTIC run-log line, the two lines post_run_check.py writes, a FinalParams file and a night.json transit window.
The fixture text copies the real tools' wording, so if their wording changes and a rule can no longer see what it
looks for, this test fails instead of the rule passing quietly. (Why this exists: in the private version, the Moon
rule read a flag nobody ever wrote, and passed every night for two days.)

Moon dates: 2026-11-09 is new moon; 2026-11-24 07:00 UT is near full moon, high over Arizona.
usage: python3 test_inclusion_check.py   (exit 0 = all cases behave)
"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DARK = [2461353.70, 2461353.75, 2461353.80]       # 2026-11-09 ~04:50-07:10 UT, new moon
MOONLIT = [2461368.75, 2461368.79, 2461368.83]    # 2026-11-24 ~06:00-07:55 UT, Moon ~full and up
MOON_DOWN = [2461369.25, 2461369.29, 2461369.33]  # 2026-11-24 ~18:00-20:00 UT, Moon ~full but below the horizon


def make(d, qc='PASS', warn='', rprs='0.165 +/- 0.004', pub=0.165, window=DARK, noise=True, finalparams=True,
         night=True, clocks=None, cov=None):
    out, frames = os.path.join(d, 'out'), os.path.join(d, 'frames')
    os.makedirs(os.path.join(out, 'Diagnostics')); os.makedirs(frames)
    inits = {'user_info': {'Directory to Save Plots': out, 'Directory with FITS files': frames},
             'planetary_parameters': {'Planet Name': 'Test-1 b', 'Ratio of Planet to Stellar Radius (Rp/Rs)': pub,
                                      'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty': 0.002}}
    json.dump(inits, open(os.path.join(d, 'inits.json'), 'w'))
    if qc:
        open(os.path.join(out, 'Diagnostics', 'EXOTIC_RunLog_1.log'), 'w').write(
            f'[11:00]                 Transit detection QC: {qc} - KTMF supports a pass-quality transit fit.\n')
    lines = cov if cov is not None else ['transit-fit comparison #1 (310,240): coverage 80/80 frames, 0 PSF rejects' + warn]
    if noise:
        lines.append('noise: 80 points, median error 0.600%, residual scatter 0.640%, chi2/dof 1.14 '
                     '(errors x1.07 to match), beta 0.95')
    open(os.path.join(out, 'TMID_BAR_CHECK.txt'), 'w').write('\n'.join(lines) + '\n')
    if finalparams:
        json.dump({'FINAL PLANETARY PARAMETERS': {'Ratio of Planet to Stellar Radius (Rp/R*)': rprs}},
                  open(os.path.join(out, 'FinalParams_1.json'), 'w'))
    if night:
        json.dump({'window_jd': window}, open(os.path.join(frames, 'night.json'), 'w'))
    args = [os.path.join(d, 'inits.json')]
    if clocks is not None:
        json.dump(clocks, open(os.path.join(d, 'clocks.json'), 'w'))
        args += ['--clocks', os.path.join(d, 'clocks.json')]
    return args


CASES = [  # name, fixture overrides, substrings of rules that must FAIL, expected exit code, text that must appear
    ('clean night counts', {}, [], 0, None),
    ('clean night with a verified clock', {'clocks': ['Test-1 b']}, [], 0, None),
    ('QC FAIL', {'qc': 'FAIL'}, ['transit QC'], 1, None),
    ('partial-coverage comparison', {'warn': ' -- WARN partial-coverage comparison; the fit lost frames.'},
     ['partial-coverage'], 1, None),
    ('saturating comparison', {'warn': ' -- WARN the comparison saturates in over a tenth of the frames.'},
     ['saturation'], 1, None),
    ('Rp/Rs far from published', {'rprs': '0.330 +/- 0.010'}, ['Rp/Rs'], 1, None),
    ('Rp/Rs asymmetric bars parsed', {'rprs': '0.164 -0.003/+0.005'}, [], 0, None),
    ('moonlit night', {'window': MOONLIT}, ['Moon'], 1, None),
    ('target without a verified clock', {'clocks': ['Other-2 b']}, ['verified ephemeris'], 1, None),
    ('bright Moon below the horizon passes', {'window': MOON_DOWN}, [], 0, None),
    ('WARN on an earlier comparison line still counts',
     {'cov': ['transit-fit comparison #1 (310,240): coverage 60/80 frames, 20 PSF rejects -- WARN partial-coverage comparison; the fit lost frames.',
              'transit-fit comparison #2 (120,90): coverage 80/80 frames, 0 PSF rejects']}, ['partial-coverage'], 1, None),
    ('coverage unreadable: UNDECIDED, never a pass',
     {'cov': ['transit-fit comparison #1 (310,240): COULD NOT READ its coverage from the log (format changed?); check it by hand before quoting Tmid.']},
     [], 2, '[????] post_run_check could not read the comparison coverage'),
    ('published Rp/Rs malformed: UNDECIDED, not a crash', {'pub': ''}, [], 2, '[????] Rp/Rs not parsed'),
    ('no QC line', {'qc': None}, [], 2, '[????] no QC line'),
    ('no noise line', {'noise': False}, [], 2, '[????] no noise line'),
    ('no FinalParams', {'finalparams': False}, [], 2, '[????] no FinalParams'),
    ('no night.json: Moon unknown, never a pass', {'night': False}, [], 2, '[????] Moon not computed'),
]

bad = 0
for name, over, must_fail, code, must_print in CASES:
    with tempfile.TemporaryDirectory() as d:
        args = make(d, **over)
        p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'inclusion_check.py')] + args,
                           capture_output=True, text=True)
    failed = [l for l in p.stdout.splitlines() if l.strip().startswith('[FAIL]')]
    problems = [r for r in must_fail if not any(r in l for l in failed)]
    if p.returncode != code:
        problems.append(f'exit {p.returncode}, expected {code}' + (f' ({p.stderr.strip()[-200:]})' if p.stderr.strip() else ''))
    if must_print and must_print not in p.stdout:
        problems.append('missing: ' + must_print)
    if code == 0 and failed:
        problems.append('unexpected FAIL: ' + '; '.join(failed))
    print(('ok   ' if not problems else 'BAD  ') + name + ('' if not problems else '  -> ' + ', '.join(problems)))
    bad += bool(problems)
sys.exit(1 if bad else 0)
