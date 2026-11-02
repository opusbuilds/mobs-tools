#!/usr/bin/env python3
"""oc_consistency.py must flag a row whose note disagrees with its column, and only that row.
usage: python3 test_oc_consistency.py   (exit 0 = behaves)"""
import json, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = [
    {'id': 1, 'target': 'Agree-1 b', 'obsDate': '2026-10-11', 'note': 'O-C +0.3 +/- 2.9 min vs ExoClock.',
     'ocMin': {'value': 0.3, 'clock': 'ExoClock'}},
    {'id': 2, 'target': 'Stale-2 b', 'obsDate': '2026-09-01', 'note': 'Scored by hand: +1.9 min.',
     'ocMin': {'value': 85.8, 'clock': 'NASA Exoplanet Archive'}},
    {'id': 3, 'target': 'Many-3 b', 'obsDate': '2026-09-02', 'note': 'Another observer got -14.7 min; refit leaving +2.2 +/- 1.1.',
     'ocMin': {'value': 2.0, 'clock': 'ExoClock'}},
    {'id': 5, 'target': 'BareBar-5 b', 'obsDate': '2026-09-04', 'note': 'Fit O-C +6.0 (+/- 2.9 min).',
     'ocMin': {'value': -2.9, 'clock': 'ExoClock'}},          # the bar's '-2.9' must not count as a stated value
    {'id': 6, 'target': 'NullCol-6 b', 'obsDate': '2026-09-05', 'note': '+1.0 min', 'ocMin': {'value': None, 'clock': 'x'}},
    {'id': 4, 'target': 'NoNote-4 b', 'obsDate': '2026-09-03', 'note': '', 'ocMin': {'value': 40.0, 'clock': 'x'}},
]
with tempfile.TemporaryDirectory() as d:
    f = os.path.join(d, 'rows.json'); json.dump(ROWS, open(f, 'w'))
    p = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'oc_consistency.py'), '--published', f],
                       capture_output=True, text=True)
problems = []
if p.returncode != 1: problems.append(f'exit {p.returncode}, expected 1')
for t in ('Stale-2 b', 'BareBar-5 b'):
    if t not in p.stdout: problems.append(f'{t} was not flagged')
for t in ('Agree-1 b', 'Many-3 b', 'NoNote-4 b', 'NullCol-6 b'):
    if t in p.stdout: problems.append(f'{t} flagged wrongly')
print('ok' if not problems else 'BAD -> ' + ', '.join(problems)); print(p.stdout.strip())
sys.exit(1 if problems else 0)
