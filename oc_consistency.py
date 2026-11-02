#!/usr/bin/env python3
"""Compare each ledger row's PROSE O-C against its computed O-C column.

A ledger where people write notes ("+1.9 min against ExoClock") next to a computed O-C column can disagree with
itself without anyone noticing. In the ledger this was written for, two rows said +1.9 and +1.7 min in their notes
for six weeks while the computed column on the same public rows said +85.8 and +86.4, because the column used an
older ephemeris. The check nobody ran was the row against itself.

Input: a JSON list of rows, each with at least
  id, target, obsDate, note (free text), ocMin: {value: minutes, clock: name of the ephemeris used}
(rows without a note or an ocMin are skipped). This is the observations.json written for
observatory.opusgarden.dev; any ledger exported in that shape works.
Exit 1 if any row's note states minute values and none is within --tol of the column.
usage: oc_consistency.py --published observations.json [--tol 3]
"""
import json, re, sys, argparse

ap = argparse.ArgumentParser()
ap.add_argument('--published', required=True, help='JSON list of ledger rows (see above)')
ap.add_argument('--tol', type=float, default=3.0, help='minutes of allowed disagreement')
a = ap.parse_args()

rows = json.load(open(a.published))
# Phrases in the notes that state an O-C in minutes.
# Extraction is deliberately GREEDY: every signed minute-quantity in the note,
# however it is phrased. The rule below only complains when NONE of them is near
# the column, so over-extracting costs nothing and under-extracting cost three
# false alarms on the first run (notes that said "+6.6 min (0.8 sigma) vs the
# archive" and "leaving +6.6 +/- 2.3", which a phrase-anchored regex missed).
PAT = [
    # (?<!/) : the '-' of a bare '+/- 2.9 min' is an error bar, not a stated value (found by review, 2026-11-02)
    re.compile(r'(?<!/)([+-]\s*\d+(?:\.\d+)?)\s*(?:\+/-\s*[\d.]+\s*)?(?:\([^)]*\)\s*)?min', re.I),
    re.compile(r'O-C\s*(?:=\s*)?([+-]\s*\d+(?:\.\d+)?)', re.I),
    re.compile(r'leaving\s+([+-]\s*\d+(?:\.\d+)?)', re.I),
]
bad = []
for r in rows:
    oc = r.get('ocMin')
    note = r.get('note') or ''
    if not oc or not note or not isinstance(oc.get('value'), (int, float)):
        continue
    stated = []
    for p in PAT:
        for m in p.finditer(note):
            try: stated.append(float(m.group(1).replace(' ', '')))
            except ValueError: pass
    if not stated:
        continue
    # A note may legitimately quote several numbers (other observers, refits).
    # Only complain when NONE of them is near the published column.
    if min(abs(s - oc['value']) for s in stated) > a.tol:
        bad.append((r['id'], r['target'], r['obsDate'], oc['value'], oc.get('clock', '?'), sorted(set(stated))))

print(f"checked {len(rows)} rows; {len(bad)} disagree by more than {a.tol} min")
for i, t, d, col, clock, stated in bad:
    print(f"  {i:>4} {t:13} {d}  column {col:+.1f} min vs {clock}")
    print(f"       note states: {', '.join(f'{s:+.1f}' for s in stated)}")
sys.exit(1 if bad else 0)
