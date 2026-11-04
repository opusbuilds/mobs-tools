#!/usr/bin/env python3
"""Every rule ledger.py enforces must be able to fail, and the export must feed oc_consistency.py.

Uses a temporary database; no telescope data. The O-C case uses a real night (HAT-P-32 b, 2026-11-03, MObs,
reduced with EXOTIC) against the ExoClock ephemeris, where the answer was computed independently: -0.21 min.
usage: python3 test_ledger.py   (exit 0 = all cases behave)
"""
import json, os, sqlite3, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CLOCKS = {'HAT-P-32 b': {'t0': 2456819.856194, 't0Err': 3.5e-5, 'period': 2.150008096, 'periodErr': 4.8e-8,
                         'source': 'ExoClock, Kokori et al. 2025'}}
bad = 0


def run(db, *args):
    return subprocess.run([sys.executable, '-I', os.path.join(HERE, 'ledger.py'), '--db', db, *args],
                          capture_output=True, text=True)


def check(name, ok, detail=''):
    global bad
    print(('ok   ' if ok else 'BAD  ') + name + ('' if ok else '  -> ' + detail))
    bad += not ok


with tempfile.TemporaryDirectory() as d:
    db = os.path.join(d, 'l.db'); out = os.path.join(d, 'out'); clocks = os.path.join(d, 'clocks.json')
    json.dump(CLOCKS, open(clocks, 'w'))
    p = run(db, 'add', '--target', 'X-1 b', '--date', '2026-11-01', '--status', 'rejected')
    check('rejected row without --reason is refused', p.returncode == 1, p.stdout + p.stderr)
    p = run(db, 'add', '--target', 'X-1 b', '--date', '2026-11-01', '--status', 'rejected', '--reason', 'fog')
    check('reason outside the seven is refused', p.returncode == 2, p.stderr[-200:])
    p = run(db, 'add', '--target', 'X-1 b', '--date', '2026-11-01', '--tmid', '2461346.5', '--tmid-err', '0.002')
    check('timed row without --bar-scale is refused', p.returncode == 1, p.stdout + p.stderr)
    p = run(db, 'add', '--target', 'X-1 b', '--date', '2026-11-01', '--notes', 'misspelt flag')
    check('unknown flag (--notes) is an error, not silence', p.returncode == 2, p.stderr[-200:])
    n0 = sqlite3.connect(db).execute('SELECT count(*) FROM reductions').fetchone()[0]
    check('nothing was written by the refused adds', n0 == 0, f'{n0} rows')

    p = run(db, 'add', '--target', 'HAT-P-32 b', '--date', '2026-11-03', '--tmid', '2461347.7731', '--tmid-err', '0.00095',
            '--bar-scale', '1.20', '--rprs', '0.146', '--scatter', '0.479')
    row = sqlite3.connect(db).execute('SELECT tmid_err, tmid_err_raw, bar_scale FROM reductions').fetchone()
    check('stored bar = printed bar x scale, both kept', p.returncode == 0 and abs(row[0] - 0.00114) < 1e-9 and row[1] == 0.00095,
          f'{p.stderr} {row}')
    p = run(db, 'add', '--target', 'Y-2 b', '--date', '2026-11-02', '--status', 'rejected', '--reason', 'cloud')
    check('rejected row with a valid reason is recorded', p.returncode == 0, p.stderr)

    p = run(db, 'append', '--id', '2', '--note', 'Cloud from ingress on.')
    note = sqlite3.connect(db).execute('SELECT note FROM reductions WHERE id = 2').fetchone()[0]
    check('append to an empty note writes the text', p.returncode == 0 and note == 'Cloud from ingress on.', repr(note))
    p = run(db, 'append', '--id', '99', '--note', 'x')
    check('append to a missing row is refused', p.returncode == 1, p.stdout)
    p = run(db, 'crosscheck', '--id', '99', '--source', 'observer A')
    check('cross-check on a missing row is refused', p.returncode == 1, p.stdout)
    p = run(db, 'crosscheck', '--id', '1', '--source', 'observer A', '--tmid', '2461347.7746', '--tmid-err', '0.0037')
    check('cross-check recorded', p.returncode == 0, p.stderr)

    p = run(db, 'export', '--out', out, '--clocks', clocks)
    rows = json.load(open(os.path.join(out, 'observations.json'))) if p.returncode == 0 else []
    hp = [r for r in rows if r['target'] == 'HAT-P-32 b']
    oc = hp[0]['ocMin'] if hp else None
    check('export computes O-C against the clock (-0.21 min, as computed by hand)',
          oc is not None and abs(oc['value'] + 0.21) < 0.01 and oc['epoch'] == 2106, f'{p.stderr} {oc}')
    # by hand: fit 0.00114 d = 1.642 min; clock hypot(3.5e-5, 2106*4.8e-8) d = 0.154 min; total 1.65 min
    check('O-C errors: clock 0.15 min, total 1.65 min', oc is not None and abs(oc['clockErr'] - 0.15) < 0.006
          and abs(oc['err'] - 1.65) < 0.006, str(oc))
    check('a rejected row gets no O-C', all(r['ocMin'] is None for r in rows if r['status'] == 'rejected'), str(rows))

    # The export feeds oc_consistency: a note that states the right O-C passes, a stale one is flagged.
    run(db, 'append', '--id', '1', '--note', 'O-C -0.2 +/- 1.6 min vs ExoClock.')
    run(db, 'export', '--out', out, '--clocks', clocks)
    q = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'oc_consistency.py'), '--published',
                        os.path.join(out, 'observations.json')], capture_output=True, text=True)
    check('oc_consistency passes a consistent row from the export', q.returncode == 0, q.stdout + q.stderr)
    sqlite3.connect(db).execute("UPDATE reductions SET note = 'Earlier hand score: +86.0 min.' WHERE id = 1").connection.commit()
    run(db, 'export', '--out', out, '--clocks', clocks)
    q = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'oc_consistency.py'), '--published',
                        os.path.join(out, 'observations.json')], capture_output=True, text=True)
    check('oc_consistency flags a stale note from the export', q.returncode == 1 and 'HAT-P-32 b' in q.stdout, q.stdout)

    p = run(db, 'append', '--id', '2', '--note', 'Second line.')
    note = sqlite3.connect(db).execute('SELECT note FROM reductions WHERE id = 2').fetchone()[0]
    check('append to an existing note joins with a space', note == 'Cloud from ingress on. Second line.', repr(note))
    p = run(db, 'reason', '--id', '2', '--reason', 'partial')
    r2 = sqlite3.connect(db).execute('SELECT reason FROM reductions WHERE id = 2').fetchone()[0]
    check('reason command updates the reason', p.returncode == 0 and r2 == 'partial', p.stderr + str(r2))
    p = run(db, 'add', '--target', 'Z-3 b', '--date', '2026-11-04', '--tmid', '2461348.5')
    check('--tmid without --tmid-err is refused', p.returncode == 1, p.stdout + p.stderr)
    p = run(db, 'add', '--target', 'Z-3 b', '--date', '2026-11-04', '--tmid', '2461348.5', '--tmid-err', '0.001', '--bar-scale', '0')
    check('--bar-scale 0 is refused', p.returncode == 1, p.stdout + p.stderr)
    lst = os.path.join(d, 'list.json'); json.dump(['HAT-P-32 b'], open(lst, 'w'))
    p = run(db, 'export', '--out', out, '--clocks', lst)
    check('a list-form clocks file is refused by the ledger with a clear message', p.returncode == 1 and 'dict form' in p.stderr, p.stderr)

    old = os.path.join(d, 'old.db'); con = sqlite3.connect(old)
    con.execute('CREATE TABLE reductions (id INTEGER PRIMARY KEY AUTOINCREMENT, target TEXT NOT NULL, obs_date TEXT NOT NULL,'
                ' status TEXT NOT NULL DEFAULT "reduced", reduced_at TEXT NOT NULL)'); con.commit(); con.close()
    p = run(old, 'add', '--target', 'Old-1 b', '--date', '2026-01-01', '--status', 'rejected', '--reason', 'cloud', '--dataset', 'x')
    check('an older database is migrated, not crashed on', p.returncode == 0, p.stderr[-200:])

    bad_out = os.path.join(d, 'badout'); os.makedirs(bad_out); open(os.path.join(bad_out, 'observations.json'), 'w').write('{not json')
    p = run(db, 'export', '--out', bad_out)
    check('an unreadable previous export is refused, not overwritten', p.returncode == 2 and 'unreadable' in p.stderr, p.stderr)

    con = sqlite3.connect(db); con.execute('DELETE FROM reductions WHERE id = 2'); con.commit()
    p = run(db, 'export', '--out', out, '--clocks', clocks)
    check('export refuses to drop a night the last export had', p.returncode == 2 and 'Y-2 b' in p.stderr, p.stderr)
sys.exit(1 if bad else 0)
