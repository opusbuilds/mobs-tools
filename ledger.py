#!/usr/bin/env python3
"""A transit-timing ledger: one SQLite row per reduced (or rejected) night, exported as JSON.

Why a ledger and not a folder of results: rejected nights are data too. Most nights a small telescope takes are
unusable, and a record that keeps only the successes cannot show how often the instrument fails or why. Every
night gets a row; a rejected row must say which of seven reasons applies.

Rules the ledger enforces, each learned from a mistake:
  - Unknown flags are errors. (A misspelt --notes once wrote eight rows of empty notes, and nothing complained.)
  - A rejected row needs --reason, one of: cloud, floor, no comparison, partial, no transit, saturated, qc fail.
  - A timed row's --tmid-err is the fit's printed bar and needs --bar-scale, the factor post_run_check.py prints
    ("errors xK to match"); the stored bar is the printed one times that factor, and both are kept.
  - append verifies the note actually grew. (Appending to a NULL note once "succeeded" while writing nothing.)
  - export refuses to drop a row that the previous export had: the export is generated, never edited by hand.
Convention (not enforced): cross-checks are other observers' reductions of the SAME frames, under anonymous
labels such as "observer A", never names.

usage:
  ledger.py [--db ledger.db] add --target "HAT-P-32 b" --date 2026-11-03 [--dataset ID] [--status reduced]
            [--tmid BJD --tmid-err D --bar-scale K [--bar-source TEXT]] [--rprs R] [--scatter PCT] [--reason R] [--note TEXT]
  ledger.py [--db ledger.db] append --id N --note TEXT
  ledger.py [--db ledger.db] crosscheck --id N --source "observer A" [--tmid ...] [--tmid-err ...] [--depth ...]
            [--depth-err ...] [--scatter ...] [--note ...]
  ledger.py [--db ledger.db] reason --id N --reason R
  ledger.py [--db ledger.db] export --out DIR [--clocks clocks.json]
  ledger.py [--db ledger.db] list

export writes DIR/observations.json (and DIR/crosschecks.json): a list of rows with id, target, obsDate, dataset,
status, reason, tmid, tmidErr, tmidErrRaw, barScale, rprs, scatterPct, note, and, when --clocks gives the
target's ephemeris, ocMin = {value, err, clockErr, clock} in minutes. That is the input oc_consistency.py reads.
clocks.json is a dict: {"HAT-P-32 b": {"t0": BJD_TDB, "t0Err": d, "period": d, "periodErr": d, "source": "..."}}.
Its keys also work as inclusion_check.py --clocks (it only asks whether the target is in it); the ledger itself needs
the dict form, since it computes from t0 and period. The clock bar propagates t0Err and n*periodErr in quadrature,
ignoring their covariance, which usually overstates it slightly at large epoch n. --db goes before the command.
"""
import argparse, json, math, os, sqlite3, sys
from datetime import datetime, timezone

REASONS = ("cloud", "floor", "no comparison", "partial", "no transit", "saturated", "qc fail")
STATUSES = ("reduced", "validation", "submitted", "rejected")
SCHEMA = """
CREATE TABLE IF NOT EXISTS reductions (
  id INTEGER PRIMARY KEY AUTOINCREMENT, target TEXT NOT NULL, obs_date TEXT NOT NULL, dataset TEXT,
  tmid REAL, tmid_err REAL, tmid_err_raw REAL, bar_scale REAL, bar_scale_source TEXT,
  rprs REAL, scatter_pct REAL, status TEXT NOT NULL DEFAULT 'reduced', reason TEXT,
  reduced_at TEXT NOT NULL, note TEXT);
CREATE TABLE IF NOT EXISTS crosschecks (
  id INTEGER PRIMARY KEY AUTOINCREMENT, reduction_id INTEGER NOT NULL REFERENCES reductions(id),
  source TEXT NOT NULL, tmid REAL, tmid_err REAL, depth_pct REAL, depth_err REAL, scatter_pct REAL,
  note TEXT, added_at TEXT NOT NULL);
"""
COLUMNS = {'dataset': 'TEXT', 'tmid': 'REAL', 'tmid_err': 'REAL', 'tmid_err_raw': 'REAL', 'bar_scale': 'REAL',
           'bar_scale_source': 'TEXT', 'rprs': 'REAL', 'scatter_pct': 'REAL', 'reason': 'TEXT', 'note': 'TEXT'}


def migrate(db):
    """Add any column an older database lacks, so an old ledger keeps working instead of raising."""
    have = {r[1] for r in db.execute('PRAGMA table_info(reductions)')}
    for col, typ in COLUMNS.items():
        if col not in have:
            db.execute(f'ALTER TABLE reductions ADD COLUMN {col} {typ}')
    db.commit()


def load_clocks(path):
    if not path:
        return {}
    try:
        c = json.load(open(path))
    except (OSError, ValueError) as e:
        die(f'clocks file {path} unreadable ({type(e).__name__}: {e})')
    if not isinstance(c, dict):
        die('--clocks for the ledger must be the dict form {"Planet b": {"t0": ..., "period": ...}}; a list of names '
            'only works for inclusion_check.py')
    for name, v in c.items():
        num = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
        if not (isinstance(v, dict) and num(v.get('t0')) and num(v.get('period')) and v['period'] > 0):
            die(f'clock for {name!r} needs numeric t0 and a positive period')
    return c


def write_json(obj, path):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def die(msg, code=1):
    print(msg, file=sys.stderr)
    sys.exit(code)


def parser():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], allow_abbrev=False)
    ap.add_argument('--db', default='ledger.db')
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('add', allow_abbrev=False)
    a.add_argument('--target', required=True); a.add_argument('--date', required=True)
    a.add_argument('--dataset'); a.add_argument('--status', default='reduced', choices=STATUSES)
    a.add_argument('--tmid', type=float); a.add_argument('--tmid-err', type=float)
    a.add_argument('--bar-scale', type=float); a.add_argument('--bar-source')
    a.add_argument('--rprs', type=float); a.add_argument('--scatter', type=float)
    a.add_argument('--reason', choices=REASONS); a.add_argument('--note')
    p = sub.add_parser('append', allow_abbrev=False)
    p.add_argument('--id', type=int, required=True); p.add_argument('--note', required=True)
    c = sub.add_parser('crosscheck', allow_abbrev=False)
    c.add_argument('--id', type=int, required=True); c.add_argument('--source', required=True)
    for f in ('--tmid', '--tmid-err', '--depth', '--depth-err', '--scatter'):
        c.add_argument(f, type=float)
    c.add_argument('--note')
    r = sub.add_parser('reason', allow_abbrev=False)
    r.add_argument('--id', type=int, required=True); r.add_argument('--reason', required=True, choices=REASONS)
    e = sub.add_parser('export', allow_abbrev=False)
    e.add_argument('--out', required=True); e.add_argument('--clocks')
    sub.add_parser('list', allow_abbrev=False)
    return ap


def oc_minutes(row, clock):
    """O-C in minutes against a linear ephemeris, with the fit bar and the propagated clock bar."""
    n = round((row['tmid'] - clock['t0']) / clock['period'])
    pred = clock['t0'] + n * clock['period']
    clock_err = math.hypot(clock.get('t0Err', 0.0), n * clock.get('periodErr', 0.0)) * 1440
    fit_err = (row['tmid_err'] or 0.0) * 1440
    return {'value': round((row['tmid'] - pred) * 1440, 2), 'err': round(math.hypot(fit_err, clock_err), 2),
            'clockErr': round(clock_err, 2), 'clock': clock.get('source', 'unnamed clock'), 'epoch': n}


def export(db, out, clocks_path):
    clocks = load_clocks(clocks_path)
    rows = [dict(r) for r in db.execute('SELECT * FROM reductions ORDER BY obs_date DESC, id DESC')]
    path = os.path.join(out, 'observations.json')
    if os.path.exists(path):   # drop guard: never make a published row disappear
        try:
            prev = json.load(open(path))
            keys = [(o.get('id'), o['target'], o['obsDate']) for o in prev]
        except (ValueError, KeyError, TypeError) as e:
            die(f'EXPORT REFUSED: the previous {path} is unreadable ({type(e).__name__}); move it aside deliberately.', 2)
        have_ids = {r['id'] for r in rows}; have = {(r['target'], r['obs_date']) for r in rows}
        orphans = [{'target': t, 'obsDate': d} for i, t, d in keys
                   if (i not in have_ids if i is not None else (t, d) not in have)]
        if orphans:
            die('EXPORT REFUSED: %d row(s) in %s are not in the ledger: %s. Record them with `add`; the export is '
                'generated and hand edits are lost.' % (len(orphans), path,
                                                        ', '.join(f"{o['target']} {o['obsDate']}" for o in orphans)), 2)
    obs = []
    for r in rows:
        o = {'id': r['id'], 'target': r['target'], 'obsDate': r['obs_date'], 'dataset': r['dataset'],
             'status': r['status'], 'reason': r['reason'], 'tmid': r['tmid'], 'tmidErr': r['tmid_err'],
             'tmidErrRaw': r['tmid_err_raw'], 'barScale': r['bar_scale'], 'rprs': r['rprs'],
             'scatterPct': r['scatter_pct'], 'note': r['note'], 'ocMin': None}
        if r['tmid'] is not None and r['status'] != 'rejected' and r['target'] in clocks:
            o['ocMin'] = oc_minutes(r, clocks[r['target']])
        obs.append(o)
    os.makedirs(out, exist_ok=True)
    write_json(obs, path)
    checks = [dict(c) for c in db.execute(
        'SELECT c.*, r.target, r.obs_date FROM crosschecks c JOIN reductions r ON r.id = c.reduction_id ORDER BY c.id')]
    write_json(checks, os.path.join(out, 'crosschecks.json'))
    print(f'exported {len(obs)} row(s) and {len(checks)} cross-check(s) to {out}')


def main(argv=None):
    a = parser().parse_args(argv)        # argparse: an unknown flag is an error, not silence
    db = sqlite3.connect(a.db); db.row_factory = sqlite3.Row; db.executescript(SCHEMA); migrate(db)
    if a.cmd == 'add':
        if a.status == 'rejected' and not a.reason:
            die('a rejected row needs --reason (one of: %s)' % ', '.join(REASONS))
        if a.status != 'rejected' and a.tmid_err is not None and not (a.bar_scale and a.bar_scale > 0):
            die("a timed row's --tmid-err needs --bar-scale (the 'errors xK to match' factor from post_run_check.py)")
        if (a.tmid is None) != (a.tmid_err is None):
            die('--tmid and --tmid-err go together')
        if a.bar_scale is not None and a.bar_scale <= 0:
            die('--bar-scale must be positive')
        err = a.tmid_err * a.bar_scale if a.tmid_err is not None and a.bar_scale else a.tmid_err
        cur = db.execute(
            'INSERT INTO reductions (target, obs_date, dataset, tmid, tmid_err, tmid_err_raw, bar_scale, bar_scale_source,'
            ' rprs, scatter_pct, status, reason, reduced_at, note) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
            (a.target, a.date, a.dataset, a.tmid, err, a.tmid_err if a.bar_scale else None, a.bar_scale,
             (a.bar_source or 'post_run_check noise line') if a.bar_scale else None, a.rprs, a.scatter, a.status,
             a.reason, now(), a.note))
        db.commit(); print(f'recorded #{cur.lastrowid} {a.target} ({a.date}): {a.status}')
    elif a.cmd == 'append':
        before = db.execute('SELECT length(note) FROM reductions WHERE id = ?', (a.id,)).fetchone()
        if before is None:
            die(f'no row with id {a.id}')
        db.execute("UPDATE reductions SET note = COALESCE(note, '') || ? WHERE id = ?",
                   ((' ' if before[0] else '') + a.note.strip(), a.id))
        after = db.execute('SELECT length(note) FROM reductions WHERE id = ?', (a.id,)).fetchone()[0]
        if not after or after <= (before[0] or 0):
            db.rollback(); die(f'APPEND FAILED: note length {before[0]} -> {after}')
        db.commit(); print(f'appended to #{a.id}: note {before[0] or 0} -> {after} chars')
    elif a.cmd == 'crosscheck':
        if not db.execute('SELECT id FROM reductions WHERE id = ?', (a.id,)).fetchone():
            die(f'no row with id {a.id}')
        db.execute('INSERT INTO crosschecks (reduction_id, source, tmid, tmid_err, depth_pct, depth_err, scatter_pct,'
                   ' note, added_at) VALUES (?,?,?,?,?,?,?,?,?)',
                   (a.id, a.source, a.tmid, a.tmid_err, a.depth, a.depth_err, a.scatter, a.note, now()))
        db.commit(); print(f'cross-check on #{a.id} ({a.source})')
    elif a.cmd == 'reason':
        row = db.execute('SELECT reason FROM reductions WHERE id = ?', (a.id,)).fetchone()
        if row is None:
            die(f'no row with id {a.id}')
        db.execute('UPDATE reductions SET reason = ? WHERE id = ?', (a.reason, a.id)); db.commit()
        print(f'#{a.id}: reason {row[0] or "(none)"} -> {a.reason}')
    elif a.cmd == 'export':
        export(db, a.out, a.clocks)
    else:
        for r in db.execute('SELECT * FROM reductions ORDER BY obs_date, id'):
            t = f"Tmid {r['tmid']} +/- {r['tmid_err']}" if r['tmid'] is not None else (r['reason'] or '')
            print(f"#{r['id']} [{r['status']}] {r['target']} {r['obs_date']} {t}")


if __name__ == '__main__':
    main()
