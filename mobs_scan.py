#!/usr/bin/env python3
"""
Scan the MicroObservatory Image Directory for one UT night and say, per
exoplanet target, whether the frames hold a transit at all, before anything is
downloaded.

    python3 tools/mobs_scan.py            # last night (today's UT date)
    python3 tools/mobs_scan.py 260911

For every target with enough frames it maps the MObs name to an archive planet
name (HATP-17 -> HAT-P-17 b, TOI2570 -> TOI-2570 b, TRES-5 -> TrES-5 b, falling
back to the archive's alias service), pulls the ephemeris from pscomppars, reads
the window from the filename timestamps, and prints one line: frames, window,
predicted ingress/mid/egress, coverage (full / ingress-only / egress-only /
none), pre and post baseline, V and depth. The line ends with the mobs_night
command for the nights worth opening.

Non-exoplanet targets (nebulae, planets, darks) resolve to nothing at the
archive once and are cached as such in tools/.mobs_scan_cache.json, so a scan
costs one listing fetch plus one archive query per NEW target name.
Nothing here downloads a frame or submits anything anywhere.
"""
import json, os, re, sys, urllib.parse, urllib.request, collections
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
LISTING = 'https://waps.cfa.harvard.edu/microobservatory/MOImageDirectory/ImageDirectory.php?SortBy=Filename&SortPos=DESC'
TAP = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?'
ALIAS = 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/Lookup/nph-aliaslookup.py?objname='
CACHE = os.path.join(HERE, '.mobs_scan_cache.json')
MIN_FRAMES = 10


def jd(yymmdd, hms):
    t = datetime.strptime(yymmdd + hms, '%y%m%d%H%M%S').replace(tzinfo=timezone.utc)
    return t.timestamp() / 86400.0 + 2440587.5


def ut(j):
    return datetime.fromtimestamp((j - 2440587.5) * 86400, tz=timezone.utc).strftime('%H:%M')


def candidates(name):
    """Archive-style spellings to try for a MObs target name."""
    n = name.strip('-')
    out = [n]
    m = re.match(r'^([A-Za-z]+)-?(\d.*)$', n)
    if m:
        pre, num = m.groups()
        fixed = {'HATP': 'HAT-P', 'HATS': 'HATS', 'TOI': 'TOI', 'TRES': 'TrES', 'COROT': 'CoRoT', 'WASP': 'WASP',
                 'QATAR': 'Qatar', 'KELT': 'KELT', 'XO': 'XO', 'GJ': 'GJ', 'KEPLER': 'Kepler', 'K2': 'K2', 'HD': 'HD',
                 'HIP': 'HIP', 'TRAPPIST': 'TRAPPIST', 'NGTS': 'NGTS', 'WD': 'WD', 'LHS': 'LHS', 'LTT': 'LTT'}
        p = fixed.get(pre.upper(), pre)
        out += [f'{p}-{num}', f'{p} {num}']
    seen, res = set(), []
    for c in out:
        for s in (c + ' b', c):
            if s not in seen:
                seen.add(s); res.append(s)
    return res


def tap(q):
    url = TAP + urllib.parse.urlencode({'query': q, 'format': 'json'})
    return json.loads(urllib.request.urlopen(url, timeout=60).read().decode() or '[]')


COLS = 'pl_name,hostname,pl_orbper,pl_orbpererr1,pl_tranmid,pl_tranmiderr1,pl_trandur,pl_ratror,pl_trandep,pl_radj,st_rad,sy_vmag'


def lookup(name, cache):
    if name in cache:
        return cache[name]
    row = None
    for c in candidates(name):
        rows = tap(f"select {COLS} from pscomppars where pl_name='{c}' or (hostname='{c}' and pl_letter='b')")
        if rows:
            row = rows[0]; break
    if row is None:
        # the archive's alias service, as EXOTIC does since #1412 (HAT-P-10 b -> WASP-11 b)
        for c in candidates(name)[:3]:
            try:
                a = json.loads(urllib.request.urlopen(ALIAS + urllib.parse.quote(c), timeout=30).read().decode())
            except Exception:
                continue
            resolved = ((a.get('manifest') or {}).get('resolved_name')) or None
            if resolved:
                rows = tap(f"select {COLS} from pscomppars where pl_name='{resolved}' or (hostname='{resolved}' and pl_letter='b')")
                if rows:
                    row = rows[0]; break
    cache[name] = row
    return row


def main():
    yymmdd = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime('%y%m%d')
    html = urllib.request.urlopen(LISTING, timeout=90).read().decode(errors='ignore')
    by = collections.defaultdict(list)
    for name, hms in set(re.findall(r'ImageDirectory/([A-Za-z0-9\-]+?)' + yymmdd + r'(\d{6})\.FITS', html)):
        by[name].append(hms)
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    print(f'MObs {yymmdd}: {len(by)} target names in the listing, {sum(len(v) for v in by.values())} frames')
    # The MObs night in Arizona runs until dawn, ~12:00-13:30 UT depending on season. A scan of TODAY's UT
    # date made before then sees a night still in progress: on 2026-09-13 the 09:07 UT scan showed 35 of
    # TOI-3693's eventual 102 frames, and the night was rejected as ingress-only on that partial listing.
    now = datetime.now(timezone.utc)
    if yymmdd == now.strftime('%y%m%d') and now.hour < 14:
        print(f'  NOTE: scanned at {now:%H:%M} UT on the same UT date; the MObs night runs until ~13:30 UT, so this '
              f'listing may still be growing. Do not reject a night on geometry from it; re-scan after 14:00 UT.')
    lines, worth = [], []
    for name, hms in sorted(by.items()):
        if len(hms) < MIN_FRAMES:
            continue
        r = lookup(name, cache)
        if not r:
            continue
        hms.sort()
        j0, j1 = jd(yymmdd, hms[0]), jd(yymmdd, hms[-1])
        P, T0 = r['pl_orbper'], r['pl_tranmid']
        if not P or not T0:
            lines.append(f'  {name:12s} {len(hms):3d} frames  {ut(j0)}-{ut(j1)}  {r["pl_name"]}: no ephemeris in pscomppars')
            continue
        n = round(((j0 + j1) / 2 - T0) / P)
        tmid = T0 + n * P
        dur_h = r['pl_trandur'] or 2.0
        half = dur_h / 48.0
        ing, egr = tmid - half, tmid + half
        pre, post = (ing - j0) * 1440, (j1 - egr) * 1440
        bar = ((r['pl_tranmiderr1'] or 1e-3) + abs(n) * (r['pl_orbpererr1'] or 1e-6)) * 1440
        if egr < j0 or ing > j1:
            cov = 'NONE'
        elif ing > j0 and egr < j1:
            cov = 'full'
        elif ing > j0:
            cov = 'ingress-only'
        else:
            cov = 'egress-only'
        rprs = r['pl_ratror']
        if rprs is None and r['pl_radj'] and r['st_rad']:
            rprs = r['pl_radj'] * 0.10045 / r['st_rad']
        depth = f'{100 * rprs ** 2:.2f}%' if rprs else (f'{r["pl_trandep"]:.2f}%' if r['pl_trandep'] else '?')
        v = f'V {r["sy_vmag"]:.2f}' if r['sy_vmag'] else 'V ?'
        last_dt = datetime.fromtimestamp((j1 - 2440587.5) * 86400, tz=timezone.utc)
        growing = yymmdd == now.strftime('%y%m%d') and now.hour < 14 and (now - last_dt).total_seconds() < 3 * 3600
        line = (f'  {name:12s} {len(hms):3d} frames  {ut(j0)}-{ut(j1)} UT  {r["pl_name"]:14s} '
                f'ingress {ut(ing)} mid {ut(tmid)} egress {ut(egr)}  {cov:12s} pre {pre:+.0f} post {post:+.0f} min  '
                f'bar {bar:.0f} min  {v}  depth {depth}' + ('  [last frame < 3 h ago: still running?]' if growing else ''))
        lines.append(line)
        if cov != 'NONE':
            worth.append(f'  venv/bin/python tools/mobs_night.py {name} {yymmdd} --planet "{r["pl_name"]}"')
    json.dump(cache, open(CACHE, 'w'), indent=0, sort_keys=True)
    print('\n'.join(lines) if lines else '  no exoplanet target with %d+ frames' % MIN_FRAMES)
    if worth:
        print('worth opening:')
        print('\n'.join(worth))
    else:
        print('nothing to open: no listed target holds a transit')


if __name__ == '__main__':
    main()
