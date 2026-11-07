"""HTTP helpers shared by mobs_night.py, mobs_scan.py and check_inits.py (one copy since 2026-11-07; before that each
tool carried its own, and only mobs_night's had the archive's browser-challenge fallback)."""
import io, os, urllib.error, urllib.request

# MObs's listing began refusing Python-urllib's default User-Agent with 403 on 2026-10-21; every fetch here identifies
# itself instead (an honest UA, not a browser disguise).
UA = 'mobs-tools/1.0 (+https://github.com/opusbuilds/mobs-tools)'
def _ua_urlopen(url, *a, **kw):
    if isinstance(url, str):
        url = urllib.request.Request(url, headers={'User-Agent': UA})
    return urllib.request.urlopen(url, *a, **kw)


def _urlopen_retry(url, timeout, tries=2):
    """The NASA Exoplanet Archive occasionally stalls for a full timeout and then
    answers the next request in two seconds (first seen 2026-09-13 from the home
    box). One retry covers it; a second failure is a real outage."""
    import time
    for i in range(tries):
        try:
            return _ua_urlopen(url, timeout=timeout)
        except urllib.error.HTTPError as e:
            # 2026-09-17: the archive put a Cloudflare browser challenge in front
            # of TAP/sync (403, cf-mitigated: challenge). A real browser passes it;
            # fetch through headless Chromium and hand back a file-like object.
            if e.code == 403 and e.headers.get('cf-mitigated') == 'challenge' and 'exoplanetarchive' in url:
                import subprocess
                # Optional: a browser-backed fetcher (see browser/tap-fetch.mjs in
                # opus-infra; set MOBS_TAP_FETCH to its path). Without one, the 403
                # is raised as before.
                fetcher = os.environ.get('MOBS_TAP_FETCH', '/opt/opus-infra/browser/tap-fetch.mjs')
                r = subprocess.run(['node', fetcher, url], capture_output=True, timeout=150) if os.path.exists(fetcher) else None
                if r is not None and r.returncode == 0 and r.stdout.strip():
                    print('  archive: API behind a browser challenge; fetched through headless Chromium instead')
                    return io.BytesIO(r.stdout)
            if i == tries - 1:
                raise
            time.sleep(3)
        except (TimeoutError, OSError) as e:
            if i == tries - 1:
                raise
            time.sleep(3)

