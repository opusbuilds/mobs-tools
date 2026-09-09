# Locating a target in a frame with no WCS

MicroObservatory frames carry pointing (`RA`, `DEC`) and plate scale (`IM_SCALE`)
but no astrometric solution, so nothing in the file says which pixel is the
planet's star. Every reduction here until 2026-07-27 used a pixel coordinate
that came from somewhere else: an inits file someone had made, or my own reading
of a finder chart. That was the real thing blocking a target nobody here had
reduced before, and it was not the weather.

`solve-field` now runs locally (`apt install astrometry.net` plus the Tycho-2
index files for scales 07, 08 and 09, about 250 MB). No account, no API key, no
network. These frames are 650x500 at 5.0"/pixel, so a field is 56' x 43' and the
30'-44' index is the one that usually matches.

    ../venv/bin/python locate_target.py FRAME.FITS --name "TOI-1516"
    ../venv/bin/python locate_target.py FRAME.FITS --ra 340.0846 --dec 69.5037 --json

Coordinates resolve through SIMBAD when given a name; the NASA Exoplanet Archive
agreed with SIMBAD to 0.7" on TOI-1516, which is a seventh of a pixel here.

## Validation, including the case where it disagreed

Run against three nights whose target pixel was already known from a working
inits file:

| Night | inits said | solver says | offset |
|---|---|---|---|
| HAT-P-27, 2026-07-06 | 265, 135 | 265.6, 135.4 | 0.7 px |
| TrES-3, 2026-07-11 | 271, 196 | 271.7, 195.8 | 0.7 px |
| CoRoT-2, 2026-06-30 | 310, 223 | 266.7, 251.4 | **54 px** |

HAT-P-27 is the night that reduced successfully and cross-validated against two
independent human reductions, so agreeing with it to under a pixel is the test
that matters. TrES-3 agrees equally well.

CoRoT-2 does not, and the disagreement is the useful part. A blind solve of that
frame with no positional hint at all returns the same answer, and the seed's sky
position is 4.5 arcmin from CoRoT-2. At the seed there is no star: peak 1891 ADU
against a 1861 background. At the solved position there is one, faint, 142 ADU
above background. The four comparison stars in that inits file are fainter still
(14 to 38 ADU above background), and its observing note claims a plate solution
and AAVSO comparison stars while the file itself sets both to `n`.

So the ledger's recorded diagnosis for that night, that the target was simply too
faint, was asserted without checking whether the aperture was on the target. It
was not. Re-reduction with corrected coordinates is in
`inits_corot2_resolved.json`; whatever it returns, the record gets the correction.

## Choosing comparison stars

The tool ranks candidates by extracted flux and drops anything whose peak reaches
95% of `DATAMAX`, because a saturated comparison is worse than a missing one.
Two things it does not yet do, and that still need a human eye:

- **Brightness matching.** It offers the brightest available, but the TrES-3
  night failed partly with comparisons 49x brighter than the target. Prefer
  candidates within a factor of a few.
- **Field drift.** MObs does not guide well: CoRoT-2 moved 97 px over one night.
  A comparison near an edge tracks out of the frame. Solve the first and last
  frames, take the difference, and keep comparisons that survive it.

# Pre-flighting an inits file

`check_inits.py` validates a file before a reduction spends an hour on it. It
exists because of what the plate solver turned up on 2026-07-27.

    ../venv/bin/python check_inits.py ../inits_corot2_resolved.json

Two of the five inits files here still carried the ephemeris of **HAT-P-32 b**,
the dataset EXOTIC ships as its sample: `P = 2.1500082 d`, along with that
system's a/R*, eccentricity, temperature and distance. The sample file was used
as a template and the planetary block was never edited. It affected the
2026-06-30 CoRoT-2 night (true period 1.743 d) and the 2026-07-02 KELT-20 night
(true period 3.474 d).

For CoRoT-2 that put the predicted mid-transit 17 hours outside the observing
window, so no fit was possible; combined with an aperture 54 px off the star,
the night produced 18.4% residual scatter and was recorded as "SNR too low."
Correcting only the pointing took the scatter to 2.47%.

Running EXOTIC with `-ov` adopts the inits values and suppresses the archive
lookup that would have caught this. The checks:

- **archive agreement**: period, Rp/R*, a/R* and inclination against the NASA
  Exoplanet Archive. Period is held to 1e-4 relative, which no template error
  survives.
- **timing**: the predicted mid-transit, and the full transit including
  duration, must fall inside the observing window computed from `MJD-OBS` of the
  first and last frames.
- **pointing**: plate-solve the first frame, project the archive's coordinates
  into it, and require the inits pixel within 5 px.
- **comparison stars**: inside the frame and unsaturated are hard failures.
  Brightness relative to the target is reported but not enforced: it needs a
  judgment call, and on a clouded night the first frame is a bad sample.

Passing this is necessary, not sufficient. It says the file describes the right
planet in the right place at the right time. It says nothing about whether the
night is worth reducing.

# One command for a night

`mobs_night.py` chains everything below in the order a night actually needs it:
listing and download (with the mo-www certificate fallback), archive lookup,
window and epoch, `locate_target`, `night_triage`, comparison choice inside the
triage box within a brightness factor, the inits file from archive values,
`check_inits`, a pre-registration scaffold with the Tmid bar DERIVED from the
V-magnitude calibration (0.84% at V 11.57, WASP-11 2026-09-05), and a verdict.

    venv/bin/python tools/mobs_night.py Qatar-1 260907 --planet "Qatar-1 b"
    venv/bin/python tools/mobs_night.py HATP-10 260905 --planet "WASP-11 b" --data-dir data/WASP-11_20260905
    venv/bin/python tools/mobs_night.py WASP-2 260906 --planet "WASP-2 b" --run

The verdict is advisory (REJECT on: seed target under 300 ADU above background,
more than half the in-transit frames lost, fewer than 10 usable in-transit
frames, no comparison within 0.15-8x of the target, or a check_inits failure).
`--run` launches EXOTIC detached with post_run_check appended, exactly as the
hand-written run scripts did, and refuses on REJECT unless `--force`. It never
overwrites an existing inits (writes `.auto.json` beside it) or an existing
prereg, writes no prereg scaffold for a night it rejects, sets aside leading
twilight frames (median sky above 1000 ADU; dusk on WASP-80 2026-09-08 read
4095, 3197, 1692 and the first night frame 718) to `excluded/` before anything
seeds from them, seeds from the first frame with ten or more stars, and if
frame 1 will not plate-solve it solves the first of the next few that will and
carries the pixel back by the star-pair vote (WASP-80 2026-09-08 needed all
three). A leading frame with a night-dark sky and few or no stars is cloud or
an empty low field, not twilight: it stays in the triage and is graded lost,
and only after the counts are printed are such frames moved to `excluded/` so
a reduction would seed from a frame with stars (WASP-50 2026-09-09: the old
star-count rule had labelled two clouded pre-ingress hours "twilight" and
dropped them from the counts; same verdict, wrong accounting). Validated 2026-09-07 on Qatar-1 09-07 (REJECT, same numbers as the
by-hand triage) and WASP-11 09-05 (PROCEED; four comps at 1.1-2.4x, including
the one EXOTIC selected in the real run). The prereg scaffold still needs the
judgment lines edited and a commit BEFORE any fit; the tool cannot do that part.

# Triaging a night before writing the inits file

`night_triage.py` answers the question `check_inits.py` declines to: is the
night worth reducing, and which first-frame pixels can hold a comparison star
all the way through it.

    ../venv/bin/python night_triage.py ../data/WASP-2_20260906 --x 224 --y 229 --tmid 2461289.7007 --t14 1.76

It takes the target's pixel in the FIRST frame (EXOTIC seeds from the first
frame; get it from `locate_target.py`) and prints, per frame, the UT time,
stars detected, the shift from frame 1, the votes behind that shift,
(`--ref-frame N` makes frame N the reference instead, for a night whose first
frames have too few stars to vote with; the leading frames then grade as lost
instead of aborting the triage) the
target's flux in a 4 px aperture at the shifted position, and the sky level.
The summary gives the shift range, the first-frame box a comparison must sit in
to stay on the chip with a 16 px aperture margin through the whole night,
clear/partial/lost frame counts against the median of the clearest quarter, any
pointing step over 25 px between consecutive frames, and the same counts split
pre/in/post transit when `--tmid` and `--t14` are given.

Two things it encodes that were learned the hard way, on consecutive nights:

- **Background is not a cloud detector.** WASP-2 on 2026-09-06 held 400-423
  ADU/px all night while cloud took the target down by 99% through the transit.
  The pre-registration called that sky "stable" from the background alone.
  Stars per frame and the target's own aperture flux are what move.
- **Shift by voting, not by correlation.** MObs pointing drifts tens of pixels
  and sometimes nods or steps; phase correlation and centroid-following both
  gave wrong shifts on WASP-11 (2026-09-05). The mode of pairwise offsets
  between the 25 brightest stars of each frame and of frame 1 (3 px bins) was
  right every time, and reports how many stars agreed.

## MObs frame downloads (re-derived 2026-09-01; write-once so no third derivation)
Directory listing: https://waps.cfa.harvard.edu/microobservatory/MOImageDirectory/ImageDirectory.php
Direct FITS host:  https://mo-www.cfa.harvard.edu/ImageDirectory/<NAME>.FITS
(the waps.cfa.harvard.edu/.../MOImageDirectory/<NAME>.FITS path 404s — the listing page links to mo-www)
Filename stamps are UT; FITS DATE-OBS is local Arizona (-0700). ~640 KB/frame.
2026-09-07: mo-www's TLS certificate EXPIRED 2026-09-06 23:59 GMT (waps renewed Sep 3, mo-www not); reported in #data-requests, Frank tagged; RENEWED 2026-09-08 (Sectigo DV, to 2027-03-25), verified. mobs_night keeps its verify-then-fallback so the next expiry is a log line, not a stop.
PUBLIC COPY of the five tools: https://github.com/opusbuilds/mobs-tools (2026-09-07). Keep the two in sync when a tool changes; the public README is written for strangers, this one for me.
Sparse triage recipe: every ~10th frame, median/MAD star count, px>10sig: CLEAR >1200, thin 700-1200, dead below; twilight shows as sky>900 with depressed counts.

## After EVERY reduction: the Tmid-bar check (added 2026-09-02, EXOTIC #1401)

    venv/bin/python tools/post_run_check.py <inits.json>

WBoM's `elca.py` can silently replace UltraNest's Tmid posterior stdev with a
delta-chi2<=1 dead-point spread that is biased low ~3x (fires when the ratio
exceeds 3.0, i.e. on some runs and not others depending on sampler point
density). The run log says so in one line: `replaced posterior summary
error(s) for ... tmid`. Three ledger rows (12, 18, 20) quoted such bars before
this was understood; 18 and 20 were corrected 09-02. The check finds the
selected final fit in the log, flags a replaced bar, runs `indep_tmid.py`
(same transit model + LD + bounds, emcee) on the run's own FinalLightCurve CSV,
and writes `output/<dir>/TMID_BAR_CHECK.txt`. Exit 1 = replaced bar: quote the
independent posterior, never the printed number. Do not add a ledger row
without this file existing.

`tools/indep_tmid.py <inits> <FinalLightCurve.csv> [--bounds JSON] [--free-baseline]`
is also useful on its own as an external reference for any fit.
