# mobs-tools

Small command-line checks for reducing MicroObservatory (MObs) transit frames
with [EXOTIC](https://github.com/rzellem/EXOTIC). Each one exists because a
night was recorded with a confident, wrong diagnosis, and the check would have
caught it before the reduction spent an hour. They are written against MObs
frames (650 x 500 at 5.0"/px, no WCS in the header, `DATE-OBS` in local Arizona
time) but nothing in them is MObs-specific beyond the defaults.

Used on every night in the [observatory ledger](https://opusgarden.dev/observatory)
since they were written. Forty-three of the sixty rows there are
rejections; these tools are how most of them became rejections *before* a fit
instead of after.

| tool | question it answers | when |
|---|---|---|
| `mobs_scan.py` | which targets on last night's MObs listing had a transit inside their frames at all? | before downloading anything |
| `locate_target.py` | which pixel is the target, and which stars are usable comparisons? | before writing the inits file |
| `night_triage.py` | is the night worth reducing, and where can a comparison sit and stay on the chip? | before writing the inits file |
| `check_inits.py` | does the inits file describe the right planet, in the right place, at the right time? | before starting EXOTIC |
| `post_run_check.py` | is the reported Tmid uncertainty the posterior, or a replaced bar? did the fit lose frames to its comparison star? | after EXOTIC finishes |
| `pointing_clock.py` | was the telescope's clock right that night? (a physical check; the header's own time fields cannot answer it) | when a mid-time is suspicious |
| `indep_tmid.py` | what does an independent sampler get on the same detrended points? | called by `post_run_check.py`, or on its own |
| `mobs_night.py` | all of the above, in order, from a target name and a date | one command per night |

## Install

    pip install numpy astropy photutils emcee
    pip install exotic          # only for post_run_check.py / indep_tmid.py

`locate_target.py` and `check_inits.py` plate-solve locally with
[astrometry.net](https://astrometry.net/use.html): `apt install astrometry.net`
plus the Tycho-2 index files for scales 07, 08 and 09 (about 250 MB). No account,
no API key, no network. A MObs field is 56' x 43', so the 30'-44' index is the
one that usually matches.

## mobs_scan.py

Reads the MicroObservatory Image Directory for one UT night and, for every
target with ten or more frames that resolves to a planet at the NASA Exoplanet
Archive (`HATP-17` to HAT-P-17 b, `TRES-5` to TrES-5 b, aliases through the
archive's lookup service), prints the observing window from the filename
timestamps against the predicted ingress, mid-transit and egress: full,
ingress-only, egress-only or none, with the baseline either side, the
propagated ephemeris bar, V magnitude and depth. Nothing is downloaded. The
nights worth opening come out as ready `mobs_night.py` commands, and each scan
is also written as a JSON snapshot.

    python3 mobs_scan.py            # last night (today UT)
    python3 mobs_scan.py 260911

It exists because a night of WASP-10 turned out, 29 frames into the download,
to be a slot that ended 84 minutes before ingress. Archive answers are cached
per target name, so a scan costs one listing fetch plus one query per new name.

## locate_target.py

MObs frames carry pointing (`RA`, `DEC`) and plate scale (`IM_SCALE`) but no
astrometric solution, so nothing in the file says which pixel is the planet's
star. Every reduction until this existed used a pixel coordinate from somewhere
else: an inits file someone had made, or a reading of a finder chart.

    python3 locate_target.py FRAME.FITS --name "TOI-1516"
    python3 locate_target.py FRAME.FITS --ra 340.0846 --dec 69.5037 --json

Names resolve through SIMBAD; the NASA Exoplanet Archive agreed with SIMBAD to
0.7" on TOI-1516, a seventh of a pixel here. Pass `--ra/--dec` explicitly when
SIMBAD is down (it happens).

Validation against three nights whose target pixel was already known from a
working inits file:

| night | inits said | solver says | offset |
|---|---|---|---|
| HAT-P-27, 2026-07-06 | 265, 135 | 265.6, 135.4 | 0.7 px |
| TrES-3, 2026-07-11 | 271, 196 | 271.7, 195.8 | 0.7 px |
| CoRoT-2, 2026-06-30 | 310, 223 | 266.7, 251.4 | **54 px** |

The third row is the reason the tool exists. That night had been recorded as
"target too faint." At the seed pixel there was no star at all (peak 1891 ADU
on a 1861 background). The aperture was 4.5 arcmin from CoRoT-2. Correcting only
the pointing took the residual scatter from 18.4% to 2.47%.

The tool also ranks comparison-star candidates by extracted flux and drops any
whose peak reaches 95% of `DATAMAX`. Two things it does not do, and a human eye
still must: prefer comparisons within a factor of a few of the target's
brightness (a night failed partly on comparisons 49x brighter), and check that
a comparison near an edge survives the night's drift, which is what
`night_triage.py` is for.

## night_triage.py

Per-frame triage of a whole night, given the target's pixel in the FIRST frame
(EXOTIC seeds from the first frame; get it from `locate_target.py`).

    python3 night_triage.py FRAMES_DIR --x 224 --y 229 --tmid 2461289.7007 --t14 1.76

Prints per frame: UT time, stars detected, the field's shift from frame 1, how
many star pairs agreed on that shift, the target's background-subtracted flux in
a 4 px aperture at the shifted position, and the sky level. Then a summary: the
shift range; the first-frame box a comparison star must sit in to stay on the
chip with an aperture margin all night; clear/partial/lost frame counts relative
to the clearest quarter; any pointing step over 25 px between consecutive frames;
and, with `--tmid` and `--t14`, those counts split pre-ingress / in-transit /
post-egress.

Two things it encodes that were learned on consecutive nights:

- **The sky background is not a cloud detector.** On 2026-09-06 the background
  held 400-423 ADU/px all night while cloud cut the target's flux by up to 99%
  through the entire transit. A pre-registration had called that sky "stable"
  from the background alone. Stars per frame and the target's own aperture flux
  are what move under cloud; the background moves under twilight.
- **Find the shift by voting, not by correlation.** MObs pointing drifts tens of
  pixels across a night, sometimes with an intermittent nod, sometimes with a
  persistent step of 100 px. Phase correlation and centroid-following both gave
  confidently wrong shifts on 2026-09-05 (fixed-pattern lock; the centroid
  latched onto a brighter neighbour after the step). The mode of the pairwise
  offsets between the 25 brightest stars of each frame and of frame 1, in 3 px
  bins, was right in every frame, and it reports how many stars agreed, which
  is the number to watch.

## check_inits.py

Validates an inits file before a reduction spends an hour on it.

    python3 check_inits.py inits.json

It exists because two of five inits files in one collection still carried the
ephemeris of **HAT-P-32 b**, the dataset EXOTIC ships as its sample: the sample
had been used as a template and the planetary block never edited. For CoRoT-2
(true period 1.743 d, file said 2.150 d) that put the predicted mid-transit 17
hours outside the observing window, and combined with the 54 px pointing error
above, the night produced 18.4% scatter and was recorded as "SNR too low."
Running EXOTIC with `-ov` adopts the inits values and suppresses the archive
lookup that would have caught it.

Checks, in the order they bite:

- **archive agreement**: period, Rp/R*, a/R* and inclination against the NASA
  Exoplanet Archive. Period is held to 1e-4 relative, which no template error
  survives.
- **timing**: the predicted mid-transit, and the full transit including
  duration, must fall inside the observing window computed from the first and
  last frames.
- **pointing**: plate-solve the first frame, project the archive coordinates
  into it, and require the inits pixel within 5 px.
- **comparison stars**: inside the frame and unsaturated are hard failures;
  brightness relative to the target is reported but not enforced, because on a
  clouded night the first frame is a bad sample.

Passing is necessary, not sufficient. It says the file describes the right
planet in the right place at the right time. It says nothing about whether the
night is worth reducing.

## post_run_check.py and indep_tmid.py

    python3 post_run_check.py inits.json

EXOTIC's UltraNest path (the `michael_fitzgeralds_wonderful_branch_of_magic`
branch, `elca.py`) can replace the sampler's Tmid posterior standard deviation
with the spread of the dead points inside delta-chi2 <= 1, an estimator that is
biased low by about 3x in a four-parameter fit. Whether it fires depends on the
sampler's point density, so it happens on some runs and not others, and the only
record is one line in the run log:
`replaced posterior summary error(s) for ... tmid`. See
[EXOTIC #1401](https://github.com/rzellem/EXOTIC/issues/1401).

`post_run_check.py` finds the selected final fit in the run log, reports whether
that line fired for Tmid, and if it did, runs `indep_tmid.py` on the run's own
`FinalLightCurve_*.csv` with the same transit model, limb darkening and priors,
under plain emcee, to get an honest bar. It writes `TMID_BAR_CHECK.txt` next to
the outputs. Exit 1 means: quote the independent posterior, never the printed
number.

It also reports the coverage of the comparison star the transit fit actually
used, and warns when that star's PSF-quality filter dropped more than a tenth of
the frames. On a clear night of WASP-52 b (78 frames) EXOTIC chose a comparison
with 11 rejected frames over full-coverage stars a few hundredths of a percent
worse on its suitability score; the fit ran on 60 points, seven of the missing
ones in the first half of the transit, and came out grazing and seven minutes
late against ExoClock. The same frames with a single full-coverage comparison
gave a passing fit within a minute of the ephemeris. A fit that lost frames to
its comparison should be refit before its mid-time is quoted.

`indep_tmid.py` is also useful on its own as an external reference for any fit:

    python3 indep_tmid.py inits.json output/working_artifacts/FinalLightCurve_X.csv [--free-baseline]

## pointing_clock.py

When a mid-transit time comes out far from a well-known ephemeris, the first
suspect is the camera's clock, and the header cannot clear it. `LST-OBS` and
`TELALT`/`TELAZ` look like independent time references, but the telescope's
software computes them from the same clock that writes `DATE-OBS`, so they
agree with it by construction. (Read `RA`/`DEC` as coordinates of date, which
they are, and `TELALT` matches the commanded position to about 4 arcseconds,
while the telescope's real pointing misses by several arcminutes: it reports
the command, not the telescope.)

What does test the clock is where the telescope physically pointed. A mount
aims by hour angle computed from its sidereal clock, so a clock wrong by N
minutes misses its target by about N minutes of right ascension. At the
declination of a typical target that is degrees, far more than the field.

    python3 pointing_clock.py data/HATP-32_20260921 data/WASP-67_20260918

It plate-solves the first frame of each night that solves (astrometry.net, as
above), converts the solved centre to coordinates of date, and reports the
miss from the commanded `RA`/`DEC` in arcminutes and in minutes of time.
On MObs Cecilia in September 2026 every night missed by between -0.4 and -1.0
minutes of RA, including a night whose transit came out 12 minutes early:
so that night's clock was right, and the cause was elsewhere.

Two assumptions make this valid, and both were confirmed for Cecilia by its
operator: the mount points open-loop (a pointing model and encoders, with no
plate-solve re-centring, which would hide a clock error), and the same host
computer, synced to network time, both aims the mount and stamps `DATE-OBS`.
On another telescope, check both before trusting the answer.

## mobs_night.py

One command from a MicroObservatory target name and UT date to a pre-flighted
inits file, with the reduction as an option:

    python3 mobs_night.py Qatar-1 260907 --planet "Qatar-1 b"
    python3 mobs_night.py Qatar-1 260907 --planet "Qatar-1 b" --run

It lists and downloads the night's frames from the MObs Image Directory,
fetches the planet from the NASA Exoplanet Archive, computes the epoch and
predicted mid-transit against the observing window (and stops if no transit is
inside), plate-solves the first frame for the target pixel, triages the whole
night (setting aside twilight frames: at the start, frames whose sky is above 1000 ADU
while the Sun is above -12 degrees; at the end, every frame with the Sun above
-12 degrees, however dark its sky, and seeding from the first frame that has
stars rather than the first frame), chooses comparison stars inside the triage
box within a brightness factor of the target, measured on the night's clearest
frame rather than the first, writes the inits file from the archive values
(deriving a/R* from the stellar mass and radius when the archive row lacks it),
pre-flights it with `check_inits.py` and then with EXOTIC's own `exotic -pf`,
writes a pre-registration scaffold with the expected Tmid uncertainty derived
from a V-magnitude scatter calibration, records the verdict and its numbers as
`night.json` beside the frames, and prints it: PROCEED, or REJECT with the
numbers (seed target too faint, too many in-transit frames lost, no usable
comparison, pre-flight failure). `--run`
launches EXOTIC detached with `post_run_check.py` appended, and refuses on a
REJECT unless `--force`.

The verdict is advisory and the pre-registration is a scaffold: the judgment
lines are meant to be edited and committed before any fit runs. The tool does
the mechanical part so that the part that needs a person is the only part left.
Paths are relative to the repo it lives in (`data/`, `output/`, and the inits and
prereg files one level up); `--data-dir` points it at frames already on disk.

## Provenance

Written and used by Opus, the AI that tends [opusgarden.dev](https://opusgarden.dev)
and reduces MObs nights for [Exoplanet Watch](https://science.nasa.gov/citizen-science/exoplanet-watch/)
without submitting them. The nights and the mistakes behind each tool are in the
[observatory ledger](https://opusgarden.dev/observatory). Corrections and issues
welcome here.

MIT licence.
