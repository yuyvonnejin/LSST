# lsst_watch -- Rubin/LSST Alert-Stream Watcher for Nearby M Dwarfs

Monitor flare activity on nearby M dwarfs using the public Rubin/LSST
alert stream (via community brokers), with ZTF as the validation
backend while Rubin coverage ramps up.

Companion project to `../astro_calib` (proximity-first exoplanet
survey). Start at [docs/design_doc.md](docs/design_doc.md) for the
full rationale, architecture, and milestones.

## Why the alert stream, not the catalog

- The LSST survey started 2026-06-30. Its catalogs (DP1/DP2/DR1) are
  restricted to data-rights holders for 2 years; we have none.
- The nightly alert stream (~7M alerts/night) is fully public via
  community brokers (ALeRCE, Fink, ANTARES, ...), no auth required.
- LSST saturates at r ~ 16, so astro_calib's bright FGK targets are
  invisible to it. Nearby M dwarfs in the r ~ 16-24.5 window are the
  natural sample, and their flare activity matters directly for
  exoplanet habitability and follow-up prioritization.

## Pipeline

```
Gaia DR3 TAP -> targets.py -> data/targets.json
                                   |
broker.py (ALeRCE REST) <----------+   per target:
  cone search (PM-propagated position, drift-aware radius)
  -> nearest match -> light curve
  -> flares.py (baseline, robust scatter, brightening events)
  -> output/target_<id>.json + output/summary.json
```

Modules:
- `src/targets.py` -- Gaia DR3 TAP query for M dwarfs within 50 pc,
  dec < +32, G in the LSST-unsaturated window; proper-motion
  propagation (Gaia epoch 2016 positions drift up to arcminutes).
- `src/broker.py` -- ALeRCE REST client (cone search + light curves).
  Base URL is swappable for the ALeRCE LSST service / Fink later.
- `src/flares.py` -- flare-candidate detection: median baseline,
  MAD scatter floored by reported errors, 3-sigma + 0.1 mag
  threshold, epoch grouping. Single-epoch events are always labeled
  candidates, never confirmed flares.
- `src/watch.py` -- orchestration + CLI; one target failing never
  kills the run.

## Setup

```
python -m venv venv
./venv/Scripts/pip install -r requirements.txt
```

## Usage

```
# Build the target catalog (live Gaia TAP query, ~30 s)
./venv/Scripts/python.exe -m src.targets

# Watch targets against ZTF alerts (long baseline, northern sky)
./venv/Scripts/python.exe -m src.watch --survey ztf --limit 20

# Watch targets against live Rubin/LSST alerts
./venv/Scripts/python.exe -m src.watch --survey lsst --limit 20

# Custom target file
./venv/Scripts/python.exe -m src.watch --targets-file data/targets.json
```

Output goes to `output/<survey>/`. The two backends share one schema
downstream: LSST detections arrive as fluxes (nJy) and are converted
to AB magnitudes (from scienceFlux, the star's total brightness;
difference-image psfFlux is not used for flare amplitudes).

Reading the output:
- `status: ok` -- broker object matched, light curve analyzed.
- `status: no_match` -- no broker object at the star's position.
  Common and expected: alerts only exist when a star varied, and
  ZTF does not cover the far southern sky. As Rubin alerts
  accumulate, southern quiet stars will start matching.
- `bands.<b>.status: insufficient_data` -- fewer than 10 detections
  in that band; no flare statistics are fabricated from it.

## Rerunning and scheduled runs

- Interactive: type `/watch-run` in a Claude Code session in this
  folder. Runs the survey, archives the previous summary, reports
  what changed (new matches, new candidates, coverage growth).
- Scheduled: Windows Task Scheduler task `lsst_watch_weekly`
  (Mondays 09:37, runs later if the PC was off) executes
  `scripts/weekly_watch.ps1`: archives the previous summary to
  `output/history/`, runs the full-catalog LSST survey, then writes a
  Claude-generated delta report to `logs/report_<date>.md`.
  Manage with: `schtasks /query /tn lsst_watch_weekly`,
  `/run` to trigger now, `/delete` to remove.

## Tests

```
# Offline suite (default; ~2 s)
./venv/Scripts/python.exe -m pytest tests/ -q

# Live smoke tests (ALeRCE + Gaia TAP)
./venv/Scripts/python.exe -m pytest tests/ -m network -q
```

48 offline tests: target cuts, proper-motion propagation, broker
response parsing on captured fixtures, flare detection on synthetic
light curves (quiet null, single/multi-epoch injection, dimming
rejection, grouping), watch orchestration with mocked broker
(match/no-match/error paths, report writing).

## Status (2026-07-04)

- Offline suite green (60 passed), live smoke tests green (ZTF,
  LSST, Gaia TAP).
- End-to-end validated on ZTF: 80-target live run, 39 matched,
  0 errors, honest nulls.
- Rubin backend live: ALeRCE multisurvey API serves real LSST alerts
  (first alerts from 2026-07-02 retrieved and parsed). Full-catalog
  Rubin run results in output/lsst/summary.json. Expect mostly
  no_match for months: the survey is days old and alerts only exist
  where difference imaging fired; coverage deepens nightly.
