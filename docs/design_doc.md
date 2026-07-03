# Design Doc: lsst_watch -- Rubin/LSST Alert-Stream Watcher for Nearby M Dwarfs

Date: 2026-07-03
Status: v1 implementation (this session)

## 1. Context and data-access reality

- 2026-06-30: Rubin Observatory started the 10-year LSST survey.
  What went online is the survey itself and its nightly alert stream
  (~7 million alerts/night), not a queryable public catalog.
- Catalog products (DP1, DP2 ~Q3 2026, DR1 ~2028) sit behind a 2-year
  proprietary period for LSST data-rights holders. We have no data
  rights, so direct Rubin Science Platform / TAP catalog access is
  NOT a viable path.
- The alert stream IS public, no proprietary period. It is served by
  community brokers with open REST APIs: ALeRCE, Fink, ANTARES,
  Lasair, and others. Fink and ALeRCE have processed Rubin alerts
  since Feb 2026.
- Practical consequence: any contribution from this machine must go
  through broker APIs. Brokers also serve ZTF alerts (7+ years of
  data), which use the same access patterns and let us validate the
  whole chain today while Rubin coverage ramps up over months.
- Verified live (2026-07-03): ALeRCE ZTF API
  (https://api.alerce.online/ztf/v1/) answers cone searches with real
  data. Fink API did not respond to probes from this network; it is a
  secondary backend to add later. ALeRCE's Rubin/LSST endpoint is not
  yet documented publicly; the client keeps the base URL configurable
  so it can be pointed at the LSST service when published.

## 2. Science case and link to astro_calib

astro_calib (../astro_calib) is a proximity-first exoplanet program:
survey the nearest stars, report honest detection limits. Its RV
targets are bright FGK dwarfs (V 1-9).

- LSST saturates at roughly r ~ 16 in a standard 30 s visit. Bright
  nearby FGK stars are unusable; do not point LSST tooling at the
  astro_calib FGK shell catalog.
- Nearby mid/late M dwarfs land in LSST's dynamic range
  (r ~ 16-24.5). They are also the highest-value exoplanet hosts for
  detection (deep transits, large RV signals) and the population
  where flare activity is decisive: flares set habitability
  constraints, contaminate transit photometry, and correlate with the
  RV jitter floor.
- Contribution ladder (mirrors astro_calib's):
  1. Working public-data pipeline: uniform flare monitoring for a
     defined nearby M dwarf sample (value even with zero events).
  2. Per-target activity summaries (flare candidate rate, amplitudes)
     usable to prioritize exoplanet follow-up.
  3. Flagging unusual activity on nearby stars for community
     follow-up (not discovery claims).
  4. Positioned for Rubin scale-up: same code path, deeper cadence,
     southern sky -- where ZTF coverage is poor and LSST adds the
     most.

## 3. Architecture

```
Gaia DR3 TAP (public)            broker REST API (public)
      |                                 |
      v                                 v
 targets.py  ---- targets.json ---> broker.py
 (select nearby M dwarfs           (cone search, object
  in LSST magnitude window,         light curves; ALeRCE now,
  propagate proper motion)          Fink/LSST later)
                                        |
                                        v
                                    flares.py
                                   (baseline, robust sigma,
                                    brightening events)
                                        |
                                        v
                                    watch.py
                                   (orchestrate per target,
                                    JSON report + summary)
```

Modules:

- `src/targets.py`
  - ADQL query to the Gaia DR3 TAP sync endpoint (plain HTTP POST,
    no astroquery dependency): parallax, magnitude, color, proper
    motion, quality cuts.
  - Default sample: parallax > 20 mas (d < 50 pc), bp_rp > 2.0
    (approx M0 and later), M_G > 8 (dwarfs, rejects giants),
    G >= 15.5 (proxy for r >~ 16.5, avoids saturation; for M dwarfs
    r is fainter than G so this is conservative), dec < +32
    (LSST footprint), parallax_over_error > 10, ruwe < 1.4.
  - Writes `data/targets.json`. Selection cuts are pure functions,
    offline-testable.
  - `propagate_position(ra, dec, pmra, pmdec, dt_years)`: Gaia
    coordinates are epoch 2016.0; nearby M dwarfs have proper motion
    up to ~10 arcsec/yr, so a 2016 position can be tens of arcsec off
    by 2026. Positions must be propagated to the observation epoch
    before cone searching.

- `src/broker.py`
  - Thin client over broker REST APIs. Functions:
    `cone_search(ra_deg, dec_deg, radius_arcsec)` -> list of objects,
    `get_lightcurve(object_id)` -> list of detections
    (mjd, band, mag, magerr).
  - Backend = ALeRCE ZTF v1 (verified working). Base URL is a module
    constant, overridable, so the ALeRCE LSST service or Fink can be
    swapped in without touching callers.
  - Response parsing is separated from HTTP so tests run offline on
    captured JSON fixtures.

- `src/flares.py`
  - Input: one band's detections (mjd, mag, magerr).
  - Quiescent baseline: median magnitude. Scatter: MAD * 1.4826,
    floored by the median reported magerr.
  - Flare candidate: detection brighter (lower mag) than baseline by
    > max(k * scatter, min_amplitude), k = 3, min_amplitude = 0.1 mag.
    Consecutive candidate epochs within `group_gap_days` merge into
    one event.
  - Per event: peak amplitude (delta mag), peak flux ratio, time
    span, n_points. Sparse-cadence honesty: a single-epoch
    brightening is a "flare candidate", never a confirmed flare;
    the report says which.
  - Guardrails: needs >= `min_epochs` (default 10) quiescent points,
    otherwise returns "insufficient data" instead of fake results.

- `src/watch.py`
  - For each target: propagate position to survey mid-epoch, cone
    search (radius = base 2 arcsec + PM drift margin over the survey
    span), pick nearest match, fetch light curve, run flare detection
    per band, write one scorecard per target plus a survey summary
    JSON to `output/`.
  - One target failing never kills the run (same rule as
    astro_calib's survey driver).
  - CLI: `python -m src.watch [--limit N] [--targets-file path]`.

## 4. What is deliberately out of scope for v1

- Kafka/streaming subscription to the live alert firehose (brokers
  expose it, but polling REST per target is enough at this sample
  size).
- Flare energy in erg (needs quiescent luminosity per star; possible
  later by joining Gaia photometry, kept out of v1).
- Fink/ANTARES/Lasair backends (interface allows adding them).
- Any claim of catalog-based science (no data rights).

## 5. Test plan

Offline (default `pytest`, no network):
- `test_targets.py`: selection cuts on synthetic rows (parallax,
  color, magnitude, dec edges); proper-motion propagation against
  hand-computed values including a Barnard's-star-scale case;
  targets.json round-trip.
- `test_broker.py`: parsing of captured ALeRCE cone-search and
  light-curve JSON fixtures; radius/parameter construction; empty
  results; malformed rows skipped, not fatal.
- `test_flares.py`: synthetic quiet light curve -> zero events;
  injected single-epoch and multi-epoch flares recovered with correct
  amplitude and grouping; insufficient-data guard; noise floor uses
  magerr when scatter underestimates.
- `test_watch.py`: full per-target flow with a mocked broker
  (monkeypatched HTTP): match found, no match, broker error ->
  target marked failed, run continues; report file written.

Network (marked `@pytest.mark.network`, run manually):
- Live ALeRCE cone search on a known ZTF object returns rows.
- Live Gaia TAP query with tight limits returns M dwarf rows.

## 6. Validation targets

- Known active M dwarf with ZTF coverage in the usable magnitude
  window (e.g. a mid-M flare star fainter than r ~ 15) -- expect
  detections retrievable and at least plausible candidate events.
- A photometrically quiet M dwarf -- expect zero events (honest
  null, same philosophy as Tau Ceti in astro_calib).

## 7. Milestones

1. v1 (this session): modules + offline tests green + live smoke of
   the ZTF path end to end on a handful of targets.
2. Point broker client at ALeRCE/Fink Rubin endpoints as they are
   published; first Rubin-based target scorecards.
3. Flare rates per target normalized by exposure (epochs, baseline
   span); compare northern (ZTF) vs southern (LSST-only) samples.
4. Optional: energy calibration, streaming ingestion, candidate
   bulletin generation.
