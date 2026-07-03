"""Watch orchestration: broker alerts + flare analysis per target.

For each target from data/targets.json:
1. Propagate the Gaia epoch-2016 position to the survey mid-epoch.
2. Cone search the broker; radius covers residual proper-motion drift.
3. Pick the nearest matched object, fetch its light curve.
4. Run flare detection per band.
5. Write one scorecard per target + a survey summary JSON.

One target failing never kills the run.
"""

import argparse
import json
import logging
import math
from pathlib import Path

from . import broker
from .flares import analyze_by_band
from .targets import (GAIA_EPOCH, TARGETS_FILE, load_targets,
                      pm_drift_arcsec, propagate_position)

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output"

# ZTF started 2018.2; use its mid-epoch for position propagation.
# For a Rubin backend this becomes ~2026.7 onward.
SURVEY_START_YEAR = 2018.2
SURVEY_END_YEAR = 2026.5

BASE_RADIUS_ARCSEC = 2.0


def search_radius_arcsec(target, mid_epoch_dt, half_span_years):
    """Base match radius plus proper-motion drift over half the survey."""
    drift = pm_drift_arcsec(target.get("pmra_mas_yr"),
                            target.get("pmdec_mas_yr"),
                            half_span_years)
    return BASE_RADIUS_ARCSEC + drift


def angular_sep_arcsec(ra1, dec1, ra2, dec2):
    """Small-angle separation, fine at arcminute scales."""
    dra = (ra1 - ra2 + 180.0) % 360.0 - 180.0
    dra *= math.cos(math.radians(0.5 * (dec1 + dec2)))
    ddec = dec1 - dec2
    return math.hypot(dra, ddec) * 3600.0


def watch_target(target, survey_start=SURVEY_START_YEAR,
                 survey_end=SURVEY_END_YEAR):
    """Run the full chain for one target. Returns a scorecard dict."""
    mid_epoch = 0.5 * (survey_start + survey_end)
    half_span = 0.5 * (survey_end - survey_start)
    dt = mid_epoch - GAIA_EPOCH

    ra, dec = propagate_position(
        target["ra"], target["dec"],
        target.get("pmra_mas_yr"), target.get("pmdec_mas_yr"), dt)
    radius = search_radius_arcsec(target, dt, half_span)

    card = {
        "source_id": target["source_id"],
        "distance_pc": target.get("distance_pc"),
        "g_mag": target.get("g_mag"),
        "search_ra": round(ra, 6),
        "search_dec": round(dec, 6),
        "search_radius_arcsec": round(radius, 2),
        "status": "no_match",
        "object_id": None,
        "match_sep_arcsec": None,
        "n_detections": 0,
        "bands": {},
        "n_flare_candidates": 0,
    }

    matches = broker.cone_search(ra, dec, radius)
    if not matches:
        return card

    best = min(matches,
               key=lambda m: angular_sep_arcsec(ra, dec, m["ra"], m["dec"]))
    card["object_id"] = best["object_id"]
    card["match_sep_arcsec"] = round(
        angular_sep_arcsec(ra, dec, best["ra"], best["dec"]), 2)

    detections = broker.get_lightcurve(best["object_id"])
    card["n_detections"] = len(detections)
    card["bands"] = analyze_by_band(detections)
    card["n_flare_candidates"] = sum(
        len(b["events"]) for b in card["bands"].values())
    card["status"] = "ok"
    return card


def run_watch(targets, output_dir=OUTPUT_DIR):
    """Watch every target; write per-target cards and a summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cards = []
    for t in targets:
        sid = t["source_id"]
        try:
            card = watch_target(t)
        except Exception:
            logger.exception("target %s failed", sid)
            card = {"source_id": sid, "status": "error"}
        cards.append(card)
        logger.info("%-22s %-10s det=%-5s flare_candidates=%s",
                    sid, card["status"], card.get("n_detections", "-"),
                    card.get("n_flare_candidates", "-"))
        with open(output_dir / f"target_{sid}.json", "w") as f:
            json.dump(card, f, indent=1)

    summary = {
        "n_targets": len(cards),
        "n_ok": sum(1 for c in cards if c["status"] == "ok"),
        "n_no_match": sum(1 for c in cards if c["status"] == "no_match"),
        "n_error": sum(1 for c in cards if c["status"] == "error"),
        "n_flare_candidates": sum(
            c.get("n_flare_candidates", 0) for c in cards),
        "targets_with_candidates": [
            c["source_id"] for c in cards
            if c.get("n_flare_candidates", 0) > 0],
    }
    with open(output_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    logger.info("summary: %s", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Broker alert watch over nearby M dwarf targets")
    parser.add_argument("--targets-file", default=str(TARGETS_FILE))
    parser.add_argument("--limit", type=int, default=None,
                        help="watch only the first N targets")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    targets = load_targets(args.targets_file)
    if args.limit:
        targets = targets[:args.limit]
    run_watch(targets)


if __name__ == "__main__":
    main()
