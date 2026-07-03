"""Alert-broker REST client.

Backend: ALeRCE (https://alerce.science). The ZTF v1 API is the
working service today; ALeRCE also processes Rubin/LSST alerts and
the same client can be pointed at that service via BASE_URL once its
public endpoint is published. Parsing is separated from HTTP so tests
run offline on captured JSON.

Alerts and broker access are fully public (no LSST data rights
required) -- see docs/design_doc.md section 1.
"""

import logging

import requests

logger = logging.getLogger(__name__)

# Swap this to the ALeRCE LSST service (or a Fink adapter) when available.
BASE_URL = "https://api.alerce.online/ztf/v1"

# ZTF filter ids -> band names. Rubin will use different ids; extend then.
BAND_NAMES = {1: "g", 2: "r", 3: "i"}

TIMEOUT = 60


def parse_cone_search(payload):
    """Extract matched objects from an ALeRCE objects/ response.

    Returns list of dicts: object_id, ra, dec, n_det, first_mjd, last_mjd.
    Malformed items are skipped, not fatal.
    """
    objects = []
    for item in payload.get("items") or []:
        try:
            objects.append({
                "object_id": item["oid"],
                "ra": float(item["meanra"]),
                "dec": float(item["meandec"]),
                "n_det": int(item.get("ndet") or 0),
                "first_mjd": item.get("firstmjd"),
                "last_mjd": item.get("lastmjd"),
            })
        except (KeyError, TypeError, ValueError):
            logger.warning("skipping malformed cone-search item: %r", item)
    return objects


def parse_lightcurve(payload):
    """Extract detections from an ALeRCE lightcurve response.

    Returns list of dicts: mjd, band, mag, magerr. Prefers corrected
    magnitudes (magpsf_corr) when present and finite, falls back to
    magpsf. Rows without a usable magnitude are skipped.
    """
    detections = []
    for det in payload.get("detections") or []:
        mag = det.get("magpsf_corr")
        magerr = det.get("sigmapsf_corr")
        if mag is None or not (0 < mag < 30):
            mag = det.get("magpsf")
            magerr = det.get("sigmapsf")
        if mag is None or det.get("mjd") is None:
            continue
        detections.append({
            "mjd": float(det["mjd"]),
            "band": BAND_NAMES.get(det.get("fid"), str(det.get("fid"))),
            "mag": float(mag),
            "magerr": float(magerr) if magerr is not None else None,
        })
    detections.sort(key=lambda d: d["mjd"])
    return detections


def cone_search(ra_deg, dec_deg, radius_arcsec, base_url=BASE_URL,
                page_size=20):
    """Objects within radius_arcsec of (ra, dec). Public API, no auth."""
    resp = requests.get(
        f"{base_url}/objects/",
        params={
            "ra": ra_deg,
            "dec": dec_deg,
            "radius": radius_arcsec,
            "page_size": page_size,
        },
        timeout=TIMEOUT,
    )
    resp.raise_for_status()
    return parse_cone_search(resp.json())


def get_lightcurve(object_id, base_url=BASE_URL):
    """All detections for one broker object, sorted by mjd."""
    resp = requests.get(
        f"{base_url}/objects/{object_id}/lightcurve",
        timeout=TIMEOUT,
    )
    resp.raise_for_status()
    return parse_lightcurve(resp.json())
