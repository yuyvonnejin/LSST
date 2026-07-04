"""Alert-broker REST client (ALeRCE), survey-parametrized.

Two backends behind one interface:
- "ztf": legacy ALeRCE ZTF v1 API. 7+ years of alerts; magnitudes.
- "lsst": ALeRCE multisurvey API serving live Rubin/LSST alerts
  (verified 2026-07-04). Detections carry fluxes in nJy, converted
  here to AB magnitudes so downstream code sees one schema:
  {mjd, band, mag, magerr}.

LSST notes:
- scienceFlux is the star's total brightness on the science image;
  psfFlux is the difference-image flux (star minus template, can be
  negative). Flare detection works on total light, so magnitudes come
  from scienceFlux.
- Rows with isNegative, solar-system associations (ssObjectId != 0)
  or non-positive flux are skipped.
- Cone-search pages can repeat an object (one row per classifier
  ranking); results are deduped by oid.

Alerts and broker access are fully public (no LSST data rights
required) -- see docs/design_doc.md section 1.
"""

import logging
import math

import requests

logger = logging.getLogger(__name__)

ZTF_BASE_URL = "https://api.alerce.online/ztf/v1"
LSST_BASE_URL = "https://api-lsst.alerce.online"

ZTF_BAND_NAMES = {1: "g", 2: "r", 3: "i"}
# fallback if an LSST payload omits band_map
LSST_BAND_NAMES = {1: "g", 2: "r", 3: "i", 4: "z", 5: "y", 6: "u"}

# AB magnitude zero point for flux in nJy: m = -2.5 log10(f) + 31.4
AB_ZP_NJY = 31.4

TIMEOUT = 60

SURVEYS = ("ztf", "lsst")


def flux_njy_to_mag(flux, flux_err=None):
    """(AB mag, mag error) from flux in nJy. None if flux <= 0."""
    if flux is None or flux <= 0:
        return None, None
    mag = -2.5 * math.log10(flux) + AB_ZP_NJY
    magerr = None
    if flux_err is not None and flux_err > 0:
        magerr = 1.0857 * flux_err / flux
    return mag, magerr


def parse_cone_search_ztf(payload):
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


def parse_cone_search_lsst(payload):
    """Same output shape as ZTF; dedupes repeated oids (one row per
    classifier ranking) and keys on n_det instead of ndet."""
    objects = []
    seen = set()
    for item in payload.get("items") or []:
        try:
            oid = item["oid"]
            if oid in seen:
                continue
            objects.append({
                "object_id": oid,
                "ra": float(item["meanra"]),
                "dec": float(item["meandec"]),
                "n_det": int(item.get("n_det") or 0),
                "first_mjd": item.get("firstmjd"),
                "last_mjd": item.get("lastmjd"),
            })
            seen.add(oid)
        except (KeyError, TypeError, ValueError):
            logger.warning("skipping malformed cone-search item: %r", item)
    return objects


def parse_lightcurve_ztf(payload):
    """Detections as {mjd, band, mag, magerr}. Prefers corrected
    magnitudes (magpsf_corr) when present and sane."""
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
            "band": ZTF_BAND_NAMES.get(det.get("fid"), str(det.get("fid"))),
            "mag": float(mag),
            "magerr": float(magerr) if magerr is not None else None,
        })
    detections.sort(key=lambda d: d["mjd"])
    return detections


def parse_lightcurve_lsst(payload):
    """Detections as {mjd, band, mag, magerr} from LSST diaSource rows.

    Total brightness from scienceFlux (nJy) -> AB magnitude. Skips
    negative detections, solar-system objects, and rows without a
    positive flux.
    """
    detections = []
    for det in payload.get("detections") or []:
        if det.get("mjd") is None:
            continue
        if det.get("isNegative"):
            continue
        if det.get("ssObjectId"):
            continue
        mag, magerr = flux_njy_to_mag(det.get("scienceFlux"),
                                      det.get("scienceFluxErr"))
        if mag is None:
            continue
        band_map = det.get("band_map") or {}
        band = det.get("band")
        band_name = (band_map.get(str(band)) or band_map.get(band)
                     or LSST_BAND_NAMES.get(band) or str(band))
        detections.append({
            "mjd": float(det["mjd"]),
            "band": band_name,
            "mag": round(mag, 4),
            "magerr": round(magerr, 4) if magerr is not None else None,
        })
    detections.sort(key=lambda d: d["mjd"])
    return detections


def _get(url, params):
    resp = requests.get(url, params=params, timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.json()


def cone_search(ra_deg, dec_deg, radius_arcsec, survey="ztf",
                page_size=20):
    """Objects within radius_arcsec of (ra, dec). Public API, no auth."""
    if survey == "ztf":
        payload = _get(f"{ZTF_BASE_URL}/objects/", {
            "ra": ra_deg, "dec": dec_deg, "radius": radius_arcsec,
            "page_size": page_size,
        })
        return parse_cone_search_ztf(payload)
    if survey == "lsst":
        payload = _get(f"{LSST_BASE_URL}/object_api/list_objects", {
            "survey": "lsst", "ra": ra_deg, "dec": dec_deg,
            "radius": radius_arcsec, "page_size": page_size,
        })
        return parse_cone_search_lsst(payload)
    raise ValueError(f"survey must be one of {SURVEYS}")


def get_lightcurve(object_id, survey="ztf"):
    """All usable detections for one broker object, sorted by mjd."""
    if survey == "ztf":
        payload = _get(f"{ZTF_BASE_URL}/objects/{object_id}/lightcurve", {})
        return parse_lightcurve_ztf(payload)
    if survey == "lsst":
        payload = _get(f"{LSST_BASE_URL}/lightcurve_api/lightcurve", {
            "survey_id": "lsst", "oid": object_id,
        })
        return parse_lightcurve_lsst(payload)
    raise ValueError(f"survey must be one of {SURVEYS}")
