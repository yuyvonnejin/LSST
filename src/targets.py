"""Target selection: nearby M dwarfs in the LSST magnitude window.

Queries the public Gaia DR3 TAP service (plain HTTP, no astroquery)
and writes data/targets.json. Selection cuts are pure functions so
they can be tested offline.

Gaia positions are epoch 2016.0. Nearby M dwarfs move fast (up to
~10 arcsec/yr), so positions must be propagated to the observation
epoch before any cone search. See propagate_position().
"""

import json
import logging
import math
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_EPOCH = 2016.0

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
TARGETS_FILE = DATA_DIR / "targets.json"

# Selection defaults (rationale in docs/design_doc.md section 3)
MIN_PARALLAX_MAS = 20.0      # d < 50 pc
MIN_BP_RP = 2.0              # approx M0V and later
MIN_ABS_G = 8.0              # dwarfs only, rejects giants
MIN_G_MAG = 15.5             # proxy for r >~ 16.5: avoid LSST saturation
MAX_G_MAG = 20.5             # keep decent alert-photometry S/N
MAX_DEC_DEG = 32.0           # LSST footprint (southern + low northern)
MIN_PLX_OVER_ERR = 10.0
MAX_RUWE = 1.4


def abs_g_mag(phot_g_mean_mag, parallax_mas):
    """Absolute G magnitude from apparent G and parallax in mas."""
    return phot_g_mean_mag + 5.0 * math.log10(parallax_mas / 100.0)


def passes_cuts(row,
                min_parallax_mas=MIN_PARALLAX_MAS,
                min_bp_rp=MIN_BP_RP,
                min_abs_g=MIN_ABS_G,
                min_g_mag=MIN_G_MAG,
                max_g_mag=MAX_G_MAG,
                max_dec_deg=MAX_DEC_DEG,
                min_plx_over_err=MIN_PLX_OVER_ERR,
                max_ruwe=MAX_RUWE):
    """Apply the M dwarf selection cuts to one Gaia row (dict).

    Returns False on any missing required field.
    """
    required = ("parallax", "parallax_error", "phot_g_mean_mag",
                "bp_rp", "dec", "ruwe")
    if any(row.get(k) is None for k in required):
        return False
    if row["parallax"] < min_parallax_mas:
        return False
    if row["parallax"] / row["parallax_error"] < min_plx_over_err:
        return False
    if row["bp_rp"] < min_bp_rp:
        return False
    if not (min_g_mag <= row["phot_g_mean_mag"] <= max_g_mag):
        return False
    if row["dec"] > max_dec_deg:
        return False
    if row["ruwe"] > max_ruwe:
        return False
    if abs_g_mag(row["phot_g_mean_mag"], row["parallax"]) < min_abs_g:
        return False
    return True


def propagate_position(ra_deg, dec_deg, pmra_mas_yr, pmdec_mas_yr,
                       dt_years):
    """Propagate a position by proper motion over dt_years.

    pmra is the Gaia convention: pmra = mu_alpha* = d(alpha)/dt * cos(dec),
    already projected on the sky. Linear propagation, fine for the
    <= arcminute drifts involved here.

    Returns (ra_deg, dec_deg) at the new epoch.
    """
    if pmra_mas_yr is None or pmdec_mas_yr is None:
        return ra_deg, dec_deg
    cos_dec = math.cos(math.radians(dec_deg))
    if abs(cos_dec) < 1e-9:
        cos_dec = 1e-9
    dra_deg = (pmra_mas_yr * dt_years / 3.6e6) / cos_dec
    ddec_deg = pmdec_mas_yr * dt_years / 3.6e6
    return (ra_deg + dra_deg) % 360.0, dec_deg + ddec_deg


def pm_drift_arcsec(pmra_mas_yr, pmdec_mas_yr, dt_years):
    """Total sky drift in arcsec over dt_years."""
    if pmra_mas_yr is None or pmdec_mas_yr is None:
        return 0.0
    pm_total = math.hypot(pmra_mas_yr, pmdec_mas_yr)
    return pm_total * abs(dt_years) / 1000.0


def build_adql(limit=2000):
    """ADQL for the M dwarf sample. Cuts that need derived quantities
    (absolute magnitude) are applied client-side in passes_cuts."""
    return f"""
        SELECT TOP {limit}
            source_id, ra, dec, parallax, parallax_error,
            pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe
        FROM gaiadr3.gaia_source
        WHERE parallax > {MIN_PARALLAX_MAS}
          AND parallax_over_error > {MIN_PLX_OVER_ERR}
          AND bp_rp > {MIN_BP_RP}
          AND phot_g_mean_mag BETWEEN {MIN_G_MAG} AND {MAX_G_MAG}
          AND dec < {MAX_DEC_DEG}
          AND ruwe < {MAX_RUWE}
        ORDER BY parallax DESC
    """


def query_gaia(adql, timeout=120):
    """Run a synchronous TAP query, return list of row dicts."""
    resp = requests.post(
        GAIA_TAP_SYNC,
        data={
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "json",
            "QUERY": adql,
        },
        timeout=timeout,
    )
    resp.raise_for_status()
    payload = resp.json()
    cols = [c["name"] for c in payload["metadata"]]
    return [dict(zip(cols, row)) for row in payload["data"]]


def build_target_catalog(limit=2000, out_file=TARGETS_FILE):
    """Query Gaia, apply cuts, write the target catalog. Returns targets."""
    rows = query_gaia(build_adql(limit=limit))
    logger.info("Gaia returned %d rows", len(rows))
    targets = [
        {
            "source_id": str(r["source_id"]),
            "ra": r["ra"],
            "dec": r["dec"],
            "parallax_mas": r["parallax"],
            "pmra_mas_yr": r.get("pmra"),
            "pmdec_mas_yr": r.get("pmdec"),
            "g_mag": r["phot_g_mean_mag"],
            "bp_rp": r["bp_rp"],
            "abs_g": round(abs_g_mag(r["phot_g_mean_mag"], r["parallax"]), 3),
            "distance_pc": round(1000.0 / r["parallax"], 2),
        }
        for r in rows
        if passes_cuts(r)
    ]
    logger.info("%d targets pass cuts", len(targets))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, "w") as f:
        json.dump({"gaia_epoch": GAIA_EPOCH, "targets": targets}, f, indent=1)
    logger.info("wrote %s", out_file)
    return targets


def load_targets(path=TARGETS_FILE):
    with open(path) as f:
        return json.load(f)["targets"]


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    build_target_catalog()
