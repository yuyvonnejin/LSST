"""Flare-candidate detection on sparse survey light curves.

Works in magnitudes. Baseline = median; scatter = MAD * 1.4826,
floored by the median reported magnitude error (sparse light curves
can produce an unrealistically small MAD). A candidate epoch is
brighter than baseline by more than max(k * scatter, min_amplitude).
Consecutive candidate epochs within group_gap_days merge into one
event.

Sparse-cadence honesty: with survey cadence a flare is usually one
epoch. Events are always "candidates"; single_epoch=True marks the
ones that can never be confirmed from this data alone.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

MIN_EPOCHS = 10
SIGMA_THRESHOLD = 3.0
MIN_AMPLITUDE_MAG = 0.1
GROUP_GAP_DAYS = 0.5


def quiescent_stats(mags, magerrs=None):
    """(baseline_mag, scatter_mag) robust to flares in the data."""
    mags = np.asarray(mags, dtype=float)
    baseline = float(np.median(mags))
    mad = float(np.median(np.abs(mags - baseline)))
    scatter = 1.4826 * mad
    if magerrs is not None:
        errs = np.asarray([e for e in magerrs if e is not None], dtype=float)
        if errs.size:
            scatter = max(scatter, float(np.median(errs)))
    if scatter <= 0:
        scatter = 0.01  # photometry is never perfect
    return baseline, scatter


def find_flares(detections,
                sigma_threshold=SIGMA_THRESHOLD,
                min_amplitude_mag=MIN_AMPLITUDE_MAG,
                group_gap_days=GROUP_GAP_DAYS,
                min_epochs=MIN_EPOCHS):
    """Find brightening events in one band's detections.

    detections: list of {mjd, mag, magerr} sorted by mjd.

    Returns dict:
      status: "ok" | "insufficient_data"
      n_epochs, baseline_mag, scatter_mag
      events: list of {mjd_start, mjd_end, peak_mjd, peak_amplitude_mag,
                       peak_flux_ratio, n_points, single_epoch}
    """
    n = len(detections)
    if n < min_epochs:
        return {"status": "insufficient_data", "n_epochs": n,
                "baseline_mag": None, "scatter_mag": None, "events": []}

    mjds = np.array([d["mjd"] for d in detections], dtype=float)
    mags = np.array([d["mag"] for d in detections], dtype=float)
    magerrs = [d.get("magerr") for d in detections]

    baseline, scatter = quiescent_stats(mags, magerrs)
    threshold = max(sigma_threshold * scatter, min_amplitude_mag)

    # brightening = magnitude below baseline
    amp = baseline - mags
    is_candidate = amp > threshold

    events = []
    i = 0
    while i < n:
        if not is_candidate[i]:
            i += 1
            continue
        j = i
        while (j + 1 < n and is_candidate[j + 1]
               and mjds[j + 1] - mjds[j] <= group_gap_days):
            j += 1
        seg = slice(i, j + 1)
        peak_idx = i + int(np.argmax(amp[seg]))
        peak_amp = float(amp[peak_idx])
        events.append({
            "mjd_start": float(mjds[i]),
            "mjd_end": float(mjds[j]),
            "peak_mjd": float(mjds[peak_idx]),
            "peak_amplitude_mag": round(peak_amp, 3),
            "peak_flux_ratio": round(10.0 ** (0.4 * peak_amp), 3),
            "n_points": j - i + 1,
            "single_epoch": j == i,
        })
        i = j + 1

    return {
        "status": "ok",
        "n_epochs": n,
        "baseline_mag": round(baseline, 3),
        "scatter_mag": round(scatter, 4),
        "events": events,
    }


def analyze_by_band(detections, **kwargs):
    """Split detections by band and run find_flares per band."""
    bands = {}
    for d in detections:
        bands.setdefault(d.get("band", "?"), []).append(d)
    return {band: find_flares(dets, **kwargs)
            for band, dets in sorted(bands.items())}
