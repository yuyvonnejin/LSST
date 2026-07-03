import numpy as np
import pytest

from src.flares import analyze_by_band, find_flares, quiescent_stats


def make_lc(mags, start_mjd=60000.0, cadence_days=1.0, magerr=0.05,
            band="r"):
    return [{"mjd": start_mjd + i * cadence_days, "mag": float(m),
             "magerr": magerr, "band": band}
            for i, m in enumerate(mags)]


def quiet_mags(n=50, level=17.0, noise=0.02, seed=42):
    rng = np.random.default_rng(seed)
    return level + rng.normal(0, noise, n)


class TestQuiescentStats:
    def test_baseline_is_median(self):
        baseline, _ = quiescent_stats([17.0, 17.1, 16.9, 17.0, 17.0])
        assert baseline == pytest.approx(17.0)

    def test_scatter_floored_by_magerr(self):
        # identical mags -> MAD 0, but reported errors are 0.05
        _, scatter = quiescent_stats([17.0] * 20, magerrs=[0.05] * 20)
        assert scatter == pytest.approx(0.05)

    def test_scatter_never_zero(self):
        _, scatter = quiescent_stats([17.0] * 20)
        assert scatter > 0

    def test_baseline_robust_to_flare(self):
        mags = [17.0] * 30 + [15.0]  # one big flare
        baseline, _ = quiescent_stats(mags)
        assert baseline == pytest.approx(17.0)


class TestFindFlares:
    def test_quiet_star_no_events(self):
        result = find_flares(make_lc(quiet_mags()))
        assert result["status"] == "ok"
        assert result["events"] == []

    def test_insufficient_data(self):
        result = find_flares(make_lc([17.0] * 5))
        assert result["status"] == "insufficient_data"
        assert result["events"] == []

    def test_single_epoch_flare_recovered(self):
        mags = quiet_mags()
        mags[25] -= 1.0  # 1 mag brightening
        result = find_flares(make_lc(mags))
        assert len(result["events"]) == 1
        ev = result["events"][0]
        assert ev["single_epoch"] is True
        assert ev["peak_amplitude_mag"] == pytest.approx(1.0, abs=0.1)
        assert ev["peak_flux_ratio"] == pytest.approx(2.51, abs=0.3)

    def test_multi_epoch_flare_grouped(self):
        mags = quiet_mags()
        mags[20] -= 0.8
        mags[21] -= 0.5  # same event, adjacent epochs...
        lc = make_lc(mags, cadence_days=0.3)  # ...within group gap
        result = find_flares(lc)
        assert len(result["events"]) == 1
        ev = result["events"][0]
        assert ev["n_points"] == 2
        assert ev["single_epoch"] is False
        assert ev["peak_amplitude_mag"] == pytest.approx(0.8, abs=0.1)

    def test_separated_flares_not_grouped(self):
        mags = quiet_mags()
        mags[10] -= 1.0
        mags[40] -= 1.0  # 30 days apart
        result = find_flares(make_lc(mags))
        assert len(result["events"]) == 2

    def test_dimming_is_not_a_flare(self):
        mags = quiet_mags()
        mags[25] += 1.5  # eclipse/dimming, not a flare
        result = find_flares(make_lc(mags))
        assert result["events"] == []

    def test_small_bump_below_min_amplitude(self):
        # tiny noise -> 3 sigma could be < 0.1 mag; min_amplitude blocks
        mags = quiet_mags(noise=0.005)
        mags[25] -= 0.06
        result = find_flares(make_lc(mags, magerr=0.005))
        assert result["events"] == []


class TestAnalyzeByBand:
    def test_bands_analyzed_separately(self):
        g = make_lc(quiet_mags(level=18.0), band="g")
        r_mags = quiet_mags(level=17.0)
        r_mags[25] -= 1.0
        r = make_lc(r_mags, band="r")
        result = analyze_by_band(g + r)
        assert set(result) == {"g", "r"}
        assert result["g"]["events"] == []
        assert len(result["r"]["events"]) == 1

    def test_sparse_band_insufficient(self):
        result = analyze_by_band(make_lc([17.0] * 3, band="i"))
        assert result["i"]["status"] == "insufficient_data"
