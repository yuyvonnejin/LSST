import json

import numpy as np
import pytest

from src import broker, watch
from src.watch import angular_sep_arcsec, run_watch, watch_target


def make_target(**overrides):
    t = {
        "source_id": "123456789",
        "ra": 280.0,
        "dec": -7.8,
        "parallax_mas": 50.0,
        "pmra_mas_yr": 100.0,
        "pmdec_mas_yr": -50.0,
        "g_mag": 16.5,
        "bp_rp": 3.0,
        "distance_pc": 20.0,
    }
    t.update(overrides)
    return t


def fake_lightcurve(n=40, flare_at=None):
    rng = np.random.default_rng(1)
    mags = 17.0 + rng.normal(0, 0.02, n)
    if flare_at is not None:
        mags[flare_at] -= 1.0
    return [{"mjd": 60000.0 + i, "band": "r", "mag": float(m),
             "magerr": 0.05} for i, m in enumerate(mags)]


class TestAngularSep:
    def test_zero(self):
        assert angular_sep_arcsec(100.0, -30.0, 100.0, -30.0) == 0.0

    def test_one_arcsec_dec(self):
        sep = angular_sep_arcsec(100.0, 0.0, 100.0, 1.0 / 3600.0)
        assert sep == pytest.approx(1.0)

    def test_ra_wrap(self):
        sep = angular_sep_arcsec(359.9999, 0.0, 0.0001, 0.0)
        assert sep == pytest.approx(0.72, abs=0.01)


class TestWatchTarget:
    def test_match_with_flare(self, monkeypatch):
        monkeypatch.setattr(broker, "cone_search", lambda ra, dec, r: [
            {"object_id": "ZTFfake1", "ra": ra, "dec": dec,
             "n_det": 40, "first_mjd": 60000.0, "last_mjd": 60040.0}])
        monkeypatch.setattr(broker, "get_lightcurve",
                            lambda oid: fake_lightcurve(flare_at=20))
        card = watch_target(make_target())
        assert card["status"] == "ok"
        assert card["object_id"] == "ZTFfake1"
        assert card["n_detections"] == 40
        assert card["n_flare_candidates"] == 1

    def test_no_match(self, monkeypatch):
        monkeypatch.setattr(broker, "cone_search", lambda ra, dec, r: [])
        card = watch_target(make_target())
        assert card["status"] == "no_match"
        assert card["object_id"] is None

    def test_picks_nearest_match(self, monkeypatch):
        def two_matches(ra, dec, r):
            return [
                {"object_id": "far", "ra": ra + 0.001, "dec": dec,
                 "n_det": 5, "first_mjd": None, "last_mjd": None},
                {"object_id": "near", "ra": ra, "dec": dec + 1e-5,
                 "n_det": 5, "first_mjd": None, "last_mjd": None},
            ]
        monkeypatch.setattr(broker, "cone_search", two_matches)
        monkeypatch.setattr(broker, "get_lightcurve",
                            lambda oid: fake_lightcurve())
        card = watch_target(make_target())
        assert card["object_id"] == "near"

    def test_position_propagated(self, monkeypatch):
        seen = {}

        def capture(ra, dec, r):
            seen["ra"], seen["dec"], seen["radius"] = ra, dec, r
            return []

        monkeypatch.setattr(broker, "cone_search", capture)
        t = make_target(pmra_mas_yr=0.0, pmdec_mas_yr=1000.0)
        watch_target(t)
        # mid-epoch 2022.35, dt from 2016.0 = 6.35 yr at 1 arcsec/yr
        assert seen["dec"] == pytest.approx(
            t["dec"] + 6.35 / 3600.0, abs=1e-5)
        # radius covers base 2" + drift over half span (4.15 yr)
        assert seen["radius"] == pytest.approx(2.0 + 4.15, abs=0.1)


class TestRunWatch:
    def test_failure_does_not_kill_run(self, monkeypatch, tmp_path):
        calls = {"n": 0}

        def flaky(ra, dec, r):
            calls["n"] += 1
            if calls["n"] == 1:
                raise RuntimeError("broker down")
            return []

        monkeypatch.setattr(broker, "cone_search", flaky)
        targets = [make_target(source_id="a"), make_target(source_id="b")]
        summary = run_watch(targets, output_dir=tmp_path)
        assert summary["n_targets"] == 2
        assert summary["n_error"] == 1
        assert summary["n_no_match"] == 1

    def test_reports_written(self, monkeypatch, tmp_path):
        monkeypatch.setattr(broker, "cone_search", lambda ra, dec, r: [
            {"object_id": "ZTFfake1", "ra": ra, "dec": dec,
             "n_det": 40, "first_mjd": None, "last_mjd": None}])
        monkeypatch.setattr(broker, "get_lightcurve",
                            lambda oid: fake_lightcurve(flare_at=10))
        summary = run_watch([make_target(source_id="xyz")],
                            output_dir=tmp_path)
        card = json.loads((tmp_path / "target_xyz.json").read_text())
        assert card["status"] == "ok"
        assert card["bands"]["r"]["status"] == "ok"
        loaded = json.loads((tmp_path / "summary.json").read_text())
        assert loaded == summary
        assert summary["targets_with_candidates"] == ["xyz"]
