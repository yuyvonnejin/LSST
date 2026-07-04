import pytest

from src.broker import (flux_njy_to_mag, parse_cone_search_lsst,
                        parse_cone_search_ztf, parse_lightcurve_lsst,
                        parse_lightcurve_ztf)

# Trimmed real ALeRCE ZTF v1 response shapes (captured 2026-07-03)
CONE_PAYLOAD = {
    "total": None, "page": None, "has_next": False,
    "items": [
        {
            "oid": "ZTF25abaaish",
            "ndet": 2,
            "meanra": 280.697399900019,
            "meandec": -7.785294249999999,
            "firstmjd": 60862.413888900075,
            "lastmjd": 61162.41991900001,
        },
        {  # malformed: missing meandec
            "oid": "ZTFbroken",
            "meanra": 280.0,
        },
    ],
}

LC_PAYLOAD = {
    "detections": [
        {"mjd": 60002.5, "fid": 2, "magpsf": 17.20, "sigmapsf": 0.05,
         "magpsf_corr": 17.10, "sigmapsf_corr": 0.04},
        {"mjd": 60001.5, "fid": 1, "magpsf": 18.00, "sigmapsf": 0.08,
         "magpsf_corr": None, "sigmapsf_corr": None},
        {"mjd": 60003.5, "fid": 2, "magpsf": 17.25, "sigmapsf": 0.05,
         "magpsf_corr": 100.0, "sigmapsf_corr": 99.0},  # sentinel junk
        {"mjd": None, "fid": 2, "magpsf": 17.0, "sigmapsf": 0.05},
        {"mjd": 60004.5, "fid": 2, "magpsf": None, "sigmapsf": None},
    ],
    "non_detections": [],
}


class TestParseConeSearch:
    def test_parses_good_item(self):
        objs = parse_cone_search_ztf(CONE_PAYLOAD)
        assert len(objs) == 1
        obj = objs[0]
        assert obj["object_id"] == "ZTF25abaaish"
        assert abs(obj["ra"] - 280.6974) < 1e-3
        assert obj["n_det"] == 2

    def test_empty_payload(self):
        assert parse_cone_search_ztf({"items": []}) == []
        assert parse_cone_search_ztf({}) == []
        assert parse_cone_search_ztf({"items": None}) == []


class TestParseLightcurve:
    def test_prefers_corrected_mag(self):
        dets = parse_lightcurve_ztf(LC_PAYLOAD)
        r_dets = [d for d in dets if d["band"] == "r"]
        assert r_dets[0]["mag"] == 17.10  # corrected, not 17.20

    def test_falls_back_when_corrected_missing_or_junk(self):
        dets = parse_lightcurve_ztf(LC_PAYLOAD)
        g = [d for d in dets if d["band"] == "g"][0]
        assert g["mag"] == 18.00  # corr was None
        junk = [d for d in dets if d["mjd"] == 60003.5][0]
        assert junk["mag"] == 17.25  # corr was sentinel 100.0

    def test_skips_unusable_rows(self):
        dets = parse_lightcurve_ztf(LC_PAYLOAD)
        # 5 raw rows, 2 unusable (no mjd, no mag)
        assert len(dets) == 3

    def test_sorted_by_mjd(self):
        dets = parse_lightcurve_ztf(LC_PAYLOAD)
        mjds = [d["mjd"] for d in dets]
        assert mjds == sorted(mjds)

    def test_empty(self):
        assert parse_lightcurve_ztf({}) == []
        assert parse_lightcurve_ztf({"detections": None}) == []


# Trimmed real ALeRCE multisurvey LSST response shapes (captured 2026-07-04)
LSST_CONE_PAYLOAD = {
    "total": 4, "current_page": 1,
    "items": [
        {"oid": 170591527609303944, "meanra": 305.5822327501884,
         "meandec": -18.7909207179724, "firstmjd": 61218.348303451654,
         "lastmjd": 61218.348303451654, "n_det": 1,
         "class_name": "SN", "ranking": 1},
        # same object again, different classifier ranking
        {"oid": 170591527609303944, "meanra": 305.5822327501884,
         "meandec": -18.7909207179724, "firstmjd": 61218.348303451654,
         "lastmjd": 61218.348303451654, "n_det": 1,
         "class_name": "VS", "ranking": 2},
        {"oid": 170591547426865168, "meanra": 347.0995208808436,
         "meandec": -5.778599497547491, "firstmjd": 61218.42718196748,
         "lastmjd": 61221.37129479152, "n_det": 2},
    ],
}

BAND_MAP = {"6": "u", "1": "g", "2": "r", "3": "i", "4": "z", "5": "y"}

LSST_LC_PAYLOAD = {
    "detections": [
        {"band_map": BAND_MAP, "band": 4, "mjd": 61218.348303451654,
         "oid": 170591527609303944, "ssObjectId": 0, "isNegative": False,
         "psfFlux": 32207.861, "psfFluxErr": 764.197,
         "scienceFlux": 268571.7, "scienceFluxErr": 723.83276},
        {"band_map": BAND_MAP, "band": 2, "mjd": 61219.5,
         "oid": 170591527609303944, "ssObjectId": 0, "isNegative": False,
         "scienceFlux": 100000.0, "scienceFluxErr": 1000.0},
        # negative difference detection: skipped
        {"band_map": BAND_MAP, "band": 2, "mjd": 61220.5,
         "ssObjectId": 0, "isNegative": True,
         "scienceFlux": 90000.0, "scienceFluxErr": 900.0},
        # solar system object: skipped
        {"band_map": BAND_MAP, "band": 2, "mjd": 61220.6,
         "ssObjectId": 998877, "isNegative": False,
         "scienceFlux": 90000.0, "scienceFluxErr": 900.0},
        # non-positive flux: no magnitude exists, skipped
        {"band_map": BAND_MAP, "band": 2, "mjd": 61220.7,
         "ssObjectId": 0, "isNegative": False,
         "scienceFlux": -50.0, "scienceFluxErr": 900.0},
    ],
}


class TestFluxToMag:
    def test_known_value(self):
        # 3631 Jy = AB 0; 1 nJy -> 31.4
        mag, _ = flux_njy_to_mag(1.0)
        assert mag == pytest.approx(31.4, abs=0.001)

    def test_science_flux_from_capture(self):
        # hand-computed: -2.5*log10(268571.7) + 31.4 = 17.827
        mag, magerr = flux_njy_to_mag(268571.7, 723.83276)
        assert mag == pytest.approx(17.827, abs=0.002)
        assert magerr == pytest.approx(1.0857 * 723.83276 / 268571.7,
                                       rel=1e-6)

    def test_nonpositive_flux(self):
        assert flux_njy_to_mag(0.0) == (None, None)
        assert flux_njy_to_mag(-10.0) == (None, None)
        assert flux_njy_to_mag(None) == (None, None)


class TestParseConeSearchLsst:
    def test_dedupes_by_oid(self):
        objs = parse_cone_search_lsst(LSST_CONE_PAYLOAD)
        assert len(objs) == 2
        assert objs[0]["object_id"] == 170591527609303944
        assert objs[0]["n_det"] == 1
        assert objs[1]["object_id"] == 170591547426865168

    def test_empty(self):
        assert parse_cone_search_lsst({}) == []
        assert parse_cone_search_lsst({"items": None}) == []


class TestParseLightcurveLsst:
    def test_converts_science_flux_to_ab_mag(self):
        dets = parse_lightcurve_lsst(LSST_LC_PAYLOAD)
        z = [d for d in dets if d["band"] == "z"][0]
        assert z["mag"] == pytest.approx(17.827, abs=0.002)
        assert z["magerr"] == pytest.approx(0.0029, abs=0.0005)

    def test_band_names_from_payload_map(self):
        dets = parse_lightcurve_lsst(LSST_LC_PAYLOAD)
        assert {d["band"] for d in dets} == {"z", "r"}

    def test_skips_negative_ss_and_bad_flux(self):
        dets = parse_lightcurve_lsst(LSST_LC_PAYLOAD)
        # 5 raw rows: isNegative, ssObjectId, flux<0 skipped
        assert len(dets) == 2

    def test_sorted_and_shape(self):
        dets = parse_lightcurve_lsst(LSST_LC_PAYLOAD)
        mjds = [d["mjd"] for d in dets]
        assert mjds == sorted(mjds)
        assert set(dets[0]) == {"mjd", "band", "mag", "magerr"}

    def test_empty(self):
        assert parse_lightcurve_lsst({}) == []
        assert parse_lightcurve_lsst({"detections": None}) == []
