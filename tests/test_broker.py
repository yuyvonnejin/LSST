from src.broker import parse_cone_search, parse_lightcurve

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
        objs = parse_cone_search(CONE_PAYLOAD)
        assert len(objs) == 1
        obj = objs[0]
        assert obj["object_id"] == "ZTF25abaaish"
        assert abs(obj["ra"] - 280.6974) < 1e-3
        assert obj["n_det"] == 2

    def test_empty_payload(self):
        assert parse_cone_search({"items": []}) == []
        assert parse_cone_search({}) == []
        assert parse_cone_search({"items": None}) == []


class TestParseLightcurve:
    def test_prefers_corrected_mag(self):
        dets = parse_lightcurve(LC_PAYLOAD)
        r_dets = [d for d in dets if d["band"] == "r"]
        assert r_dets[0]["mag"] == 17.10  # corrected, not 17.20

    def test_falls_back_when_corrected_missing_or_junk(self):
        dets = parse_lightcurve(LC_PAYLOAD)
        g = [d for d in dets if d["band"] == "g"][0]
        assert g["mag"] == 18.00  # corr was None
        junk = [d for d in dets if d["mjd"] == 60003.5][0]
        assert junk["mag"] == 17.25  # corr was sentinel 100.0

    def test_skips_unusable_rows(self):
        dets = parse_lightcurve(LC_PAYLOAD)
        # 5 raw rows, 2 unusable (no mjd, no mag)
        assert len(dets) == 3

    def test_sorted_by_mjd(self):
        dets = parse_lightcurve(LC_PAYLOAD)
        mjds = [d["mjd"] for d in dets]
        assert mjds == sorted(mjds)

    def test_empty(self):
        assert parse_lightcurve({}) == []
        assert parse_lightcurve({"detections": None}) == []
