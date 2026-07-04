"""Live smoke tests. Run manually: pytest -m network"""

import pytest

from src.broker import cone_search, get_lightcurve
from src.targets import build_adql, passes_cuts, query_gaia

pytestmark = pytest.mark.network


def test_alerce_cone_search_live():
    # position verified to return a ZTF object (2026-07-03)
    objs = cone_search(280.694, -7.783, 30.0)
    assert objs
    assert objs[0]["object_id"].startswith("ZTF")


def test_alerce_lightcurve_live():
    objs = cone_search(280.694, -7.783, 30.0)
    dets = get_lightcurve(objs[0]["object_id"])
    assert isinstance(dets, list)
    for d in dets:
        assert d["mjd"] > 0
        assert 5 < d["mag"] < 30


def test_gaia_tap_live():
    rows = query_gaia(build_adql(limit=20))
    assert rows
    assert any(passes_cuts(r) for r in rows)


def test_alerce_lsst_cone_search_live():
    # position of a real Rubin alert captured 2026-07-04
    objs = cone_search(305.5822, -18.7909, 30.0, survey="lsst")
    assert objs
    assert isinstance(objs[0]["object_id"], int)


def test_alerce_lsst_lightcurve_live():
    objs = cone_search(305.5822, -18.7909, 30.0, survey="lsst")
    dets = get_lightcurve(objs[0]["object_id"], survey="lsst")
    assert isinstance(dets, list)
    for d in dets:
        assert d["mjd"] > 61000
        assert 5 < d["mag"] < 30
        assert d["band"] in ("u", "g", "r", "i", "z", "y")
