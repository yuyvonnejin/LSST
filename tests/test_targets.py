import math

import pytest

from src.targets import (abs_g_mag, build_adql, passes_cuts,
                         pm_drift_arcsec, propagate_position)


def good_row(**overrides):
    row = {
        "parallax": 50.0,          # 20 pc
        "parallax_error": 0.05,
        "phot_g_mean_mag": 16.5,
        "bp_rp": 3.0,
        "dec": -30.0,
        "ruwe": 1.0,
    }
    row.update(overrides)
    return row


class TestAbsGMag:
    def test_100pc_star(self):
        # at 100 pc (10 mas): M_G = G - 5
        assert abs_g_mag(15.0, 10.0) == pytest.approx(10.0)

    def test_10pc_star_equals_apparent(self):
        # at 10 pc (100 mas): M_G = G
        assert abs_g_mag(15.0, 100.0) == pytest.approx(15.0)

    def test_20pc(self):
        # d = 20 pc: M_G = G - 5*log10(20/10) = G - 1.505
        assert abs_g_mag(16.5, 50.0) == pytest.approx(16.5 - 5 * math.log10(2))


class TestPassesCuts:
    def test_good_m_dwarf_passes(self):
        assert passes_cuts(good_row())

    def test_too_far(self):
        assert not passes_cuts(good_row(parallax=10.0))

    def test_too_blue(self):
        assert not passes_cuts(good_row(bp_rp=1.0))

    def test_too_bright_saturates(self):
        assert not passes_cuts(good_row(phot_g_mean_mag=12.0))

    def test_too_faint(self):
        assert not passes_cuts(good_row(phot_g_mean_mag=21.5))

    def test_north_of_footprint(self):
        assert not passes_cuts(good_row(dec=45.0))

    def test_bad_astrometry_ruwe(self):
        assert not passes_cuts(good_row(ruwe=2.0))

    def test_low_parallax_snr(self):
        assert not passes_cuts(good_row(parallax_error=10.0))

    def test_giant_rejected(self):
        # bright absolute magnitude: G=15.5 at 50 pc -> M_G = 12,
        # push distance out via small parallax is blocked by distance
        # cut, so fake a giant with abs G < 8: G=16 at 20pc -> M_G=14.5
        # need M_G < 8: G=15.5, parallax=20 -> M_G=12. Can't make a
        # giant inside 50pc in the allowed G window; verify the cut
        # directly instead.
        assert not passes_cuts(good_row(), min_abs_g=15.0)

    def test_missing_field_rejected(self):
        assert not passes_cuts(good_row(bp_rp=None))


class TestProperMotion:
    def test_no_pm_returns_input(self):
        assert propagate_position(100.0, -20.0, None, None, 10.0) == \
            (100.0, -20.0)

    def test_dec_only(self):
        # 1000 mas/yr for 10 yr = 10 arcsec in dec
        ra, dec = propagate_position(100.0, 0.0, 0.0, 1000.0, 10.0)
        assert ra == pytest.approx(100.0)
        assert dec == pytest.approx(10.0 / 3600.0)

    def test_ra_cos_dec_projection(self):
        # pmra is already sky-projected; at dec=60 the RA coordinate
        # change is doubled (cos 60 = 0.5)
        ra, dec = propagate_position(100.0, 60.0, 1000.0, 0.0, 10.0)
        expected_dra = (10.0 / 3600.0) / math.cos(math.radians(60.0))
        assert ra == pytest.approx(100.0 + expected_dra)
        assert dec == pytest.approx(60.0)

    def test_barnard_scale_drift(self):
        # ~10.3 arcsec/yr over 8.15 yr (2016 -> 2024.15) ~ 84 arcsec
        drift = pm_drift_arcsec(-802.8, 10362.5, 8.15)
        assert drift == pytest.approx(84.7, abs=0.5)

    def test_ra_wraps(self):
        ra, _ = propagate_position(359.9999, 0.0, 3600.0 * 1000.0, 0.0, 1.0)
        assert 0.0 <= ra < 360.0


class TestAdql:
    def test_contains_cuts_and_limit(self):
        q = build_adql(limit=123)
        assert "TOP 123" in q
        assert "gaiadr3.gaia_source" in q
        assert "parallax > 20.0" in q
        assert "pmra" in q and "pmdec" in q
