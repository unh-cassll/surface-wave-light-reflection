"""Spectrum-level checks: shape, convergence, slope statistics, drag."""

import numpy as np
import pytest

from seapol import spectrum as sp


def test_omni_positive_and_finite():
    k = np.geomspace(1e-3, 2e4, 500)
    S = sp.elfouhaily_omni(k, 7.0)
    assert np.all(np.isfinite(S))
    assert np.all(S >= 0)


def test_omni_peak_near_kp():
    U10 = 7.0
    Omega_c = sp.inverse_wave_age(U10)
    kp = sp.GRAVITY * (Omega_c / U10) ** 2
    k = np.geomspace(1e-3, 1e3, 4000)
    S = sp.elfouhaily_omni(k, U10)
    k_peak = k[np.argmax(S)]
    assert 0.5 * kp < k_peak < 2.0 * kp


def test_variance_integral_converges_with_taper():
    """The low-k capillary taper keeps Hs physical; without it the
    elevation variance diverges toward low k."""
    k = np.geomspace(1e-5, 2e4, 40000)
    for U10 in [5.0, 7.0, 10.0]:
        S = sp.elfouhaily_omni(k, U10, lowk_capillary_taper=True)
        Hs = 4.0 * np.sqrt(np.trapezoid(S, k))
        Hs_PM = 0.21 * U10**2 / sp.GRAVITY
        assert 0.5 * Hs_PM < Hs < 1.6 * Hs_PM, (U10, Hs, Hs_PM)
    S_raw = sp.elfouhaily_omni(k, 5.0, lowk_capillary_taper=False)
    Hs_raw = 4.0 * np.sqrt(np.trapezoid(S_raw, k))
    assert Hs_raw > 10.0  # documents the artifact


def test_total_mss_vs_cox_munk():
    for U10 in [3.0, 5.0, 7.0, 10.0, 13.0]:
        mss = sp.total_mean_square_slope(U10, drag_model="elfouhaily")
        mss_cm = 0.003 + 5.12e-3 * U10
        assert 0.6 * mss_cm < mss < 1.4 * mss_cm, (U10, mss, mss_cm)


def test_cutoff_slope_variances():
    sa_lo, sc_lo = sp.cutoff_slope_variances(7.0, k_cutoff=10.0)
    sa_hi, sc_hi = sp.cutoff_slope_variances(7.0, k_cutoff=100.0)
    assert sa_lo > sa_hi > 0
    assert sc_lo > sc_hi > 0
    assert sa_lo > sc_lo  # along-wind exceeds cross-wind
    assert sp.cutoff_slope_variances(7.0, k_cutoff=3e4) == (0.0, 0.0)


def test_spread_normalization():
    k = np.array([0.1, 1.0, 100.0])
    phi = np.linspace(-np.pi, np.pi, 720, endpoint=False)
    dphi = 2 * np.pi / phi.size
    for one_sided in (False, True):
        D = sp.directional_spread(k[:, None], phi[None, :], 7.0,
                                  one_sided=one_sided)
        np.testing.assert_allclose(D.sum(axis=1) * dphi, 1.0, rtol=2e-3)


def test_drag_models():
    for model in ("logistic", "elfouhaily"):
        Cd, us = sp.drag_coefficient(np.array([3.0, 10.0, 25.0]), model)
        assert np.all(np.diff(us) > 0)
        assert np.all(Cd > 0)
    Cd, _ = sp.drag_coefficient(np.array([1.0, 60.0]), "logistic")
    assert Cd[0] > 0.8e-3 and Cd[1] < 3.2e-3
    with pytest.raises(ValueError):
        sp.drag_coefficient(5.0, "nope")
