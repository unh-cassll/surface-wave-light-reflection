"""FM98 solver checks: convergence, harmonic structure, gauge fixing."""

import numpy as np

from seapol import fm98


def test_linear_phase_speed_minimum():
    k = np.geomspace(10, 3000, 2000)
    c = np.array([fm98.linear_phase_speed(ki) for ki in k])
    k_min = k[np.argmin(c)]
    k_m = np.sqrt(fm98.G_EARTH / fm98.T_OVER_RHO)
    assert abs(k_min / k_m - 1.0) < 0.05
    assert abs(c.min() - 0.23) < 0.01


def test_continuation_solve_5cm():
    sol = fm98.solve_fm98_continuation(0.05, 0.20, p_target=1e-3,
                                       M=24, n_steps=6)
    assert sol.residual_norm < 1e-3
    mags = np.abs(sol.a)
    assert mags[0] > mags[1] > mags[2]          # Stokes-like decay
    assert 0.10 < mags[1] / mags[0] < 0.40      # strong bound harmonic
    assert 0.95 < sol.c / sol.c0 < 1.10
    # amplitude constraint satisfied (residual-norm-level tolerance)
    assert abs(0.5 * (sol.Y.max() - sol.Y.min()) * sol.k - 0.20) < 5e-3


def test_gauge_fix_invariants():
    rng = np.random.default_rng(0)
    a = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    a_fixed = fm98.gauge_fix(a)
    m = np.arange(1, 6)
    np.testing.assert_allclose(np.abs(a_fixed), np.abs(a), atol=1e-12)
    assert abs(np.angle(a_fixed[0]) - np.pi) < 1e-12 or \
        abs(np.angle(a_fixed[0]) + np.pi) < 1e-12
    inv_before = m * np.angle(a[0]) - np.angle(a)
    inv_after = m * np.angle(a_fixed[0]) - np.angle(a_fixed)
    np.testing.assert_allclose(np.exp(1j * inv_before),
                               np.exp(1j * inv_after), atol=1e-12)


def test_plausibility_filter():
    assert not fm98.is_physically_plausible(
        np.array([0.2 + 0j, 0.3 + 0j, 0.1 + 0j]), 0.2)   # non-monotone
    assert not fm98.is_physically_plausible(
        np.array([0.01 + 0j, 0.005 + 0j, 0.001 + 0j]), 0.2)  # |a1| != ak
    assert fm98.is_physically_plausible(
        np.array([-0.2 + 0j, 0.04 - 0.01j, 0.01 + 0j]), 0.2)


def test_table_interp_linear_limit():
    tbl = fm98.FM98Table(k_grid=np.array([100.0, 200.0]),
                         ak_grid=np.array([0.1, 0.2]),
                         a_coeffs=np.full((2, 2, 3), 0.1 + 0j),
                         c_over_c0=np.ones((2, 2)),
                         converged=np.ones((2, 2), bool),
                         M_keep=3, p=1e-3)
    a, c = tbl.interp(150.0, 0.01)           # far below tabulated ak
    assert a[0] == 0.01 + 0j and np.all(a[1:] == 0) and c == 1.0
    a, _ = tbl.interp(150.0, 0.15)
    assert abs(abs(a[0]) - 0.15) < 1e-12     # |a_1| rescaled to ak
