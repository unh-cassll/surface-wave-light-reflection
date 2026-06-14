"""Monte Carlo tracer checks: flat-surface exactness, energy closure,
periodic boundaries."""

import numpy as np

from seapol import (build_cell_tables, effective_mueller_for_incident,
                    generate_sea_surface, heightmap_intersect, trace_forward)
from seapol.polarization import fresnel_mueller


def test_flat_surface_mueller_exact():
    """On a flat surface every up-escaping ray carries exactly
    M_R(theta)/p, so the bin-mean Mueller is deterministic and the
    per-launched sums match the analytic Fresnel matrix."""
    eta = np.zeros((16, 16))
    dx = 0.5
    th_deg = 40.0
    res = effective_mueller_for_incident(eta, dx, th_deg, n_rays=40_000,
                                         rng=np.random.default_rng(0))
    M_ref, M_T_ref, _ = fresnel_mueller(
        np.array(np.cos(np.deg2rad(th_deg))), 1.34)
    p = M_ref[0, 0]

    # exact energy bookkeeping: every ray classified, unit weights
    np.testing.assert_allclose(res["R_total"] + res["T_total"], 1.0,
                               atol=1e-9)
    # branch statistics: R_total ~ Binomial(p)/n
    sigma = np.sqrt(p * (1 - p) / res["launched"])
    assert abs(res["R_total"] - p) < 4 * sigma

    # specular bin-mean is exactly M_R / p
    counts = res["counts"]
    ti, pj = np.unravel_index(np.argmax(counts), counts.shape)
    assert counts[ti, pj] == res["n_up"]  # all up rays in one bin
    np.testing.assert_allclose(res["M_bin_mean"][ti, pj], M_ref / p,
                               atol=1e-9)
    # per-launched sums approach the analytic matrix
    M_eff = res["M_eff"].sum(axis=(0, 1))
    np.testing.assert_allclose(M_eff, M_ref, atol=5 * sigma * np.max(
        np.abs(M_ref)) / p)


def test_flat_surface_single_bounce():
    eta = np.zeros((8, 8))
    rng = np.random.default_rng(1)
    n = 2000
    origins = np.stack([rng.uniform(0, 4, n), rng.uniform(0, 4, n),
                        np.full(n, 1.0)], axis=1)
    th = np.deg2rad(30.0)
    d = np.array([np.sin(th), 0.0, -np.cos(th)])
    res = trace_forward(eta, 0.5, origins, np.broadcast_to(d, (n, 3)),
                        rng=rng)
    assert np.all(res["bounces"] == 1)
    assert np.all(res["escaped_up"] | res["escaped_down"])
    # transmitted rays are flagged in-water
    assert np.array_equal(res["in_water"], res["escaped_down"])


def test_energy_closure_rough_surface():
    rng = np.random.default_rng(2)
    surf = generate_sea_surface(L=16.0, N=32, U10=7.0, rng=rng)
    res = effective_mueller_for_incident(surf.eta, surf.info["dx"], 45.0,
                                         n_rays=4000, rng=rng)
    np.testing.assert_allclose(res["R_total"] + res["T_total"], 1.0,
                               atol=1e-9)
    assert res["mean_bounces"] >= 1.0


def test_reflectance_grows_with_incidence():
    rng = np.random.default_rng(3)
    surf = generate_sea_surface(L=16.0, N=32, U10=5.0, rng=rng)
    R = []
    for th in [20.0, 60.0, 75.0]:
        res = effective_mueller_for_incident(surf.eta, surf.info["dx"], th,
                                             n_rays=4000,
                                             rng=np.random.default_rng(4))
        R.append(res["R_total"])
    assert R[0] < R[1] < R[2]


def test_rough_close_to_flat_at_moderate_angle():
    rng = np.random.default_rng(5)
    surf = generate_sea_surface(L=16.0, N=32, U10=4.0, rng=rng)
    res = effective_mueller_for_incident(surf.eta, surf.info["dx"], 40.0,
                                         n_rays=8000, rng=rng)
    M_ref, _, _ = fresnel_mueller(np.array(np.cos(np.deg2rad(40.0))), 1.34)
    R_flat = M_ref[0, 0]
    assert 0.5 * R_flat < res["R_total"] < 2.0 * R_flat


def test_periodic_intersection():
    """Rays that start outside [0, L) or cross the boundary still hit."""
    eta = np.zeros((16, 16))
    dx = 0.5  # L = 8
    tbl = build_cell_tables(eta, dx)
    th = np.deg2rad(80.0)  # grazing: long horizontal run, crosses boundary
    origins = np.array([[7.9, 4.0, 0.3],
                        [-2.0, 4.0, 0.5],
                        [25.0, 4.0, 0.5]])
    d = np.array([np.sin(th), 0.0, -np.cos(th)])
    dirs = np.broadcast_to(d, (3, 3))
    t_hit, hit_pos, n_hat = heightmap_intersect(origins, dirs, tbl)
    assert np.all(np.isfinite(t_hit))
    np.testing.assert_allclose(hit_pos[:, 2], 0.0, atol=1e-9)
    np.testing.assert_allclose(t_hit[0], 0.3 / np.cos(th), rtol=1e-9)
    np.testing.assert_allclose(n_hat, [[0, 0, 1]] * 3, atol=1e-12)
    # wrapped hit coordinates stay inside the period
    assert np.all((hit_pos[:, :2] >= 0) & (hit_pos[:, :2] < 8.0))
