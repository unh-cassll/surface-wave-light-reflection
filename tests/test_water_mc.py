"""In-water polarized scattering: phase matrix limits, sampler
statistics, energy closure, and the single-scattering analytic limit."""

import numpy as np

from seapol import (WaterOptics, effective_mueller_for_incident,
                    generate_sea_surface, rayleigh_phase_mueller,
                    sample_rayleigh_scattering, trace_forward)
from seapol.polarization import fresnel_mueller


def _delta(depol):
    return (1.0 - depol) / (1.0 + 0.5 * depol)


def test_water_optics_properties():
    w = WaterOptics(absorption=0.1, scattering=0.05)
    np.testing.assert_allclose(w.attenuation, 0.15)
    np.testing.assert_allclose(w.albedo, 1.0 / 3.0)
    assert WaterOptics(absorption=0.0, scattering=0.0).albedo == 0.0


def test_phase_matrix_normalization_and_limits():
    mu = np.linspace(-1.0, 1.0, 20001)
    for depol in (0.0, 0.039, 0.1):
        P = rayleigh_phase_mueller(mu, depol)
        # <P00> over the sphere = 1
        norm = np.trapezoid(P[:, 0, 0], mu) / 2.0
        np.testing.assert_allclose(norm, 1.0, atol=1e-6)
        # symmetric off-diagonals, physical bound |P01| <= P00
        np.testing.assert_allclose(P[:, 0, 1], P[:, 1, 0])
        assert np.all(np.abs(P[:, 0, 1]) <= P[:, 0, 0] + 1e-12)

    # pure Rayleigh: full polarization at 90 deg
    P90 = rayleigh_phase_mueller(np.array(0.0), 0.0)
    np.testing.assert_allclose(-P90[0, 1] / P90[0, 0], 1.0, rtol=1e-12)

    # depolarized: DoLP at 90 deg = (1 - d) / (1 + d)
    d = 0.039
    P90 = rayleigh_phase_mueller(np.array(0.0), d)
    np.testing.assert_allclose(-P90[0, 1] / P90[0, 0],
                               (1.0 - d) / (1.0 + d), rtol=1e-12)


def test_sampler_matches_phase_function():
    depol = 0.039
    rng = np.random.default_rng(11)
    mu = sample_rayleigh_scattering(300_000, rng, depol)
    Delta = _delta(depol)

    # moments of the mixture pdf
    np.testing.assert_allclose(mu.mean(), 0.0, atol=0.005)
    np.testing.assert_allclose(mu.var(), Delta * 0.4 + (1 - Delta) / 3.0,
                               atol=0.005)

    # Kolmogorov distance against the analytic CDF
    s = np.sort(mu)
    F = (Delta * (3.0 / 8.0) * (s + s**3 / 3.0 + 4.0 / 3.0)
         + (1 - Delta) * 0.5 * (s + 1.0))
    emp = np.arange(1, s.size + 1) / s.size
    assert np.max(np.abs(emp - F)) < 0.004


def test_pure_absorber_closure():
    """b = 0: every transmitted ray is absorbed at its first collision,
    so R is the Fresnel reflectance and R + A = 1 exactly."""
    eta = np.zeros((8, 8))
    water = WaterOptics(absorption=1.0, scattering=0.0)
    res = effective_mueller_for_incident(eta, 0.5, 40.0, n_rays=20_000,
                                         water=water,
                                         rng=np.random.default_rng(2))
    assert res["T_total"] == 0.0
    assert res["W_total"] == 0.0
    assert res["R_water"] == 0.0
    np.testing.assert_allclose(res["R_total"] + res["A_total"], 1.0,
                               atol=1e-9)
    p = fresnel_mueller(np.array(np.cos(np.deg2rad(40.0))), 1.34)[0][0, 0]
    sigma = np.sqrt(p * (1 - p) / res["launched"])
    assert abs(res["R_total"] - p) < 4 * sigma


def test_single_scattering_analytic_limit():
    """Index-matched surface (n = 1), normal-incidence beam, omega << 1:
    the up-escape fraction approaches the analytic single-scattering
    albedo of a semi-infinite Rayleigh medium,

        R1 = omega [Delta (3/8)(11/6 - 2 ln 2) + (1 - Delta)(1 - ln 2)/2],

    from int p(Theta) mu/(1 + mu) dOmega over the upward hemisphere."""
    omega = 0.05
    depol = 0.039
    Delta = _delta(depol)
    R1 = omega * (Delta * (3.0 / 8.0) * (11.0 / 6.0 - 2.0 * np.log(2.0))
                  + (1 - Delta) * 0.5 * (1.0 - np.log(2.0)))

    water = WaterOptics(absorption=1.0 - omega, scattering=omega,
                        depolarization=depol)
    res = effective_mueller_for_incident(np.zeros((8, 8)), 0.5, 0.0,
                                         n_rays=200_000, n_water=1.0,
                                         max_bounces=30, water=water,
                                         rng=np.random.default_rng(3))
    np.testing.assert_allclose(
        res["R_total"] + res["T_total"] + res["A_total"] + res["W_total"],
        1.0, atol=1e-9)
    # index-matched: no glint, all up-escapes are water-leaving
    assert res["R_glint"] == 0.0
    # multiple scattering only adds O(omega^2): R in [R1, 1.1 R1] up to noise
    sigma = np.sqrt(R1 * (1 - R1) / res["launched"])
    assert R1 - 4 * sigma < res["R_total"] < 1.1 * R1 + 4 * sigma


def test_conservative_scattering_closure():
    """a = 0 in semi-infinite water: nothing is absorbed or transmitted;
    photons diffuse back out, R + W = 1 exactly and R -> 1."""
    rng = np.random.default_rng(4)
    surf = generate_sea_surface(L=16.0, N=32, U10=5.0, rng=rng)
    water = WaterOptics(absorption=0.0, scattering=0.5)
    res = effective_mueller_for_incident(surf.eta, surf.info["dx"], 40.0,
                                         n_rays=1500, max_bounces=300,
                                         water=water, rng=rng)
    assert res["A_total"] == 0.0
    assert res["T_total"] == 0.0
    np.testing.assert_allclose(res["R_total"] + res["W_total"], 1.0,
                               atol=1e-9)
    assert res["R_total"] > 0.85


def test_single_scatter_dolp_pattern():
    """Index-matched surface, pure Rayleigh scattering (depol = 0):
    up-escaping rays with exactly one volume scattering carry the exact
    single-scattering DoLP (1 - mu^2) / (1 + mu^2) for unpolarized
    input, ray by ray."""
    rng = np.random.default_rng(5)
    n = 60_000
    water = WaterOptics(absorption=0.0, scattering=1.0, depolarization=0.0)
    origins = np.stack([rng.uniform(0, 4, n), rng.uniform(0, 4, n),
                        np.full(n, 0.5)], axis=1)
    dirs = np.broadcast_to([0.0, 0.0, -1.0], (n, 3))
    res = trace_forward(np.zeros((8, 8)), 0.5, origins, dirs, n_water=1.0,
                        max_bounces=8, rng=rng, water=water)

    sel = res["escaped_up"] & (res["scatters"] == 1)
    assert sel.sum() > 1000
    S = res["mueller"][sel][:, :, 0]           # response to (1, 0, 0, 0)
    mu = -res["dir"][sel][:, 2]                # cos(Theta) from the -z beam
    dolp = np.hypot(S[:, 1], S[:, 2]) / S[:, 0]
    np.testing.assert_allclose(dolp, (1 - mu**2) / (1 + mu**2), atol=1e-9)


def test_rough_surface_energy_closure_with_water():
    """Full physics on a rough surface: exact four-way energy closure
    and a nonzero water-leaving share."""
    rng = np.random.default_rng(6)
    surf = generate_sea_surface(L=16.0, N=32, U10=7.0, rng=rng)
    water = WaterOptics(absorption=0.1, scattering=0.05)
    res = effective_mueller_for_incident(surf.eta, surf.info["dx"], 45.0,
                                         n_rays=4000, max_bounces=60,
                                         water=water, rng=rng)
    np.testing.assert_allclose(
        res["R_total"] + res["T_total"] + res["A_total"] + res["W_total"],
        1.0, atol=1e-9)
    np.testing.assert_allclose(res["R_total"],
                               res["R_glint"] + res["R_water"], atol=1e-12)
    assert res["R_water"] > 0.0
    assert res["A_total"] > 0.5
