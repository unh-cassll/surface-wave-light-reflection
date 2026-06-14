"""
Near-surface scattering tests: energy closures and analytic limits of
the plane-parallel water-column Monte Carlo, frame consistency of the
table lookup against the first-order model, and renderer integration.
"""

import numpy as np
import pytest

from seapol.render import (CameraGeometry, SubpixelSlopes, make_clear_sky,
                           make_unpolarized_sky, render_facet_stokes)
from seapol.scattering import (UpwellingRadianceTable, build_upwelling_table,
                               load_table, save_table,
                               water_leaving_from_table)
from seapol.polarization import stokes_dolp
from seapol.surface import generate_sea_surface
from seapol.water import WATER_TYPES, WaterBody, WaterColumn, WaterType, \
    water_leaving_stokes


def _black_sky(dirs):
    return np.zeros(np.asarray(dirs).shape[:-1] + (4,))


def test_energy_budget_closes():
    """Every photon ends absorbed, exited, or unresolved; the budget is
    exact by construction."""
    col = WATER_TYPES["case1"].column()
    tab = build_upwelling_table(make_unpolarized_sky(1.0), col,
                                sun=(40.0, 0.0, 3.0), n_photons=20000,
                                rng=np.random.default_rng(0))
    i = tab.info
    total = (i["absorbed"] + i["exited"] + i["unresolved"]) / i["E_down"]
    np.testing.assert_allclose(total, 1.0, atol=1e-9)


def test_conservative_index_matched_closure():
    """omega_0 = 1, n = 1: everything injected eventually crosses back
    up; E_u + unresolved accounts for E_down."""
    col = WaterColumn(absorption=0.0, rayleigh_scattering=0.5,
                      particulate_scattering=0.0, n_water=1.0)
    tab = build_upwelling_table(make_unpolarized_sky(1.0), col,
                                n_photons=20000, max_events=2000,
                                rng=np.random.default_rng(1))
    ratio = (tab.info["E_u"] + tab.info["unresolved"]) / tab.info["E_down"]
    np.testing.assert_allclose(ratio, 1.0, atol=0.02)
    assert tab.info["E_u"] / tab.info["E_down"] > 0.9


def test_single_scattering_analytic_limit():
    """Normal beam, n = 1, omega << 1: E_u/E_down approaches the
    analytic single-scattering albedo of a semi-infinite Rayleigh
    medium (same benchmark as the surface Monte Carlo)."""
    omega, depol = 0.05, 0.039
    Delta = (1.0 - depol) / (1.0 + 0.5 * depol)
    R1 = omega * (Delta * (3.0 / 8.0) * (11.0 / 6.0 - 2.0 * np.log(2.0))
                  + (1 - Delta) * 0.5 * (1.0 - np.log(2.0)))
    col = WaterColumn(absorption=1.0 - omega, rayleigh_scattering=omega,
                      particulate_scattering=0.0, n_water=1.0)
    tab = build_upwelling_table(_black_sky, col, sun=(0.0, 0.0, 1.0),
                                n_photons=200000, max_events=60,
                                rng=np.random.default_rng(2))
    R = tab.info["E_u"] / tab.info["E_down"]
    sigma = np.sqrt(R1 / tab.info["n_photons"])
    assert R1 - 4 * sigma < R < 1.1 * R1 + 4 * sigma


def test_pure_absorber_dark():
    col = WaterColumn(absorption=0.5, rayleigh_scattering=0.0,
                      particulate_scattering=0.0, n_water=1.0)
    tab = build_upwelling_table(make_unpolarized_sky(1.0), col,
                                n_photons=5000,
                                rng=np.random.default_rng(3))
    assert tab.info["E_u"] == 0.0


def test_bubble_layer_brightens():
    base = WaterColumn(absorption=0.06, rayleigh_scattering=0.002,
                       particulate_scattering=0.05)
    bub = WaterColumn(absorption=0.06, rayleigh_scattering=0.002,
                      particulate_scattering=0.05,
                      bubble_scattering=2.0, bubble_efold_m=0.3)
    sky = make_unpolarized_sky(1.0)
    t0 = build_upwelling_table(sky, base, n_photons=30000,
                               rng=np.random.default_rng(4))
    t1 = build_upwelling_table(sky, bub, n_photons=30000,
                               rng=np.random.default_rng(4))
    assert t1.info["E_u"] > 2.0 * t0.info["E_u"]


def test_antisolar_brightening():
    """Upwelling radiance propagating away from the sun azimuth exceeds
    the sunward side (large-angle scattering of the refracted beam vs
    near-backscatter of a forward-peaked phase function)."""
    sky = make_clear_sky(50.0, 0.0, 1.0)
    col = WaterColumn(absorption=0.06, rayleigh_scattering=0.002,
                      particulate_scattering=0.2)
    tab = build_upwelling_table(sky, col, sun=(50.0, 0.0, 6.0),
                                n_photons=80000,
                                rng=np.random.default_rng(5))
    S = np.asarray(tab.S)
    phim = 0.5 * (np.asarray(tab.phi_edges[:-1])
                  + np.asarray(tab.phi_edges[1:]))
    sunward = S[:, np.abs(phim) < np.pi / 3, 0].mean()
    anti = S[:, np.abs(phim) > 2 * np.pi / 3, 0].mean()
    assert anti > 1.2 * sunward


def test_spectral_blue_vs_red_clear_water():
    """Clear ocean water: the 450 nm band returns far more light than
    the 650 nm band (a_w(650) ~ 50x a_w(450))."""
    wt = WATER_TYPES["clear"]
    col = wt.column(np.array([450.0, 650.0]))
    sky = make_unpolarized_sky(1.0)
    t_blue = build_upwelling_table(sky, col.at_band(0), n_photons=40000,
                                   max_events=600,
                                   rng=np.random.default_rng(6))
    t_red = build_upwelling_table(sky, col.at_band(1), n_photons=40000,
                                  max_events=600,
                                  rng=np.random.default_rng(6))
    assert t_blue.info["E_u"] > 5.0 * t_red.info["E_u"]


def _isotropic_table(L_u: float, n_water: float = 1.34, n_mu=16, n_phi=24):
    S = np.zeros((n_mu, n_phi, 4))
    S[..., 0] = L_u
    return UpwellingRadianceTable(S=S,
                                  mu_edges=np.linspace(0, 1, n_mu + 1),
                                  phi_edges=np.linspace(-np.pi, np.pi,
                                                        n_phi + 1),
                                  n_water=n_water)


def test_isotropic_table_matches_first_order():
    """A uniform unpolarized table reduces to the first-order WaterBody
    model: identical intensity transfer and identical DoLP (the Q/U
    split differs only by the meridian rotation the first-order model
    omits)."""
    R_w, E_d = 0.02, np.pi          # L_u = R_w E_d / pi = 0.02
    tab = _isotropic_table(R_w * E_d / np.pi)
    body = WaterBody(reflectance=R_w)
    rng = np.random.default_rng(7)
    th = np.deg2rad(rng.uniform(5.0, 70.0, 200))
    ph = rng.uniform(-np.pi, np.pi, 200)
    d_out = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph),
                      np.cos(th)], axis=-1)
    sx = rng.normal(0, 0.15, 200)
    sy = rng.normal(0, 0.15, 200)
    n_hat = np.stack([-sx, -sy, np.ones(200)], axis=-1)
    n_hat /= np.linalg.norm(n_hat, axis=-1, keepdims=True)

    S_tab = water_leaving_from_table(d_out, n_hat, tab)
    S_fo = water_leaving_stokes(d_out, n_hat, body, E_d)
    np.testing.assert_allclose(S_tab[..., 0], S_fo[..., 0], rtol=1e-10)
    np.testing.assert_allclose(stokes_dolp(S_tab), stokes_dolp(S_fo),
                               atol=1e-10)
    # vertical view through a level facet: rotation degenerate, exact
    S_t0 = water_leaving_from_table(np.array([0.0, 0.0, 1.0]),
                                    np.array([0.0, 0.0, 1.0]), tab)
    S_f0 = water_leaving_stokes(np.array([0.0, 0.0, 1.0]),
                                np.array([0.0, 0.0, 1.0]), body, E_d)
    np.testing.assert_allclose(S_t0, S_f0, atol=1e-12)


def test_table_lookup_interpolates():
    n_mu, n_phi = 8, 12
    tab = _isotropic_table(1.0, n_mu=n_mu, n_phi=n_phi)
    S = np.asarray(tab.S).copy()
    mu_c = 0.5 * (tab.mu_edges[:-1] + tab.mu_edges[1:])
    S[..., 0] = mu_c[:, None]            # linear in mu
    tab.S = S
    from seapol.scattering import _table_lookup
    from seapol.backend import NUMPY_XP
    out = _table_lookup(tab, np.asarray(mu_c[2:6]),
                        np.zeros(4), NUMPY_XP)
    np.testing.assert_allclose(out[..., 0], mu_c[2:6], rtol=1e-12)


def test_save_load_roundtrip(tmp_path):
    col = WATER_TYPES["clear"].column()
    tab = build_upwelling_table(make_unpolarized_sky(1.0), col,
                                n_photons=5000,
                                rng=np.random.default_rng(8))
    p = tmp_path / "tab.npz"
    save_table(p, tab)
    tab2 = load_table(p)
    np.testing.assert_allclose(tab2.S, np.asarray(tab.S))
    assert tab2.n_water == pytest.approx(tab.n_water)


def test_renderer_integration_breaks_mirror():
    """The scattering table brightens S0, reduces its relative contrast,
    and damps DoLP versus the specular-only render."""
    surf = generate_sea_surface(20.0, 48, 7.0, rng=np.random.default_rng(9))
    sky = make_clear_sky(50.0, 0.0, 1.0, turbidity=0.1)
    col = WATER_TYPES["coastal_case2"].column()
    tab = build_upwelling_table(sky, col, sun=(50.0, 0.0, 6.0),
                                n_photons=40000,
                                rng=np.random.default_rng(10))
    cam = CameraGeometry(45.0, 180.0, 50.0)
    sp = SubpixelSlopes.from_cox_munk(7.0, 100.0)
    kw = dict(sky=sky, subpixel=sp, n_subpixel=6)
    S0 = render_facet_stokes(surf.eta, surf.info["dx"], cam,
                             rng=np.random.default_rng(11), **kw)
    S1 = render_facet_stokes(surf.eta, surf.info["dx"], cam, water=tab,
                             rng=np.random.default_rng(11), **kw)
    i0, i1 = S0[..., 0], S1[..., 0]
    assert np.nanmean(i1) > 1.5 * np.nanmean(i0)
    cv0 = np.nanstd(i0) / np.nanmean(i0)
    cv1 = np.nanstd(i1) / np.nanmean(i1)
    assert cv1 < cv0
    assert np.nanmean(stokes_dolp(S1)) < 0.7 * np.nanmean(stokes_dolp(S0))


def test_water_type_presets_sane():
    for name, wt in WATER_TYPES.items():
        col = wt.column()
        assert col.absorption > 0 and col.n_bands == 1
    coastal = WATER_TYPES["coastal_case2"].column(np.array([440.0]))
    clear = WATER_TYPES["clear"].column(np.array([440.0]))
    assert float(np.asarray(coastal.absorption)[0]) \
        > 3.0 * float(np.asarray(clear.absorption)[0])


def test_spectral_column_band_mismatch_raises():
    col = WATER_TYPES["clear"].column(np.array([450.0, 550.0]))
    with pytest.raises(ValueError):
        build_upwelling_table(make_unpolarized_sky(1.0), col,
                              n_photons=100)


# ---------------------------------------------------------------------------
# Fournier-Forand particulate phase function
# ---------------------------------------------------------------------------

def _sphere_avg(P00):
    """<P00> over the sphere on a forward-resolved mu grid (= 1 for a
    normalized phase function)."""
    lin = np.linspace(0.0, 1.0, 1_000_000)
    mu = 1.0 - 2.0 * (1.0 - lin) ** 2.5
    return np.trapezoid(P00(mu), mu) / 2.0, mu


def test_ff_phase_normalized_and_backscatter():
    from seapol.water import fournier_forand_phase as ff
    avg, mu = _sphere_avg(lambda m: ff(m, 1.05, 3.5))
    assert abs(avg - 1.0) < 0.05            # <P00> = 1 (forward-trunc < 5%)

    def bbb(n, muj):
        P = ff(mu, n, muj)
        return float(np.trapezoid(np.where(mu < 0, P, 0.0), mu)
                     / np.trapezoid(P, mu))
    # backscatter fraction is in the measured ocean range and rises with
    # both the particle index and the Junge slope
    assert 0.003 < bbb(1.05, 3.5) < 0.05
    assert bbb(1.10, 3.5) > bbb(1.05, 3.5)
    assert bbb(1.05, 4.0) > bbb(1.05, 3.5)
    assert (ff(mu, 1.05, 3.5) >= -1e-9).all()


def test_ff_requires_valid_junge_slope():
    from seapol.water import fournier_forand_phase
    with pytest.raises(ValueError):
        fournier_forand_phase(np.array([0.0]), 1.05, 3.0)


def test_ff_sampler_matches_phase():
    from seapol.water import fournier_forand_phase, sample_ff_scattering
    s = sample_ff_scattering(2_000_000, np.random.default_rng(0), 1.05, 3.5)
    assert s.min() >= -1.0 and s.max() <= 1.0
    assert s.mean() > 0.9                   # strongly forward-peaked
    hist, edges = np.histogram(s, bins=50, range=(-1, 1), density=True)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    pdf = fournier_forand_phase(ctr, 1.05, 3.5) / 2.0
    sel = slice(8, 46)                      # away from the forward spike
    rel = np.abs(hist[sel] - pdf[sel]) / np.maximum(pdf[sel], 1e-9)
    assert np.median(rel) < 0.1


def test_ff_phase_mueller_polarization():
    """ff_phase_mueller carries the depolarized-Rayleigh polarization:
    P01/P00 is negative and peaks near 90 deg, like hg_phase_mueller."""
    from seapol.water import ff_phase_mueller
    mu = np.cos(np.deg2rad(np.array([20.0, 90.0, 160.0])))
    P = ff_phase_mueller(mu, 1.05, 3.5, depol=0.2)
    ratio = P[..., 0, 1] / P[..., 0, 0]
    assert ratio[1] < ratio[0] and ratio[1] < ratio[2]   # most negative at 90
    assert abs(ratio[1]) > 0.5


def test_ff_monte_carlo_closes_and_backscatters():
    """FF particulate scattering closes the photon budget and raises the
    water-leaving reflectance vs HG at the same scattering coefficient
    (its larger backscatter sends more light back up)."""
    sky = make_unpolarized_sky(1.0)
    base = dict(absorption=0.15, rayleigh_scattering=0.002,
                particulate_scattering=0.8, n_water=1.34)
    out = {}
    for tag, kw in [("hg", dict(particulate_phase="hg",
                                particulate_g=0.924)),
                    ("ff", dict(particulate_phase="ff", ff_n=1.10,
                                ff_mu_junge=3.6))]:
        col = WaterColumn(**base, **kw)
        t = build_upwelling_table(sky, col, n_photons=60000,
                                  max_events=400,
                                  rng=np.random.default_rng(1))
        i = t.info
        np.testing.assert_allclose(
            (i["absorbed"] + i["exited"] + i["unresolved"]) / i["E_down"],
            1.0, atol=1e-9)
        out[tag] = i["E_u"] / i["E_down"]
    assert out["ff"] > out["hg"]
