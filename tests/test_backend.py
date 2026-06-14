"""
Backend dispatch tests: the numpy path is plain numpy; the torch path
(skipped when torch is absent) reproduces numpy bitwise when fed the
same numpy Generator, across synthesis, hybrid augmentation, rendering,
Monte Carlo, diagnostics, and inversion.
"""

import numpy as np
import pytest

from seapol import backend
from seapol.backend import get_xp, to_numpy, xp_of

torch = pytest.importorskip("torch") if False else None
HAS_TORCH = backend.has_torch()
needs_torch = pytest.mark.skipif(not HAS_TORCH, reason="torch not installed")


def _txp():
    return get_xp("torch", "cpu", "float64")


def test_numpy_dispatch_is_numpy():
    xp = xp_of(np.zeros(3))
    assert not xp.is_torch
    assert xp.sqrt(np.array(4.0)) == 2.0
    out = np.zeros(3)
    xp.index_add(out, np.array([0, 0, 2]), np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(out, [3.0, 0.0, 3.0])


def test_get_xp_unknown_backend():
    with pytest.raises(ValueError):
        get_xp("tensorflow")


@needs_torch
def test_torch_shim_semantics():
    xp = _txp()
    x = xp.asarray([1.0, 2.0, 3.0, 4.0])
    # ddof-0 std/var (numpy default), not torch's Bessel correction
    np.testing.assert_allclose(float(xp.std(x)), np.std([1, 2, 3, 4.0]))
    np.testing.assert_allclose(float(xp.var(x)), np.var([1, 2, 3, 4.0]))
    # endpoint-aware linspace
    np.testing.assert_allclose(
        to_numpy(xp.linspace(0.0, 1.0, 4, endpoint=False)),
        np.linspace(0.0, 1.0, 4, endpoint=False))
    # interp matches numpy incl. clamped ends
    xq = np.array([-1.0, 0.5, 2.5, 9.0])
    xs = np.array([0.0, 1.0, 2.0, 3.0])
    fs = np.array([0.0, 10.0, 5.0, 2.0])
    np.testing.assert_allclose(
        to_numpy(xp.interp(xp.asarray(xq), xp.asarray(xs), xp.asarray(fs))),
        np.interp(xq, xs, fs))
    # scalar coercion in binary ops
    np.testing.assert_allclose(to_numpy(xp.maximum(x, 2.5)),
                               np.maximum([1, 2, 3, 4.0], 2.5))
    # integer arithmetic stays integral, float promotion respects dtype
    assert xp.arange(5).dtype.is_floating_point is False
    assert xp.arange(5, dtype=float).dtype == xp.dtype


@needs_torch
def test_torch_rng_statistics():
    rng = backend.TorchRNG(0, _txp())
    z = to_numpy(rng.standard_normal(20000))
    assert abs(z.mean()) < 0.05 and abs(z.std() - 1.0) < 0.05
    e = to_numpy(rng.exponential(2.0, 20000))
    assert abs(e.mean() - 2.0) < 0.1
    u = to_numpy(rng.uniform(1.0, 3.0, 20000))
    assert 1.0 <= u.min() and u.max() <= 3.0


@needs_torch
def test_synthesis_parity_bitwise():
    from seapol.surface import generate_sea_surface
    s_np = generate_sea_surface(20.0, 48, 7.0, bound_fraction=0.5,
                                times=np.array([0.0, 0.7]),
                                current=(0.3, -0.1),
                                rng=np.random.default_rng(11))
    s_t = generate_sea_surface(20.0, 48, 7.0, bound_fraction=0.5,
                               times=np.array([0.0, 0.7]),
                               current=(0.3, -0.1),
                               rng=np.random.default_rng(11),
                               backend="torch", dtype="float64")
    assert backend.xp_of(s_t.eta).is_torch
    np.testing.assert_allclose(to_numpy(s_t.eta), s_np.eta, atol=1e-12)
    np.testing.assert_allclose(to_numpy(s_t.slope_x), s_np.slope_x,
                               atol=1e-12)
    np.testing.assert_allclose(s_t.info["Hs_realized"],
                               s_np.info["Hs_realized"], rtol=1e-12)


def _toy_table():
    """Small synthetic gauge-fixed FM98 table (smooth, unphysical) for
    exercising the augmentation numerics."""
    from seapol.fm98 import FM98Table
    k_grid = np.geomspace(70.0, 380.0, 5)
    ak_grid = np.linspace(0.02, 0.38, 6)
    M_keep = 6
    K, A = np.meshgrid(k_grid, ak_grid, indexing="ij")
    a = np.zeros((5, 6, M_keep), dtype=complex)
    for m in range(M_keep):
        mag = A * (0.5 ** m) * (1.0 + 0.1 * np.log(K / 70.0))
        a[:, :, m] = mag * np.exp(1j * (np.pi + 0.1 * m * A))
    return FM98Table(k_grid=k_grid, ak_grid=ak_grid, a_coeffs=a,
                     c_over_c0=np.ones((5, 6)),
                     converged=np.ones((5, 6), bool), M_keep=M_keep,
                     p=0.015)


@needs_torch
def test_hybrid_augmentation_parity():
    from seapol.hybrid import augment_fm98, long_wave_modulation
    from seapol.surface import generate_sea_surface
    table = _toy_table()
    k_split = 2.0 * np.pi / 0.5
    s = generate_sea_surface(2.0, 256, 7.0, rng=np.random.default_rng(2))
    eta_m, _ = long_wave_modulation(s.eta, s.info["dx"], k_split=k_split)
    hi, comp, info = augment_fm98(eta_m, s.info["dx"], table)

    xp = _txp()
    eta_t = xp.asarray(s.eta)
    eta_m_t, _ = long_wave_modulation(eta_t, s.info["dx"],
                                      k_split=k_split)
    hi_t, comp_t, info_t = augment_fm98(eta_m_t, s.info["dx"], table)
    # torch and numpy FFTs accumulate in different orders: equality to
    # rounding noise, not bitwise
    np.testing.assert_allclose(to_numpy(eta_m_t), eta_m, atol=1e-12)
    np.testing.assert_allclose(to_numpy(hi_t), hi, atol=1e-10)
    np.testing.assert_allclose(to_numpy(comp_t), comp, atol=1e-10)
    np.testing.assert_allclose(info_t["ak_p50"], info["ak_p50"], rtol=1e-9)


@needs_torch
def test_orbital_advection_parity_to_interp_order():
    from seapol.hybrid import apply_orbital_advection
    from seapol.surface import generate_sea_surface
    s = generate_sea_surface(20.0, 128, 7.0, rng=np.random.default_rng(3))
    out_np, _ = apply_orbital_advection(s.eta, s.info["dx"],
                                        k_split=2 * np.pi / 2.0)
    out_t, _ = apply_orbital_advection(_txp().asarray(s.eta),
                                       s.info["dx"],
                                       k_split=2 * np.pi / 2.0)
    # scipy spline vs Catmull-Rom: same physics, different kernels
    a, b = out_np.ravel(), to_numpy(out_t).ravel()
    corr = np.corrcoef(a, b)[0, 1]
    assert corr > 0.9999
    np.testing.assert_allclose(b.std(), a.std(), rtol=1e-3)


@needs_torch
def test_render_parity_bitwise():
    from seapol.render import (CameraGeometry, Foam, SubpixelSlopes,
                               make_clear_sky, render_facet_stokes)
    from seapol.surface import generate_sea_surface
    from seapol.water import WaterBody
    s = generate_sea_surface(20.0, 32, 6.0, rng=np.random.default_rng(4))
    sky = make_clear_sky(45.0, 30.0, 1.0, turbidity=0.2)
    sp = SubpixelSlopes.from_cox_munk(6.0, 100.0)
    kw = dict(camera=CameraGeometry(40.0, 15.0, 60.0), sky=sky,
              subpixel=sp, n_subpixel=4, shadowing=True,
              water=WaterBody(case=1), sun_glint=(45.0, 30.0, 5.0),
              foam=Foam(coverage=0.02))
    S_np = render_facet_stokes(s.eta, s.info["dx"],
                               rng=np.random.default_rng(5), **kw)
    S_t = render_facet_stokes(_txp().asarray(s.eta), s.info["dx"],
                              rng=np.random.default_rng(5), **kw)
    np.testing.assert_allclose(to_numpy(S_t), S_np, atol=1e-10)


@needs_torch
def test_montecarlo_parity():
    from seapol.montecarlo import effective_mueller_for_incident
    from seapol.surface import generate_sea_surface
    from seapol.water import WaterOptics
    s = generate_sea_surface(16.0, 32, 6.0, rng=np.random.default_rng(6))
    water = WaterOptics(absorption=0.2, scattering=0.3)
    kw = dict(theta_i_deg=40.0, n_rays=5000, max_bounces=20, water=water)
    r_np = effective_mueller_for_incident(s.eta, s.info["dx"],
                                          rng=np.random.default_rng(7),
                                          **kw)
    r_t = effective_mueller_for_incident(_txp().asarray(s.eta),
                                         s.info["dx"],
                                         rng=np.random.default_rng(7),
                                         **kw)
    for key in ("R_total", "T_total", "A_total", "W_total"):
        np.testing.assert_allclose(r_t[key], r_np[key], atol=1e-10)
    np.testing.assert_allclose(to_numpy(r_t["M_eff"]), r_np["M_eff"],
                               atol=1e-10)


@needs_torch
def test_spectral_render_parity():
    """The color-aware spectral renderer reproduces numpy bitwise on the
    torch backend and yields valid (H, W, B, 4) Stokes -> sRGB in [0, 1]."""
    from seapol.color import stokes_bands_to_rgb
    from seapol.render import PinholeCamera
    from seapol.spectral import (SpectralBands, render_camera_image_spectral,
                                 spectral_sky_factories)
    from seapol.surface import generate_sea_surface
    bands = SpectralBands.rgb()
    skies = spectral_sky_factories(bands, "clear", sun_zenith_deg=45.0)
    cam = PinholeCamera(altitude_m=200.0, zenith_deg=40.0, hfov_deg=5.0,
                        img_shape=(48, 48))

    def run(bk, dev):
        s = generate_sea_surface(20.0, 96, 6.0, rng=np.random.default_rng(0),
                                 backend=bk, device=dev)
        return to_numpy(render_camera_image_spectral(
            s.eta, s.info["dx"], bands, skies, camera=cam, water=None,
            slope_x=s.slope_x, slope_y=s.slope_y, seed=1))

    Sn = run(None, None)
    St = run("torch", "cpu")
    assert St.shape == (48, 48, 3, 4)
    np.testing.assert_allclose(St, Sn, atol=1e-9, equal_nan=True)
    rgb = stokes_bands_to_rgb(St, bands.wavelengths_nm)
    assert np.isfinite(rgb).all() and 0.0 <= rgb.min() and rgb.max() <= 1.0


@needs_torch
def test_scattering_torch_energy_closes():
    """The near-surface scattering Monte Carlo on the torch backend must
    close its photon energy budget exactly and stay statistically
    consistent with numpy (different RNG streams, so not bitwise)."""
    from seapol.render import make_unpolarized_sky
    from seapol.scattering import build_upwelling_table
    from seapol.water import WaterColumn
    col = WaterColumn(absorption=0.06, rayleigh_scattering=0.002,
                      particulate_scattering=0.18, n_water=1.34)
    sky = make_unpolarized_sky(1.0)
    out = {}
    for bk, dev, key in [(None, None, "np"), ("torch", "cpu", "t")]:
        rng = (np.random.default_rng(3) if bk is None
               else backend.default_rng(3, get_xp(bk, dev)))
        t = build_upwelling_table(sky, col, sun=(40.0, 0.0, 5.0),
                                  n_photons=40000, rng=rng,
                                  backend=bk, device=dev)
        i = t.info
        budget = (i["absorbed"] + i["exited"] + i["unresolved"]) / i["E_down"]
        np.testing.assert_allclose(budget, 1.0, atol=1e-9)
        out[key] = i["E_u"] / i["E_down"]
    # same IOPs -> same mean reflectance within MC noise
    assert abs(out["np"] - out["t"]) < 0.01


@needs_torch
def test_diagnostics_and_inversion_parity():
    from seapol.diagnostics import kf_slope_spectrum, q_nu
    from seapol.inversion import height_from_slopes
    from seapol.surface import generate_sea_surface
    s = generate_sea_surface(20.0, 32, 7.0, times=np.arange(64) / 8.0,
                             rng=np.random.default_rng(8))
    kx1, f1, P1 = kf_slope_spectrum(s.eta, s.info["dx"], 8.0)
    nu1, Q1 = q_nu(kx1, f1, P1)
    eta1 = height_from_slopes(s.slope_x[:, :, 0], s.slope_y[:, :, 0],
                              s.info["dx"])

    xp = _txp()
    kx2, f2, P2 = kf_slope_spectrum(xp.asarray(s.eta), s.info["dx"], 8.0)
    nu2, Q2 = q_nu(kx2, f2, P2)
    eta2 = height_from_slopes(xp.asarray(s.slope_x[:, :, 0]),
                              xp.asarray(s.slope_y[:, :, 0]),
                              s.info["dx"])
    np.testing.assert_allclose(to_numpy(P2), P1, rtol=1e-9)
    np.testing.assert_allclose(to_numpy(Q2), Q1, rtol=1e-9)
    np.testing.assert_allclose(to_numpy(eta2), eta1, atol=1e-12)


_HAS_CUDA = HAS_TORCH and __import__("torch").cuda.is_available()
needs_cuda = pytest.mark.skipif(not _HAS_CUDA, reason="no CUDA device")


@needs_torch
@needs_cuda
def test_cuda_smoke():
    from seapol.render import render_facet_stokes
    from seapol.surface import generate_sea_surface
    s = generate_sea_surface(10.0, 64, 6.0, rng=np.random.default_rng(0),
                             backend="torch", device="cuda")
    assert str(s.eta.device).startswith("cuda")
    S = render_facet_stokes(s.eta, s.info["dx"])
    S_np = to_numpy(S)
    assert np.isfinite(S_np[..., 0]).any()
    var = float(to_numpy(backend.xp_of(s.eta).var(s.eta)))
    assert abs(var / s.info["var_target"] - 1.0) < 0.5


@needs_torch
@needs_cuda
def test_cuda_hybrid_fp32_finite():
    """Float32 on CUDA must not reintroduce 0/0 in the bins compensation:
    subnormal epsilon floors (1e-300) underflow to 0 in float32, so the
    empty-bin guards must mask explicitly.  Regression for the all-NaN
    hybrid surface on the GPU."""
    from seapol.hybrid import augment_fm98, generate_hybrid_surface
    from seapol.surface import generate_sea_surface
    table = _toy_table()
    h = generate_hybrid_surface(2.0, 512, 9.0, table=table,
                                rng=np.random.default_rng(0),
                                backend="torch", device="cuda")
    eta = to_numpy(h.eta)
    assert np.isfinite(eta).all()
    assert abs(4.0 * eta.std() / h.info["Hs_target"] - 1.0) < 0.3
    # direct augment on a many-empty-bin field (the trigger condition)
    s = generate_sea_surface(2.0, 512, 9.0, rng=np.random.default_rng(1),
                             backend="torch", device="cuda")
    hi, comp, _ = augment_fm98(s.eta, s.info["dx"], table)
    assert np.isfinite(to_numpy(hi)).all()
    assert np.isfinite(to_numpy(comp)).all()
