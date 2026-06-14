"""Hybrid surface checks: FM98 skewness signatures, variance closure,
time evolution.  The FM98 table build dominates runtime, so it is shared
module-wide."""

import numpy as np
import pytest

from seapol import (augment_fm98, build_fm98_table, generate_hybrid_surface,
                    long_wave_modulation)
from seapol.hybrid import _analytic_signal_along_wind


@pytest.fixture(scope="module")
def tiny_table():
    # Carrier band 2-8 cm: k in [78.5, 314] rad/m
    k_grid = np.geomspace(70.0, 350.0, 3)
    ak_grid = np.array([0.05, 0.15, 0.25])
    return build_fm98_table(k_grid, ak_grid, M_keep=4, M_solve=16,
                            n_steps=6)


def test_table_mostly_converged(tiny_table):
    """k = 350 (lam = 1.8 cm) never converges (the capillary resonance
    sits at the fundamental, outside the Class-1 framing), and
    k = 157, ak = 0.25 sits near the short-carrier branch fold where
    the small solve budget lands collapsed states that the c/c0 gate
    rejects to the linear fallback."""
    assert tiny_table.converged.sum() >= 5  # of 9 points


def test_skewness_signatures(tiny_table):
    """FM98 bound harmonics break Gaussianity with an along-wind
    asymmetry.  On the well-forced (published) branch the SIGNS of the
    bound-field-only statistics are regime mixtures — crest sharpening
    (sk_eta > 0) dominates for long carriers whose resonance sits
    beyond the kept harmonics, while near-resonant short steep
    carriers invert it (exact m <= 4 truncations: +0.31 at
    lam = 9 cm, ak = 0.25 vs -0.26 at lam = 4 cm, ak = 0.25) — so the
    robust assertions are non-Gaussianity and along-wind dominance,
    not fixed signs."""
    surf = generate_hybrid_surface(L=1.0, N=512, U10=7.0, table=tiny_table,
                                   rng=np.random.default_rng(0))
    dx = surf.info["dx"]
    eh = surf.eta_high - surf.eta_high.mean()
    sk_eta = float(np.mean(eh**3) / eh.std() ** 3)
    ex = np.gradient(surf.eta_high, dx, axis=1)
    ey = np.gradient(surf.eta_high, dx, axis=0)
    sk_x = float(np.mean((ex - ex.mean()) ** 3) / ex.std() ** 3)
    sk_y = float(np.mean((ey - ey.mean()) ** 3) / ey.std() ** 3)
    assert abs(sk_eta) > 0.05, sk_eta
    assert abs(sk_x) > 0.05, sk_x
    assert abs(sk_y) < abs(sk_x), (sk_y, sk_x)


def test_variance_closure(tiny_table):
    """The shrink correction keeps the hybrid variance at the linear
    (spectral-target) level."""
    surf = generate_hybrid_surface(L=1.0, N=512, U10=7.0, table=tiny_table,
                                   rng=np.random.default_rng(1))
    ratio = surf.eta.var() / surf.eta_lin.var()
    assert abs(ratio - 1.0) < 0.05, ratio


def test_augment_rejects_coarse_grid(tiny_table):
    eta = np.zeros((64, 64))
    with pytest.raises(ValueError):
        augment_fm98(eta, dx=0.05, table=tiny_table)  # dx = 5 cm


def _long_plus_short(N=512, dx=0.004, lam_L=1.0, a_L=0.02, seed=0):
    """Monochromatic long wave (+x) plus broadband short field."""
    rng = np.random.default_rng(seed)
    x = np.arange(N) * dx
    k_L = 2 * np.pi / lam_L
    eta_L = a_L * np.cos(k_L * x)[None, :] * np.ones((N, 1))
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    noise = np.fft.fft2(rng.standard_normal((N, N)))
    eta_S = np.fft.ifft2(np.where((K > 60) & (K < 400), noise, 0)).real
    eta_S *= 0.001 / eta_S.std()
    return eta_L + eta_S, k_L


def test_long_wave_modulation_phase_locking():
    """Short-wave envelope energy must peak at the prescribed long-wave
    phase and short variance must be preserved."""
    eta, k_L = _long_plus_short()
    dx = 0.004
    theta = 90.0  # enhancement mid forward face
    eta_mod, info = long_wave_modulation(eta, dx, k_split=20.0, mtf=8.0,
                                         mtf_phase_deg=theta)
    kx = 2 * np.pi * np.fft.fftfreq(512, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    Z = np.fft.fft2(eta_mod)
    eta_S_mod = np.fft.ifft2(np.where(K > 20.0, Z, 0)).real
    eta_L = np.fft.ifft2(np.where((K > 0) & (K <= 20.0), np.fft.fft2(eta),
                                  0)).real
    # envelope energy (group-scale smoothed) vs long-wave phase
    from scipy.ndimage import uniform_filter
    env2 = np.abs(_analytic_signal_along_wind(eta_S_mod, kx, kx, 0.0)) ** 2
    env2 = uniform_filter(env2, size=25)        # ~10 cm: kills speckle
    phi_L = np.angle(_analytic_signal_along_wind(eta_L, kx, kx, 0.0))
    target = np.cos(phi_L - np.deg2rad(theta))
    c = np.corrcoef(env2.ravel(), target.ravel())[0, 1]
    assert c > 0.6, c
    # short-wave variance preserved
    eta_S = eta - eta_L
    assert abs(eta_S_mod.var() / eta_S.var() - 1.0) < 0.02
    assert info["M_E_max"] > 1.5  # modulation actually strong


def test_long_wave_modulation_domain_guard():
    eta = np.zeros((64, 64))
    with pytest.raises(ValueError):
        long_wave_modulation(eta, dx=0.002, k_split=2 * np.pi / 0.5)


def test_hybrid_capillaries_bound_to_long_waves(tiny_table):
    """With the MTF stage on, bound-harmonic energy concentrates at the
    prescribed long-wave phase."""
    from scipy.ndimage import uniform_filter
    surf = generate_hybrid_surface(L=2.0, N=1024, U10=8.0,
                                   table=tiny_table, long_wave_mtf=8.0,
                                   mtf_phase_deg=45.0,
                                   k_split=2 * np.pi / 0.5,
                                   rng=np.random.default_rng(5))
    dx = surf.info["dx"]
    N = 1024
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    Z = np.fft.fft2(surf.eta_lin)
    eta_L = np.fft.ifft2(np.where((K > 0) & (K <= 2 * np.pi / 0.5), Z,
                                  0)).real
    phi_L = np.angle(_analytic_signal_along_wind(eta_L, kx, kx, 0.0))
    energy = uniform_filter(surf.eta_high**2, size=16)
    target = np.cos(phi_L - np.deg2rad(45.0))
    c = np.corrcoef(energy.ravel(), target.ravel())[0, 1]
    assert c > 0.1, c
    # decomposition still closes with modulation enabled
    np.testing.assert_allclose(
        surf.eta, surf.eta_lin + surf.eta_shrink + surf.eta_high,
        atol=1e-12)


def test_mss_closure_high_wind(tiny_table):
    """Elfouhaily is an empirical total (tuned to Cox-Munk): the FM98
    augmentation must redistribute slope variance, not add it, even when
    the harmonic tables saturate at high wind."""
    from seapol import generate_sea_surface
    kw = dict(L=2.0, N=1024, U10=10.0)
    lin = generate_sea_surface(rng=np.random.default_rng(9), **kw)
    hyb = generate_hybrid_surface(table=tiny_table, long_wave_mtf=8.0,
                                  k_split=2 * np.pi / 0.5,
                                  rng=np.random.default_rng(9), **kw)
    mss_lin = lin.slope_x.var() + lin.slope_y.var()
    mss_hyb = hyb.slope_x.var() + hyb.slope_y.var()
    assert abs(mss_hyb / mss_lin - 1.0) < 0.1, (mss_hyb, mss_lin)


def _mono_carrier(N, dx, mult_x, mult_y=0, ak=0.20):
    """Periodic monochromatic carrier on grid harmonics (mult_x, mult_y);
    returns (eta, kx1, ky1)."""
    dk = 2 * np.pi / (N * dx)
    kx1, ky1 = mult_x * dk, mult_y * dk
    k1 = np.hypot(kx1, ky1)
    x = np.arange(N) * dx
    eta = (ak / k1) * np.cos(kx1 * x[None, :] + ky1 * x[:, None])
    return eta, kx1, ky1


def test_local_carrier_matches_representative_at_true_k(tiny_table):
    """For a monochromatic carrier the local phase gradient recovers the
    true wavenumber exactly, so carrier_k='local' must reproduce the
    representative scheme evaluated at k_rep = k_true -- at any carrier
    scale, which the band-averaged k_rep cannot do."""
    N, dx = 512, 0.002
    for mult in (16, 40):                       # ~98 and ~245 rad/m
        eta, kx1, _ = _mono_carrier(N, dx, mult)
        hi_loc, _, info = augment_fm98(eta, dx, tiny_table,
                                       compensation="none",
                                       carrier_k="local")
        hi_rep, _, _ = augment_fm98(eta, dx, tiny_table,
                                    compensation="none", k_rep=float(kx1))
        assert abs(info["k_loc_p50"] - kx1) < 0.01 * kx1
        c = np.corrcoef(hi_loc.ravel(), hi_rep.ravel())[0, 1]
        assert c > 0.999, (mult, c)
        np.testing.assert_allclose(hi_loc.var(), hi_rep.var(), rtol=1e-3)


def test_local_harmonics_follow_carrier_direction(tiny_table):
    """An oblique carrier's bound harmonics must sit at exact integer
    multiples of its wavevector: ripple crests parallel to the local
    carrier, not to the wind."""
    N, dx = 512, 0.002
    eta, kx1, ky1 = _mono_carrier(N, dx, 14, 8)   # ~30 deg off-wind
    hi, _, _ = augment_fm98(eta, dx, tiny_table, compensation="none",
                            carrier_k="local")
    F = np.abs(np.fft.fft2(hi))
    kg = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.hypot(kg[None, :], kg[:, None])
    k1 = np.hypot(kx1, ky1)
    F[K < 1.5 * k1] = 0.0
    iy, ix = np.unravel_index(np.argmax(F), F.shape)
    # peak at +-m * (kx1, ky1) for integer m >= 2 (real field: the
    # conjugate twin at -k carries equal magnitude)
    m_est = abs(kg[ix] / kx1)
    assert abs(m_est - round(m_est)) < 0.05 and round(m_est) >= 2, m_est
    np.testing.assert_allclose(kg[iy] / kg[ix], ky1 / kx1, rtol=1e-6)


def test_time_evolution_shapes(tiny_table):
    surf = generate_hybrid_surface(L=0.5, N=256, U10=6.0, table=tiny_table,
                                   times=np.array([0.0, 0.1]),
                                   rng=np.random.default_rng(2))
    assert surf.eta.shape == (256, 256, 2)
    assert np.isfinite(surf.eta).all()
    assert np.isfinite(surf.slope_x).all()
    # frames evolve, components decompose
    assert np.abs(surf.eta[:, :, 1] - surf.eta[:, :, 0]).max() > 1e-6
    np.testing.assert_allclose(
        surf.eta, surf.eta_lin + surf.eta_shrink + surf.eta_high,
        atol=1e-12)
