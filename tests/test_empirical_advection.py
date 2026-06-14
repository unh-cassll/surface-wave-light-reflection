"""Empirical beta(k) fitting and orbital advection."""

from pathlib import Path

import numpy as np
import pytest

from seapol import (apply_orbital_advection, bound_fraction_from_kf_reduction,
                    generate_sea_surface, smooth_beta_curve)

ASIT_NPZ = Path(__file__).parent.parent / "demos/output/asit_kf_reduced.npz"


def test_smooth_beta_curve_fills_gaps():
    k = np.geomspace(5, 1000, 30)
    beta = np.full(30, 0.8)
    beta[10:14] = np.nan
    out = smooth_beta_curve(k, beta)
    assert np.all(np.isfinite(out))
    assert np.all((out >= 0) & (out <= 1))
    np.testing.assert_allclose(out, 0.8, atol=1e-6)


@pytest.mark.skipif(not ASIT_NPZ.exists(), reason="ASIT reduction not built")
def test_bound_fraction_from_asit():
    beta_fn = bound_fraction_from_kf_reduction(ASIT_NPZ)
    K = np.array([1.0, 5.0, 50.0, 300.0, 2000.0])
    b = beta_fn(K)
    assert np.all((b >= 0) & (b <= 0.99))
    assert b[0] < 0.1            # dominant waves stay free
    assert b[2] > 0.5            # short-gravity range mostly bound
    # works as a synthesis input
    surf = generate_sea_surface(L=8.0, N=128, U10=7.0,
                                bound_fraction=beta_fn, bound_speed=2.0,
                                rng=np.random.default_rng(0))
    assert np.isfinite(surf.eta).all()


def test_bound_speed_spectrum_mode():
    """bound_speed='spectrum' derives carrier speeds from the resolved
    long-wave band and preserves total power."""
    free = generate_sea_surface(L=8.0, N=256, U10=7.0,
                                rng=np.random.default_rng(4))
    surf = generate_sea_surface(L=8.0, N=256, U10=7.0, bound_fraction=0.9,
                                bound_speed="spectrum",
                                rng=np.random.default_rng(4))
    mss_f = free.slope_x.var() + free.slope_y.var()
    mss_s = surf.slope_x.var() + surf.slope_y.var()
    assert abs(mss_s / mss_f - 1.0) < 0.05
    with pytest.raises(ValueError):
        generate_sea_surface(L=8.0, N=64, U10=7.0, bound_fraction=0.5,
                             bound_speed="nope",
                             rng=np.random.default_rng(0))


def test_orbital_advection_preserves_short_variance():
    rng = np.random.default_rng(1)
    N, dx = 512, 0.004
    x = np.arange(N) * dx
    eta_L = 0.02 * np.cos(2 * np.pi * x / 1.0)[None, :] * np.ones((N, 1))
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    noise = np.fft.fft2(rng.standard_normal((N, N)))
    eta_S = np.fft.ifft2(np.where((K > 60) & (K < 400), noise, 0)).real
    eta_S *= 0.001 / eta_S.std()

    eta_w, info = apply_orbital_advection(eta_L + eta_S, dx, k_split=20.0)
    eta_S_w = eta_w - eta_L
    assert abs(eta_S_w.var() / eta_S.var() - 1.0) < 0.05
    assert info["D_rms"] > 0.005   # displacements actually applied


def test_orbital_advection_broadens_free_ridge():
    """Per-frame advection by propagating long waves Doppler-broadens the
    short-wave frequency content at fixed k."""
    fs, T = 24.0, 96
    times = np.arange(T) / fs
    surf = generate_sea_surface(L=2.0, N=128, U10=8.0, times=times,
                                compute_slopes=False,
                                rng=np.random.default_rng(3))
    dx = surf.info["dx"]

    def f_width(stack):
        w_t = np.hanning(T)[None, None, :]
        F = np.fft.fft(np.fft.fft(stack * w_t, axis=1), axis=2)
        P = (np.abs(F) ** 2).mean(axis=0)
        kx = 2 * np.pi * np.fft.fftfreq(128, d=dx)
        f = np.fft.fftfreq(T, d=1.0 / fs)
        i_k = np.argmin(np.abs(kx - 60.0))     # short-gravity column
        col = P[i_k, f < 0]
        fneg = -f[f < 0]
        c = (col * fneg).sum() / col.sum()
        return np.sqrt((col * (fneg - c) ** 2).sum() / col.sum())

    plain = surf.eta
    warped = np.empty_like(plain)
    for it in range(T):
        warped[:, :, it], _ = apply_orbital_advection(plain[:, :, it], dx,
                                                      k_split=20.0)
    w0, w1 = f_width(plain), f_width(warped)
    assert w1 > 1.2 * w0, (w0, w1)
