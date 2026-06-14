"""k-f diagnostics: window leakage control and Q(nu) bookkeeping."""

import numpy as np
import pytest

from seapol.diagnostics import kf_slope_spectrum, q_nu, slow_nu_fraction


def _advected_wave(N=64, T=256, dx=0.01, fs=40.0, n_cyc=12, c=0.75):
    """Grid-periodic monochromatic wave advected rigidly at speed c (+x).
    The wavenumber sits exactly on a grid bin so the only leakage is
    temporal (window sidelobes)."""
    k0 = n_cyc * 2 * np.pi / (N * dx)
    x = np.arange(N) * dx
    t = np.arange(T) / fs
    eta = np.cos(k0 * (x[None, :, None] - c * t[None, None, :]))
    return np.broadcast_to(eta, (N, N, T)).copy()


def test_blackmanharris_kills_leakage_floor():
    """A pure constant-c ridge must not leak into slow phase speeds; the
    Blackman-Harris floor should sit orders of magnitude below Hann's."""
    stack = _advected_wave()
    fracs = {}
    for win in ("hann", "blackmanharris"):
        kx, f, P = kf_slope_spectrum(stack, 0.01, 40.0, window=win)
        nu, Q = q_nu(kx, f, P)
        fracs[win] = slow_nu_fraction(nu, Q, nu_split=2.0)
    assert fracs["blackmanharris"] < 0.05 * fracs["hann"], fracs
    assert fracs["blackmanharris"] < 1e-4


def test_ridge_at_right_speed():
    c = 1.5
    stack = _advected_wave(c=c, n_cyc=6)
    kx, f, P = kf_slope_spectrum(stack, 0.01, 40.0)
    i, j = np.unravel_index(np.argmax(P), P.shape)
    c_est = 2 * np.pi * f[j] / kx[i]
    assert abs(c_est / c - 1.0) < 0.1


def test_welch_segments():
    stack = _advected_wave(T=512)
    kx, f1, P1 = kf_slope_spectrum(stack, 0.01, 40.0, n_seg=1)
    kx, f3, P3 = kf_slope_spectrum(stack, 0.01, 40.0, n_seg=3)
    assert f3.size < f1.size              # shorter segments
    assert np.isfinite(P3).all() and P3.sum() > 0
    with pytest.raises(ValueError):
        kf_slope_spectrum(stack, 0.01, 40.0, n_seg=200)
