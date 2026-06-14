"""
Wavenumber-frequency and inverse-phase-speed diagnostics for evolving
surface stacks.

The spatial directions are periodic (FFT synthesis), so no spatial window
is needed; the time axis is finite and the strong non-dispersive bound
ridges leak across frequency through the window sidelobes.  A Hann window
(-31 dB sidelobes) leaves a leakage floor that contaminates slow-phase-
speed statistics; the default Blackman-Harris window (-92 dB) suppresses
it, with optional Welch segment averaging for variance reduction.
"""

from __future__ import annotations

import numpy as np

from .backend import xp_of

__all__ = ["kf_slope_spectrum", "q_nu", "slow_nu_fraction"]


def _window(name: str, n: int) -> np.ndarray:
    if name == "hann":
        return np.hanning(n)
    if name == "blackmanharris":
        from scipy.signal.windows import blackmanharris
        return blackmanharris(n)
    raise ValueError(f"unknown window: {name}")


def kf_slope_spectrum(stack, dx: float, fs: float,
                      window: str = "blackmanharris",
                      n_seg: int = 1, overlap: float = 0.5,
                      slope_weight: bool = True):
    """Along-wind (kx > 0, f > 0) spectrum of an (ny, nx, T) stack,
    averaged over rows and over Welch time segments.

    slope_weight multiplies by kx^2 (along-wind slope spectrum).
    Returns (kx, f, P)."""
    xp = xp_of(stack)
    Ny, Nx, T = stack.shape
    if n_seg < 1:
        raise ValueError("n_seg >= 1")
    seg_len = T if n_seg == 1 else int(T / (1 + (n_seg - 1) * (1 - overlap)))
    if seg_len < 16:
        raise ValueError("record too short for n_seg segments")
    step = max(1, int(seg_len * (1 - overlap)))
    starts = [i * step for i in range(n_seg)]
    if starts[-1] + seg_len > T:
        raise ValueError("segments exceed record; reduce n_seg")

    w_t = xp.asarray(_window(window, seg_len), dtype=float)[None, None, :]
    P = None
    for s0 in starts:
        seg = stack[:, :, s0:s0 + seg_len]
        seg = seg - xp.mean(seg, axis=2, keepdims=True)
        F = xp.fft.fft(xp.fft.fft(seg * w_t, axis=1), axis=2)
        Pi = xp.mean(xp.abs(F) ** 2, axis=0)
        P = Pi if P is None else P + Pi
    P = P / len(starts)

    kx = 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    f = xp.fft.fftfreq(seg_len, d=1.0 / fs)
    ik = kx > 0
    jf = f < 0          # fold: downwind-propagating power lands at f > 0
    P = xp.flip(P[ik][:, jf], axis=1)
    if slope_weight:
        P = P * kx[ik][:, None] ** 2
    return kx[ik], xp.flip(-f[jf], axis=0), P


def q_nu(kx, f, P, n_nu: int = 60, nu_max: float = 4.0,
         f_min: float = 0.3):
    """Inverse-phase-speed spectrum Q(nu), nu = kx / (2 pi f) [s/m]."""
    xp = xp_of(kx, f, P)
    KX, F = xp.meshgrid(kx, f, indexing="ij")
    nu = KX / (2 * np.pi * xp.maximum(F, 1e-9))
    edges = xp.linspace(0.05, nu_max, n_nu + 1)
    Q = xp.zeros(n_nu)
    for i in range(n_nu):
        m = (nu >= edges[i]) & (nu < edges[i + 1]) & (F > f_min)
        Q[i] = xp.sum(P[m]) / (edges[i + 1] - edges[i])
    return 0.5 * (edges[:-1] + edges[1:]), Q


def slow_nu_fraction(nu, Q, nu_split: float = 2.0) -> float:
    """Fraction of Q at phase speeds slower than 1/nu_split."""
    xp = xp_of(nu, Q)
    dq = xp.gradient(nu)
    sel = nu > nu_split
    return float(xp.nansum((Q * dq)[sel]) / xp.nansum(Q * dq))
