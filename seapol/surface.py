"""
FFT synthesis of a random sea surface from the Elfouhaily directional
spectrum, with optional time evolution under the gravity-capillary
dispersion relation (random-phase / Tessendorf approach).

Normalization: the realized grid variance matches the spectral target,
    var(eta) = sum Psi(kx, ky) dkx dky,
via the Hermitian construction
    Z(k) = [W(k) + conj(W(-k))] / sqrt(2),   E|W(k)|^2 = N^4 Psi dkx dky,
which preserves variance for asymmetric (one-sided) spectra as well.

Slopes are computed spectrally (i k Z) rather than by finite differences.

backend="torch" (with device="cuda" and float32 dtype for consumer
GPUs) runs the synthesis on the GPU; all downstream seapol stages
dispatch on the produced arrays, so the pipeline stays on-device.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

from .backend import adapt_rng, get_xp, xp_of
from .spectrum import (angular_frequency, cutoff_slope_variances,
                       directional_spectrum)


@dataclass
class SeaSurface:
    """Synthesized surface fields.  eta, slope_x, slope_y are (N, N) for a
    static surface or (N, N, T) for an evolving one; row index is y."""
    x: np.ndarray
    y: np.ndarray
    eta: np.ndarray
    slope_x: np.ndarray
    slope_y: np.ndarray
    t: np.ndarray | None = None
    info: dict = field(default_factory=dict)


def default_bound_ramp(k, beta_max: float,
                       k_on: float = 15.0, k_full: float = 60.0):
    """Smooth bound-fraction ramp: 0 below k_on, beta_max above k_full."""
    xp = xp_of(k)
    s = xp.clip((k - k_on) / max(k_full - k_on, 1e-9), 0.0, 1.0)
    return beta_max * s * s * (3.0 - 2.0 * s)


def generate_sea_surface(L: float, N: int, U10: float,
                         wind_dir_rad: float = 0.0,
                         times: np.ndarray | None = None,
                         fetch_m: float = 1.0e6,
                         one_sided: bool = True,
                         drag_model: str = "logistic",
                         lowk_capillary_taper: bool = True,
                         compute_slopes: bool = True,
                         bound_fraction=None,
                         bound_speed=None,
                         psi_override=None,
                         current: tuple[float, float] | None = None,
                         frame_callback=None,
                         rng=None,
                         backend: str | None = None,
                         device=None,
                         dtype=None
                         ) -> SeaSurface:
    """Generate a random sea surface eta(x, y[, t]) on an N x N grid of side
    length L [m].  `times` (seconds) triggers time evolution with
    omega^2 = g k + (sigma/rho) k^3.

    current = (Ux, Uy) [m/s] adds a uniform surface current: every
    component (free and bound) is Doppler-shifted to the absolute
    frequency omega + k . U, i.e. the whole pattern additionally
    advects at U.  Relevant for comparing synthetic k-f spectra with
    measurements at tidal sites.

    bound_fraction partitions the spectrum into free and bound parts.
    Field observations (e.g. ASIT 2019 Q(nu, theta) spectra) show that
    short-wave variance travels at the dominant-wave phase speed rather
    than on the free dispersion shell: the bound fraction beta(k) of each
    bin's power evolves non-dispersively, advected downwind at
    bound_speed (default: phase speed of the resolved spectral peak),
    while bin powers -- and hence Hs, MSS, Cox-Munk statistics -- are
    unchanged.  Pass a float beta_max (smooth ramp between 15 and 60
    rad/m) or a callable beta(K) returning values in [0, 1].

    backend/device/dtype select the array backend ("numpy" default, or
    "torch" with e.g. device="cuda"); rng may be a numpy Generator
    (drawn on the CPU and copied -- bitwise-reproducible across
    backends) or a backend.TorchRNG for device-native draws.
    """
    xp = get_xp(backend, device, dtype)
    rng = adapt_rng(rng, xp)

    dx = L / N
    dk = 2.0 * np.pi / L
    kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
    KX, KY = xp.meshgrid(kx, kx, indexing="xy")
    K = xp.hypot(KX, KY)

    # psi_override (callable Psi(KX, KY) or array, e.g. interpolated from
    # measured spectra via empirical.psi_from_asit) replaces the
    # Elfouhaily model; U10 then only parameterizes the sub-grid slope
    # extrapolation and the bound-speed default.
    if psi_override is None:
        Psi = directional_spectrum(KX, KY, U10, wind_dir_rad, fetch_m,
                                   drag_model, one_sided,
                                   lowk_capillary_taper)
    elif callable(psi_override):
        Psi = xp.asarray(psi_override(KX, KY), dtype=float)
    else:
        Psi = xp.asarray(psi_override, dtype=float)

    # Complex spectral draw: E|W|^2 = N^4 Psi dkx dky
    amp = (N * N) * xp.sqrt(Psi * dk * dk)
    W = amp * (rng.standard_normal((N, N))
               + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    W[0, 0] = 0.0

    # Free/bound power partition.  The bound part is an independent draw
    # (two distinct wave populations); same-draw splitting would add
    # coherently at t = 0 and inflate the variance.
    if bound_fraction is not None:
        if callable(bound_fraction):
            beta = xp.clip(xp.asarray(bound_fraction(K), dtype=float),
                           0.0, 1.0)
        else:
            beta = default_bound_ramp(K, float(bound_fraction))
        if bound_speed is None:
            # phase speed at the resolved spectral peak: 2-D argmax (a
            # kx-marginal would miss the peak for wind off the x axis)
            i_pk = int(xp.argmax(Psi))
            k_pk = float(K.reshape(-1)[i_pk])
            from .spectrum import phase_speed
            bound_speed = float(phase_speed(max(k_pk, dk)))
        # bound waves may ride a spectrum of carriers: bound_speed can be
        # a single speed, (speeds, weights), or "spectrum" -- carrier
        # speeds c(k_j) over the resolved long-wave/short-gravity band
        # weighted by its ring-integrated steepness k^2 Psi (modulation
        # strength scales with carrier slope)
        if isinstance(bound_speed, str):
            if bound_speed != "spectrum":
                raise ValueError(f"unknown bound_speed: {bound_speed!r}")
            from .spectrum import phase_speed
            n_sp = 8
            k_car_lo = max(dk, 0.5)
            k_car_hi = 80.0
            edges = np.geomspace(k_car_lo, k_car_hi, n_sp + 1)
            speeds = np.empty(n_sp)
            wts = np.empty(n_sp)
            for j in range(n_sp):
                ring = (K >= edges[j]) & (K < edges[j + 1])
                wts[j] = float(xp.sum((K[ring] ** 2) * Psi[ring]))
                speeds[j] = phase_speed(np.sqrt(edges[j] * edges[j + 1]))
            keep = wts > 0
            speeds, wts = speeds[keep], wts[keep]
            wts = wts / wts.sum()
        elif np.ndim(bound_speed) == 0:
            speeds = np.array([float(bound_speed)])
            wts = np.array([1.0])
        else:
            speeds, wts = (np.asarray(s, dtype=float) for s in bound_speed)
            wts = wts / wts.sum()
        W_b = []
        for w_j in wts:
            Wj = amp * xp.sqrt(beta * float(w_j)) \
                * (rng.standard_normal((N, N))
                   + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
            Wj[0, 0] = 0.0
            W_b.append(Wj)
        W = W * xp.sqrt(1.0 - beta)
        if psi_override is not None and wind_dir_rad != 0.0:
            warnings.warn(
                "psi_override is used as-is (not rotated by wind_dir_rad), "
                "but bound advection follows wind_dir_rad; supply a spectrum "
                "already rotated to wind_dir_rad or keep wind_dir_rad = 0",
                stacklevel=2)
        k_par = KX * np.cos(wind_dir_rad) + KY * np.sin(wind_dir_rad)
        # self-conjugate Nyquist bins cannot carry a complex advection
        # phase; leave them un-advected to preserve Hermitian symmetry
        if N % 2 == 0:
            k_par[N // 2, :] = 0.0
            k_par[:, N // 2] = 0.0
    else:
        beta = None
        W_b = None

    # Index map k -> -k for the fft layout
    u = xp.arange(N)
    neg = (-u) % N
    Wm = xp.conj(W[neg[:, None], neg[None, :]])
    if W_b is not None:
        Zb0 = [(Wj + xp.conj(Wj[neg[:, None], neg[None, :]])) / np.sqrt(2.0)
               for Wj in W_b]

    omega = angular_frequency(K)
    if current is not None:
        k_dot_U = KX * current[0] + KY * current[1]
        # self-conjugate Nyquist bins (k = -k there) cannot carry a
        # complex advection phase; leaving them un-advected preserves
        # Hermitian symmetry (and hence variance) exactly
        if N % 2 == 0:
            k_dot_U[N // 2, :] = 0.0
            k_dot_U[:, N // 2] = 0.0
    else:
        k_dot_U = None

    def synthesize(t: float):
        if t == 0.0:
            Z = (W + Wm) / np.sqrt(2.0)
            if W_b is not None:
                for Zb in Zb0:
                    Z = Z + Zb
        else:
            ph = xp.exp(-1j * omega * t)
            Z = (W * ph + Wm * xp.conj(ph)) / np.sqrt(2.0)
            if W_b is not None:
                # frozen-pattern advection: Z_b(k, t) = Z_b(k, 0)
                # exp(-i c_j k_par t) (signed exponent keeps Hermitian
                # symmetry); one term per carrier speed
                for c_j, Zb in zip(speeds, Zb0):
                    Z = Z + Zb * xp.exp(-1j * float(c_j) * k_par * t)
            if k_dot_U is not None:
                # uniform-current Doppler: exp(-i k.U t) is odd in k, so
                # Hermitian symmetry is preserved and the full pattern
                # (free + bound) advects at U
                Z = Z * xp.exp(-1j * k_dot_U * t)
        eta = xp.fft.ifft2(Z).real
        if compute_slopes:
            sx = xp.fft.ifft2(1j * KX * Z).real
            sy = xp.fft.ifft2(1j * KY * Z).real
        else:
            sx = sy = None
        return eta, sx, sy

    if times is None:
        eta, sx, sy = synthesize(0.0)
        t_out = None
    elif frame_callback is not None:
        # streaming mode: frames go to the callback, no stacks are kept
        # (long capillary-resolving records); the first frame is retained
        # for the info statistics
        times = np.atleast_1d(np.asarray(times, dtype=float))
        eta = sx = sy = None
        for it, t in enumerate(times):
            e, gx, gy = synthesize(float(t))
            if eta is None:
                eta, sx, sy = e, gx, gy
            frame_callback(it, float(t), e, gx, gy)
        t_out = times
    else:
        times = np.atleast_1d(np.asarray(times, dtype=float))
        T = times.size
        eta = xp.empty((N, N, T))
        sx = xp.empty((N, N, T)) if compute_slopes else None
        sy = xp.empty((N, N, T)) if compute_slopes else None
        for it, t in enumerate(times):
            e, gx, gy = synthesize(float(t))
            eta[:, :, it] = e
            if compute_slopes:
                sx[:, :, it] = gx
                sy[:, :, it] = gy
        t_out = times

    k_cutoff = np.pi / dx
    sa2, sc2 = cutoff_slope_variances(U10, k_cutoff, fetch_m=fetch_m,
                                      drag_model=drag_model,
                                      lowk_capillary_taper=lowk_capillary_taper)
    # The FFT corners resolve part of the k > pi/dx annulus; remove that
    # resolved portion from the sub-grid tail to avoid double counting
    corner = K > k_cutoff
    if bool(xp.any(corner)):
        ka_c = KX * np.cos(wind_dir_rad) + KY * np.sin(wind_dir_rad)
        kc_c = -KX * np.sin(wind_dir_rad) + KY * np.cos(wind_dir_rad)
        sa2 = max(sa2 - float(xp.sum((ka_c**2 * Psi)[corner])) * dk * dk, 0.0)
        sc2 = max(sc2 - float(xp.sum((kc_c**2 * Psi)[corner])) * dk * dk, 0.0)
    var_target = float(xp.sum(Psi)) * dk * dk
    mss_resolved = float(xp.sum(K**2 * Psi)) * dk * dk
    eta0 = eta if eta.ndim == 2 else eta[:, :, 0]

    info = dict(dx=dx, L=L, N=N, U10=U10, wind_dir_rad=wind_dir_rad,
                k_cutoff=k_cutoff,
                sigma_a2_cut=sa2, sigma_c2_cut=sc2,
                var_target=var_target,
                mss_resolved=mss_resolved,
                Hs_target=4.0 * np.sqrt(var_target),
                Hs_realized=4.0 * float(xp.std(eta0)),
                current=(None if current is None
                         else (float(current[0]), float(current[1]))))

    coords = xp.arange(N, dtype=float) * dx
    if not compute_slopes:
        sx = sy = xp.empty(0)
    return SeaSurface(x=coords, y=xp.copy(coords), eta=eta,
                      slope_x=sx, slope_y=sy, t=t_out, info=info)
