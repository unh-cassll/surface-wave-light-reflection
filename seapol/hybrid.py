"""
Phase-locked hybrid sea surface: linear Elfouhaily synthesis plus
Fedorov-Melville (1998) bound harmonics injected by envelope demodulation
(Melville & Fedorov 2015 picture).

On a capillary-resolving grid (dx <~ 2 mm) the Elfouhaily spectrum already
carries the parasitic-capillary variance; what linear random-phase
synthesis cannot supply is the phase coherence between short-gravity
carriers and their bound harmonics.  The augmentation:

  1. band-passes the field to a short-gravity carrier band
     (default 2-8 cm wavelength),
  2. takes the 2D analytic signal along the wind to get the local
     envelope amplitude A(r) and carrier phase phi(r),
  3. estimates the local carrier wavevector k_loc(r) = |grad phi| from
     the analytic-signal phase gradient (carrier_k="local", default),
     so every wave group is treated at its own scale and direction,
  4. forms the pointwise slope map ak(r) = A(r) k_loc(r),
  5. injects bound harmonics with FM98 Class-1 coefficients,
         Y_m(r) = (|a_m|/k_loc) cos(m phi(r) + m phi_1 - phi_m),
     interpolated from a gauge-fixed (k, ak) table at the local
     (k_loc, ak),
  6. shrinks the carrier fundamental to conserve band variance.

The phase formula (m phi + m phi_1 - phi_m) places the bound harmonics on
the forward face of the carrier, producing the two FM98 signatures absent
from Gaussian synthesis: positive elevation skewness (sharp crests) and
negative along-wind slope skewness (steep forward faces).  Because the
harmonic wavevector is m grad(phi), the parasitic ripple trains align
crest-parallel with their local carrier and inherit its propagation
direction; with carrier_k="local" their wavelength also follows the
local carrier scale (the FM98 capillary bump sits at the harmonic
number m* where c(m* k) ~ c(k), so m* shifts with k_loc) rather than a
single band-averaged k_rep.

The synthesis closes variance on the spectral target, and one-sided
spreading puts the whole carrier band in the downwind half-plane
captured by the analytic signal, so the envelope slope map operates at
the physically correct amplitude (parasitic capillaries ignite only for
ak >~ 0.2).

The FM98 coefficient tables are solved on the CPU (tiny dense Newton
systems); their application to the field dispatches on the input arrays,
so torch-backed frames are augmented on-device.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .backend import xp_of
from .fm98 import FM98Table, build_fm98_table
from .surface import generate_sea_surface

__all__ = ["HybridSeaSurface", "apply_orbital_advection", "augment_fm98",
           "default_fm98_table", "generate_hybrid_surface",
           "long_wave_modulation"]


@dataclass
class HybridSeaSurface:
    """Hybrid surface fields; (N, N) or (N, N, T), row index is y.
    eta = eta_lin + eta_shrink + eta_high."""
    x: np.ndarray
    y: np.ndarray
    eta: np.ndarray
    eta_lin: np.ndarray
    eta_high: np.ndarray
    eta_shrink: np.ndarray
    slope_x: np.ndarray
    slope_y: np.ndarray
    t: np.ndarray | None = None
    info: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Band-pass and analytic signal
# ---------------------------------------------------------------------------

def _bandpass(Z, kx, ky, k_lo: float, k_hi: float):
    xp = xp_of(Z)
    K2 = kx[None, :] ** 2 + ky[:, None] ** 2
    mask = (K2 >= k_lo**2) & (K2 <= k_hi**2)
    Z_band = xp.where(mask, Z, xp.zeros_like(Z))
    return xp.fft.ifft2(Z_band).real, Z_band


def _analytic_signal_along_wind(eta_band, kx, ky, wind_dir_rad: float):
    """2D analytic signal keeping the downwind half-plane (k . w > 0)."""
    xp = xp_of(eta_band)
    wx, wy = np.cos(wind_dir_rad), np.sin(wind_dir_rad)
    F = xp.fft.fft2(eta_band)
    k_dot_w = kx[None, :] * wx + ky[:, None] * wy
    F = F * xp.where(k_dot_w > 0, 2.0, xp.where(k_dot_w < 0, 0.0, 1.0))
    return xp.fft.ifft2(F)


# ---------------------------------------------------------------------------
# Orbital advection (geometric long-wave / short-wave coupling)
# ---------------------------------------------------------------------------

def _periodic_warp(f, fi, fj, xp):
    """Sample f at fractional indices (fi rows, fj cols) with periodic
    wrap.  numpy uses the scipy cubic-spline path; torch uses periodic
    Catmull-Rom bicubic (no spline prefilter, so the two backends agree
    to interpolation order, not bitwise)."""
    if not xp.is_torch:
        from scipy.ndimage import map_coordinates
        return map_coordinates(f, [fi, fj], order=3, mode="grid-wrap")
    N = f.shape[0]
    i0 = xp.astype(xp.floor(fi), int)
    j0 = xp.astype(xp.floor(fj), int)
    ti = fi - i0
    tj = fj - j0

    def cr_weights(t):
        t2 = t * t
        t3 = t2 * t
        return (-0.5 * t3 + t2 - 0.5 * t,
                1.5 * t3 - 2.5 * t2 + 1.0,
                -1.5 * t3 + 2.0 * t2 + 0.5 * t,
                0.5 * t3 - 0.5 * t2)

    wi = cr_weights(ti)
    wj = cr_weights(tj)
    out = xp.zeros_like(f)
    for a in range(4):
        ia = (i0 + (a - 1)) % N
        for b in range(4):
            jb = (j0 + (b - 1)) % N
            out = out + wi[a] * wj[b] * f[ia, jb]
    return out


def apply_orbital_advection(eta, dx: float,
                            k_split: float | None = None,
                            preserve_short_variance: bool = True
                            ) -> tuple[np.ndarray, dict]:
    """Warp the short-wave field by the long-wave horizontal orbital
    displacement: eta'(r) = eta_L(r) + eta_S(r - D_L(r)) with the
    deep-water parcel displacement D_L = ifft2(i K_vec/|K| Z_L).

    Short waves strain-converge toward the long-wave crests (div D < 0
    there), and because D oscillates with the propagating long waves,
    per-frame application Doppler-broadens the free dispersion ridge --
    the smearing seen in measured k-f spectra that delta-thin linear
    synthesis lacks."""
    xp = xp_of(eta)
    eta = xp.asarray(eta, dtype=float)
    N = eta.shape[0]
    if k_split is None:
        k_split = 2.0 * np.pi / 0.5

    kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
    KX, KY = kx[None, :], kx[:, None]
    K = xp.hypot(KX, KY)
    Z = xp.fft.fft2(eta)
    long_mask = (K > 0) & (K <= k_split)
    Z_L = xp.where(long_mask, Z, xp.zeros_like(Z))
    eta_L = xp.fft.ifft2(Z_L).real
    eta_S = xp.fft.ifft2(xp.where(K > k_split, Z, xp.zeros_like(Z))).real

    K_safe = xp.where(K > 0, K, xp.ones_like(K))
    D_x = xp.fft.ifft2(1j * (KX / K_safe) * Z_L).real
    D_y = xp.fft.ifft2(1j * (KY / K_safe) * Z_L).real

    idx = xp.arange(N, dtype=float)
    I, J = xp.meshgrid(idx, idx, indexing="ij")
    eta_S_w = _periodic_warp(eta_S, I - D_y / dx, J - D_x / dx, xp)
    if preserve_short_variance and float(xp.std(eta_S_w)) > 0:
        # spline resampling mildly low-passes near-Nyquist content
        eta_S_w = eta_S_w * (xp.std(eta_S) / xp.std(eta_S_w))

    info = dict(k_split=float(k_split),
                D_rms=float(xp.sqrt(xp.mean(D_x**2 + D_y**2))))
    return eta_L + eta_S_w, info


# ---------------------------------------------------------------------------
# Long-wave hydrodynamic modulation (two-scale binding)
# ---------------------------------------------------------------------------

def long_wave_modulation(eta, dx: float,
                         k_split: float | None = None,
                         mtf: float = 6.5,
                         mtf_phase_deg: float = 30.0,
                         wind_dir_rad: float = 0.0,
                         preserve_short_variance: bool = True
                         ) -> tuple[np.ndarray, dict]:
    """Bind the short-wave field to the long waves by hydrodynamic
    modulation, the first stage of the binding hierarchy
    long waves -> short-gravity carriers -> FM98 parasitic capillaries.

    The field is split at k_split; the short part (k > k_split) is
    amplitude-modulated phase-locked to the long-wave profile with the
    standard modulation-transfer-function form for the energy,

        E'(r)/E = 1 + mtf * eps_L(r) * cos(phi_L(r) - theta),

    where (A_L, phi_L) come from the long-wave analytic signal along the
    wind, eps_L = A_L k_L is the local long-wave steepness, mtf ~ 4-12 is
    the empirical MTF magnitude, and theta = mtf_phase_deg sets where the
    enhancement peaks: 0 on the long-wave crest, +90 deg mid forward
    (downwind) face, -90 deg rear face.

    Because the modulation is recomputed per frame from the instantaneous
    long-wave phase, modulated short waves (and their bound capillaries)
    travel with the long wave: in (k, omega) space the modulation
    sidebands sit at the long-wave phase speed, the bound-wave signature
    seen in inverse-phase-speed spectra Q(nu, theta).

    Returns (eta_modulated, info).
    """
    xp = xp_of(eta)
    eta = xp.asarray(eta, dtype=float)
    N = eta.shape[0]
    L = N * dx
    if k_split is None:
        k_split = 2.0 * np.pi / 0.5
    if 2.0 * np.pi / k_split > L / 3.0:
        raise ValueError(
            f"domain L = {L:.2f} m too small for a long-wave band above "
            f"{2 * np.pi / k_split:.2f} m wavelength")

    kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
    K = xp.hypot(kx[None, :], kx[:, None])
    Z = xp.fft.fft2(eta)
    long_mask = (K > 0) & (K <= k_split)
    eta_L = xp.fft.ifft2(xp.where(long_mask, Z, xp.zeros_like(Z))).real
    eta_S = xp.fft.ifft2(xp.where(K > k_split, Z, xp.zeros_like(Z))).real

    Z_a = _analytic_signal_along_wind(eta_L, kx, kx, wind_dir_rad)
    A_L = xp.abs(Z_a)
    phi_L = xp.angle(Z_a)

    w = xp.where(long_mask, Z.real**2 + Z.imag**2,
                 xp.zeros_like(K))
    w_sum = float(xp.sum(w))
    k_L = float(xp.sum(K * w)) / w_sum if w_sum > 0 else k_split / 2.0

    eps_L = A_L * k_L
    theta = np.deg2rad(mtf_phase_deg)
    M_E = xp.clip(1.0 + mtf * eps_L * xp.cos(phi_L - theta), 0.05, None)
    eta_S_mod = xp.sqrt(M_E) * eta_S
    if preserve_short_variance and float(xp.std(eta_S_mod)) > 0:
        eta_S_mod = eta_S_mod * (xp.std(eta_S) / xp.std(eta_S_mod))

    info = dict(k_split=float(k_split), k_L=k_L,
                eps_L_rms=float(xp.sqrt(xp.mean(eps_L**2))),
                M_E_min=float(xp.min(M_E)), M_E_max=float(xp.max(M_E)),
                mtf=mtf, mtf_phase_deg=mtf_phase_deg)
    return eta_L + eta_S_mod, info


# ---------------------------------------------------------------------------
# FM98 envelope augmentation of a single frame
# ---------------------------------------------------------------------------

def augment_fm98(eta, dx: float, table: FM98Table,
                 wind_dir_rad: float = 0.0,
                 carrier_band: tuple[float, float] | None = None,
                 ak_clip: tuple[float, float] = (1e-3, 0.42),
                 k_rep: float | None = None,
                 compensation: str = "bins",
                 carrier_k: str = "local",
                 f_nyq: float | None = None
                 ) -> tuple[np.ndarray, np.ndarray, dict]:
    """Bound-harmonic field for one elevation frame.

    f_nyq [Hz]: optional temporal Nyquist cap for time-evolving
    records — harmonic m of a carrier k oscillates at m f(k) (it rides
    the carrier), so harmonics with m f(k_loc) > f_nyq alias when the
    field is sampled in time; pass the record's FS/2 to exclude them.

    carrier_k:
        "local"          : per-pixel carrier wavenumber from the
                           analytic-signal phase gradient, with the FM98
                           coefficients bilinearly interpolated at the
                           local (k_loc, ak); the steepness map
                           ak = A k_loc and the harmonic Nyquist
                           truncation also follow the local carrier, so
                           each wave group carries the capillary train
                           of its own scale and direction
        "representative" : single energy-weighted band wavenumber for
                           the whole frame; also selected implicitly by
                           passing k_rep

    Harmonics are added for m = 2 .. M_keep but truncated per-m at the
    grid Nyquist (per-pixel m k_loc <= k_nyq in local mode, the global
    band-top rule m k_hi <= k_nyq in representative mode), so deep
    tables (M_keep ~ 20) safely degrade on coarser grids.

    compensation:
        "bins"    : rescale the linear spectrum bin-by-bin so the total
                    power spectrum stays on the Elfouhaily target while
                    the bound-harmonic phase coherence is gained
        "carrier" : uniform variance-conserving shrink of the carrier
                    band
        "none"    : add harmonics on top (small spectral overshoot)

    Returns (eta_high, eta_comp, info); the augmented surface is
    eta + eta_comp + eta_high.
    """
    xp = xp_of(eta)
    eta = xp.asarray(eta, dtype=float)
    N = eta.shape[0]
    if eta.shape[0] != eta.shape[1]:
        raise ValueError("eta must be square")
    k_nyq = np.pi / dx
    if carrier_band is None:
        carrier_band = default_carrier_band(k_nyq)
    k_lo, k_hi = carrier_band
    if int(np.floor(k_nyq / k_hi)) < 2:
        raise ValueError(
            f"grid too coarse: second harmonic of carrier k={k_hi:.0f} "
            f"rad/m exceeds Nyquist {k_nyq:.0f} rad/m")

    kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
    ky = kx
    Z = xp.fft.fft2(eta)

    eta_band, Z_band = _bandpass(Z, kx, ky, k_lo, k_hi)
    wx, wy = np.cos(wind_dir_rad), np.sin(wind_dir_rad)
    k_dot_w = kx[None, :] * wx + ky[:, None] * wy
    F_a = Z_band * xp.where(k_dot_w > 0, 2.0,
                            xp.where(k_dot_w < 0, 0.0, 1.0))
    Z_a = xp.fft.ifft2(F_a)
    A = xp.abs(Z_a)
    phi_loc = xp.angle(Z_a)

    # Energy-weighted carrier wavenumber (representative mode and info)
    if k_rep is None:
        w = Z_band.real**2 + Z_band.imag**2
        w_sum = float(xp.sum(w))
        if w_sum > 0:
            K = xp.hypot(kx[None, :], ky[:, None])
            k_rep_eff = float(xp.sum(K * w)) / w_sum
        else:
            k_rep_eff = float(np.sqrt(k_lo * k_hi))
    else:
        k_rep_eff = float(k_rep)
        carrier_k = "representative"
    k_rep_eff = float(np.clip(k_rep_eff, table.k_grid[0],
                              table.k_grid[-1]))
    del Z_band

    if carrier_k == "local":
        # Local carrier wavevector grad(phi) = Im(conj(Z_a) grad Z_a)
        # / |Z_a|^2; where the envelope vanishes k_loc is arbitrary but
        # ak ~ 0 there, so no harmonics are placed.
        dZx = xp.fft.ifft2(1j * kx[None, :] * F_a)
        dZy = xp.fft.ifft2(1j * ky[:, None] * F_a)
        A2 = xp.maximum(A**2, 1e-30)
        k_loc = xp.hypot((xp.conj(Z_a) * dZx).imag / A2,
                         (xp.conj(Z_a) * dZy).imag / A2)
        del dZx, dZy
        k_loc = xp.clip(k_loc, k_lo, k_hi)
    elif carrier_k == "representative":
        k_loc = xp.full((N, N), k_rep_eff)
    else:
        raise ValueError(f"unknown carrier_k mode: {carrier_k!r}")
    del Z_a, F_a

    ak_map = xp.clip(A * k_loc, ak_clip[0], ak_clip[1])
    kg = xp.asarray(table.k_grid, dtype=float)
    ag = xp.asarray(table.ak_grid, dtype=float)
    a_coeffs = xp.asarray(table.a_coeffs)
    n_k = table.k_grid.size
    n_ak = table.ak_grid.size

    # Per-pixel bilinear table lookup in (log k, ak): complex
    # coefficients (gauge-fixed, so complex interpolation is
    # well-defined), magnitudes rescaled so |a_1| = ak as in
    # FM98Table.interp
    log_kg = xp.log(kg)
    kq = xp.clip(k_loc.ravel(), float(table.k_grid[0]),
                 float(table.k_grid[-1]))
    akq = xp.clip(ak_map.ravel(), float(table.ak_grid[0]),
                  float(table.ak_grid[-1]))
    # akq only indexes the table; harmonic amplitudes scale with the
    # true envelope steepness so they vanish where the envelope does
    ak_amp = ak_map.ravel()
    iq = xp.clip(xp.searchsorted(log_kg, xp.log(kq)) - 1, 0, n_k - 2)
    jq = xp.clip(xp.searchsorted(ag, akq) - 1, 0, n_ak - 2)
    tk = (xp.log(kq) - log_kg[iq]) / (log_kg[iq + 1] - log_kg[iq])
    ta = (akq - ag[jq]) / (ag[jq + 1] - ag[jq])
    w00 = (1 - tk) * (1 - ta)
    w10 = tk * (1 - ta)
    w01 = (1 - tk) * ta
    w11 = tk * ta

    def _coeff(m_idx):
        Am = a_coeffs[..., m_idx]
        return (w00 * Am[iq, jq] + w10 * Am[iq + 1, jq]
                + w01 * Am[iq, jq + 1] + w11 * Am[iq + 1, jq + 1])

    a1 = _coeff(0)
    a1_abs = xp.maximum(xp.abs(a1), 1e-30)
    phi_1 = xp.angle(a1)
    phi_q = phi_loc.ravel()

    # Temporal Nyquist: harmonic m rides its carrier at m f(k_loc) Hz
    if f_nyq is not None:
        from .spectrum import angular_frequency
        f_loc_q = angular_frequency(kq) / (2.0 * np.pi)
    else:
        f_loc_q = None

    # Harmonic Nyquist truncation: per-pixel m k_loc <= k_nyq in local
    # mode; the original global band-top rule in representative mode
    if carrier_k == "local":
        m_top = min(table.M_keep, int(np.floor(k_nyq / k_lo)))
    else:
        m_top = min(table.M_keep, int(np.floor(k_nyq / k_hi)))
    m_max = m_top

    acc = xp.zeros(N * N)
    r2 = xp.zeros(N * N)
    for m_idx in range(1, table.M_keep):
        m = m_idx + 1
        am_c = _coeff(m_idx)
        am_abs = xp.abs(am_c)
        if float(xp.max(am_abs)) < 1e-14:
            continue
        ratio = am_abs / a1_abs
        r2 += ratio**2                  # full-table shrink, as before
        if m_idx >= m_top:
            continue
        amp = ak_amp * ratio            # |a_m| with |a_1| -> ak
        if carrier_k == "local":
            amp = xp.where(m * kq <= k_nyq, amp, xp.zeros_like(amp))
        if f_loc_q is not None:
            amp = xp.where(m * f_loc_q <= f_nyq, amp, xp.zeros_like(amp))
        acc += amp * xp.cos(m * phi_q + m * phi_1 - xp.angle(am_c))

    # Variance-conserving shrink s = 1/sqrt(1 + sum |a_m/a_1|^2)
    shrink = (1.0 / xp.sqrt(1.0 + r2)).reshape(N, N)
    if compensation == "carrier":
        acc = acc * shrink.ravel()
    eta_high = (acc / kq).reshape(N, N)

    if compensation == "carrier":
        eta_comp = (float(xp.mean(shrink)) - 1.0) * eta_band
    elif compensation == "bins":
        # The Elfouhaily spectrum is an empirical TOTAL (free + bound,
        # tuned to Cox-Munk statistics), so the augmentation must
        # redistribute variance, never add it.  Bin-wise: cap the
        # bound-harmonic power at the available empirical power (bin
        # phases kept, so the phase-locking survives), then remove
        # exactly the capped power from the random-phase field.  Total
        # bin power -- and hence MSS -- is preserved by construction.
        F_high = xp.fft.fft2(eta_high)
        P_high = F_high.real**2 + F_high.imag**2
        P_lin = Z.real**2 + Z.imag**2
        # Zero-safe denominators: empty bins (P == 0) divide by 1 and the
        # numerator there is 0, so cap/factor are unaffected.  A subnormal
        # epsilon floor would underflow to 0 in float32 and reintroduce
        # 0/0, so mask explicitly instead.
        safe_high = xp.where(P_high > 0, P_high, xp.ones_like(P_high))
        cap = xp.sqrt(xp.minimum(1.0, P_lin / safe_high))
        F_high = F_high * cap
        eta_high = xp.fft.ifft2(F_high).real
        P_high_capped = P_high * cap**2
        del F_high, P_high, safe_high
        safe_lin = xp.where(P_lin > 0, P_lin, xp.ones_like(P_lin))
        factor = xp.sqrt(xp.clip(1.0 - P_high_capped / safe_lin, 0.0, 1.0))
        eta_comp = xp.fft.ifft2((factor - 1.0) * Z).real
        del P_high_capped, P_lin, factor, safe_lin
    elif compensation == "none":
        eta_comp = xp.zeros((N, N))
    else:
        raise ValueError(f"unknown compensation mode: {compensation!r}")

    # Envelope-weighted local-k percentiles (low-envelope pixels carry
    # no harmonics, so weight by A^2)
    A2w = (A**2).ravel()
    order = xp.argsort(k_loc.ravel())
    cw = xp.cumsum(A2w[order])
    cw = cw / max(float(cw[-1]), 1e-300)
    k_loc_sorted = k_loc.ravel()[order]
    k_loc_p50 = float(xp.interp(0.5, cw, k_loc_sorted))
    k_loc_p90 = float(xp.interp(0.9, cw, k_loc_sorted))

    info = dict(k_rep=k_rep_eff, carrier_k=carrier_k,
                carrier_band=carrier_band, m_max=m_max,
                k_loc_p50=k_loc_p50, k_loc_p90=k_loc_p90,
                A_rms=float(xp.sqrt(xp.mean(A**2))),
                ak_p50=float(xp.percentile(ak_map, 50)),
                ak_p99=float(xp.percentile(ak_map, 99)),
                shrink_mean=float(xp.mean(shrink)),
                var_band=float(xp.var(eta_band)),
                var_high=float(xp.var(eta_high)))
    return eta_high, eta_comp, info


def default_carrier_band(k_nyq: float, m_min: int = 4,
                         lam_longest: float = 0.08) -> tuple[float, float]:
    """Short-gravity carrier band: lam_longest down to the shortest
    wavelength whose m_min-th harmonic stays below Nyquist (>= 2 cm).
    Higher harmonics are truncated per-m at synthesis time."""
    lam_lo = max(0.02, m_min * 2.0 * np.pi / k_nyq)
    if lam_lo >= lam_longest:
        raise ValueError(
            f"grid too coarse for FM98 augmentation: need "
            f"dx <= {lam_longest * 1e3 / (2 * m_min):.1f} mm")
    return (2.0 * np.pi / lam_longest, 2.0 * np.pi / lam_lo)


def default_fm98_table(k_nyq: float, m_keep: int = 12,
                       M_solve: int = 40, n_steps: int = 10,
                       verbose: bool = False) -> FM98Table:
    """Table spanning the default carrier band on an 8 x 9 (k, ak) grid.
    m_keep ~ 12-28 is needed to render the parasitic capillary train
    itself (the ripples live at harmonic numbers ~ 10+ of a 5-8 cm
    carrier); small m only sharpens crests.  The top ak nodes target
    the Melville & Fedorov (2015) steepness range; the solver realizes
    what the branch supports (~0.32 at 5 cm to ~0.36 at 10 cm)."""
    k_lo, k_hi = default_carrier_band(k_nyq)
    k_grid = np.geomspace(max(0.8 * k_lo, 1.0), 1.2 * k_hi, 8)
    ak_grid = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                        0.34, 0.38])
    return build_fm98_table(k_grid, ak_grid, M_keep=m_keep,
                            M_solve=M_solve, n_steps=n_steps,
                            verbose=verbose)


# ---------------------------------------------------------------------------
# Top-level generator
# ---------------------------------------------------------------------------

def generate_hybrid_surface(L: float, N: int, U10: float,
                            wind_dir_rad: float = 0.0,
                            times: np.ndarray | None = None,
                            table: FM98Table | None = None,
                            carrier_band: tuple[float, float] | None = None,
                            compensation: str = "bins",
                            carrier_k: str = "local",
                            f_nyq: float | None = None,
                            long_wave_mtf: float | None = None,
                            mtf_phase_deg: float = 30.0,
                            k_split: float | None = None,
                            orbital_advection: bool = False,
                            frame_callback=None,
                            rng=None,
                            verbose: bool = False,
                            **surface_kwargs) -> HybridSeaSurface:
    """Capillary-resolving hybrid surface with phase-locked FM98 bound
    harmonics; supports the same time evolution as generate_sea_surface
    (each frame is augmented independently, so the harmonics ride their
    carriers).

    carrier_k="local" (default) places each wave group's capillary train
    at its own scale and direction from the local carrier wavevector;
    "representative" uses a single band-averaged carrier for the whole
    frame (see augment_fm98).

    long_wave_mtf (e.g. 4-12) enables the binding hierarchy: short waves
    are first hydrodynamically modulated by the long-wave field
    (long_wave_modulation), so the FM98 envelope -- and hence the
    parasitic capillaries -- concentrates on the long-wave faces selected
    by mtf_phase_deg.  Extra keyword arguments (including backend/
    device/dtype for the torch path) pass through to
    generate_sea_surface."""
    dx = L / N
    k_nyq = np.pi / dx
    if table is None:
        if verbose:
            print("building FM98 table (one-time cost)...")
        table = default_fm98_table(k_nyq)
    if carrier_band is None:
        carrier_band = default_carrier_band(k_nyq)

    def _augment_frame(eta_frame):
        xp = xp_of(eta_frame)
        kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
        if orbital_advection:
            eta_frame, _ = apply_orbital_advection(eta_frame, dx,
                                                   k_split=k_split)
        if long_wave_mtf is not None:
            eta_frame, mod_info = long_wave_modulation(
                eta_frame, dx, k_split=k_split, mtf=long_wave_mtf,
                mtf_phase_deg=mtf_phase_deg, wind_dir_rad=wind_dir_rad)
        else:
            mod_info = {}
        eta_high, eta_shrink, env = augment_fm98(
            eta_frame, dx, table, wind_dir_rad=wind_dir_rad,
            carrier_band=carrier_band, compensation=compensation,
            carrier_k=carrier_k, f_nyq=f_nyq)
        env.update(mod_info)
        total = eta_frame + eta_shrink + eta_high
        Zt = xp.fft.fft2(total)
        sx = xp.fft.ifft2(1j * kx[None, :] * Zt).real
        sy = xp.fft.ifft2(1j * kx[:, None] * Zt).real
        return total, eta_frame, eta_high, eta_shrink, sx, sy, env

    if frame_callback is not None and times is not None:
        # streaming mode: each linear frame is augmented in flight and
        # handed to frame_callback(it, t, eta_total, slope_x, slope_y);
        # the returned object carries frame 0 only.
        first = {}
        env = None

        def lin_cb(it, t, e, gx, gy):
            nonlocal env
            out = _augment_frame(e)
            if not first:
                first["fields"] = out
                first["lin"] = e
                env = out[6]
            frame_callback(it, t, out[0], out[4], out[5])

        lin = generate_sea_surface(L, N, U10, wind_dir_rad=wind_dir_rad,
                                   times=times, compute_slopes=False,
                                   frame_callback=lin_cb, rng=rng,
                                   **surface_kwargs)
        (eta, eta_lin, eta_high, eta_shrink, sx, sy,
         _) = first["fields"]
        xp = xp_of(eta)
        info = dict(lin.info)
        info.update(env)
        info.update(var_lin=float(xp.var(first["lin"])),
                    var_high=float(xp.var(eta_high)),
                    Hs_realized=4.0 * float(xp.std(eta)),
                    k_nyq=k_nyq, streamed=True)
        return HybridSeaSurface(x=lin.x, y=lin.y, eta=eta,
                                eta_lin=eta_lin, eta_high=eta_high,
                                eta_shrink=eta_shrink, slope_x=sx,
                                slope_y=sy, t=np.atleast_1d(times),
                                info=info)

    lin = generate_sea_surface(L, N, U10, wind_dir_rad=wind_dir_rad,
                               times=times, compute_slopes=False,
                               rng=rng, **surface_kwargs)
    xp = xp_of(lin.eta)

    # eta_lin holds the linear field after long-wave modulation, so the
    # decomposition eta = eta_lin + eta_shrink + eta_high always closes.
    if lin.eta.ndim == 2:
        (eta, eta_lin, eta_high, eta_shrink, sx, sy,
         env) = _augment_frame(lin.eta)
    else:
        T = lin.eta.shape[2]
        eta = xp.empty_like(lin.eta)
        eta_lin = xp.empty_like(lin.eta)
        eta_high = xp.empty_like(lin.eta)
        eta_shrink = xp.empty_like(lin.eta)
        sx = xp.empty_like(lin.eta)
        sy = xp.empty_like(lin.eta)
        env = None
        for it in range(T):
            (eta[:, :, it], eta_lin[:, :, it], eta_high[:, :, it],
             eta_shrink[:, :, it], sx[:, :, it], sy[:, :, it],
             env_it) = _augment_frame(lin.eta[:, :, it])
            env = env or env_it

    eta0 = eta if eta.ndim == 2 else eta[:, :, 0]
    info = dict(lin.info)
    info.update(env)
    info.update(var_lin=float(xp.var(lin.eta)),
                var_high=float(xp.var(eta_high)),
                Hs_realized=4.0 * float(xp.std(eta0)),
                k_nyq=k_nyq)
    if verbose:
        print(f"hybrid: k_rep={env['k_rep']:.0f} rad/m, "
              f"ak p50/p99 = {env['ak_p50']:.3f}/{env['ak_p99']:.3f}, "
              f"var(high)/var(lin) = {info['var_high'] / info['var_lin']:.2e}")

    return HybridSeaSurface(x=lin.x, y=lin.y, eta=eta, eta_lin=eta_lin,
                            eta_high=eta_high, eta_shrink=eta_shrink,
                            slope_x=sx, slope_y=sy, t=lin.t, info=info)
