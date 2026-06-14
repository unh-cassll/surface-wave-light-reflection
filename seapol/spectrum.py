"""
Elfouhaily et al. (1997) unified directional wave spectrum.

References:
    Elfouhaily, T., Chapron, B., Katsaros, K., Vandemark, D. (1997).
    "A unified directional spectrum for long and short wind-driven waves."
    J. Geophys. Res. 102(C7), 15781-15796.

    Mobley, C. D. (2016). "Modeling Sea Surfaces: A Tutorial on Fourier
    Transform Methods." Sequoia Scientific.

Conventions:
    S(k)        omnidirectional elevation spectrum [m^3/rad];
                <eta^2> = integral S(k) dk
    D(k, phi)   angular spreading; integral over phi in [-pi, pi] = 1
    Psi(kx, ky) = S(k) D(k, phi) / k  [m^4/rad^2];
                <eta^2> = integral Psi dkx dky over the full plane
"""

from __future__ import annotations

import numpy as np

from .backend import xp_of

GRAVITY = 9.81             # m/s^2
SURFACE_TENSION = 0.072    # N/m
WATER_DENSITY = 1030.0     # kg/m^3
KM = 370.0                 # rad/m, wavenumber of minimum phase speed
CM = 0.23                  # m/s, minimum phase speed
X0 = 2.2e4                 # dimensionless fetch constant


def drag_coefficient(U10: float, model: str = "logistic"):
    """10-m neutral drag coefficient Cd and friction velocity u* [m/s].

    'logistic'   : logistic fit to Edson et al. (2013) for U10 < 25 m/s and
                   Curcic & Haus (2020) above (port of logistic_fit_drag.m).
    'elfouhaily' : Cd = 1e-3 (0.8 + 0.065 U10), as used in the ECKV paper.
    """
    U10 = np.asarray(U10, dtype=float)
    if model == "logistic":
        Cd = (3.2e-3 - 0.8e-3) / (1.0 + np.exp(-0.2 * (U10 - 17.0))) + 0.8e-3
    elif model == "elfouhaily":
        Cd = 1.0e-3 * (0.8 + 0.065 * U10)
    else:
        raise ValueError(f"unknown drag model: {model!r}")
    u_star = np.sqrt(Cd) * U10
    return Cd, u_star


def inverse_wave_age(U10: float, fetch_m: float = 1.0e6) -> float:
    """Fetch-dependent inverse wave age Omega_c (ECKV Eq. 37).
    Omega_c -> 0.84 for fully developed seas (large fetch)."""
    X = GRAVITY * fetch_m / U10**2
    return float(0.84 * np.tanh((X / X0) ** 0.4) ** (-0.75))


def angular_frequency(k):
    """Gravity-capillary dispersion omega(k) [rad/s], deep water."""
    xp = xp_of(k)
    k = xp.asarray(k, dtype=float)
    return xp.sqrt(GRAVITY * k + (SURFACE_TENSION / WATER_DENSITY) * k**3)


def phase_speed(k):
    """Gravity-capillary phase speed c(k) [m/s], deep water."""
    xp = xp_of(k)
    k = xp.asarray(k, dtype=float)
    return xp.sqrt(GRAVITY / k + (SURFACE_TENSION / WATER_DENSITY) * k)


def elfouhaily_omni(k: np.ndarray, U10: float,
                    fetch_m: float = 1.0e6,
                    drag_model: str = "logistic",
                    lowk_capillary_taper: bool = True) -> np.ndarray:
    """Omnidirectional ECKV elevation spectrum S(k) [m^3/rad].

    lowk_capillary_taper multiplies the capillary branch B_h by the
    Pierson-Moskowitz cutoff exp(-1.25 (kp/k)^2).  Taken literally, the
    published B_h ~ sqrt(k) leaves S_h ~ k^(-5/2) at low k, which makes the
    elevation variance integral diverge; the taper removes that artifact
    without touching the equilibrium/capillary range (k >> kp).  Equivalent
    in intent to the low-k trim in the original Elfouhaily_omni.m.
    """
    xp = xp_of(k)
    k = xp.asarray(k, dtype=float)
    _, u_star = drag_coefficient(U10, drag_model)
    u_star = float(u_star)

    Omega_c = inverse_wave_age(U10, fetch_m)
    kp = GRAVITY * (Omega_c / U10) ** 2
    cp = np.sqrt(GRAVITY / kp)
    c = phase_speed(xp.maximum(k, 1e-300))

    # Long-wave (gravity) curvature branch B_l
    Lpm = xp.exp(-1.25 * (kp / k) ** 2)
    gamma_j = 1.7 if Omega_c <= 1.0 else 1.7 + 6.0 * np.log10(Omega_c)
    sigma_j = 0.08 * (1.0 + 4.0 * Omega_c ** (-3))
    Gamma = xp.exp(-0.5 * (xp.sqrt(k / kp) - 1.0) ** 2 / sigma_j**2)
    Jp = gamma_j**Gamma
    Fp = Lpm * Jp * xp.exp(-(Omega_c / np.sqrt(10.0)) * (xp.sqrt(k / kp) - 1.0))
    alpha_p = 6.0e-3 * np.sqrt(Omega_c)
    B_l = 0.5 * alpha_p * cp / c * Fp

    # Short-wave (capillary) curvature branch B_h
    if u_star <= CM:
        alpha_m = 1.0e-2 * (1.0 + np.log(max(u_star / CM, 1e-12)))
    else:
        alpha_m = 1.0e-2 * (1.0 + 3.0 * np.log(u_star / CM))
    alpha_m = max(alpha_m, 0.0)
    Fm = xp.exp(-0.25 * (k / KM - 1.0) ** 2)
    B_h = 0.5 * alpha_m * CM / c * Fm
    if lowk_capillary_taper:
        B_h = B_h * Lpm

    with np.errstate(divide="ignore", invalid="ignore"):
        S = (B_l + B_h) / k**3
    S = xp.where(k > 0, S, 0.0)
    return xp.maximum(S, 0.0)


def elfouhaily_delta(k: np.ndarray, U10: float,
                     fetch_m: float = 1.0e6,
                     drag_model: str = "logistic") -> np.ndarray:
    """Spreading amplitude Delta(k) of the ECKV angular function (Eq. 57)."""
    xp = xp_of(k)
    k = xp.asarray(k, dtype=float)
    _, u_star = drag_coefficient(U10, drag_model)
    Omega_c = inverse_wave_age(U10, fetch_m)
    kp = GRAVITY * (Omega_c / U10) ** 2
    cp = np.sqrt(GRAVITY / kp)
    c = phase_speed(xp.maximum(k, 1e-300))
    a0 = np.log(2.0) / 4.0
    ap = 4.0
    am = 0.13 * float(u_star) / CM
    return xp.tanh(a0 + ap * (c / cp) ** 2.5 + am * (CM / c) ** 2.5)


def directional_spread(k: np.ndarray, phi: np.ndarray, U10: float,
                       fetch_m: float = 1.0e6,
                       drag_model: str = "logistic",
                       one_sided: bool = False) -> np.ndarray:
    """Angular spreading D(k, phi); integral over phi in [-pi, pi] = 1.

    one_sided doubles the downwind half (|phi| < pi/2) and zeroes the upwind
    half.  Because cos(2 phi) has period pi, this preserves the total
    variance and all even slope moments while making the synthesized field
    propagate downwind under time evolution.
    """
    xp = xp_of(k, phi)
    Delta = elfouhaily_delta(k, U10, fetch_m, drag_model)
    D = (1.0 + Delta * xp.cos(2.0 * phi)) / (2.0 * np.pi)
    if one_sided:
        c = xp.cos(phi)
        # half weight exactly crosswind so discrete grids stay normalized
        w = xp.where(xp.abs(c) < 1e-12, xp.ones_like(c),
                     2.0 * xp.astype(c > 0.0, float))
        D = D * w
    return D


def directional_spectrum(kx: np.ndarray, ky: np.ndarray, U10: float,
                         wind_dir_rad: float = 0.0,
                         fetch_m: float = 1.0e6,
                         drag_model: str = "logistic",
                         one_sided: bool = False,
                         lowk_capillary_taper: bool = True) -> np.ndarray:
    """Cartesian directional spectrum Psi(kx, ky) [m^4/rad^2].
    <eta^2> = sum Psi dkx dky over the full (kx, ky) plane."""
    xp = xp_of(kx, ky)
    K = xp.hypot(kx, ky)
    phi = xp.arctan2(ky, kx) - wind_dir_rad
    good = K > 0
    S = xp.zeros_like(K)
    S[good] = elfouhaily_omni(K[good], U10, fetch_m, drag_model,
                              lowk_capillary_taper)
    D = directional_spread(xp.maximum(K, 1e-12), phi, U10,
                           fetch_m, drag_model, one_sided)
    return xp.where(good, S / xp.maximum(K, 1e-12) * D,
                    xp.zeros_like(K))


def cutoff_slope_variances(U10: float, k_cutoff: float,
                           k_max: float = 2.0e4,
                           n_k: int = 600, n_phi: int = 180,
                           fetch_m: float = 1.0e6,
                           drag_model: str = "logistic",
                           lowk_capillary_taper: bool = True
                           ) -> tuple[float, float]:
    """Sub-grid slope variances (along-wind, cross-wind) from k_cutoff up.

    sigma_a2 = integral k^2 cos^2(phi) Psi k dk dphi over k > k_cutoff;
    sigma_c2 with sin^2(phi).  Used as the unresolved (subpixel) facet
    slope statistics in the reflection models.
    """
    if k_cutoff >= k_max:
        return 0.0, 0.0
    k = np.geomspace(k_cutoff, k_max, n_k)
    phi = np.linspace(-np.pi, np.pi, n_phi, endpoint=False)
    dphi = 2.0 * np.pi / n_phi
    S = elfouhaily_omni(k, U10, fetch_m, drag_model, lowk_capillary_taper)
    Delta = elfouhaily_delta(k, U10, fetch_m, drag_model)
    # Psi k = S(k) D(phi) ; slope integrand k^2 * that
    base = k**2 * S
    cos2_int = np.sum((1.0 + Delta[:, None] * np.cos(2.0 * phi)[None, :])
                      / (2.0 * np.pi) * np.cos(phi)[None, :] ** 2,
                      axis=1) * dphi
    sin2_int = np.sum((1.0 + Delta[:, None] * np.cos(2.0 * phi)[None, :])
                      / (2.0 * np.pi) * np.sin(phi)[None, :] ** 2,
                      axis=1) * dphi
    sa2 = np.trapezoid(base * cos2_int, k)
    sc2 = np.trapezoid(base * sin2_int, k)
    return float(sa2), float(sc2)


def total_mean_square_slope(U10: float,
                            k_min: float = 1.0e-3, k_max: float = 2.0e4,
                            n_k: int = 2000,
                            fetch_m: float = 1.0e6,
                            drag_model: str = "logistic",
                            lowk_capillary_taper: bool = True) -> float:
    """Total mean-square slope: integral k^2 S(k) dk (for Cox-Munk checks)."""
    k = np.geomspace(k_min, k_max, n_k)
    S = elfouhaily_omni(k, U10, fetch_m, drag_model, lowk_capillary_taper)
    return float(np.trapezoid(k**2 * S, k))
