"""
Water-body optics: first-order water-leaving radiance, spectral inherent
optical properties, and phase Mueller matrices for in-water polarized
Monte Carlo propagation.

First-order model (WaterBody): the water body is represented by a
band-integrated irradiance reflectance R_w just below the surface:
upwelling radiance is taken isotropic, L_u = R_w E_d / pi, refracted
through the facet with the n^2 radiance law.  This is the lowest-order
ocean-color term; seapol.scattering replaces it with a directional,
polarized upwelling radiance field computed by Monte Carlo.

Spectral model (WaterType -> WaterColumn): bulk IOPs assembled from
nominal open-literature components --
    * pure-water absorption a_w(lambda), Pope & Fry (1997);
    * pure-seawater Rayleigh scattering b_w(lambda) ~ lambda^-4.32
      (Morel 1974) with depolarization 0.039;
    * phytoplankton absorption a_ph = 0.06 Chl^0.65 A*(lambda)
      (Bricaud-type Case 1) with a CDOM tail tied to Chl;
    * explicit CDOM a_g(440) exp(-0.014 (lambda - 440)) and sediment
      scattering for Case 2;
    * particulate scattering b_p(550) = 0.416 Chl^0.766
      (Loisel & Morel 1998), Henyey-Greenstein phase with g = 0.924
      (Petzold-like) and a depolarized Rayleigh-shaped polarization
      ratio (peak single-scatter DoLP ~ 0.65, Voss & Fry 1984);
    * an optional near-surface bubble layer: exponentially decaying
      scattering b_bub(z) = b_bub0 exp(z / z_e), strongly forward
      (g ~ 0.85) and nearly unpolarizing.
Refractive index n(lambda) from Quan & Fry (1995).

These are nominal magnitudes for scene realism, not a bio-optical
retrieval model.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .backend import NUMPY_XP, xp_of
from .polarization import (frame_rotation_angle, fresnel_mueller,
                           meridian_frame, mueller_rotation, normalize)

__all__ = ["WaterBody", "WaterOptics", "WaterColumn", "WaterType",
           "WATER_TYPES", "water_leaving_stokes",
           "rayleigh_phase_mueller", "hg_phase_mueller",
           "fournier_forand_phase", "ff_phase_mueller",
           "sample_rayleigh_scattering", "sample_hg_scattering",
           "sample_ff_scattering", "polarized_scatter_event",
           "pure_water_absorption", "pure_water_scattering",
           "water_refractive_index", "RAYLEIGH_DEPOL"]

_DEFAULT_R = {1: 0.02, 2: 0.06}

# Pure-seawater Rayleigh depolarization ratio (Morel 1974)
RAYLEIGH_DEPOL = 0.039


@dataclass
class WaterBody:
    """Band-integrated water-body description (first-order model)."""
    case: int = 1
    reflectance: float | None = None     # overrides the case default

    @property
    def R_w(self) -> float:
        if self.reflectance is not None:
            return float(self.reflectance)
        try:
            return _DEFAULT_R[self.case]
        except KeyError:
            raise ValueError(f"water case must be 1 or 2, got {self.case}")


@dataclass
class WaterOptics:
    """Bulk inherent optical properties for in-water Monte Carlo
    propagation (band-integrated visible values).

    absorption     : a [1/m]
    scattering     : b [1/m]
    depolarization : Rayleigh depolarization ratio delta

    Defaults are rough green-band clear-ocean values: pure-seawater
    absorption plus a modest particulate scattering contribution.
    """
    absorption: float = 0.10
    scattering: float = 0.05
    depolarization: float = RAYLEIGH_DEPOL

    @property
    def attenuation(self) -> float:
        """Beam attenuation c = a + b [1/m]."""
        return self.absorption + self.scattering

    @property
    def albedo(self) -> float:
        """Single-scattering albedo omega_0 = b / c."""
        c = self.attenuation
        return self.scattering / c if c > 0 else 0.0


# ---------------------------------------------------------------------------
# Spectral inherent optical properties
# ---------------------------------------------------------------------------

# Pope & Fry (1997) pure-water absorption [1/m], 10 nm grid 380-730 nm
_AW_WL_NM = np.arange(380.0, 731.0, 10.0)
_AW = np.array([
    0.01137, 0.00941, 0.00663, 0.00473, 0.00454, 0.00495, 0.00635,
    0.00922, 0.00979, 0.01060, 0.01270, 0.01500, 0.02040, 0.03250,
    0.04090, 0.04340, 0.04740, 0.05650, 0.06190, 0.06950, 0.08960,
    0.13510, 0.22240, 0.26440, 0.27550, 0.29160, 0.31080, 0.34000,
    0.41000, 0.43900, 0.46500, 0.51600, 0.62400, 0.82700, 1.23100,
    1.79900])

# Normalized phytoplankton absorption shape A*(lambda), A*(440) = 1
# (coarse Bricaud-type Case 1 shape with the 675 nm red peak)
_APH_WL_NM = np.array([380., 400., 420., 440., 460., 480., 500., 520.,
                       540., 560., 580., 600., 620., 640., 660., 675.,
                       690., 710., 730.])
_APH_REL = np.array([0.65, 0.69, 0.83, 1.00, 0.91, 0.79, 0.61, 0.45,
                     0.35, 0.28, 0.26, 0.24, 0.27, 0.29, 0.38, 0.50,
                     0.33, 0.07, 0.02])


def pure_water_absorption(wavelength_nm):
    """Pure-water absorption a_w(lambda) [1/m] (Pope & Fry 1997)."""
    wl = np.asarray(wavelength_nm, dtype=float)
    return np.interp(wl, _AW_WL_NM, _AW)


def pure_water_scattering(wavelength_nm):
    """Pure-seawater scattering b_w(lambda) = 0.0031 (450/lambda)^4.32
    [1/m] (Morel 1974)."""
    wl = np.asarray(wavelength_nm, dtype=float)
    return 0.0031 * (450.0 / wl) ** 4.32


def water_refractive_index(wavelength_nm, salinity_psu: float = 35.0,
                           temperature_c: float = 20.0):
    """Seawater refractive index n(S, T, lambda) (Quan & Fry 1995)."""
    lam = np.asarray(wavelength_nm, dtype=float)
    S, T = salinity_psu, temperature_c
    return (1.31405 + (1.779e-4 + -1.05e-6 * T + 1.6e-8 * T**2) * S
            + -2.02e-6 * T**2
            + (15.868 + 0.01155 * S + -0.00423 * T) / lam
            + -4382.0 / lam**2 + 1.1455e6 / lam**3)


@dataclass
class WaterColumn:
    """Bulk IOPs of a horizontally homogeneous, semi-infinite water
    column for the seapol.scattering Monte Carlo.  Fields are scalars
    for a single band or 1-D arrays over wavelength bands.

    absorption             : a [1/m]
    rayleigh_scattering    : b_w [1/m], molecular (Rayleigh) component
    particulate_scattering : b_p [1/m], Henyey-Greenstein component
    particulate_g          : HG asymmetry parameter (Petzold-like 0.924)
    particulate_depol      : depolarization of the particulate phase
                             matrix (0.2 -> peak DoLP ~ 0.65)
    bubble_scattering      : b_bub at the surface [1/m]; decays as
                             exp(z / bubble_efold_m) with depth
    bubble_efold_m         : bubble-layer e-folding depth [m]
    bubble_g               : bubble HG asymmetry
    depolarization         : Rayleigh depolarization (molecular)
    n_water                : refractive index (per band when spectral)
    wavelengths_nm         : band wavelengths (None for band-integrated)
    """
    absorption: float | np.ndarray = 0.10
    rayleigh_scattering: float | np.ndarray = 0.0022
    particulate_scattering: float | np.ndarray = 0.05
    particulate_g: float = 0.924
    particulate_depol: float = 0.2
    particulate_phase: str = "hg"           # "hg" or "ff" (Fournier-Forand)
    ff_n: float = 1.05                       # FF particle index (rel. water)
    ff_mu_junge: float = 3.5                 # FF Junge slope (> 3)
    bubble_scattering: float | np.ndarray = 0.0
    bubble_efold_m: float = 0.3
    bubble_g: float = 0.85
    depolarization: float = RAYLEIGH_DEPOL
    n_water: float | np.ndarray = 1.34
    wavelengths_nm: np.ndarray | None = None

    @property
    def n_bands(self) -> int:
        return 1 if self.wavelengths_nm is None else len(self.wavelengths_nm)

    def at_band(self, i: int) -> "WaterColumn":
        """Scalar-IOP column for band i."""
        def pick(v):
            return float(np.asarray(v).ravel()[i]) if np.ndim(v) > 0 \
                else float(v)
        if self.wavelengths_nm is None and i != 0:
            raise IndexError("band-integrated column has a single band")
        return replace(self, absorption=pick(self.absorption),
                       rayleigh_scattering=pick(self.rayleigh_scattering),
                       particulate_scattering=pick(
                           self.particulate_scattering),
                       bubble_scattering=pick(self.bubble_scattering),
                       n_water=pick(self.n_water),
                       wavelengths_nm=None)

    @property
    def attenuation_max(self):
        """Majorant beam attenuation (bubble layer at full strength)."""
        return (self.absorption + self.rayleigh_scattering
                + self.particulate_scattering + self.bubble_scattering)


@dataclass
class WaterType:
    """Nominal bio-optical water description generating spectral IOPs.

    chlorophyll_mg_m3 : Case 1 pigment load (drives a_ph, b_p and a
                        proportional CDOM tail)
    cdom_a440         : additional explicit CDOM absorption at 440 nm
                        [1/m] (Case 2)
    sediment_b550     : additional sediment scattering at 550 nm [1/m]
                        (Case 2, lambda^-0.5, g = 0.94 folded into the
                        particulate component)
    bubble_scattering : near-surface bubble-layer scattering at z = 0
                        [1/m] (spectrally flat)
    """
    chlorophyll_mg_m3: float = 0.1
    cdom_a440: float = 0.0
    sediment_b550: float = 0.0
    bubble_scattering: float = 0.0
    bubble_efold_m: float = 0.3
    salinity_psu: float = 35.0
    temperature_c: float = 20.0
    particulate_phase: str = "hg"            # "hg" or "ff" (Fournier-Forand)
    ff_n: float = 1.05                       # FF particle index (rel. water)
    ff_mu_junge: float = 3.5                 # FF Junge slope (> 3)

    def column(self, wavelengths_nm=None) -> WaterColumn:
        """Spectral (or 550 nm band-integrated) WaterColumn."""
        wl = np.atleast_1d(np.asarray(
            550.0 if wavelengths_nm is None else wavelengths_nm,
            dtype=float))
        chl = max(self.chlorophyll_mg_m3, 0.0)
        a_ph440 = 0.06 * chl ** 0.65 if chl > 0 else 0.0
        a_ph = a_ph440 * np.interp(wl, _APH_WL_NM, _APH_REL)
        a_g440 = 0.2 * a_ph440 + self.cdom_a440
        a_g = a_g440 * np.exp(-0.014 * (wl - 440.0))
        a = pure_water_absorption(wl) + a_ph + a_g

        b_w = pure_water_scattering(wl)
        b_p = (0.416 * chl ** 0.766 if chl > 0 else 0.0) \
            * (550.0 / wl) ** 0.5
        b_s = self.sediment_b550 * (550.0 / wl) ** 0.5
        # sediment folded into the particulate HG component with a
        # scattering-weighted asymmetry
        b_part = b_p + b_s
        g_eff = 0.924 if b_part.max() <= 0 else \
            float(((0.924 * b_p + 0.94 * b_s) / np.maximum(b_part, 1e-30)
                   ).mean())

        n = water_refractive_index(wl, self.salinity_psu,
                                   self.temperature_c)
        squeeze = wavelengths_nm is None
        def out(v):
            return float(v[0]) if squeeze else v
        return WaterColumn(
            absorption=out(a), rayleigh_scattering=out(b_w),
            particulate_scattering=out(b_part), particulate_g=g_eff,
            particulate_phase=self.particulate_phase,
            ff_n=self.ff_n, ff_mu_junge=self.ff_mu_junge,
            bubble_scattering=self.bubble_scattering,
            bubble_efold_m=self.bubble_efold_m,
            n_water=out(n),
            wavelengths_nm=None if squeeze else wl)


WATER_TYPES = {
    "clear": WaterType(chlorophyll_mg_m3=0.03),
    "case1": WaterType(chlorophyll_mg_m3=0.3),
    "productive_case1": WaterType(chlorophyll_mg_m3=3.0),
    "coastal_case2": WaterType(chlorophyll_mg_m3=1.0, cdom_a440=0.25,
                               sediment_b550=0.6),
    # turbid sediment-laden water with the Fournier-Forand particulate
    # phase function (realistic backscatter that a single HG under-states)
    "turbid_ff": WaterType(chlorophyll_mg_m3=2.0, cdom_a440=0.15,
                           sediment_b550=1.2, particulate_phase="ff",
                           ff_n=1.10, ff_mu_junge=3.6),
}


# ---------------------------------------------------------------------------
# Phase Mueller matrices and scattering-angle sampling
# ---------------------------------------------------------------------------

def rayleigh_phase_mueller(cos_theta, depol: float = RAYLEIGH_DEPOL):
    """Rayleigh scattering phase Mueller matrix with depolarization
    (Hansen & Travis 1974), shape cos_theta.shape + (4, 4), acting in
    the (p, s) scattering-plane basis.

    Normalized so the (0, 0) element averages to 1 over the sphere:
    int P[0, 0] dOmega / (4 pi) = 1.  delta = 0 recovers pure Rayleigh
    (DoLP = 1 at 90 deg); finite delta caps it at (1 - d) / (1 + d).
    """
    xp = xp_of(cos_theta)
    mu = xp.asarray(cos_theta, dtype=float)
    Delta = (1.0 - depol) / (1.0 + 0.5 * depol)
    Delta_p = (1.0 - 2.0 * depol) / (1.0 - depol)

    P = xp.zeros(mu.shape + (4, 4))
    ray = 0.75 * (1.0 + mu**2)
    P[..., 0, 0] = Delta * ray + (1.0 - Delta)
    P[..., 0, 1] = -Delta * 0.75 * (1.0 - mu**2)
    P[..., 1, 0] = P[..., 0, 1]
    P[..., 1, 1] = Delta * ray
    P[..., 2, 2] = Delta * 1.5 * mu
    P[..., 3, 3] = Delta * Delta_p * 1.5 * mu
    return P


def hg_phase_mueller(cos_theta, g: float, depol: float = 0.2):
    """Particulate phase Mueller matrix: Henyey-Greenstein intensity
    with the polarization structure of a depolarized Rayleigh matrix,

        P(mu) = P_HG(mu; g) * R(mu; depol) / R00(mu; depol),

    so P[0, 0] is exactly HG (<P00> = 1 over the sphere) and the
    polarization ratio P01/P00 peaks at ~0.65 for depol = 0.2 near
    90 deg, the Voss & Fry (1984) average-ocean magnitude."""
    xp = xp_of(cos_theta)
    mu = xp.asarray(cos_theta, dtype=float)
    hg = (1.0 - g * g) / (1.0 + g * g - 2.0 * g * mu) ** 1.5
    R = rayleigh_phase_mueller(mu, depol)
    return R * (hg / R[..., 0, 0])[..., None, None]


def sample_rayleigh_scattering(n: int, rng, depol: float = RAYLEIGH_DEPOL,
                               xp=None):
    """n scattering-angle cosines drawn from the unpolarized phase
    function P[0, 0] / (4 pi): a Delta-weighted mixture of the
    (1 + mu^2) Rayleigh kernel (closed-form cubic inversion) and an
    isotropic remainder.  Azimuth is uniform and sampled by the caller.
    """
    Delta = (1.0 - depol) / (1.0 + 0.5 * depol)
    u = rng.random(n)
    xp = xp_of(u) if xp is None else xp
    use_ray = rng.random(n) < Delta

    mu = 2.0 * u - 1.0
    # Rayleigh branch: invert (3/8)(mu + mu^3/3 + 4/3) = u via the
    # depressed cubic mu^3 + 3 mu = 8u - 4
    m = 8.0 * u[use_ray] - 4.0
    w = xp.cbrt(0.5 * m + xp.sqrt(0.25 * m**2 + 1.0))
    mu[use_ray] = w - 1.0 / w
    return xp.clip(mu, -1.0, 1.0)


def sample_hg_scattering(n: int, rng, g: float, xp=None):
    """n scattering-angle cosines from the Henyey-Greenstein phase
    function (closed-form inversion)."""
    u = rng.random(n)
    xp = xp_of(u) if xp is None else xp
    if abs(g) < 1e-6:
        return 2.0 * u - 1.0
    frac = (1.0 - g * g) / (1.0 - g + 2.0 * g * u)
    return xp.clip((1.0 + g * g - frac * frac) / (2.0 * g), -1.0, 1.0)


# ---------------------------------------------------------------------------
# Fournier-Forand particulate phase function (turbid Case 2 water)
# ---------------------------------------------------------------------------

def _ff_beta(mu, n_particle: float, mu_junge: float, xp):
    """Raw (Mobley-normalized) Fournier-Forand phase function beta(mu),
    int beta dOmega = 1.  mu = cos(scattering angle); the forward
    singularity at mu = 1 is clipped (measure zero in the MC)."""
    mu = xp.clip(xp.asarray(mu, dtype=float), -1.0, 1.0 - 1e-9)
    s2 = xp.maximum((1.0 - mu) / 2.0, 1e-12)        # sin^2(theta/2)
    nu = (3.0 - mu_junge) / 2.0
    c = 4.0 / (3.0 * (n_particle - 1.0) ** 2)
    d = xp.clip(c * s2, 1e-12, None)
    d180 = c                                        # delta at theta = 180
    dnu = d ** nu
    d180nu = d180 ** nu
    t1 = (1.0 / (4.0 * np.pi * (1.0 - d) ** 2 * dnu)) * (
        nu * (1.0 - d) - (1.0 - dnu)
        + (d * (1.0 - dnu) - nu * (1.0 - d)) / s2)
    t2 = ((1.0 - d180nu) / (16.0 * np.pi * (d180 - 1.0) * d180nu)
          * (3.0 * mu**2 - 1.0))
    return t1 + t2


_FF_CACHE: dict = {}


def _ff_model(n_particle: float, mu_junge: float):
    """Cached (norm, mu_grid, cdf) for a Fournier-Forand model.

    The Mobley beta integrates to 1 over the sphere analytically, so
    <P00> = 1 requires norm = 4 pi (the forward peak is too singular to
    normalize reliably by quadrature).  The sampling CDF over mu is built
    on a grid refined toward the forward peak mu = 1, where the phase
    function -- and hence the CDF -- changes fastest."""
    if mu_junge <= 3.0:
        raise ValueError("Fournier-Forand needs mu_junge > 3 (the Junge "
                         "slope; nu = (3 - mu_junge)/2 must be < 0)")
    key = (round(float(n_particle), 5), round(float(mu_junge), 5))
    if key not in _FF_CACHE:
        # mu in [-1, 1], increasing, dense toward the forward peak mu = 1
        lin = np.linspace(0.0, 1.0, 200001)
        mu = 1.0 - 2.0 * (1.0 - lin) ** 2.5
        pdf = np.clip(_ff_beta(mu, n_particle, mu_junge, NUMPY_XP), 0.0, None)
        cdf = np.concatenate([[0.0], np.cumsum(
            0.5 * (pdf[1:] + pdf[:-1]) * np.diff(mu))])
        cdf = cdf / cdf[-1]
        _FF_CACHE[key] = (4.0 * np.pi, mu, cdf)
    return _FF_CACHE[key]


def fournier_forand_phase(cos_theta, n_particle: float = 1.05,
                          mu_junge: float = 3.5):
    """Fournier-Forand phase function P00(cos_theta), normalized so its
    average over the sphere is 1 (<P00> = 1), matching the Rayleigh/HG
    convention.

    n_particle : real refractive index of the particles relative to
                 water (~1.02-1.20)
    mu_junge   : Junge (hyperbolic) slope of the particle size
                 distribution (> 3; ocean particulates ~3.5-4.5)

    The backscatter fraction bb/b grows with both parameters, spanning
    the measured ocean range (~0.005 clear to ~0.03 turbid) -- the
    realistic backscattering that a single Henyey-Greenstein term
    under-represents.  Fournier & Forand (1994); Mobley (2002)."""
    xp = xp_of(cos_theta)
    norm, _, _ = _ff_model(n_particle, mu_junge)
    return norm * _ff_beta(cos_theta, n_particle, mu_junge, xp)


def ff_phase_mueller(cos_theta, n_particle: float = 1.05,
                     mu_junge: float = 3.5, depol: float = 0.2):
    """Particulate phase Mueller matrix with the Fournier-Forand
    intensity (realistic backscatter) and a depolarized-Rayleigh
    polarization structure, exactly as hg_phase_mueller does for HG."""
    xp = xp_of(cos_theta)
    mu = xp.asarray(cos_theta, dtype=float)
    P00 = fournier_forand_phase(mu, n_particle, mu_junge)
    R = rayleigh_phase_mueller(mu, depol)
    return R * (P00 / R[..., 0, 0])[..., None, None]


def sample_ff_scattering(n: int, rng, n_particle: float = 1.05,
                         mu_junge: float = 3.5, xp=None):
    """n scattering-angle cosines drawn from the Fournier-Forand phase
    function by inverse-CDF sampling (numerically tabulated)."""
    u = rng.random(n)
    xp = xp_of(u) if xp is None else xp
    _, mu_grid, cdf = _ff_model(n_particle, mu_junge)
    return xp.clip(xp.interp(u, xp.asarray(cdf), xp.asarray(mu_grid)),
                   -1.0, 1.0)


def polarized_scatter_event(d, S_path, mu, p00, p01, P, rng):
    """Generic polarized volume-scattering event for K rays sharing one
    set of sampled scattering cosines mu and phase Mueller matrices P
    (with intensity/polarized kernel coefficients p00 = P[..., 0, 0],
    p01 = P[..., 0, 1]).

    The azimuth is rejection-sampled from the polarized intensity
    kernel N(phi) = p00 I' + p01 (cos 2a Q' + sin 2a U') of the path
    Stokes S_path (meridian frame of d), and the event matrix is scaled
    by I'/N so the I-component path weight is exactly preserved.

    Returns (d_new, M_step) with M_step including both frame rotations.
    """
    xp = xp_of(d, S_path)
    K = d.shape[0]
    sn = xp.sqrt(xp.clip(1.0 - mu**2, 0.0, 1.0))

    v_in, h_in = meridian_frame(d)
    I_p = xp.maximum(S_path[:, 0], 0.0)
    env = p00 * I_p + xp.abs(p01) * xp.hypot(S_path[:, 1], S_path[:, 2])

    d_new = xp.empty_like(d)
    s_axis = xp.empty_like(d)
    a_in = xp.empty(K)
    N_val = xp.empty(K)
    rem = xp.arange(K)
    for it in range(200):
        n_rem = rem.shape[0]
        ph = rng.uniform(0.0, 2.0 * np.pi, n_rem)
        dn = normalize(mu[rem, None] * d[rem]
                       + sn[rem, None] * (xp.cos(ph)[:, None] * v_in[rem]
                                          + xp.sin(ph)[:, None] * h_in[rem]))
        s_ax = xp.cross(d[rem], dn)
        s_n = xp.linalg.norm(s_ax, axis=-1, keepdims=True)
        s_ax = xp.where(s_n > 1e-9, s_ax / xp.maximum(s_n, 1e-300),
                        h_in[rem])
        p_in = xp.cross(s_ax, d[rem])
        a = frame_rotation_angle(d[rem], v_in[rem], p_in)
        N = (p00[rem] * I_p[rem]
             + p01[rem] * (xp.cos(2.0 * a) * S_path[rem, 1]
                           + xp.sin(2.0 * a) * S_path[rem, 2]))
        # zero-weight degenerate paths fall back to uniform azimuth
        acc = ((rng.random(n_rem) * env[rem] <= N)
               | (env[rem] <= 1e-300) | (it == 199))
        tgt = rem[acc]
        d_new[tgt] = dn[acc]
        s_axis[tgt] = s_ax[acc]
        a_in[tgt] = a[acc]
        N_val[tgt] = N[acc]
        rem = rem[~acc]
        if rem.shape[0] == 0:
            break

    v_out, _ = meridian_frame(d_new)
    p_out = xp.cross(s_axis, d_new)
    a_out = frame_rotation_angle(d_new, p_out, v_out)
    scale = xp.where(N_val > 1e-300, I_p / xp.maximum(N_val, 1e-300),
                     1.0 / p00)
    M_step = (mueller_rotation(a_out) @ (P * scale[:, None, None])
              @ mueller_rotation(a_in))
    return d_new, M_step


def water_leaving_stokes(d_out, n_hat, water: WaterBody, E_d=np.pi,
                         n_water: float = 1.34):
    """Water-leaving Stokes contribution (..., 4) along d_out (surface ->
    camera) through facets with upward normals n_hat (first-order
    isotropic model).

    E_d is the downwelling irradiance entering the water; for a uniform
    unpolarized sky of radiance I_sky, E_d = pi I_sky."""
    xp = xp_of(d_out, n_hat)
    d_out = normalize(xp.asarray(d_out, dtype=float))
    n_hat = normalize(xp.asarray(n_hat, dtype=float))

    cos_air = xp.clip(xp.sum(d_out * n_hat, axis=-1), 0.0, 1.0)
    sin_w = xp.sqrt(xp.clip(1.0 - cos_air**2, 0.0, 1.0)) / n_water
    cos_w = xp.sqrt(xp.clip(1.0 - sin_w**2, 0.0, 1.0))

    _, M_T, _ = fresnel_mueller(cos_w, 1.0 / n_water)
    L_u = water.R_w * E_d / np.pi
    S_u = xp.zeros(cos_w.shape + (4,))
    S_u[..., 0] = L_u
    S = xp.einsum("...ij,...j->...i", M_T, S_u) / n_water**2
    return xp.where((cos_air > 0.0)[..., None], S, xp.zeros_like(S))
