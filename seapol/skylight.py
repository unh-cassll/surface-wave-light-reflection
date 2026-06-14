"""
Rayleigh single-scattering polarized skylight, after Xue et al. (2021)
Sec. 3.1 and Coulson (1988), plus nominal spectral radiance scalings for
the color-aware mode.

The polarization geometry is computed with vectors rather than the arccos
angle formulas, which preserves the sign of U across the sun's meridian.
Stokes vectors are returned in the meridian frame of the (downward)
propagation direction; +Q = vertical (parallel to the meridian plane).

Spectral model (single-scattering, optically thin): sky radiance scales
as E_sun(lambda) tau_R(lambda) for the molecular sky, with an aerosol
component tau_a ~ lambda^-alpha mixed in by turbidity; the direct beam
is Beer-attenuated through the Rayleigh + aerosol column, which reddens
the low sun.  All spectra are normalized to 1 at 550 nm so the scalar
(band-integrated) API is recovered at a single mid-visible band.
"""

from __future__ import annotations

import numpy as np

from .backend import xp_of
from .polarization import meridian_frame, normalize

# Empirical maximum sky DoLP vs solar zenith angle (Coulson 1988)
_SZA_TABLE_DEG = np.array([0.0, 30.0, 60.0, 90.0])
_DOLP_MAX_TABLE = np.array([0.56, 0.63, 0.70, 0.77])


def dolp_max(sun_zenith_deg: float) -> float:
    """Interpolated empirical maximum sky DoLP."""
    return float(np.interp(sun_zenith_deg, _SZA_TABLE_DEG, _DOLP_MAX_TABLE))


def direction_from_angles(zenith, azimuth):
    """Unit vector(s) (..., 3) pointing from the surface toward a sky point
    at the given zenith/azimuth [rad].  Azimuth from +x, counter-clockwise."""
    xp = xp_of(zenith, azimuth)
    zenith = xp.asarray(zenith, dtype=float)
    azimuth = xp.asarray(azimuth, dtype=float)
    return xp.stack([xp.sin(zenith) * xp.cos(azimuth),
                     xp.sin(zenith) * xp.sin(azimuth),
                     xp.cos(zenith) * xp.ones_like(azimuth)], axis=-1)


def rayleigh_sky_stokes(sky_dirs, sun_dir, I_sky=1.0,
                        sun_zenith_deg: float | None = None):
    """Stokes vector (..., 4) of skylight arriving from sky_dirs.

    sky_dirs : (..., 3) unit vectors surface -> sky point
    sun_dir  : (3,) unit vector surface -> sun
    I_sky    : scalar (or broadcastable) sky radiance

    The Stokes frame is the meridian frame of the downward propagation
    direction d = -sky_dirs.  Rayleigh E-field lies perpendicular to the
    scattering plane spanned by sky_dirs and sun_dir.
    """
    xp = xp_of(sky_dirs)
    sky_dirs = normalize(xp.asarray(sky_dirs, dtype=float))
    sun_dir = normalize(xp.asarray(sun_dir, dtype=float))
    if sun_zenith_deg is None:
        sun_zenith_deg = float(np.degrees(np.arccos(
            np.clip(float(sun_dir[2]), -1.0, 1.0))))

    sun_b = xp.broadcast_to(sun_dir, sky_dirs.shape)
    cos_g = xp.clip(xp.sum(sky_dirs * sun_b, axis=-1), -1.0, 1.0)
    Dmax = dolp_max(sun_zenith_deg)
    DoLP = Dmax * (1.0 - cos_g**2) / (1.0 + cos_g**2)

    d = -sky_dirs
    e_pol = xp.cross(sky_dirs, sun_b)
    n = xp.linalg.norm(e_pol, axis=-1, keepdims=True)
    # Degenerate at/opposite the sun: DoLP -> 0, axis arbitrary
    x_hat = xp.broadcast_to(xp.asarray([1.0, 0.0, 0.0]), sky_dirs.shape)
    e_pol = xp.where(n > 1e-12, e_pol / xp.maximum(n, 1e-300), x_hat)

    v, h = meridian_frame(d)
    psi = xp.arctan2(xp.sum(e_pol * h, axis=-1), xp.sum(e_pol * v, axis=-1))

    I = xp.broadcast_to(xp.asarray(I_sky, dtype=float), DoLP.shape)
    S = xp.zeros(DoLP.shape + (4,))
    S[..., 0] = I
    S[..., 1] = I * DoLP * xp.cos(2.0 * psi)
    S[..., 2] = I * DoLP * xp.sin(2.0 * psi)
    return S


def rayleigh_sky_stokes_angles(sky_zenith, sky_azimuth,
                               sun_zenith: float, sun_azimuth: float,
                               I_sky=1.0):
    """Angle-based wrapper; all angles in radians."""
    sky_dirs = direction_from_angles(sky_zenith, sky_azimuth)
    xp = xp_of(sky_dirs)
    sun_dir = xp.asarray([np.sin(sun_zenith) * np.cos(sun_azimuth),
                          np.sin(sun_zenith) * np.sin(sun_azimuth),
                          np.cos(sun_zenith)])
    return rayleigh_sky_stokes(sky_dirs, sun_dir, I_sky,
                               sun_zenith_deg=float(np.degrees(sun_zenith)))


def overcast_sky_stokes(sky_dirs, I_zenith=1.0):
    """Unpolarized overcast sky with the Moon-Spencer luminance gradation
    L(theta) = I_zenith (1 + 2 cos theta) / 3."""
    xp = xp_of(sky_dirs)
    sky_dirs = normalize(xp.asarray(sky_dirs, dtype=float))
    cos_z = xp.clip(sky_dirs[..., 2], 0.0, 1.0)
    S = xp.zeros(cos_z.shape + (4,))
    S[..., 0] = I_zenith * (1.0 + 2.0 * cos_z) / 3.0
    return S


def sun_disk_stokes(sky_dirs, sun_dir, I_sun=100.0,
                    halfwidth_deg: float = 1.5):
    """Unpolarized direct-solar term: Gaussian glow of the given angular
    half-width around the sun direction (a tractable stand-in for the
    0.27-deg disk plus circumsolar aureole)."""
    xp = xp_of(sky_dirs)
    sky_dirs = normalize(xp.asarray(sky_dirs, dtype=float))
    sun_dir = normalize(xp.asarray(sun_dir, dtype=float))
    sun_b = xp.broadcast_to(sun_dir, sky_dirs.shape)
    cos_g = xp.clip(xp.sum(sky_dirs * sun_b, axis=-1), -1.0, 1.0)
    gamma = xp.arccos(cos_g)
    hw = np.deg2rad(halfwidth_deg)
    S = xp.zeros(gamma.shape + (4,))
    S[..., 0] = I_sun * xp.exp(-0.5 * (gamma / hw) ** 2)
    return S


# ---------------------------------------------------------------------------
# Nominal spectral radiance scalings (color-aware mode)
# ---------------------------------------------------------------------------

# Smoothed extraterrestrial-to-surface solar spectral irradiance shape
# (ASTM G173 global tilt, heavily smoothed), relative to 550 nm, 10 nm
# grid 380-730 nm.  Captures the broad solar curve without Fraunhofer
# line detail -- intended for nominal scene color, not radiometry.
_SOLAR_WL_NM = np.arange(380.0, 731.0, 10.0)
_SOLAR_REL = np.array([
    0.58, 0.70, 0.83, 0.94, 1.00, 1.02, 1.04, 1.06, 1.07, 1.07,
    1.06, 1.05, 1.04, 1.03, 1.02, 1.01, 1.00, 0.99, 0.98, 0.97,
    0.96, 0.95, 0.93, 0.92, 0.90, 0.89, 0.87, 0.86, 0.84, 0.83,
    0.81, 0.80, 0.78, 0.77, 0.75, 0.74])


def solar_spectrum(wavelength_nm):
    """Relative surface solar spectral irradiance, normalized to 1 at
    550 nm (smoothed ASTM G173 shape)."""
    wl = np.asarray(wavelength_nm, dtype=float)
    ref = np.interp(550.0, _SOLAR_WL_NM, _SOLAR_REL)
    return np.interp(wl, _SOLAR_WL_NM, _SOLAR_REL) / ref


def rayleigh_optical_depth(wavelength_nm):
    """Sea-level molecular (Rayleigh) optical depth tau_R(lambda)
    (Bodhaine et al. 1999 fit): ~0.097 at 550 nm, ~0.36 at 400 nm."""
    lam_um = np.asarray(wavelength_nm, dtype=float) / 1000.0
    return 0.008569 * lam_um**-4 * (1.0 + 0.0113 * lam_um**-2
                                    + 0.00013 * lam_um**-4)


def aerosol_optical_depth(wavelength_nm, turbidity: float = 0.0,
                          angstrom_alpha: float = 1.3):
    """Nominal aerosol optical depth: tau_a(550) = 0.25 * turbidity with
    an Angstrom-law spectral slope (turbidity is the same 0-1 haze knob
    as in make_clear_sky)."""
    wl = np.asarray(wavelength_nm, dtype=float)
    return 0.25 * float(turbidity) * (wl / 550.0) ** (-angstrom_alpha)


def sky_radiance_spectrum(wavelength_nm, turbidity: float = 0.0):
    """Relative diffuse-sky radiance spectrum (normalized to 1 at
    550 nm): single-scattering mix of the blue Rayleigh component
    (E_sun tau_R) and a spectrally flat aerosol haze component weighted
    by their optical depths."""
    wl = np.asarray(wavelength_nm, dtype=float)
    E = solar_spectrum(wl)
    tau_r = rayleigh_optical_depth(wl)
    tau_a = aerosol_optical_depth(wl, turbidity)
    num = E * (tau_r + tau_a)
    den = (solar_spectrum(550.0)
           * (rayleigh_optical_depth(550.0)
              + aerosol_optical_depth(550.0, turbidity)))
    return num / den


def sun_beam_spectrum(wavelength_nm, sun_zenith_deg: float,
                      turbidity: float = 0.0):
    """Relative direct-beam spectral irradiance at the surface,
    normalized to 1 at 550 nm: solar spectrum Beer-attenuated through
    the Rayleigh + aerosol column at the solar airmass (reddens the low
    sun)."""
    wl = np.asarray(wavelength_nm, dtype=float)
    mu = max(np.cos(np.deg2rad(sun_zenith_deg)), 0.05)
    tau = rayleigh_optical_depth(wl) + aerosol_optical_depth(wl, turbidity)
    tau550 = (rayleigh_optical_depth(550.0)
              + aerosol_optical_depth(550.0, turbidity))
    return (solar_spectrum(wl) * np.exp(-tau / mu)
            / (solar_spectrum(550.0) * np.exp(-tau550 / mu)))


def overcast_spectrum(wavelength_nm):
    """Relative overcast-sky radiance spectrum (normalized at 550 nm):
    cloud-transmitted sunlight, spectrally neutral scattering."""
    wl = np.asarray(wavelength_nm, dtype=float)
    return solar_spectrum(wl) / solar_spectrum(550.0)
