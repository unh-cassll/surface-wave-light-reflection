"""
Single-bounce polarized reflection renderers for a synthesized sea surface.

Two viewing models:
    * render_facet_stokes : one Stokes vector per surface grid cell, finite
      camera position (improved port of
      produce_simulated_sea_surface_modeled_reflection.m).
    * render_camera_image : pinhole camera with per-pixel rays
      (improved port of sea_polarization_image.py).

Both share the same Mueller chain (polarization.reflection_chain):
sky meridian frame -> scattering plane -> Fresnel reflection -> camera
meridian frame, with optional subpixel slope ensembles and Smith/Saunders
bistatic shadowing.

The `water` argument of both renderers accepts either a water.WaterBody
(first-order isotropic water-leaving radiance) or a
scattering.UpwellingRadianceTable (directional polarized sub-surface
light field from the in-water Monte Carlo), which is what breaks the
mirror-only look of the scene.

All array work dispatches on the input arrays (numpy or torch).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .backend import adapt_rng, xp_of
from .polarization import apply_mueller, normalize, reflection_chain
from .skylight import (direction_from_angles, overcast_sky_stokes,
                       rayleigh_sky_stokes, sun_disk_stokes)
from .spectrum import cutoff_slope_variances


# ---------------------------------------------------------------------------
# Sky models
# ---------------------------------------------------------------------------

def make_rayleigh_sky(sun_zenith_deg: float, sun_azimuth_deg: float,
                      I_sky: float = 1.0):
    """Callable sky_dirs (..., 3) -> Stokes (..., 4) for a Rayleigh sky."""
    sun_dir = direction_from_angles(np.deg2rad(sun_zenith_deg),
                                    np.deg2rad(sun_azimuth_deg))

    def sky(sky_dirs):
        return rayleigh_sky_stokes(sky_dirs, sun_dir, I_sky,
                                   sun_zenith_deg=sun_zenith_deg)
    return sky


def make_unpolarized_sky(I_sky: float = 1.0):
    """Uniform unpolarized sky."""
    def sky(sky_dirs):
        xp = xp_of(sky_dirs)
        d = xp.asarray(sky_dirs)
        S = xp.zeros(d.shape[:-1] + (4,))
        S[..., 0] = I_sky
        return S
    return sky


def compose_skies(*sky_fns):
    """Additive composition of sky models (Stokes vectors add)."""
    def sky(sky_dirs):
        return sum(fn(sky_dirs) for fn in sky_fns)
    return sky


def make_clear_sky(sun_zenith_deg: float, sun_azimuth_deg: float,
                   I_sky: float = 1.0, turbidity: float = 0.0,
                   I_sun: float = 0.0, sun_halfwidth_deg: float = 1.5):
    """Clear sky: Rayleigh polarization pattern with an aerosol turbidity
    knob (0 = pure Rayleigh; 1 = fully depolarized haze) and an optional
    unpolarized direct-solar glow (I_sun > 0)."""
    sun_dir = direction_from_angles(np.deg2rad(sun_zenith_deg),
                                    np.deg2rad(sun_azimuth_deg))
    scale = float(np.clip(1.0 - turbidity, 0.0, 1.0))

    def sky(sky_dirs):
        S = rayleigh_sky_stokes(sky_dirs, sun_dir, I_sky,
                                sun_zenith_deg=sun_zenith_deg)
        S[..., 1:3] *= scale
        if I_sun > 0.0:
            S = S + sun_disk_stokes(sky_dirs, sun_dir, I_sun,
                                    sun_halfwidth_deg)
        return S
    return sky


def make_overcast_sky(I_zenith: float = 1.0):
    """Unpolarized overcast sky (Moon-Spencer gradation)."""
    def sky(sky_dirs):
        return overcast_sky_stokes(sky_dirs, I_zenith)
    return sky


def make_partly_cloudy_sky(sun_zenith_deg: float, sun_azimuth_deg: float,
                           I_sky: float = 1.0,
                           cloud_fraction: float = 0.4,
                           cloud_brightness: float = 3.0,
                           n_grid: int = 96,
                           rng: np.random.Generator | None = None):
    """Clear Rayleigh sky with a procedural broken-cloud field: smooth
    random texture on a direction grid, thresholded at cloud_fraction;
    cloudy directions are replaced by bright unpolarized radiance.  If a
    cloud covers the sun, the polarized clear-sky pattern persists (single
    scattering from the remaining clear paths) but real occlusion effects
    are not modeled."""
    if rng is None:
        rng = np.random.default_rng(0)
    from scipy.ndimage import gaussian_filter
    sun_dir = direction_from_angles(np.deg2rad(sun_zenith_deg),
                                    np.deg2rad(sun_azimuth_deg))
    # cloud texture on a (zenith, azimuth) grid, periodic in azimuth
    tex = gaussian_filter(rng.standard_normal((n_grid, 2 * n_grid)),
                          sigma=n_grid / 12, mode="wrap")
    thresh = np.quantile(tex, 1.0 - cloud_fraction)
    cloudy = tex > thresh

    def sky(sky_dirs):
        xp = xp_of(sky_dirs)
        S = rayleigh_sky_stokes(sky_dirs, sun_dir, I_sky,
                                sun_zenith_deg=sun_zenith_deg)
        d = xp.asarray(sky_dirs, dtype=float)
        zen = xp.arccos(xp.clip(d[..., 2], -1.0, 1.0))
        az = xp.mod(xp.arctan2(d[..., 1], d[..., 0]), 2.0 * np.pi)
        iz = xp.clip(xp.astype(zen / (np.pi / 2) * (n_grid - 1), int),
                     0, n_grid - 1)
        ia = xp.clip(xp.astype(az / (2 * np.pi) * (2 * n_grid - 1), int),
                     0, 2 * n_grid - 1)
        mask = xp.asarray(cloudy)[iz, ia]
        S_cloud = xp.zeros_like(S)
        S_cloud[..., 0] = I_sky * cloud_brightness
        return xp.where(mask[..., None], S_cloud, S)
    return sky


# ---------------------------------------------------------------------------
# Smith / Saunders shadowing
# ---------------------------------------------------------------------------

def smith_lambda(theta, sigma_slope: float):
    """Smith (1967) shadowing integral Lambda(theta) for an isotropic
    Gaussian surface with total slope std sigma_slope."""
    xp = xp_of(theta)
    theta = xp.asarray(theta, dtype=float)
    if sigma_slope <= 0.0:
        return xp.zeros_like(theta)
    mu = 1.0 / (xp.tan(xp.clip(theta, 1e-9, np.pi / 2 - 1e-9))
                * np.sqrt(2.0) * sigma_slope)
    lam = 0.5 * (xp.exp(-mu**2) / (np.sqrt(np.pi) * mu) - xp.erfc(mu))
    return xp.maximum(lam, 0.0)


def saunders_bistatic(theta_i, theta_r, sigma_slope: float):
    """Bistatic shadowing factor S = 1 / (1 + Lambda_i + Lambda_r)."""
    return 1.0 / (1.0 + smith_lambda(theta_i, sigma_slope)
                  + smith_lambda(theta_r, sigma_slope))


# ---------------------------------------------------------------------------
# Geometry configuration
# ---------------------------------------------------------------------------

@dataclass
class CameraGeometry:
    """Finite camera for the facet renderer.  The camera sits at height_m,
    displaced horizontally by height_m * tan(incidence) from the scene
    center along the azimuth direction (azimuth from +x, CCW)."""
    incidence_deg: float = 30.0
    azimuth_deg: float = 0.0
    height_m: float = 100.0

    def position(self, center_xy: tuple[float, float]) -> np.ndarray:
        th = np.deg2rad(self.incidence_deg)
        ph = np.deg2rad(self.azimuth_deg)
        r = self.height_m * np.tan(th)
        return np.array([center_xy[0] + r * np.cos(ph),
                         center_xy[1] + r * np.sin(ph),
                         self.height_m])


@dataclass
class PinholeCamera:
    """Perspective camera for the image renderer, looking at the scene
    center from the given zenith/azimuth at altitude_m."""
    altitude_m: float = 300.0
    zenith_deg: float = 30.0
    azimuth_deg: float = 0.0
    hfov_deg: float = 8.0
    img_shape: tuple[int, int] = (256, 256)


@dataclass
class SubpixelSlopes:
    """Statistics of unresolved (sub-grid) facet slopes: Gaussian widths
    plus optional Gram-Charlier shape coefficients in the Cox & Munk
    (1954) form, mapped to direct moments of the standardized along-
    (a) and cross-wind (c) slopes:

        skew(a) = c03,  E[c^2 a] = c21,
        excess kurtosis: kurt(c) = c40, E[c^2 a^2] - 1 = c22/...,
        kurt(a) = c04.

    With c03 < 0 (the field value and the emergent FM98 statistic) the
    along-wind slope PDF is skewed toward steep forward faces, which is
    what produces the observed upwind/downwind glitter asymmetry.  All
    coefficients default to 0 (pure Gaussian)."""
    sigma_a2: float                 # along-wind slope variance
    sigma_c2: float                 # cross-wind slope variance
    wind_dir_rad: float = 0.0
    c21: float = 0.0
    c03: float = 0.0
    c40: float = 0.0
    c22: float = 0.0
    c04: float = 0.0

    @classmethod
    def from_spectrum(cls, U10: float, k_cutoff: float,
                      wind_dir_rad: float = 0.0, **kwargs):
        sa2, sc2 = cutoff_slope_variances(U10, k_cutoff, **kwargs)
        return cls(sigma_a2=sa2, sigma_c2=sc2, wind_dir_rad=wind_dir_rad)

    @classmethod
    def from_cox_munk(cls, U10: float, k_cutoff: float,
                      wind_dir_rad: float = 0.0, **kwargs):
        """Sub-grid widths from the spectral tail plus the Cox-Munk
        (1954) clean-surface Gram-Charlier shape coefficients.  The
        non-Gaussianity of sea-surface slopes lives at the parasitic-
        capillary scales, which at rendering grid resolutions are
        exactly the unresolved scales, so the full-sea shape
        coefficients are the right leading model for the sub-pixel
        ensemble."""
        sa2, sc2 = cutoff_slope_variances(U10, k_cutoff, **kwargs)
        return cls(sigma_a2=sa2, sigma_c2=sc2, wind_dir_rad=wind_dir_rad,
                   c21=0.01 - 0.0086 * U10, c03=0.04 - 0.033 * U10,
                   c40=0.40, c22=0.12, c04=0.23)

    @property
    def sigma_total(self) -> float:
        return float(np.sqrt(self.sigma_a2 + self.sigma_c2))

    @property
    def is_gaussian(self) -> bool:
        return not (self.c21 or self.c03 or self.c40 or self.c22
                    or self.c04)

    def gc_factor(self, a, c):
        """Gram-Charlier correction factor G(a, c) >= 0 multiplying the
        Gaussian PDF; a, c are along/cross-wind slope deviations (not
        standardized).  G integrates to 1 against the Gaussian exactly
        (Hermite orthogonality); clipping the small negative far-tail
        lobes adds a bias that is negligible for field-size
        coefficients.  The correction is smoothly tapered beyond ~4
        standard deviations (retaining the moment-carrying core): the
        GC series is asymptotic and its Hermite polynomials dominate
        the Gaussian in the far tails, producing spurious bright wings
        in log-scale glint transects (Cox & Munk limited their fit to
        ~2.5 sigma)."""
        xp = xp_of(a, c)
        at = xp.asarray(a, dtype=float) / np.sqrt(max(self.sigma_a2, 1e-12))
        ct = xp.asarray(c, dtype=float) / np.sqrt(max(self.sigma_c2, 1e-12))
        corr = (0.5 * self.c21 * (ct**2 - 1.0) * at
                + (self.c03 / 6.0) * (at**3 - 3.0 * at)
                + (self.c40 / 24.0) * (ct**4 - 6.0 * ct**2 + 3.0)
                + 0.25 * self.c22 * (ct**2 - 1.0) * (at**2 - 1.0)
                + (self.c04 / 24.0) * (at**4 - 6.0 * at**2 + 3.0))
        r2 = at**2 + ct**2
        corr = corr * xp.exp(-((r2 / 16.0) ** 4))
        return xp.maximum(1.0 + corr, 0.0)

    def sample(self, shape, rng):
        """Earth-frame slope perturbations (xi_x, xi_y) of given shape."""
        xa = rng.normal(0.0, np.sqrt(max(self.sigma_a2, 0.0)), shape)
        xc = rng.normal(0.0, np.sqrt(max(self.sigma_c2, 0.0)), shape)
        cw, sw = np.cos(self.wind_dir_rad), np.sin(self.wind_dir_rad)
        return cw * xa - sw * xc, sw * xa + cw * xc

    def sample_weighted(self, shape, rng):
        """Gaussian draws plus the Gram-Charlier importance weight per
        draw: ensemble averages weighted by w follow the non-Gaussian
        slope PDF.  Returns (xi_x, xi_y, w)."""
        xa = rng.normal(0.0, np.sqrt(max(self.sigma_a2, 0.0)), shape)
        xc = rng.normal(0.0, np.sqrt(max(self.sigma_c2, 0.0)), shape)
        xp = xp_of(xa)
        w = xp.ones(shape) if self.is_gaussian else self.gc_factor(xa, xc)
        cw, sw = np.cos(self.wind_dir_rad), np.sin(self.wind_dir_rad)
        return cw * xa - sw * xc, sw * xa + cw * xc, w


# ---------------------------------------------------------------------------
# Whitecap foam
# ---------------------------------------------------------------------------

@dataclass
class Foam:
    """Whitecap foam as an unpolarized Lambertian cover on the steepest
    facets.  coverage is the areal fraction; the facets selected are
    those with the largest resolved slope magnitude (a breaking-crest
    proxy, so foam rides the steep faces the hybrid synthesis sharpens
    rather than being sprinkled at random).  albedo is a band-integrated
    effective foam reflectance (~0.4-0.55 in the visible for fresh
    foam, decaying with age)."""
    coverage: float
    albedo: float = 0.45


def monahan_coverage(U10: float) -> float:
    """Monahan & O'Muircheartaigh (1980) whitecap fraction
    W = 3.84e-6 U10^3.41."""
    return float(3.84e-6 * max(U10, 0.0) ** 3.41)


def _apply_foam(S, slope_x, slope_y, foam: Foam, sky, sun_glint):
    """Replace the Stokes vectors of the steepest-slope facets with
    unpolarized Lambertian foam radiance albedo * E_d / pi."""
    xp = xp_of(S)
    cov = float(np.clip(foam.coverage, 0.0, 1.0))
    if cov <= 0.0:
        return S
    mag = xp.hypot(slope_x, slope_y)
    valid = xp.isfinite(S[..., 0])
    if not bool(xp.any(valid)):
        return S
    thresh = xp.nanquantile(xp.where(valid, mag, xp.nan * mag), 1.0 - cov)
    mask = (mag >= thresh) & valid
    E_d = downwelling_irradiance(sky)
    if sun_glint is not None:
        zen_deg, _, E_sun = sun_glint
        E_d += E_sun * max(np.cos(np.deg2rad(zen_deg)), 0.0)
    S_foam = xp.zeros_like(S)
    S_foam[..., 0] = foam.albedo * E_d / np.pi
    return xp.where(mask[..., None], S_foam, S)


# ---------------------------------------------------------------------------
# Shared single-sample Stokes evaluation
# ---------------------------------------------------------------------------

def sun_glint_stokes(d_out, sx0, sy0, subpixel: SubpixelSlopes,
                     sun_zenith_deg: float, sun_azimuth_deg: float,
                     E_sun: float, n_water: float = 1.34,
                     shadow_sigma: float | None = None):
    """Analytic Cox-Munk sun glint: the sub-pixel Gaussian slope PDF
    evaluated at the specular slope, with the full Fresnel Mueller chain.

        S = E_sun p(z*) / (4 cos(theta_v) cos^4(theta_n))
            R(out) M_R(omega) R(in) [1, 0, 0, 0]

    where z* is the facet slope that mirrors the sun into d_out and p is
    the subpixel slope PDF centered on the resolved macro slope
    (sx0, sy0).  Replaces Monte Carlo sun-disk sampling, which
    under-resolves the glint probability and produces speckle.  Do not
    also include a sun disk in the sky model (double counting)."""
    xp = xp_of(d_out, sx0, sy0)
    sun_dir = direction_from_angles(np.deg2rad(sun_zenith_deg),
                                    np.deg2rad(sun_azimuth_deg))
    d_out = normalize(xp.asarray(d_out, dtype=float))
    sun_b = xp.broadcast_to(xp.asarray(sun_dir, dtype=float), d_out.shape)
    H = normalize(sun_b + d_out)
    Hz = H[..., 2]
    valid = (Hz > 1e-3) & (d_out[..., 2] > 1e-3) & (float(sun_dir[2]) > 0)
    Hz_s = xp.clip(Hz, 1e-3, 1.0)

    # specular facet slope and the subpixel PDF there
    zx = -H[..., 0] / Hz_s
    zy = -H[..., 1] / Hz_s
    da_x = zx - sx0
    da_y = zy - sy0
    cw = np.cos(subpixel.wind_dir_rad)
    sw = np.sin(subpixel.wind_dir_rad)
    da = cw * da_x + sw * da_y
    dc = -sw * da_x + cw * da_y
    sa2 = max(subpixel.sigma_a2, 1e-12)
    sc2 = max(subpixel.sigma_c2, 1e-12)
    p = xp.exp(-0.5 * (da**2 / sa2 + dc**2 / sc2)) \
        / (2.0 * np.pi * np.sqrt(sa2 * sc2))
    if not subpixel.is_gaussian:
        # Gram-Charlier slope PDF: skewed toward steep forward faces
        # (c03 < 0), reproducing the upwind/downwind glitter asymmetry
        p = p * subpixel.gc_factor(da, dc)

    M, _, chain_ok = reflection_chain(d_out, H, n_water)
    geom = p / (4.0 * xp.clip(d_out[..., 2], 1e-3, 1.0) * Hz_s**4)
    S = geom[..., None] * (E_sun * M[..., :, 0])
    if shadow_sigma is not None and shadow_sigma > 0.0:
        th_i = np.arccos(np.clip(sun_dir[2], 1e-9, 1.0))
        th_r = xp.arccos(xp.clip(d_out[..., 2], 1e-9, 1.0))
        S = S * saunders_bistatic(th_i, th_r, shadow_sigma)[..., None]
    return xp.where((valid & chain_ok)[..., None], S, xp.zeros_like(S))


def downwelling_irradiance(sky_fn, n_zen: int = 16, n_az: int = 32) -> float:
    """Plane downwelling irradiance E_d = integral I cos(theta) dOmega of a
    sky model, by simple hemisphere quadrature."""
    zen = (np.arange(n_zen) + 0.5) * (np.pi / 2) / n_zen
    az = (np.arange(n_az) + 0.5) * 2 * np.pi / n_az
    ZEN, AZ = np.meshgrid(zen, az, indexing="ij")
    dirs = direction_from_angles(ZEN, AZ)
    I = np.asarray(sky_fn(dirs))[..., 0]
    d_omega = np.sin(ZEN) * (np.pi / 2 / n_zen) * (2 * np.pi / n_az)
    return float(np.sum(I * np.cos(ZEN) * d_omega))


def _water_leaving_term(d_out, n_hat, water, E_d, n_water):
    """Dispatch the water-leaving radiance term on the water model type:
    water.WaterBody (first-order isotropic) or
    scattering.UpwellingRadianceTable (directional, polarized)."""
    from .water import WaterBody, water_leaving_stokes
    if isinstance(water, WaterBody):
        return water_leaving_stokes(d_out, n_hat, water, E_d, n_water)
    from .scattering import water_leaving_from_table
    return water_leaving_from_table(d_out, n_hat, water, n_water=n_water)


def _stokes_one_sample(d_out, sx, sy, sky, n_water, shadow_sigma,
                       water=None, E_d=0.0):
    """Reflected (+ water-leaving) Stokes (..., 4) and validity (...,)."""
    xp = xp_of(d_out, sx, sy)
    n_hat = normalize(xp.stack([-sx, -sy, xp.ones_like(sx)], axis=-1))
    M, d_in, valid = reflection_chain(d_out, n_hat, n_water)
    S_sky = sky(-d_in)
    S = apply_mueller(M, S_sky)
    if shadow_sigma is not None and shadow_sigma > 0.0:
        th_i = xp.arccos(xp.clip(-d_in[..., 2], 1e-9, 1.0))
        th_r = xp.arccos(xp.clip(d_out[..., 2], 1e-9, 1.0))
        S = S * saunders_bistatic(th_i, th_r, shadow_sigma)[..., None]
    if water is not None:
        S = S + _water_leaving_term(d_out, n_hat, water, E_d, n_water)
    return S, valid


def _accumulate_subpixel(d_out, sx0, sy0, sky, n_water, subpixel,
                         n_subpixel, shadowing, rng, water=None):
    """Ensemble average over subpixel slope samples (one slice at a time
    to bound memory)."""
    xp = xp_of(d_out, sx0, sy0)
    shadow_sigma = subpixel.sigma_total if (shadowing and subpixel) else None
    E_d = downwelling_irradiance(sky) if water is not None else 0.0
    if subpixel is None or n_subpixel < 1:
        S, valid = _stokes_one_sample(d_out, sx0, sy0, sky, n_water,
                                      shadow_sigma, water, E_d)
        return xp.where(valid[..., None], S, xp.nan * S)

    S_sum = xp.zeros(sx0.shape + (4,))
    count = xp.zeros(sx0.shape)
    for _ in range(n_subpixel):
        xx, xy, w = subpixel.sample_weighted(sx0.shape, rng)
        S, valid = _stokes_one_sample(d_out, sx0 + xx, sy0 + xy, sky,
                                      n_water, shadow_sigma, water, E_d)
        vw = xp.where(valid, w, xp.zeros_like(w))
        S_sum += vw[..., None] * xp.where(valid[..., None], S,
                                          xp.zeros_like(S))
        count += vw
    with np.errstate(invalid="ignore", divide="ignore"):
        S_avg = S_sum / count[..., None]
    return xp.where(count[..., None] > 0, S_avg, xp.nan * S_avg)


# ---------------------------------------------------------------------------
# Facet renderer
# ---------------------------------------------------------------------------

def render_facet_stokes(eta, dx: float,
                        camera: CameraGeometry = CameraGeometry(),
                        sky=None,
                        slope_x=None,
                        slope_y=None,
                        n_water: float = 1.34,
                        subpixel: SubpixelSlopes | None = None,
                        n_subpixel: int = 64,
                        shadowing: bool = False,
                        water=None,
                        sun_glint: tuple | None = None,
                        foam: Foam | None = None,
                        rng=None):
    """Reflected Stokes vector (H, W, 4) for every facet of eta(y, x).

    Each facet reflects the single sky direction that mirrors into the
    camera; with `subpixel`, results are averaged over an ensemble of
    unresolved slope perturbations (Gram-Charlier weighted when the
    subpixel statistics carry shape coefficients).  `water` (a
    water.WaterBody or scattering.UpwellingRadianceTable) adds the
    water-leaving radiance.  `sun_glint = (sun_zenith_deg,
    sun_azimuth_deg, E_sun)` adds the analytic Cox-Munk glint term
    (requires `subpixel`; leave the sun disk out of the sky model to
    avoid double counting).  `foam` (a render.Foam) covers the steepest
    facets with unpolarized Lambertian whitecap radiance.  Facets with
    no valid sky ray are NaN."""
    xp = xp_of(eta, slope_x, slope_y)
    if sky is None:
        sky = make_unpolarized_sky(1.0)
    rng = adapt_rng(rng, xp)
    eta = xp.asarray(eta, dtype=float)
    H, W = eta.shape

    if slope_x is None or slope_y is None:
        gy, gx = xp.gradient(eta, dx)
        slope_x = gx if slope_x is None else slope_x
        slope_y = gy if slope_y is None else slope_y

    xg = xp.arange(W, dtype=float) * dx
    yg = xp.arange(H, dtype=float) * dx
    X, Y = xp.meshgrid(xg, yg, indexing="xy")
    P = xp.stack([X, Y, eta], axis=-1)

    center = ((W - 1) * dx / 2.0, (H - 1) * dx / 2.0)
    cam = xp.asarray(camera.position(center), dtype=float)
    d_out = normalize(cam[None, None, :] - P)

    S = _accumulate_subpixel(d_out, slope_x, slope_y, sky, n_water,
                             subpixel, n_subpixel, shadowing, rng,
                             water=water)
    if sun_glint is not None:
        if subpixel is None:
            raise ValueError("sun_glint requires subpixel slope statistics")
        S = S + sun_glint_stokes(
            d_out, slope_x, slope_y, subpixel, *sun_glint,
            n_water=n_water,
            shadow_sigma=subpixel.sigma_total if shadowing else None)
    if foam is not None:
        S = _apply_foam(S, slope_x, slope_y, foam, sky, sun_glint)
    return S


def render_facet_stokes_stack(eta, dx: float,
                              slope_x=None,
                              slope_y=None,
                              **kwargs):
    """Per-frame facet rendering of an (H, W, T) elevation stack.
    Returns (H, W, 4, T)."""
    xp = xp_of(eta)
    eta = xp.asarray(eta, dtype=float)
    H, W, T = eta.shape
    out = xp.empty((H, W, 4, T))
    for it in range(T):
        sx = None if slope_x is None else slope_x[:, :, it]
        sy = None if slope_y is None else slope_y[:, :, it]
        out[:, :, :, it] = render_facet_stokes(eta[:, :, it], dx,
                                               slope_x=sx, slope_y=sy,
                                               **kwargs)
    return out


# ---------------------------------------------------------------------------
# Pinhole-camera image renderer
# ---------------------------------------------------------------------------

def _bilinear(field, dx: float, x, y):
    """Bilinear sample of field(y, x); returns (values, in_bounds)."""
    xp = xp_of(field, x, y)
    H, W = field.shape
    fi = y / dx
    fj = x / dx
    i0 = xp.astype(xp.floor(fi), int)
    j0 = xp.astype(xp.floor(fj), int)
    ai = fi - i0
    aj = fj - j0
    ok = (i0 >= 0) & (i0 < H - 1) & (j0 >= 0) & (j0 < W - 1)
    i0c = xp.clip(i0, 0, H - 2)
    j0c = xp.clip(j0, 0, W - 2)
    v = ((1 - ai) * (1 - aj) * field[i0c, j0c]
         + (1 - ai) * aj * field[i0c, j0c + 1]
         + ai * (1 - aj) * field[i0c + 1, j0c]
         + ai * aj * field[i0c + 1, j0c + 1])
    return v, ok


def _camera_rays(cam: PinholeCamera, L: float, xp):
    """Pixel ray origin (3,) numpy and directions (H, W, 3) in xp."""
    H, W = cam.img_shape
    th = np.deg2rad(cam.zenith_deg)
    ph = np.deg2rad(cam.azimuth_deg)
    z0 = cam.altitude_m
    center = np.array([L / 2.0, L / 2.0, 0.0])
    origin = center + np.array([np.tan(th) * z0 * np.cos(ph),
                                np.tan(th) * z0 * np.sin(ph),
                                z0])
    look = (center - origin) / np.linalg.norm(center - origin)
    world_up = np.array([0.0, 0.0, 1.0])
    right = np.cross(look, world_up)
    right = right / np.linalg.norm(right)
    up = np.cross(right, look)

    half = np.tan(np.deg2rad(cam.hfov_deg))
    ys = np.linspace(-half, half, H) * (H / max(H, W))
    xs = np.linspace(-half, half, W) * (W / max(H, W))
    XS, YS = np.meshgrid(xs, -ys)
    dirs = (look[None, None, :]
            + right[None, None, :] * XS[..., None]
            + up[None, None, :] * YS[..., None])
    dirs = dirs / np.linalg.norm(dirs, axis=-1, keepdims=True)
    return origin, xp.asarray(dirs, dtype=float)


def render_camera_image(eta, dx: float,
                        camera: PinholeCamera = PinholeCamera(),
                        sky=None,
                        slope_x=None,
                        slope_y=None,
                        n_water: float = 1.34,
                        subpixel: SubpixelSlopes | None = None,
                        n_subpixel: int = 0,
                        shadowing: bool = False,
                        water=None,
                        sun_glint: tuple | None = None,
                        foam: Foam | None = None,
                        rng=None):
    """Stokes image (H, W, 4) seen by a pinhole camera looking at the
    surface patch.  `water` (a water.WaterBody or
    scattering.UpwellingRadianceTable) adds the water-leaving radiance;
    `sun_glint = (zen_deg, az_deg, E_sun)` adds the analytic Cox-Munk
    glint (requires `subpixel`); `foam` (a render.Foam) covers the
    steepest-slope pixels with unpolarized Lambertian whitecap radiance.
    Pixels that miss the patch or have no valid sky ray are NaN."""
    xp = xp_of(eta, slope_x, slope_y)
    if sky is None:
        sky = make_unpolarized_sky(1.0)
    rng = adapt_rng(rng, xp)
    eta = xp.asarray(eta, dtype=float)
    Hs, Ws = eta.shape
    L = Ws * dx

    if slope_x is None or slope_y is None:
        gy, gx = xp.gradient(eta, dx)
        slope_x = gx if slope_x is None else slope_x
        slope_y = gy if slope_y is None else slope_y

    origin, dirs = _camera_rays(camera, L, xp)

    # Intersect z = 0, then one fixed-point refinement against eta
    dz = dirs[..., 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (0.0 - origin[2]) / dz
    hit = (t > 0) & xp.isfinite(t)
    X = origin[0] + dirs[..., 0] * t
    Y = origin[1] + dirs[..., 1] * t
    for _ in range(2):
        h, _ = _bilinear(eta, dx, X, Y)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (h - origin[2]) / dz
        X = origin[0] + dirs[..., 0] * t
        Y = origin[1] + dirs[..., 1] * t

    sx0, ok_x = _bilinear(slope_x, dx, X, Y)
    sy0, _ = _bilinear(slope_y, dx, X, Y)
    in_scene = hit & ok_x

    d_out = -dirs
    S = _accumulate_subpixel(d_out, sx0, sy0, sky, n_water,
                             subpixel, n_subpixel, shadowing, rng,
                             water=water)
    if sun_glint is not None:
        if subpixel is None:
            raise ValueError("sun_glint requires subpixel slope statistics")
        S = S + sun_glint_stokes(
            d_out, sx0, sy0, subpixel, *sun_glint, n_water=n_water,
            shadow_sigma=subpixel.sigma_total if shadowing else None)
    S = xp.where(in_scene[..., None], S, xp.nan * S)
    if foam is not None:
        S = _apply_foam(S, sx0, sy0, foam, sky, sun_glint)
    return S
