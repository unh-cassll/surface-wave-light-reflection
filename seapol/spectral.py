"""
Color-aware (multi-band) rendering: per-band skies, water IOPs, and
refractive index orchestrated over the scalar renderers.

The scalar pipeline renders one band-integrated Stokes image; this
module loops wavelength bands with consistent physics per band --

    * sky radiance spectra (blue Rayleigh lambda^-4 + aerosol haze,
      neutral overcast, Beer-reddened low sun; seapol.skylight),
    * water-type spectral IOPs and per-band upwelling-radiance tables
      (seapol.water / seapol.scattering),
    * n_water(lambda) dispersion in both the Fresnel chain and the
      water-leaving transmission (Quan & Fry 1995),
    * the same surface realization, sub-pixel ensemble draws, and
      cloud field in every band (band loops re-seed the RNG so spectral
      noise stays luminance noise, not chroma noise) --

and stacks (H, W, B, 4) spectral Stokes images, displayable through
seapol.color.  Cost scales linearly with the number of bands; the
per-band scattering tables are the expensive part and should be built
once and cached (build_spectral_tables + scattering.save_table).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .render import (PinholeCamera, CameraGeometry, make_clear_sky,
                     make_overcast_sky, make_partly_cloudy_sky,
                     make_rayleigh_sky, render_camera_image,
                     render_facet_stokes)
from .scattering import UpwellingRadianceTable, build_upwelling_table
from .skylight import (overcast_spectrum, sky_radiance_spectrum,
                       sun_beam_spectrum)
from .water import WaterBody, WaterColumn, WaterType

__all__ = ["SpectralBands", "spectral_sky_factories",
           "spectral_sun_irradiance", "build_spectral_tables",
           "render_facet_stokes_spectral", "render_camera_image_spectral"]


@dataclass
class SpectralBands:
    """Wavelength bands for color-aware rendering."""
    wavelengths_nm: np.ndarray

    @classmethod
    def rgb(cls) -> "SpectralBands":
        """Nominal three-band mode (450, 550, 650 nm)."""
        return cls(np.array([450.0, 550.0, 650.0]))

    @classmethod
    def visible(cls, n: int = 13) -> "SpectralBands":
        """n bands spanning 400-700 nm (quantitative CIE mode)."""
        return cls(np.linspace(400.0, 700.0, n))

    @property
    def n_bands(self) -> int:
        return len(self.wavelengths_nm)


def spectral_sun_irradiance(bands: SpectralBands, sun_zenith_deg: float,
                            E_sun_550: float, turbidity: float = 0.0
                            ) -> np.ndarray:
    """Per-band direct-beam irradiance: E_sun_550 scaled by the
    Beer-attenuated solar spectrum (reddens the low sun)."""
    return E_sun_550 * sun_beam_spectrum(bands.wavelengths_nm,
                                         sun_zenith_deg, turbidity)


def spectral_sky_factories(bands: SpectralBands, kind: str = "clear",
                           sun_zenith_deg: float = 50.0,
                           sun_azimuth_deg: float = 0.0,
                           I_sky_550: float = 1.0,
                           turbidity: float = 0.0,
                           cloud_fraction: float = 0.4,
                           cloud_brightness: float = 3.0,
                           cloud_seed: int = 0,
                           I_sun_550: float = 0.0) -> list:
    """Per-band sky callables with nominal sky colors.

    kind:
        "clear" / "rayleigh" : blue molecular sky (radiance ~
            E_sun tau_R ~ lambda^-4) whitened by the aerosol component
            as turbidity rises; "clear" also accepts an unpolarized
            sun glow I_sun_550 (Beer-reddened per band)
        "overcast"           : spectrally neutral cloud-transmitted sun
        "partly_cloudy"      : blue clear-sky pattern with white clouds
            (identical cloud field in every band)
    """
    wl = bands.wavelengths_nm
    skies = []
    if kind in ("clear", "rayleigh"):
        I_b = I_sky_550 * sky_radiance_spectrum(wl, turbidity)
        E_b = I_sun_550 * sun_beam_spectrum(wl, sun_zenith_deg, turbidity)
        for i in range(bands.n_bands):
            if kind == "rayleigh":
                skies.append(make_rayleigh_sky(sun_zenith_deg,
                                               sun_azimuth_deg,
                                               float(I_b[i])))
            else:
                skies.append(make_clear_sky(sun_zenith_deg, sun_azimuth_deg,
                                            float(I_b[i]), turbidity,
                                            I_sun=float(E_b[i])))
    elif kind == "overcast":
        I_b = I_sky_550 * overcast_spectrum(wl)
        skies = [make_overcast_sky(float(I)) for I in I_b]
    elif kind == "partly_cloudy":
        I_b = I_sky_550 * sky_radiance_spectrum(wl, turbidity)
        # clouds are white: absolute cloud radiance follows the neutral
        # solar spectrum, so the per-band brightness ratio compensates
        # for the blue clear-sky normalization
        C_b = I_sky_550 * cloud_brightness * overcast_spectrum(wl)
        for i in range(bands.n_bands):
            skies.append(make_partly_cloudy_sky(
                sun_zenith_deg, sun_azimuth_deg, float(I_b[i]),
                cloud_fraction=cloud_fraction,
                cloud_brightness=float(C_b[i] / max(I_b[i], 1e-30)),
                rng=np.random.default_rng(cloud_seed)))
    else:
        raise ValueError(f"unknown sky kind: {kind!r}")
    return skies


def build_spectral_tables(skies: list, water, bands: SpectralBands,
                          sun: tuple | None = None,
                          turbidity: float = 0.0,
                          n_photons: int = 200_000,
                          rng_seed: int = 0,
                          verbose: bool = False,
                          **table_kwargs) -> list[UpwellingRadianceTable]:
    """Per-band upwelling-radiance tables for a WaterType or a spectral
    WaterColumn.  sun = (zenith_deg, azimuth_deg, E_sun_550) is scaled
    per band by the attenuated solar spectrum.  This is the expensive
    color-mode step (one Monte Carlo per band) -- cache the result
    (scattering.save_table / load_table) when iterating on a scene."""
    if isinstance(water, WaterType):
        column = water.column(bands.wavelengths_nm)
    elif isinstance(water, WaterColumn):
        column = water
    else:
        raise TypeError("water must be a WaterType or spectral WaterColumn")
    if column.n_bands != bands.n_bands:
        raise ValueError("water column bands do not match render bands")

    tables = []
    for i in range(bands.n_bands):
        sun_i = None
        if sun is not None:
            zen, az, E550 = sun
            E_i = float(E550 * sun_beam_spectrum(
                float(bands.wavelengths_nm[i]), zen, turbidity))
            sun_i = (zen, az, E_i)
        if verbose:
            print(f"  table band {bands.wavelengths_nm[i]:.0f} nm ...")
        tables.append(build_upwelling_table(
            skies[i], column.at_band(i), sun=sun_i,
            n_photons=n_photons,
            rng=np.random.default_rng(rng_seed + i), **table_kwargs))
    return tables


def _spectral_loop(render_fn, eta, dx, bands: SpectralBands, skies: list,
                   water, n_water, sun_glint, turbidity, seed,
                   **render_kwargs):
    from .backend import xp_of
    xp = xp_of(eta)
    if isinstance(water, (WaterType, WaterColumn)):
        raise TypeError("pass per-band tables from build_spectral_tables "
                        "(or a WaterBody) as `water`, not raw IOPs -- "
                        "table building is a separate, cacheable step")
    if n_water is None:
        n_water = [t.n_water if isinstance(t, UpwellingRadianceTable)
                   else 1.34 for t in water] if isinstance(water, list) \
            else [1.34] * bands.n_bands
    n_water = np.broadcast_to(np.asarray(n_water, dtype=float),
                              (bands.n_bands,))

    out = None
    for i in range(bands.n_bands):
        water_i = water[i] if isinstance(water, list) else water
        glint_i = None
        if sun_glint is not None:
            zen, az, E550 = sun_glint
            E_i = float(E550 * sun_beam_spectrum(
                float(bands.wavelengths_nm[i]), zen, turbidity))
            glint_i = (zen, az, E_i)
        S_i = render_fn(eta, dx, sky=skies[i], water=water_i,
                        n_water=float(n_water[i]), sun_glint=glint_i,
                        rng=np.random.default_rng(seed),
                        **render_kwargs)
        if out is None:
            out = xp.empty(S_i.shape[:-1] + (bands.n_bands, 4))
        out[..., i, :] = S_i
    return out


def render_facet_stokes_spectral(eta, dx: float, bands: SpectralBands,
                                 skies: list,
                                 camera: CameraGeometry = CameraGeometry(),
                                 water=None,
                                 n_water=None,
                                 sun_glint: tuple | None = None,
                                 turbidity: float = 0.0,
                                 seed: int = 0,
                                 **render_kwargs):
    """Spectral facet rendering: (H, W, B, 4).

    skies     : per-band sky callables (spectral_sky_factories)
    water     : None, a WaterBody, or per-band tables
                (build_spectral_tables)
    n_water   : per-band indices; defaults to the tables' values
    sun_glint : (zenith_deg, azimuth_deg, E_sun_550); the band scaling
                (low-sun reddening) is applied here
    seed      : sub-pixel ensemble seed, identical across bands
    """
    def fn(eta_, dx_, **kw):
        return render_facet_stokes(eta_, dx_, camera=camera, **kw)
    return _spectral_loop(fn, eta, dx, bands, skies, water, n_water,
                          sun_glint, turbidity, seed, **render_kwargs)


def render_camera_image_spectral(eta, dx: float, bands: SpectralBands,
                                 skies: list,
                                 camera: PinholeCamera = PinholeCamera(),
                                 water=None,
                                 n_water=None,
                                 sun_glint: tuple | None = None,
                                 turbidity: float = 0.0,
                                 seed: int = 0,
                                 **render_kwargs):
    """Spectral pinhole-camera rendering: (H, W, B, 4).  Arguments as in
    render_facet_stokes_spectral."""
    def fn(eta_, dx_, **kw):
        return render_camera_image(eta_, dx_, camera=camera, **kw)
    return _spectral_loop(fn, eta, dx, bands, skies, water, n_water,
                          sun_glint, turbidity, seed, **render_kwargs)
