"""
Color-aware mode tests: nominal sky/solar spectra, water-type spectral
IOPs, refractive-index dispersion, CIE conversion, and the spectral
renderers (band consistency, water color signatures).
"""

import numpy as np
import pytest

from seapol.color import cie_cmf, srgb_encode, stokes_bands_to_rgb
from seapol.render import CameraGeometry, SubpixelSlopes
from seapol.skylight import (overcast_spectrum, rayleigh_optical_depth,
                             sky_radiance_spectrum, solar_spectrum,
                             sun_beam_spectrum)
from seapol.spectral import (SpectralBands, build_spectral_tables,
                             render_facet_stokes_spectral,
                             spectral_sky_factories,
                             spectral_sun_irradiance)
from seapol.surface import generate_sea_surface
from seapol.water import (WATER_TYPES, pure_water_absorption,
                          water_refractive_index)


def test_sky_spectra_orderings():
    wl = np.array([450.0, 550.0, 650.0])
    blue = sky_radiance_spectrum(wl)
    assert blue[0] > 2.0 and blue[2] < 0.6            # blue clear sky
    hazy = sky_radiance_spectrum(wl, turbidity=1.0)
    assert hazy[0] < blue[0]                          # haze whitens
    low_sun = sun_beam_spectrum(wl, 80.0)
    assert low_sun[0] < 0.8 < 1.0 <= low_sun[2] + 0.1  # reddened beam
    high_sun = sun_beam_spectrum(wl, 10.0)
    assert high_sun[0] > low_sun[0]
    oc = overcast_spectrum(wl)
    assert 0.8 < oc.min() and oc.max() < 1.2          # neutral
    np.testing.assert_allclose(solar_spectrum(550.0), 1.0)
    assert rayleigh_optical_depth(400.0) > 3.0 * rayleigh_optical_depth(600.0)


def test_water_spectral_iops():
    assert pure_water_absorption(650.0) > 30.0 * pure_water_absorption(450.0)
    wl = np.array([450.0, 550.0, 650.0])
    clear = WATER_TYPES["clear"].column(wl)
    prod = WATER_TYPES["productive_case1"].column(wl)
    coastal = WATER_TYPES["coastal_case2"].column(wl)
    a_clear = np.asarray(clear.absorption)
    a_coastal = np.asarray(coastal.absorption)
    # CDOM eats blue: coastal blue/green absorption ratio exceeds clear's
    assert (a_coastal[0] / a_coastal[1]) > (a_clear[0] / a_clear[1])
    assert np.asarray(prod.particulate_scattering)[1] \
        > 10.0 * np.asarray(clear.particulate_scattering)[1]


def test_refractive_index_dispersion():
    n450 = water_refractive_index(450.0)
    n650 = water_refractive_index(650.0)
    assert n450 > n650
    assert 1.33 < n650 < n450 < 1.35
    # fresh cold water differs from warm salt
    assert water_refractive_index(550.0, 0.0, 5.0) \
        < water_refractive_index(550.0, 35.0, 5.0)


def test_cie_conversion_basics():
    wl = np.linspace(400.0, 700.0, 16)
    cmf = cie_cmf(wl)
    assert cmf.shape == (16, 3) and cmf.min() >= 0.0
    # flat spectrum -> near-neutral sRGB (equal-energy white)
    S = np.zeros((2, 2, 16, 4))
    S[..., 0] = 1.0
    rgb = stokes_bands_to_rgb(S, wl, exposure=None)
    m = rgb.reshape(-1, 3).mean(axis=0)
    assert np.all(np.abs(m / m.mean() - 1.0) < 0.3)
    # blue-weighted spectrum -> blue-dominant pixel
    S2 = np.zeros((1, 1, 16, 4))
    S2[..., 0] = np.exp(-0.5 * ((wl - 450.0) / 30.0) ** 2)
    rgb2 = stokes_bands_to_rgb(S2, wl)
    assert rgb2[0, 0, 2] > rgb2[0, 0, 0]


def test_direct_mode_band_order():
    wl = np.array([450.0, 550.0, 650.0])
    S = np.zeros((1, 1, 3, 4))
    S[0, 0, 2, 0] = 1.0          # only the 650 nm band
    rgb = stokes_bands_to_rgb(S, wl, mode="direct", exposure=1.0)
    assert rgb[0, 0, 0] > 0.9 and rgb[0, 0, 1] == 0.0 and rgb[0, 0, 2] == 0.0
    with pytest.raises(ValueError):
        stokes_bands_to_rgb(np.zeros((1, 1, 4, 4)),
                            np.array([400, 500, 600, 700.0]),
                            mode="direct")


def test_srgb_encode_range():
    x = np.linspace(-0.2, 1.5, 50)
    y = srgb_encode(x)
    assert y.min() >= 0.0 and y.max() <= 1.0


def test_spectral_render_band_consistency():
    """With identical per-band skies, no water, and a shared seed, every
    band image is identical: the band loop adds no spurious chroma."""
    from seapol.render import make_overcast_sky
    surf = generate_sea_surface(20.0, 32, 6.0, rng=np.random.default_rng(0))
    bands = SpectralBands(np.array([450.0, 550.0, 650.0]))
    skies = [make_overcast_sky(1.0)] * 3
    sp = SubpixelSlopes.from_cox_munk(6.0, 100.0)
    S = render_facet_stokes_spectral(surf.eta, surf.info["dx"], bands,
                                     skies,
                                     camera=CameraGeometry(40.0, 0.0, 50.0),
                                     subpixel=sp, n_subpixel=4, seed=3)
    np.testing.assert_allclose(S[..., 0, :], S[..., 1, :], atol=1e-12)
    np.testing.assert_allclose(S[..., 0, :], S[..., 2, :], atol=1e-12)


def test_spectral_render_green_window():
    """Overcast sky over coastal water: the water-leaving S0 peaks in
    the green band (CDOM absorbs blue, water absorbs red)."""
    surf = generate_sea_surface(20.0, 32, 6.0, rng=np.random.default_rng(1))
    bands = SpectralBands.rgb()
    skies = spectral_sky_factories(bands, "overcast")
    tabs = build_spectral_tables(skies, WATER_TYPES["coastal_case2"],
                                 bands, n_photons=30000, rng_seed=2)
    sp = SubpixelSlopes.from_cox_munk(6.0, 100.0)
    S = render_facet_stokes_spectral(
        surf.eta, surf.info["dx"], bands, skies,
        camera=CameraGeometry(30.0, 180.0, 50.0), water=tabs,
        subpixel=sp, n_subpixel=4, seed=4)
    assert S.shape == (32, 32, 3, 4)
    # water-leaving part: subtract a no-water render of the same seed
    S_dry = render_facet_stokes_spectral(
        surf.eta, surf.info["dx"], bands, skies,
        camera=CameraGeometry(30.0, 180.0, 50.0),
        subpixel=sp, n_subpixel=4, seed=4)
    wl_part = np.nanmean(S[..., 0] - S_dry[..., 0], axis=(0, 1))
    assert wl_part[1] > wl_part[2]            # green > red
    assert wl_part[1] > 0.5 * wl_part[0]      # not blue-dominated either


def test_partly_cloudy_same_cloud_field_across_bands():
    bands = SpectralBands.rgb()
    pc = spectral_sky_factories(bands, "partly_cloudy",
                                sun_zenith_deg=50.0, cloud_fraction=0.3,
                                cloud_seed=5)
    clear = spectral_sky_factories(bands, "clear", sun_zenith_deg=50.0)
    zen = np.linspace(0.05, np.pi / 2 - 0.05, 40)
    az = np.linspace(-np.pi, np.pi, 80, endpoint=False)
    ZEN, AZ = np.meshgrid(zen, az, indexing="ij")
    from seapol.skylight import direction_from_angles
    dirs = direction_from_angles(ZEN, AZ)
    masks = []
    for i in range(3):
        ratio = pc[i](dirs)[..., 0] / np.maximum(clear[i](dirs)[..., 0],
                                                 1e-12)
        masks.append(ratio > 1.2)
    assert np.array_equal(masks[0], masks[1])
    assert np.array_equal(masks[0], masks[2])


def test_spectral_sun_irradiance_scaling():
    bands = SpectralBands.rgb()
    E = spectral_sun_irradiance(bands, 75.0, 10.0)
    assert E[0] < E[1] < 10.0 * 1.1
    assert E[2] > E[0]


def test_spectral_loop_rejects_raw_iops():
    surf = generate_sea_surface(10.0, 16, 5.0, rng=np.random.default_rng(2))
    bands = SpectralBands.rgb()
    skies = spectral_sky_factories(bands, "overcast")
    with pytest.raises(TypeError):
        render_facet_stokes_spectral(surf.eta, surf.info["dx"], bands,
                                     skies, water=WATER_TYPES["clear"])
