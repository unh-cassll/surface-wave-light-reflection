"""Sky dome models and first-order water-leaving radiance."""

import numpy as np

from seapol import (CameraGeometry, WaterBody, direction_from_angles,
                    downwelling_irradiance, make_clear_sky,
                    make_overcast_sky, make_partly_cloudy_sky,
                    make_unpolarized_sky, render_facet_stokes,
                    water_leaving_stokes)
from seapol.polarization import stokes_dolp


def _hemisphere(n=24):
    zen = np.linspace(0.05, np.pi / 2 - 0.05, n)
    az = np.linspace(-np.pi, np.pi, 2 * n, endpoint=False)
    Z, A = np.meshgrid(zen, az, indexing="ij")
    return direction_from_angles(Z, A)


def test_overcast_unpolarized_gradation():
    sky = make_overcast_sky(I_zenith=3.0)
    dirs = _hemisphere()
    S = sky(dirs)
    assert np.all(S[..., 1:] == 0.0)
    zenith = sky(np.array([0.0, 0.0, 1.0]))
    horizon = sky(np.array([1.0, 0.0, 1e-9]))
    assert abs(zenith[0] - 3.0) < 1e-9
    assert abs(horizon[0] / zenith[0] - 1.0 / 3.0) < 1e-6


def test_clear_sky_turbidity_scales_dolp():
    dirs = _hemisphere()
    S0 = make_clear_sky(40.0, 0.0, turbidity=0.0)(dirs)
    S5 = make_clear_sky(40.0, 0.0, turbidity=0.5)(dirs)
    d0 = np.nanmax(stokes_dolp(S0))
    d5 = np.nanmax(stokes_dolp(S5))
    assert abs(d5 / d0 - 0.5) < 1e-6
    # sun disk adds unpolarized radiance near the sun only
    S_sun = make_clear_sky(40.0, 0.0, I_sun=50.0,
                           sun_halfwidth_deg=2.0)(dirs)
    excess = S_sun[..., 0] - S0[..., 0]
    sun_dir = direction_from_angles(np.deg2rad(40.0), 0.0)
    near = np.sum(dirs * sun_dir, axis=-1) > np.cos(np.deg2rad(5.0))
    assert excess[near].max() > 10.0
    assert excess[~near].max() < 2.0


def test_partly_cloudy_fraction():
    sky = make_partly_cloudy_sky(40.0, 0.0, cloud_fraction=0.4,
                                 cloud_brightness=4.0,
                                 rng=np.random.default_rng(0))
    dirs = _hemisphere(48)
    S = sky(dirs)
    cloudy = S[..., 0] > 2.0  # bright unpolarized patches
    assert 0.2 < cloudy.mean() < 0.6
    assert np.all(np.abs(S[cloudy][:, 1:3]) < 1e-12)
    assert np.any(np.abs(S[~cloudy][:, 1]) > 1e-3)  # clear part polarized


def test_downwelling_irradiance_uniform_sky():
    # E_d = pi * I for a uniform unpolarized sky
    E = downwelling_irradiance(make_unpolarized_sky(2.0))
    assert abs(E / (2.0 * np.pi) - 1.0) < 5e-3


def test_water_leaving_magnitude_and_polarization():
    d_out = np.array([[0.0, 0.0, 1.0]])           # nadir view
    n_hat = np.array([[0.0, 0.0, 1.0]])
    water = WaterBody(case=1)
    S = water_leaving_stokes(d_out, n_hat, water, E_d=np.pi, n_water=1.34)
    # nadir: L_w = R_w * T(0) / n^2, T(0) = 1 - ((n-1)/(n+1))^2
    T0 = 1.0 - ((1.34 - 1.0) / (1.34 + 1.0)) ** 2
    np.testing.assert_allclose(S[0, 0], 0.02 * T0 / 1.34**2, rtol=1e-9)
    np.testing.assert_allclose(S[0, 1:], 0.0, atol=1e-12)  # nadir unpolarized
    # oblique: slight transmission polarization, nonzero I
    th = np.deg2rad(60.0)
    d_obl = np.array([[np.sin(th), 0.0, np.cos(th)]])
    S_obl = water_leaving_stokes(d_obl, n_hat, water, E_d=np.pi)
    assert S_obl[0, 0] > 0
    assert abs(S_obl[0, 1]) > 0  # p/s transmission asymmetry
    assert WaterBody(case=2).R_w > WaterBody(case=1).R_w


def test_render_with_water_lowers_dolp():
    """Adding water-leaving radiance raises I and lowers reflected DoLP."""
    eta = np.zeros((9, 9))
    z = np.zeros((9, 9))
    cam = CameraGeometry(incidence_deg=53.0, height_m=500.0)
    kw = dict(camera=cam, slope_x=z, slope_y=z)
    S_dry = render_facet_stokes(eta, 0.01, **kw)
    S_wet = render_facet_stokes(eta, 0.01, water=WaterBody(case=2), **kw)
    c = 4
    assert S_wet[c, c, 0] > S_dry[c, c, 0]
    assert stokes_dolp(S_wet[c, c]) < stokes_dolp(S_dry[c, c])