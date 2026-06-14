"""Analytic Cox-Munk sun glint."""

import numpy as np
import pytest

from seapol import (CameraGeometry, SubpixelSlopes, make_unpolarized_sky,
                    render_facet_stokes, sun_glint_stokes)
from seapol.polarization import fresnel_mueller


SUB = SubpixelSlopes(sigma_a2=0.01, sigma_c2=0.008)


def _dout(theta_deg, az_deg=0.0):
    th = np.deg2rad(theta_deg)
    az = np.deg2rad(az_deg)
    return np.array([[np.sin(th) * np.cos(az), np.sin(th) * np.sin(az),
                      np.cos(th)]])


def test_glint_peak_matches_closed_form():
    """Specular geometry on a flat facet: glint radiance equals
    E p(0) M_R(theta) / (4 cos(theta_v))."""
    th = 40.0
    z = np.zeros(1)
    E = 100.0
    S = sun_glint_stokes(_dout(th), z, z, SUB, th, 180.0, E)
    p0 = 1.0 / (2 * np.pi * np.sqrt(SUB.sigma_a2 * SUB.sigma_c2))
    M_R, _, _ = fresnel_mueller(np.array(np.cos(np.deg2rad(th))), 1.34)
    expected = E * p0 * M_R[0, 0] / (4 * np.cos(np.deg2rad(th)))
    np.testing.assert_allclose(S[0, 0], expected, rtol=1e-6)
    assert S[0, 1] < 0          # glint is horizontally polarized


def test_glint_decays_off_specular():
    th = 40.0
    z = np.zeros(1)
    S_peak = sun_glint_stokes(_dout(th), z, z, SUB, th, 180.0, 1.0)
    far = np.array([0.6])       # ~6 sigma from specular
    S_far = sun_glint_stokes(_dout(th), far, z, SUB, th, 180.0, 1.0)
    assert S_far[0, 0] < 1e-4 * S_peak[0, 0]


def test_glint_zero_below_horizon():
    z = np.zeros(1)
    S = sun_glint_stokes(_dout(40.0), z, z, SUB, 120.0, 180.0, 1.0)
    assert np.all(S == 0.0)


def test_renderer_glint_integration():
    eta = np.zeros((9, 9))
    z = np.zeros((9, 9))
    cam = CameraGeometry(incidence_deg=40.0, azimuth_deg=0.0,
                        height_m=500.0)
    kw = dict(camera=cam, sky=make_unpolarized_sky(0.0), slope_x=z,
              slope_y=z, subpixel=SUB, n_subpixel=1,
              rng=np.random.default_rng(0))
    # camera sits at azimuth 0 looking back: specular sun azimuth ~180
    S = render_facet_stokes(eta, 0.01, sun_glint=(40.0, 180.0, 10.0), **kw)
    c = 4
    assert S[c, c, 0] > 1.0      # strong glint at the specular point
    with pytest.raises(ValueError):
        render_facet_stokes(eta, 0.01, camera=cam,
                            sky=make_unpolarized_sky(0.0),
                            slope_x=z, slope_y=z,
                            sun_glint=(40.0, 180.0, 10.0))
