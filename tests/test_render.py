"""Renderer checks against analytic flat-surface expectations."""

import numpy as np

from seapol import (CameraGeometry, PinholeCamera, SubpixelSlopes,
                    generate_sea_surface, make_rayleigh_sky,
                    make_unpolarized_sky, render_camera_image,
                    render_facet_stokes, render_facet_stokes_stack)
from seapol.polarization import brewster_angle, fresnel_mueller, stokes_dolp


def _flat(N=9, dx=0.01):
    eta = np.zeros((N, N))
    z = np.zeros((N, N))
    return eta, z, dx


def test_facet_flat_brewster():
    """Flat water at Brewster incidence under an unpolarized sky: the
    center facet is fully horizontally polarized (Q/I = -1)."""
    eta, z, dx = _flat()
    cam = CameraGeometry(incidence_deg=np.degrees(brewster_angle(1.34)),
                         azimuth_deg=0.0, height_m=500.0)
    S = render_facet_stokes(eta, dx, camera=cam, sky=make_unpolarized_sky(),
                            slope_x=z, slope_y=z)
    c = S.shape[0] // 2
    assert abs(S[c, c, 1] / S[c, c, 0] + 1.0) < 1e-6
    assert abs(S[c, c, 2]) < 1e-9


def test_facet_flat_matches_fresnel():
    eta, z, dx = _flat()
    th_deg = 40.0
    cam = CameraGeometry(incidence_deg=th_deg, azimuth_deg=0.0,
                         height_m=1000.0)
    S = render_facet_stokes(eta, dx, camera=cam, sky=make_unpolarized_sky(),
                            slope_x=z, slope_y=z)
    c = S.shape[0] // 2
    M_ref, _, _ = fresnel_mueller(np.array(np.cos(np.deg2rad(th_deg))), 1.34)
    np.testing.assert_allclose(S[c, c], M_ref @ [1, 0, 0, 0], atol=1e-6)


def test_facet_sun_in_view_plane_gives_zero_u():
    """Sun and camera in the same vertical plane: symmetry forces U = 0
    at the center facet."""
    eta, z, dx = _flat()
    cam = CameraGeometry(incidence_deg=40.0, azimuth_deg=0.0, height_m=500.0)
    sky = make_rayleigh_sky(sun_zenith_deg=35.0, sun_azimuth_deg=0.0)
    S = render_facet_stokes(eta, dx, camera=cam, sky=sky,
                            slope_x=z, slope_y=z)
    c = S.shape[0] // 2
    assert abs(S[c, c, 2]) < 1e-9
    assert np.isfinite(S[c, c, 0]) and S[c, c, 0] > 0


def test_facet_subpixel_ensemble_runs():
    rng = np.random.default_rng(0)
    surf = generate_sea_surface(L=16.0, N=32, U10=5.0, rng=rng)
    sub = SubpixelSlopes(sigma_a2=0.01, sigma_c2=0.008)
    S = render_facet_stokes(surf.eta, surf.info["dx"],
                            camera=CameraGeometry(incidence_deg=30.0,
                                                  height_m=200.0),
                            sky=make_rayleigh_sky(45.0, 90.0),
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=sub, n_subpixel=16,
                            shadowing=True, rng=rng)
    finite = np.isfinite(S[..., 0])
    assert finite.mean() > 0.9
    dolp = stokes_dolp(S[finite])
    assert np.nanmax(dolp) <= 1.0 + 1e-9


def test_facet_stack_shapes():
    rng = np.random.default_rng(1)
    surf = generate_sea_surface(L=8.0, N=16, U10=4.0,
                                times=np.array([0.0, 0.25]), rng=rng)
    S = render_facet_stokes_stack(surf.eta, surf.info["dx"],
                                  slope_x=surf.slope_x,
                                  slope_y=surf.slope_y,
                                  sky=make_unpolarized_sky())
    assert S.shape == (16, 16, 4, 2)
    assert np.isfinite(S[8, 8, 0, :]).all()


def test_camera_image_flat_center():
    eta, z, dx = _flat(N=64, dx=0.5)
    cam = PinholeCamera(altitude_m=300.0, zenith_deg=40.0, azimuth_deg=0.0,
                        hfov_deg=2.0, img_shape=(32, 32))
    S = render_camera_image(eta, dx, camera=cam, sky=make_unpolarized_sky(),
                            slope_x=z, slope_y=z)
    c = 16
    assert np.isfinite(S[c, c]).all()
    M_ref, _, _ = fresnel_mueller(np.array(np.cos(np.deg2rad(40.0))), 1.34)
    np.testing.assert_allclose(S[c, c, 0], (M_ref @ [1, 0, 0, 0])[0],
                               rtol=0.02)
    assert S[c, c, 1] < 0
