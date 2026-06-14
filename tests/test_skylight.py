"""Rayleigh skylight checks: DoLP pattern and U-sign symmetry."""

import numpy as np

from seapol import skylight as sky
from seapol.polarization import stokes_dolp


def test_dolp_zero_at_sun_max_at_90():
    sun_zen = np.deg2rad(30.0)
    sun_dir = sky.direction_from_angles(sun_zen, 0.0)
    S_at_sun = sky.rayleigh_sky_stokes(sun_dir, sun_dir)
    assert stokes_dolp(S_at_sun) < 1e-12

    # 90 deg scattering: zenith 60 deg, azimuth 90 deg from a 30-deg sun
    v = sky.direction_from_angles(np.deg2rad(60.0), np.deg2rad(90.0))
    assert abs(np.dot(v, sun_dir)) < 0.51  # sanity: large scattering angle
    S = sky.rayleigh_sky_stokes(v, sun_dir)
    g = np.arccos(np.clip(np.dot(v, sun_dir), -1, 1))
    expected = sky.dolp_max(30.0) * np.sin(g) ** 2 / (1 + np.cos(g) ** 2)
    np.testing.assert_allclose(stokes_dolp(S), expected, rtol=1e-9)


def test_intensity_and_v():
    zen = np.deg2rad(np.linspace(5, 85, 9))
    az = np.deg2rad(np.linspace(0, 350, 9))
    S = sky.rayleigh_sky_stokes_angles(zen, az, np.deg2rad(40.0), 0.0,
                                       I_sky=2.5)
    np.testing.assert_allclose(S[..., 0], 2.5)
    np.testing.assert_allclose(S[..., 3], 0.0)
    assert np.all(stokes_dolp(S) <= sky.dolp_max(40.0) + 1e-12)


def test_u_antisymmetric_about_sun_meridian():
    """Mirroring the sky point across the sun's vertical plane flips U
    and preserves Q."""
    sun_zen, sun_az = np.deg2rad(35.0), 0.0
    zen = np.deg2rad(50.0)
    for az_deg in [20.0, 60.0, 130.0]:
        Sp = sky.rayleigh_sky_stokes_angles(zen, np.deg2rad(az_deg),
                                            sun_zen, sun_az)
        Sm = sky.rayleigh_sky_stokes_angles(zen, np.deg2rad(-az_deg),
                                            sun_zen, sun_az)
        np.testing.assert_allclose(Sp[..., 1], Sm[..., 1], atol=1e-12)
        np.testing.assert_allclose(Sp[..., 2], -Sm[..., 2], atol=1e-12)
        assert abs(Sp[..., 2]) > 1e-6  # U is actually nonzero off-plane
