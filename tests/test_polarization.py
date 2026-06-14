"""Mueller/Stokes machinery checks: rotations, Fresnel, frame chain."""

import numpy as np

from seapol import polarization as pol


def _unpol():
    return np.array([1.0, 0.0, 0.0, 0.0])


def test_rotation_composition():
    a = np.array(0.3)
    b = np.array(-1.1)
    Rab = pol.mueller_rotation(a + b)
    np.testing.assert_allclose(pol.mueller_rotation(a) @ pol.mueller_rotation(b),
                               Rab, atol=1e-12)
    np.testing.assert_allclose(pol.mueller_rotation(np.array(0.0)), np.eye(4),
                               atol=1e-15)
    np.testing.assert_allclose(pol.mueller_rotation(np.array(np.pi)), np.eye(4),
                               atol=1e-12)


def test_meridian_frame_orthonormal():
    rng = np.random.default_rng(0)
    d = pol.normalize(rng.standard_normal((100, 3)))
    v, h = pol.meridian_frame(d)
    np.testing.assert_allclose(np.sum(v * h, axis=-1), 0.0, atol=1e-12)
    np.testing.assert_allclose(np.sum(v * d, axis=-1), 0.0, atol=1e-12)
    np.testing.assert_allclose(np.sum(h * d, axis=-1), 0.0, atol=1e-12)
    # right-handed: v x h = d
    np.testing.assert_allclose(np.cross(v, h), d, atol=1e-12)


def test_fresnel_energy_conservation():
    cos_i = np.cos(np.deg2rad(np.linspace(0.0, 89.0, 50)))
    for n_rel in (1.34, 1.0 / 1.34):
        M_R, M_T, tir = pol.fresnel_mueller(cos_i, n_rel)
        tot = M_R[..., 0, 0] + M_T[..., 0, 0]
        np.testing.assert_allclose(tot[~tir], 1.0, atol=1e-12)
        np.testing.assert_allclose(M_R[tir][..., 0, 0], 1.0, atol=1e-12)
        assert np.all(M_T[tir] == 0.0)


def test_brewster_full_polarization():
    n = 1.34
    cos_b = np.cos(pol.brewster_angle(n))
    M_R, _, _ = pol.fresnel_mueller(np.array(cos_b), n)
    S = M_R @ _unpol()
    assert abs(S[1] / S[0] + 1.0) < 1e-9  # Q/I = -1, fully s-polarized


def test_normal_incidence_mirror():
    n = 1.34
    M_R, _, _ = pol.fresnel_mueller(np.array(1.0), n)
    r0 = ((n - 1.0) / (n + 1.0)) ** 2
    np.testing.assert_allclose(M_R, r0 * np.diag([1.0, 1.0, -1.0, -1.0]),
                               atol=1e-9)


def test_tir_preserves_intensity():
    n_rel = 1.0 / 1.34
    cos_i = np.cos(np.deg2rad(np.array([60.0, 75.0])))  # beyond ~48.3 deg
    M_R, M_T, tir = pol.fresnel_mueller(cos_i, n_rel)
    assert np.all(tir)
    np.testing.assert_allclose(M_R[:, 0, :], [[1, 0, 0, 0]] * 2, atol=1e-12)
    np.testing.assert_allclose(np.abs(np.linalg.det(M_R)), 1.0, atol=1e-9)


def test_reflection_chain_flat_surface_matches_fresnel():
    """For a flat surface the meridian and scattering planes coincide and
    the chained Mueller matrix must equal the bare Fresnel matrix."""
    n_water = 1.34
    for th_deg in [10.0, 40.0, 53.27, 70.0]:
        th = np.deg2rad(th_deg)
        d_out = np.array([np.sin(th), 0.0, np.cos(th)])[None, :]
        n_hat = np.array([0.0, 0.0, 1.0])[None, :]
        M, d_in, valid = pol.reflection_chain(d_out, n_hat, n_water)
        assert valid[0]
        np.testing.assert_allclose(d_in[0], [np.sin(th), 0.0, -np.cos(th)],
                                   atol=1e-12)
        M_ref, _, _ = pol.fresnel_mueller(np.array(np.cos(th)), n_water)
        np.testing.assert_allclose(M[0], M_ref, atol=1e-9)
    # unpolarized input reflects horizontally polarized (Q < 0)
    S = pol.apply_mueller(M[0], _unpol())
    assert S[1] < 0


def test_reflection_chain_azimuth_invariance():
    """Rotating the whole geometry about z must leave the Mueller matrix
    unchanged (meridian frames co-rotate)."""
    n_hat = np.array([[0.0, 0.0, 1.0]])
    th = np.deg2rad(40.0)
    M0, _, _ = pol.reflection_chain(
        np.array([[np.sin(th), 0.0, np.cos(th)]]), n_hat)
    for az in [0.7, 2.0, -1.2]:
        d = np.array([[np.sin(th) * np.cos(az), np.sin(th) * np.sin(az),
                       np.cos(th)]])
        M, _, _ = pol.reflection_chain(d, n_hat)
        np.testing.assert_allclose(M, M0, atol=1e-9)


def test_dolp_aop():
    S = np.array([2.0, 0.5, -0.5, 0.0])
    assert abs(pol.stokes_dolp(S) - np.hypot(0.5, 0.5) / 2.0) < 1e-12
    assert abs(pol.stokes_aop(S) - 0.5 * np.arctan2(-0.5, 0.5)) < 1e-12
