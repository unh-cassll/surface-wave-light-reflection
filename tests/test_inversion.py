"""Polarimetric slope/height reconstruction round trip."""

import numpy as np

from seapol import (CameraGeometry, generate_sea_surface,
                    height_from_slopes, make_overcast_sky,
                    render_facet_stokes, slopes_from_stokes,
                    slopes_from_stokes_polarized)
from seapol.polarization import normalize


def _render_and_invert(surf, incidence_deg=30.0):
    dx = surf.info["dx"]
    cam = CameraGeometry(incidence_deg=incidence_deg, azimuth_deg=0.0,
                         height_m=100.0)
    S = render_facet_stokes(surf.eta, dx, camera=cam,
                            sky=make_overcast_sky(1.0),
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=None, shadowing=False)
    # per-facet view directions (z = 0 approximation, as in the field)
    N = surf.eta.shape[0]
    xg = np.arange(N) * dx
    X, Y = np.meshgrid(xg, xg)
    P = np.stack([X, Y, np.zeros_like(X)], axis=-1)
    cpos = cam.position(((N - 1) * dx / 2.0, (N - 1) * dx / 2.0))
    d_out = normalize(cpos[None, None, :] - P)
    sx_r, sy_r, valid = slopes_from_stokes(S, d_out)
    return sx_r, sy_r, valid


def test_slope_reconstruction_round_trip():
    """Render with the exact single-facet chain, invert, and recover the
    slope fields nearly pixel-for-pixel (sub-Brewster geometry)."""
    surf = generate_sea_surface(L=32.0, N=64, U10=5.0,
                                rng=np.random.default_rng(0))
    sx_r, sy_r, valid = _render_and_invert(surf)
    assert valid.mean() > 0.95
    for rec, tru in ((sx_r, surf.slope_x), (sy_r, surf.slope_y)):
        m = valid & np.isfinite(rec)
        c = np.corrcoef(rec[m], tru[m])[0, 1]
        assert c > 0.99, c
    # MSS closes
    mss_t = surf.slope_x[valid].var() + surf.slope_y[valid].var()
    mss_r = np.nanvar(sx_r[valid]) + np.nanvar(sy_r[valid])
    np.testing.assert_allclose(mss_r, mss_t, rtol=0.05)


def test_height_reconstruction_statistics():
    """Spectral integration of reconstructed slopes recovers the height
    field: Hs within a few percent, high pointwise correlation."""
    surf = generate_sea_surface(L=32.0, N=64, U10=5.0,
                                rng=np.random.default_rng(1))
    sx_r, sy_r, _ = _render_and_invert(surf)
    eta_r = height_from_slopes(sx_r, sy_r, surf.info["dx"])
    eta_t = surf.eta - surf.eta.mean()
    c = np.corrcoef(eta_r.ravel(), eta_t.ravel())[0, 1]
    assert c > 0.99, c
    np.testing.assert_allclose(4.0 * eta_r.std(), 4.0 * eta_t.std(),
                               rtol=0.05)


def _polarized_setup(N=96, U10=6.0, seed=0, incidence_deg=30.0):
    from seapol import make_clear_sky
    surf = generate_sea_surface(8.0, N, U10, rng=np.random.default_rng(seed))
    dx = surf.info["dx"]
    cam = CameraGeometry(incidence_deg=incidence_deg, azimuth_deg=0.0,
                         height_m=80.0)
    sky = make_clear_sky(50.0, 90.0, 1.0, turbidity=0.05)
    S = render_facet_stokes(surf.eta, dx, camera=cam, sky=sky,
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=None)
    xg = np.arange(N) * dx
    X, Y = np.meshgrid(xg, xg)
    P = np.stack([X, Y, surf.eta], axis=-1)
    cpos = cam.position(((N - 1) * dx / 2.0, (N - 1) * dx / 2.0))
    d_out = normalize(cpos[None, None, :] - P)
    return surf, S, d_out, sky


def test_polarized_inversion_recovers_slopes():
    """Under a clear POLARIZED sky at 30 deg incidence -- where the
    DoLP-based unpolarized inversion is biased -- the full-Mueller
    Gauss-Newton inversion recovers the slope fields (it inverts the
    exact forward model), and beats the unpolarized model."""
    surf, S, d_out, sky = _polarized_setup()
    sxp, syp, vp = slopes_from_stokes_polarized(S, d_out, sky)
    mp = vp & np.isfinite(sxp)
    assert mp.mean() > 0.97
    for rec, tru in ((sxp, surf.slope_x), (syp, surf.slope_y)):
        assert np.corrcoef(rec[mp], tru[mp])[0, 1] > 0.999
    # strictly better than the unpolarized-assumption inversion here
    sxu, _, vu = slopes_from_stokes(S, d_out)
    mu = vu & np.isfinite(sxu)
    c_pol = np.corrcoef(sxp[mp], surf.slope_x[mp])[0, 1]
    c_unp = np.corrcoef(sxu[mu], surf.slope_x[mu])[0, 1]
    assert c_pol > c_unp


def test_polarized_inversion_noise_robust():
    """Graceful degradation under per-channel sensor noise (a few percent
    of mean radiance): slope RMSE stays small, correlation high."""
    surf, S, d_out, sky = _polarized_setup(seed=3)
    rng = np.random.default_rng(7)
    S_noisy = S + rng.normal(0.0, 0.01 * np.nanmean(S[..., 0]), S.shape)
    sxp, syp, vp = slopes_from_stokes_polarized(S_noisy, d_out, sky)
    m = vp & np.isfinite(sxp)
    assert m.mean() > 0.9
    assert np.corrcoef(sxp[m], surf.slope_x[m])[0, 1] > 0.99
    assert np.sqrt(np.nanmean((sxp - surf.slope_x)[m] ** 2)) < 0.02


def test_fresnel_dolp_curve_monotone():
    """The sub-Brewster DoLP branch is strictly monotone (invertible),
    spans [0, ~1], and peaks at the Brewster angle."""
    from seapol import fresnel_dolp_curve
    from seapol.polarization import brewster_angle
    th, dolp = fresnel_dolp_curve()
    assert np.all(np.diff(dolp) > 0)
    assert dolp[0] == 0.0
    np.testing.assert_allclose(dolp[-1], 1.0, atol=1e-6)
    np.testing.assert_allclose(th[-1], brewster_angle(1.34), rtol=1e-12)


def test_height_from_slopes_guards():
    import pytest
    from seapol import height_from_slopes
    with pytest.raises(ValueError, match="square"):
        height_from_slopes(np.zeros((8, 10)), np.zeros((8, 10)), 0.1)
    with pytest.raises(ValueError, match="all-NaN"):
        height_from_slopes(np.full((8, 8), np.nan),
                           np.full((8, 8), np.nan), 0.1)


def test_height_from_slopes_on_truth():
    """Integration of the exact spectral slope fields returns eta up to
    the self-conjugate Nyquist bins, whose slope content survives the
    .real projection only as cosine components (sub-mm here)."""
    surf = generate_sea_surface(L=16.0, N=64, U10=6.0,
                                rng=np.random.default_rng(2))
    eta_r = height_from_slopes(surf.slope_x, surf.slope_y,
                               surf.info["dx"])
    eta_t = surf.eta - surf.eta.mean()
    np.testing.assert_allclose(eta_r, eta_t, atol=2e-3)
    assert np.corrcoef(eta_r.ravel(), eta_t.ravel())[0, 1] > 0.9999
