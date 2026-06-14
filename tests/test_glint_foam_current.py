"""Gram-Charlier glint statistics, whitecap foam, and current Doppler."""

import numpy as np

from seapol import (Foam, SubpixelSlopes, generate_sea_surface,
                    make_unpolarized_sky, monahan_coverage,
                    render_facet_stokes, sun_glint_stokes)
from seapol.polarization import stokes_dolp


# ---------------------------------------------------------------------------
# Gram-Charlier subpixel statistics
# ---------------------------------------------------------------------------

def test_gc_factor_moments():
    """The GC-weighted Gaussian ensemble realizes the requested moment
    structure.  Enforcing positivity (clip at 0) plus the far-tail
    taper shrinks the realized moments to ~70-80% of the nominal
    coefficients at field-size values — the series coefficients are
    shape parameters, not exact moments, once the density is made
    valid — so the assertions are sign plus a magnitude band."""
    sub = SubpixelSlopes(sigma_a2=0.02, sigma_c2=0.015,
                         c21=-0.05, c03=-0.25, c04=0.3)
    rng = np.random.default_rng(0)
    n = 2_000_000
    xa = rng.normal(0.0, np.sqrt(sub.sigma_a2), n)
    xc = rng.normal(0.0, np.sqrt(sub.sigma_c2), n)
    w = sub.gc_factor(xa, xc)
    at = xa / np.sqrt(sub.sigma_a2)
    ct = xc / np.sqrt(sub.sigma_c2)
    wm = w.mean()
    np.testing.assert_allclose(wm, 1.0, atol=5e-3)
    skew = np.mean(w * at**3) / wm
    c21_r = np.mean(w * ct**2 * at) / wm
    kurt = np.mean(w * at**4) / wm - 3.0
    assert 0.55 * abs(sub.c03) < abs(skew) < 1.1 * abs(sub.c03)
    assert np.sign(skew) == np.sign(sub.c03)
    assert 0.4 * abs(sub.c21) < abs(c21_r) < 1.2 * abs(sub.c21)
    assert np.sign(c21_r) == np.sign(sub.c21)
    assert 0.3 * sub.c04 < kurt < 1.2 * sub.c04


def test_gc_zero_coefficients_is_gaussian():
    sub = SubpixelSlopes(sigma_a2=0.02, sigma_c2=0.015)
    assert sub.is_gaussian
    rng = np.random.default_rng(1)
    _, _, w = sub.sample_weighted((100,), rng)
    np.testing.assert_array_equal(w, 1.0)


def test_glint_upwind_downwind_asymmetry():
    """Sun at zenith, two views mirrored about it along the wind axis:
    the required specular slopes are +-s, so the Gaussian glint is
    exactly symmetric and only the c03 skewness breaks the tie."""
    common = dict(sigma_a2=0.02, sigma_c2=0.015)
    gauss = SubpixelSlopes(**common)
    skewed = SubpixelSlopes(**common, c03=-0.25)

    th_v = np.deg2rad(25.0)
    d_up = np.array([[np.sin(th_v), 0.0, np.cos(th_v)]])    # toward +x
    d_dn = np.array([[-np.sin(th_v), 0.0, np.cos(th_v)]])   # toward -x
    sx0 = np.zeros(1)

    def glint_I(sub, d_out):
        S = sun_glint_stokes(d_out, sx0, sx0, sub, 0.0, 0.0, 1.0)
        return float(S[0, 0])

    r_gauss = glint_I(gauss, d_up) / glint_I(gauss, d_dn)
    r_skew = glint_I(skewed, d_up) / glint_I(skewed, d_dn)
    np.testing.assert_allclose(r_gauss, 1.0, rtol=1e-9)
    assert abs(r_skew - 1.0) > 0.05, r_skew

    # flipping the skew flips the asymmetry
    flipped = SubpixelSlopes(**common, c03=+0.25)
    r_flip = glint_I(flipped, d_up) / glint_I(flipped, d_dn)
    np.testing.assert_allclose(r_skew * r_flip, 1.0, rtol=1e-6)


def test_from_cox_munk_coefficients():
    sub = SubpixelSlopes.from_cox_munk(10.0, k_cutoff=100.0)
    np.testing.assert_allclose(sub.c03, 0.04 - 0.33, atol=1e-12)
    np.testing.assert_allclose(sub.c21, 0.01 - 0.086, atol=1e-12)
    assert sub.sigma_a2 > 0 and sub.sigma_c2 > 0
    assert not sub.is_gaussian


# ---------------------------------------------------------------------------
# Whitecap foam
# ---------------------------------------------------------------------------

def test_foam_full_coverage_saturates():
    """coverage >= 1 covers every valid facet with unpolarized foam."""
    rng = np.random.default_rng(7)
    surf = generate_sea_surface(L=16.0, N=32, U10=8.0, rng=rng)
    S = render_facet_stokes(surf.eta, surf.info["dx"],
                            sky=make_unpolarized_sky(1.0),
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            foam=Foam(coverage=1.5),
                            rng=np.random.default_rng(8))
    valid = np.isfinite(S[..., 0])
    assert valid.any()
    np.testing.assert_allclose(stokes_dolp(S)[valid], 0.0, atol=1e-12)


def test_monahan_coverage():
    np.testing.assert_allclose(monahan_coverage(10.0),
                               3.84e-6 * 10.0**3.41, rtol=1e-12)
    assert monahan_coverage(5.0) < monahan_coverage(12.0) < 0.1


def test_foam_coverage_and_depolarization():
    """Foam covers the requested fraction of the steepest facets and
    lowers the scene DoLP."""
    rng = np.random.default_rng(2)
    surf = generate_sea_surface(L=32.0, N=128, U10=10.0, rng=rng)
    sub = SubpixelSlopes(surf.info["sigma_a2_cut"],
                         surf.info["sigma_c2_cut"])
    kw = dict(sky=make_unpolarized_sky(1.0), slope_x=surf.slope_x,
              slope_y=surf.slope_y, subpixel=sub, n_subpixel=4)
    S0 = render_facet_stokes(surf.eta, surf.info["dx"],
                             rng=np.random.default_rng(3), **kw)
    S1 = render_facet_stokes(surf.eta, surf.info["dx"],
                             foam=Foam(coverage=0.05),
                             rng=np.random.default_rng(3), **kw)

    mag = np.hypot(surf.slope_x, surf.slope_y)
    changed = np.zeros(S0.shape[:2], dtype=bool)
    both = np.isfinite(S0[..., 0]) & np.isfinite(S1[..., 0])
    changed[both] = ~np.isclose(S0[..., 0][both], S1[..., 0][both])
    frac = changed.sum() / both.sum()
    np.testing.assert_allclose(frac, 0.05, atol=0.01)
    # foam sits on the steepest facets
    assert mag[changed].mean() > 2.0 * mag[both & ~changed].mean()
    # foam is unpolarized: scene-median DoLP drops
    d0 = np.nanmedian(stokes_dolp(S0))
    d1 = np.nanmedian(stokes_dolp(S1))
    assert d1 <= d0
    np.testing.assert_allclose(np.nanmax(stokes_dolp(S1)[changed]), 0.0,
                               atol=1e-12)


# ---------------------------------------------------------------------------
# Current Doppler
# ---------------------------------------------------------------------------

def test_current_doppler_phase():
    """A single spectral component evolves at omega + k.U: measure the
    mode's complex phase rotation between two frames."""
    N, L = 64, 32.0
    dk = 2 * np.pi / L
    j0 = 4                                   # mode (kx = 4 dk, ky = 0)
    kx0 = j0 * dk
    Psi = np.zeros((N, N))
    Psi[0, j0] = 1.0                         # row 0 = ky 0 in fft layout

    from seapol.spectrum import angular_frequency
    t1 = 0.05
    phases = {}
    for label, cur in (("no", None), ("with", (0.8, 0.0))):
        surf = generate_sea_surface(L=L, N=N, U10=5.0, psi_override=Psi,
                                    times=np.array([0.0, t1]),
                                    compute_slopes=False, current=cur,
                                    rng=np.random.default_rng(4))
        Z0 = np.fft.fft2(surf.eta[:, :, 0])[0, j0]
        Z1 = np.fft.fft2(surf.eta[:, :, 1])[0, j0]
        phases[label] = np.angle(Z1 / Z0)

    omega = angular_frequency(np.array(kx0))
    np.testing.assert_allclose(phases["no"], -omega * t1, atol=1e-6)
    np.testing.assert_allclose(phases["with"], -(omega + kx0 * 0.8) * t1,
                               atol=1e-6)


def test_current_preserves_variance():
    """The Doppler factor is pure phase: each frame's variance matches
    the no-current run with the same draw exactly."""
    times = np.array([0.0, 1.0])
    s_cur = generate_sea_surface(L=32.0, N=64, U10=7.0, times=times,
                                 current=(0.5, -0.3),
                                 rng=np.random.default_rng(5))
    s_no = generate_sea_surface(L=32.0, N=64, U10=7.0, times=times,
                                rng=np.random.default_rng(5))
    np.testing.assert_allclose(s_cur.eta[:, :, 1].var(),
                               s_no.eta[:, :, 1].var(), rtol=1e-9)
    assert s_cur.info["current"] == (0.5, -0.3)
