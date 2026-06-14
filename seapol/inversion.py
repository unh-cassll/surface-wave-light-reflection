"""
Polarimetric slope sensing: reconstruct facet slopes and wave height
from single-reflection polarized imagery.

Under unpolarized illumination (overcast sky) the reflected Stokes
vector of a facet encodes its orientation:

    * DoLP depends only on the incidence angle through the Fresnel
      reflection matrix, |M10/M00|(theta_i) -- monotone from 0 at
      normal incidence to 1 at Brewster (the working branch);
    * AoP gives the orientation of the scattering plane: reflected
      light is s-dominated, so the polarization direction IS the
      s-axis, and the plane of incidence is perpendicular to it.

The facet normal is then d_out rotated by theta_i within the plane of
incidence (the upward candidate of the two), and slopes follow from
n = (-eta_x, -eta_y, 1)/|.|.  Height is recovered (to its unobservable
mean) by least-squares spectral integration of the slope fields.

Limits: sub-Brewster branch only (DoLP is two-valued across Brewster;
keep the camera incidence moderate), single reflection, unpolarized
sky, no foam/water-leaving contamination.
"""

from __future__ import annotations

import numpy as np

from .backend import xp_of
from .polarization import (apply_mueller, brewster_angle, fresnel_mueller,
                           meridian_frame, normalize, reflection_chain,
                           stokes_aop, stokes_dolp)

__all__ = ["fresnel_dolp_curve", "slopes_from_stokes",
           "slopes_from_stokes_polarized", "height_from_slopes"]


def fresnel_dolp_curve(n_water: float = 1.34, n_pts: int = 2048
                       ) -> tuple[np.ndarray, np.ndarray]:
    """(theta_i, DoLP) samples of the unpolarized-illumination Fresnel
    reflection DoLP on the monotone sub-Brewster branch."""
    th = np.linspace(0.0, brewster_angle(n_water), n_pts)
    M_R, _, _ = fresnel_mueller(np.cos(th), n_water)
    dolp = np.abs(M_R[:, 1, 0]) / np.maximum(M_R[:, 0, 0], 1e-300)
    return th, dolp


def slopes_from_stokes(S, d_out, n_water: float = 1.34):
    """Facet slopes (sx, sy) from per-pixel reflected Stokes vectors S
    (..., 4) in the meridian frame of the view directions d_out
    (..., 3, surface -> camera), assuming a single reflection of
    unpolarized illumination.

    Returns (sx, sy, valid).  DoLP beyond the Brewster maximum clamps
    to the Brewster angle; the normal-sign ambiguity (AoP is defined
    mod pi) is resolved by choosing the upward candidate."""
    xp = xp_of(S, d_out)
    S = xp.asarray(S, dtype=float)
    d_out = normalize(xp.asarray(d_out, dtype=float))
    th_grid_np, dolp_grid_np = fresnel_dolp_curve(n_water)
    th_grid = xp.asarray(th_grid_np, dtype=float)
    dolp_grid = xp.asarray(dolp_grid_np, dtype=float)

    dolp = stokes_dolp(S)
    alpha = stokes_aop(S)
    th_i = xp.interp(xp.clip(dolp, 0.0, float(dolp_grid_np[-1])),
                     dolp_grid, th_grid)

    # polarization direction = s-axis; plane of incidence is normal to it
    v_out, h_out = meridian_frame(d_out)
    s_dir = xp.cos(alpha)[..., None] * v_out + xp.sin(alpha)[..., None] * h_out
    p_dir = normalize(xp.cross(s_dir, d_out))

    ct = xp.cos(th_i)[..., None]
    st = xp.sin(th_i)[..., None]
    n1 = ct * d_out + st * p_dir
    n2 = ct * d_out - st * p_dir
    n = xp.where((n1[..., 2:3] >= n2[..., 2:3]), n1, n2)

    nz = n[..., 2]
    valid = xp.isfinite(dolp) & xp.isfinite(alpha) & (nz > 1e-6)
    nz_s = xp.where(valid, nz, xp.ones_like(nz))
    sx = -n[..., 0] / nz_s
    sy = -n[..., 1] / nz_s
    nan = xp.nan * nz_s
    return xp.where(valid, sx, nan), xp.where(valid, sy, nan), valid


def slopes_from_stokes_polarized(S, d_out, sky, n_water: float = 1.34,
                                 sx0=None, sy0=None, n_iter: int = 20,
                                 damping: float = 1e-3,
                                 step_clip: float = 0.4,
                                 resid_tol: float = 0.02):
    """Facet slopes from reflected Stokes under a KNOWN, possibly
    polarized sky, by per-pixel damped Gauss-Newton (Levenberg-Marquardt)
    inversion of the single-reflection forward model.

    `slopes_from_stokes` reads the incidence angle from DoLP and so
    assumes *unpolarized* illumination on the sub-Brewster branch.  Under
    a clear or partly-cloudy sky the incident light is already polarized,
    and the reflected Q/U mix the sky polarization with the Fresnel
    rotation, so DoLP no longer maps to the facet angle.  This routine
    instead inverts the full Mueller chain against the supplied sky model
    -- the same forward model as `render_facet_stokes` -- which makes it
    valid for polarized skies and off-Brewster viewing (e.g. a DoFP
    camera at 30 deg incidence under a clear sky).

    S       : (..., 4) observed reflected Stokes
    d_out   : (..., 3) view directions, surface -> camera
    sky     : sky model callable (dirs -> Stokes), as used by the renderers
    sx0, sy0: optional initial slope (default flat); a macro-slope or the
              `slopes_from_stokes` estimate make good warm starts

    The fit matches the observed Q/I and U/I (intensity-normalized, so
    insensitive to the absolute sky radiance).  Returns (sx, sy, valid);
    pixels that do not converge below `resid_tol` are NaN.  Two facet
    orientations can occasionally reproduce the same (Q/I, U/I); a
    warm start from the resolved macro-slope resolves the branch.
    """
    xp = xp_of(S, d_out)
    S = xp.asarray(S, dtype=float)
    d_out = normalize(xp.asarray(d_out, dtype=float))
    I = S[..., 0]
    good_I = I > 0
    I_s = xp.where(good_I, I, xp.ones_like(I))
    q_obs = xp.where(good_I, S[..., 1] / I_s, xp.nan)
    u_obs = xp.where(good_I, S[..., 2] / I_s, xp.nan)

    shape = d_out.shape[:-1]
    ones = xp.ones(shape)
    sx = xp.zeros(shape) if sx0 is None else xp.asarray(sx0, dtype=float) * ones
    sy = xp.zeros(shape) if sy0 is None else xp.asarray(sy0, dtype=float) * ones

    def predict(sxx, syy):
        n_hat = normalize(xp.stack([-sxx, -syy, xp.ones_like(sxx)], axis=-1))
        M, d_in, ok = reflection_chain(d_out, n_hat, n_water)
        Sp = apply_mueller(M, sky(-d_in))
        Ip = Sp[..., 0]
        good = ok & (Ip > 0)
        Ip_s = xp.where(good, Ip, xp.ones_like(Ip))
        return Sp[..., 1] / Ip_s, Sp[..., 2] / Ip_s, good

    h = 1e-3
    for _ in range(n_iter):
        q0, u0, _ = predict(sx, sy)
        qx, ux, _ = predict(sx + h, sy)
        qy, uy, _ = predict(sx, sy + h)
        # Jacobian rows (dq, du) / d(sx, sy)
        ja, jb = (qx - q0) / h, (qy - q0) / h
        jc, jd = (ux - u0) / h, (uy - u0) / h
        rq, ru = q0 - q_obs, u0 - u_obs
        # LM normal equations (J^T J + lambda I) step = -J^T r, 2x2 closed form
        a00 = ja * ja + jc * jc + damping
        a01 = ja * jb + jc * jd
        a11 = jb * jb + jd * jd + damping
        g0 = ja * rq + jc * ru
        g1 = jb * rq + jd * ru
        det = a00 * a11 - a01 * a01
        det = xp.where(xp.abs(det) > 1e-30, det, xp.ones_like(det))
        dsx = -(a11 * g0 - a01 * g1) / det
        dsy = -(-a01 * g0 + a00 * g1) / det
        dsx = xp.clip(xp.where(xp.isfinite(dsx), dsx, 0.0), -step_clip,
                      step_clip)
        dsy = xp.clip(xp.where(xp.isfinite(dsy), dsy, 0.0), -step_clip,
                      step_clip)
        sx = sx + dsx
        sy = sy + dsy

    qf, uf, okf = predict(sx, sy)
    resid = xp.hypot(qf - q_obs, uf - u_obs)
    valid = good_I & okf & xp.isfinite(resid) & (resid < resid_tol)
    nan = xp.nan * ones
    return xp.where(valid, sx, nan), xp.where(valid, sy, nan), valid


def height_from_slopes(sx, sy, dx: float):
    """Least-squares spectral integration of periodic slope fields:
    eta_hat = -i (kx sx_hat + ky sy_hat) / k^2 (mean fixed to zero).
    NaNs are filled with the field mean before the transform."""
    xp = xp_of(sx, sy)
    sx = xp.asarray(sx, dtype=float)
    sy = xp.asarray(sy, dtype=float)
    if sx.shape != sy.shape or sx.ndim != 2 or sx.shape[0] != sx.shape[1]:
        raise ValueError("height_from_slopes requires square (N, N) "
                         f"slope fields, got {tuple(sx.shape)} and "
                         f"{tuple(sy.shape)}")
    if not (bool(xp.any(xp.isfinite(sx))) and bool(xp.any(xp.isfinite(sy)))):
        raise ValueError("height_from_slopes: slope fields are all-NaN "
                         "(no valid pixels from the inversion?)")
    sx = xp.where(xp.isfinite(sx), sx, xp.nanmean(sx))
    sy = xp.where(xp.isfinite(sy), sy, xp.nanmean(sy))
    N = sx.shape[0]
    kx = 2.0 * np.pi * xp.fft.fftfreq(N, d=dx)
    KX, KY = xp.meshgrid(kx, kx, indexing="xy")
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0
    Sx = xp.fft.fft2(sx)
    Sy = xp.fft.fft2(sy)
    eta_hat = -1j * (KX * Sx + KY * Sy) / K2
    eta_hat[0, 0] = 0.0
    return xp.fft.ifft2(eta_hat).real
