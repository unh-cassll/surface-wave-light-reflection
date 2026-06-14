"""
Stokes/Mueller machinery for polarized reflection at the air-water interface.

Conventions (Mobley 2015):
    * Stokes vectors are 4-component (I, Q, U, V), shape (..., 4).
    * A beam with unit propagation direction d carries a reference frame
      (e_par, e_perp, d), right-handed: e_par x e_perp = d.
    * The meridian frame of d uses e_perp = unit(z x d) ("horizontal") and
      e_par = e_perp x d ("vertical", in the plane containing d and zenith).
      +Q means E parallel to the meridian plane (vertical polarization).
    * Mueller matrices act on Stokes vectors in the (p, s) basis of the
      scattering plane: e_par = p (in-plane), e_perp = s (normal to plane).
    * Vectors have trailing axis 3; Mueller matrices trailing axes (4, 4).

All functions dispatch on their array inputs (numpy or torch tensors)
via seapol.backend.
"""

from __future__ import annotations

import numpy as np

from .backend import xp_of

Z_HAT = np.array([0.0, 0.0, 1.0])
Y_HAT = np.array([0.0, 1.0, 0.0])


def normalize(v, eps: float = 1e-300):
    """Unit vectors along the trailing axis; zero vectors stay zero."""
    xp = xp_of(v)
    n = xp.linalg.norm(v, axis=-1, keepdims=True)
    return v / xp.maximum(n, eps)


def meridian_frame(d):
    """(e_par, e_perp) of the meridian frame for propagation direction d.
    Vertical d falls back to e_perp = +y."""
    xp = xp_of(d)
    d = xp.asarray(d, dtype=float)
    z_hat = xp.broadcast_to(xp.asarray(Z_HAT, dtype=float), d.shape)
    y_hat = xp.broadcast_to(xp.asarray(Y_HAT, dtype=float), d.shape)
    h = xp.cross(z_hat, d)
    n = xp.linalg.norm(h, axis=-1, keepdims=True)
    h = xp.where(n > 1e-12, h / xp.maximum(n, 1e-300), y_hat)
    v = xp.cross(h, d)
    return v, h


def frame_rotation_angle(d, e_from, e_to):
    """Right-handed rotation angle about d taking basis vector e_from to
    e_to (both perpendicular to d)."""
    xp = xp_of(d, e_from, e_to)
    c = xp.sum(e_from * e_to, axis=-1)
    s = xp.sum(xp.cross(e_from, e_to) * d, axis=-1)
    return xp.arctan2(s, xp.clip(c, -1.0, 1.0))


def mueller_rotation(alpha):
    """Frame-rotation Mueller matrix R(alpha), shape (..., 4, 4):
    components of S in a basis rotated by +alpha (right-handed about d)."""
    xp = xp_of(alpha)
    alpha = xp.asarray(alpha, dtype=float)
    c = xp.cos(2.0 * alpha)
    s = xp.sin(2.0 * alpha)
    R = xp.zeros(alpha.shape + (4, 4))
    R[..., 0, 0] = 1.0
    R[..., 1, 1] = c
    R[..., 1, 2] = s
    R[..., 2, 1] = -s
    R[..., 2, 2] = c
    R[..., 3, 3] = 1.0
    return R


def fresnel_mueller(cos_i, n_rel):
    """Fresnel reflection and transmission Mueller matrices in the (p, s)
    scattering-plane basis.

    cos_i : cosine of incidence angle w.r.t. the facet normal (> 0)
    n_rel : refractive index of the transmitting medium relative to the
            incident one (water/air = n for air-side incidence, 1/n for
            water-side incidence); broadcastable against cos_i.

    Returns (M_R, M_T, tir) with shapes (..., 4, 4), (..., 4, 4), (...,).
    For total internal reflection M_R is a pure retarder with the TIR
    relative phase delta = delta_p - delta_s and M_T = 0.
    """
    xp = xp_of(cos_i, n_rel)
    cos_i = xp.clip(xp.asarray(cos_i, dtype=float), 1e-9, 1.0)
    n_rel = xp.broadcast_to(xp.asarray(n_rel, dtype=float), cos_i.shape)
    sin_i2 = 1.0 - cos_i**2
    sin_t2 = sin_i2 / n_rel**2
    tir = sin_t2 > 1.0
    cos_t = xp.sqrt(xp.clip(1.0 - sin_t2, 0.0, None))

    rp = (n_rel * cos_i - cos_t) / (n_rel * cos_i + cos_t)
    rs = (cos_i - n_rel * cos_t) / (cos_i + n_rel * cos_t)
    tp = 2.0 * cos_i / (n_rel * cos_i + cos_t)
    ts = 2.0 * cos_i / (cos_i + n_rel * cos_t)

    Rp, Rs = rp**2, rs**2
    f = n_rel * cos_t / cos_i
    Tp, Ts = f * tp**2, f * ts**2

    M_R = xp.zeros(cos_i.shape + (4, 4))
    M_R[..., 0, 0] = 0.5 * (Rp + Rs)
    M_R[..., 0, 1] = 0.5 * (Rp - Rs)
    M_R[..., 1, 0] = 0.5 * (Rp - Rs)
    M_R[..., 1, 1] = 0.5 * (Rp + Rs)
    M_R[..., 2, 2] = rp * rs
    M_R[..., 3, 3] = rp * rs

    M_T = xp.zeros_like(M_R)
    M_T[..., 0, 0] = 0.5 * (Tp + Ts)
    M_T[..., 0, 1] = 0.5 * (Tp - Ts)
    M_T[..., 1, 0] = 0.5 * (Tp - Ts)
    M_T[..., 1, 1] = 0.5 * (Tp + Ts)
    M_T[..., 2, 2] = f * tp * ts
    M_T[..., 3, 3] = f * tp * ts

    if xp.any(tir):
        # TIR relative phase between p and s reflections
        root = xp.sqrt(xp.clip(sin_i2 - n_rel**2, 0.0, None))
        with np.errstate(divide="ignore", invalid="ignore"):
            d_s = 2.0 * xp.arctan2(root, cos_i)
            d_p = 2.0 * xp.arctan2(root, n_rel**2 * cos_i)
        delta = d_p - d_s
        cd, sd = xp.cos(delta), xp.sin(delta)
        M_tir = xp.zeros_like(M_R)
        M_tir[..., 0, 0] = 1.0
        M_tir[..., 1, 1] = 1.0
        M_tir[..., 2, 2] = cd
        M_tir[..., 2, 3] = sd
        M_tir[..., 3, 2] = -sd
        M_tir[..., 3, 3] = cd
        sel = tir[..., None, None]
        M_R = xp.where(sel, M_tir, M_R)
        M_T = xp.where(sel, xp.zeros_like(M_T), M_T)

    return M_R, M_T, tir


def reflection_chain(d_out, n_hat, n_water: float = 1.34):
    """Mueller matrix mapping sky Stokes (meridian frame of the incident
    beam) to camera Stokes (meridian frame of the outgoing beam) for one
    specular reflection off a facet.

    d_out : (..., 3) unit propagation direction surface -> camera
    n_hat : (..., 3) unit upward facet normal

    Returns (M, d_in, valid):
        M     : (..., 4, 4) total Mueller matrix R(a_out) M_R R(a_in)
        d_in  : (..., 3) incident propagation direction (sky -> surface)
        valid : (...,) True where the incident ray comes from above and
                the facet faces the camera
    """
    xp = xp_of(d_out, n_hat)
    d_out = normalize(xp.asarray(d_out, dtype=float))
    n_hat = normalize(xp.asarray(n_hat, dtype=float))

    r_dot_n = xp.sum(d_out * n_hat, axis=-1)
    d_in = d_out - 2.0 * r_dot_n[..., None] * n_hat
    valid = (d_in[..., 2] < 0.0) & (r_dot_n > 0.0)

    # Scattering-plane axis s = unit(d_in x n); degenerate (normal
    # incidence) falls back to the incident meridian e_perp (alpha = 0).
    v_in, h_in = meridian_frame(d_in)
    s_axis = xp.cross(d_in, n_hat)
    s_norm = xp.linalg.norm(s_axis, axis=-1, keepdims=True)
    s_axis = xp.where(s_norm > 1e-9, s_axis / xp.maximum(s_norm, 1e-300),
                      h_in)

    p_in = xp.cross(s_axis, d_in)
    a_in = frame_rotation_angle(d_in, v_in, p_in)

    cos_i = xp.clip(-xp.sum(d_in * n_hat, axis=-1), 0.0, 1.0)
    M_R, _, _ = fresnel_mueller(cos_i, n_water)

    v_out, _ = meridian_frame(d_out)
    p_out = xp.cross(s_axis, d_out)
    a_out = frame_rotation_angle(d_out, p_out, v_out)

    M = mueller_rotation(a_out) @ M_R @ mueller_rotation(a_in)
    return M, d_in, valid


def transmission_chain(d_out, n_hat, n_water: float = 1.34):
    """Mueller matrix mapping in-water upwelling Stokes (meridian frame
    of the underwater beam) to camera Stokes (meridian frame of the
    outgoing beam) for refraction through one facet, including the n^2
    radiance law:

        S_cam = (1 / n^2) R(a_out) M_T(cos_w, 1/n) R(a_in) S_water.

    d_out : (..., 3) unit propagation direction surface -> camera
    n_hat : (..., 3) unit upward facet normal

    Returns (M, d_water, valid):
        M       : (..., 4, 4) total Mueller matrix (n^2 law included)
        d_water : (..., 3) underwater propagation direction (up toward
                  the facet) whose refraction exits along d_out
        valid   : (...,) True where the facet faces the camera (the
                  underwater direction is always inside the Snell cone)
    """
    xp = xp_of(d_out, n_hat)
    d_out = normalize(xp.asarray(d_out, dtype=float))
    n_hat = normalize(xp.asarray(n_hat, dtype=float))

    cos_air = xp.sum(d_out * n_hat, axis=-1)
    valid = cos_air > 1e-6
    cos_air_s = xp.clip(cos_air, 1e-6, 1.0)

    # Underwater direction: tangential component / n, Snell cosine up
    tang = d_out - cos_air_s[..., None] * n_hat
    sin_w2 = xp.clip((1.0 - cos_air_s**2) / n_water**2, 0.0, 1.0)
    cos_w = xp.sqrt(1.0 - sin_w2)
    d_water = normalize(tang / n_water + cos_w[..., None] * n_hat)

    v_in, h_in = meridian_frame(d_water)
    s_axis = xp.cross(d_water, n_hat)
    s_norm = xp.linalg.norm(s_axis, axis=-1, keepdims=True)
    s_axis = xp.where(s_norm > 1e-9, s_axis / xp.maximum(s_norm, 1e-300),
                      h_in)
    p_in = xp.cross(s_axis, d_water)
    a_in = frame_rotation_angle(d_water, v_in, p_in)

    _, M_T, _ = fresnel_mueller(cos_w, 1.0 / n_water)

    v_out, _ = meridian_frame(d_out)
    p_out = xp.cross(s_axis, d_out)
    a_out = frame_rotation_angle(d_out, p_out, v_out)

    M = (mueller_rotation(a_out) @ M_T @ mueller_rotation(a_in)) \
        / n_water**2
    return M, d_water, valid


def apply_mueller(M, S):
    """S_out = M @ S for trailing-axis Stokes vectors."""
    xp = xp_of(M, S)
    return xp.einsum("...ij,...j->...i", M, S)


def brewster_angle(n_rel: float) -> float:
    """Brewster angle [rad] for relative index n_rel."""
    return float(np.arctan(n_rel))


def stokes_dolp(S):
    """Degree of linear polarization sqrt(Q^2 + U^2) / I."""
    xp = xp_of(S)
    with np.errstate(divide="ignore", invalid="ignore"):
        return xp.hypot(S[..., 1], S[..., 2]) / xp.where(S[..., 0] > 0,
                                                         S[..., 0], xp.nan)


def stokes_aop(S):
    """Angle of polarization 0.5 atan2(U, Q) [rad], measured from the
    meridian (vertical) toward horizontal."""
    xp = xp_of(S)
    return 0.5 * xp.arctan2(S[..., 2], S[..., 1])
