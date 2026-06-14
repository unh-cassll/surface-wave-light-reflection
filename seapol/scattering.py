"""
Near-surface scattering: a plane-parallel polarized Monte Carlo for the
water column under the actual sky, producing the directional upwelling
Stokes radiance field just beneath the surface.

This is what keeps the rendered scene from looking like a bare mirror:
instead of the first-order isotropic L_u = R_w E_d / pi (water.WaterBody),
each facet sees a sub-surface light field L_u(mu_w, phi) that varies
with view direction, water type, sun position, and wavelength band --
brighter toward the sun's in-water beam, darker toward grazing, with the
polarization of in-water scattering attached.

Model and bookkeeping
---------------------
* Horizontally homogeneous, semi-infinite column below a FLAT mean
  interface (the standard ocean-color decoupling: facet-scale tilts act
  at the exit refraction, applied per pixel at render time through
  polarization.transmission_chain).
* Sources: the sky radiance pattern sampled on a quadrature grid
  (weights I cos(theta) dOmega, full polarized Stokes), plus an optional
  direct solar beam sun = (zenith_deg, azimuth_deg, E_sun).  The
  transmitted-and-scattered sun beam is a different light path from the
  surface sun glint, so both may be on together without double counting;
  likewise the sky-reflection term of the renderers handles all photons
  this tracer rejects at the air-side interface.
* Analog polarized transport with the same exact-I importance scheme as
  seapol.montecarlo: interface branches and scattering azimuths are
  importance-sampled with the photon's actual Stokes vector, so every
  photon keeps unit I-weight; energy closes by construction.
* Woodcock (null-collision) tracking through the depth-dependent
  attenuation of the bubble layer; at real collisions photons absorb
  with a(z)/c(z), else scatter off molecules / particulates / bubbles
  in proportion to their local scattering coefficients.
* Every upward crossing of z = 0 deposits the photon Stokes (meridian
  frame of its direction) into the (mu_w, phi) table -- the estimator is
  the steady-state upwelling radiance just below the surface, including
  all orders of internal reflection.  After depositing, the photon
  reflects back down with the polarized Fresnel probability or
  terminates (its transmitted share is exactly what the renderer's
  transmission chain applies when consuming the table).

The tracer dispatches on the backend: pass backend="torch",
device="cuda" to build tables on the GPU.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .backend import adapt_rng, get_xp, to_numpy, xp_of
from .polarization import (frame_rotation_angle, fresnel_mueller,
                           meridian_frame, mueller_rotation, normalize,
                           transmission_chain)
from .skylight import direction_from_angles
from .water import (WaterColumn, ff_phase_mueller, hg_phase_mueller,
                    polarized_scatter_event, rayleigh_phase_mueller,
                    sample_ff_scattering, sample_hg_scattering,
                    sample_rayleigh_scattering)

__all__ = ["UpwellingRadianceTable", "build_upwelling_table",
           "water_leaving_from_table", "save_table", "load_table"]


@dataclass
class UpwellingRadianceTable:
    """Sub-surface upwelling Stokes radiance L_u(mu_w, phi) [per sr,
    in the source's radiance units], meridian frame of the upward
    in-water direction; mu_w = cos(in-water zenith of the upward
    direction), phi = earth-frame azimuth (same convention as the sky
    models).

    S          : (n_mu, n_phi, 4)
    mu_edges   : (n_mu + 1,) in [0, 1]
    phi_edges  : (n_phi + 1,) in [-pi, pi]
    n_water    : refractive index the table was built with
    info       : energy budget and build parameters
    """
    S: np.ndarray
    mu_edges: np.ndarray
    phi_edges: np.ndarray
    n_water: float = 1.34
    info: dict = field(default_factory=dict)

    @property
    def E_u(self) -> float:
        """Upwelling plane irradiance just below the surface."""
        xp = xp_of(self.S)
        mu_c = 0.5 * (self.mu_edges[:-1] + self.mu_edges[1:])
        d_mu = self.mu_edges[1:] - self.mu_edges[:-1]
        d_phi = self.phi_edges[1:] - self.phi_edges[:-1]
        return float(xp.sum(self.S[..., 0] * (mu_c * d_mu)[:, None]
                            * d_phi[None, :]))


# ---------------------------------------------------------------------------
# Source sampling
# ---------------------------------------------------------------------------

def _sky_source_cells(sky_fn, n_zen: int, n_az: int):
    """Quadrature cells of the sky hemisphere: downward propagation
    directions, polarized Stokes (meridian frame), and irradiance
    weights I cos(theta) dOmega."""
    zen = (np.arange(n_zen) + 0.5) * (np.pi / 2) / n_zen
    az = -np.pi + (np.arange(n_az) + 0.5) * 2.0 * np.pi / n_az
    ZEN, AZ = np.meshgrid(zen, az, indexing="ij")
    sky_dirs = direction_from_angles(ZEN, AZ)
    S = to_numpy(sky_fn(sky_dirs)).reshape(-1, 4)
    d_omega = (np.sin(ZEN) * (np.pi / 2 / n_zen)
               * (2.0 * np.pi / n_az)).ravel()
    w = np.maximum(S[:, 0], 0.0) * np.cos(ZEN).ravel() * d_omega
    d_down = -sky_dirs.reshape(-1, 3)
    # normalize Stokes to unit I (weights carry the energy)
    I_safe = np.maximum(S[:, 0], 1e-300)
    S_unit = S / I_safe[:, None]
    return d_down, S_unit, w


def _flat_interface_event(d, S, going_up: bool, n_water: float, rng, xp):
    """Polarized Fresnel branch at the flat z = 0 interface for photons
    of direction d and unit-I Stokes S (meridian frame of d).

    going_up=False: air-side photons -- transmit into the water with the
    polarized probability T_pol (Stokes updated through the full
    rotation/Mueller chain, I renormalized to 1) or terminate (the
    sky-reflection path is the renderers' job).
    going_up=True: water-side photons -- reflect back down with R_pol or
    terminate (transmitted share applied by the table consumer).

    Returns (keep, d_new, S_new).
    """
    n_hat = xp.broadcast_to(xp.asarray([0.0, 0.0, 1.0]), d.shape)
    n_or = xp.where((xp.sum(n_hat * d, axis=-1) > 0)[..., None],
                    -n_hat, n_hat)
    cos_i = xp.clip(-xp.sum(d * n_or, axis=-1), 1e-9, 1.0)
    n_rel = (1.0 / n_water) if going_up else n_water
    M_R, M_T, tir = fresnel_mueller(cos_i, n_rel)

    v_in, h_in = meridian_frame(d)
    s_axis = xp.cross(d, n_or)
    s_n = xp.linalg.norm(s_axis, axis=-1, keepdims=True)
    s_axis = xp.where(s_n > 1e-9, s_axis / xp.maximum(s_n, 1e-300), h_in)
    p_in = xp.cross(s_axis, d)
    a_in = frame_rotation_angle(d, v_in, p_in)
    S_rot = xp.einsum("...ij,...j->...i", mueller_rotation(a_in), S)

    R_pol = (M_R[..., 0, 0] * S_rot[..., 0] + M_R[..., 0, 1] * S_rot[..., 1])
    R_pol = xp.clip(R_pol, 0.0, 1.0)
    R_pol = xp.where(tir, xp.ones_like(R_pol), R_pol)

    if going_up:
        keep = rng.random(d.shape[0]) < R_pol
        p_evt = xp.maximum(R_pol, 1e-12)
        M_evt = M_R / p_evt[..., None, None]
        d_new = d - 2.0 * xp.sum(d * n_or, axis=-1, keepdims=True) * n_or
    else:
        T_pol = 1.0 - R_pol
        keep = rng.random(d.shape[0]) < T_pol
        p_evt = xp.maximum(T_pol, 1e-12)
        M_evt = M_T / p_evt[..., None, None]
        sin_t2 = xp.clip((1.0 - cos_i**2) / n_rel**2, 0.0, 1.0)
        cos_t = xp.sqrt(1.0 - sin_t2)
        d_new = normalize(d / n_rel
                          + (cos_i / n_rel - cos_t)[..., None] * n_or)
    d_new = normalize(d_new)

    v_out, _ = meridian_frame(d_new)
    p_out = xp.cross(s_axis, d_new)
    a_out = frame_rotation_angle(d_new, p_out, v_out)
    S_new = xp.einsum("...ij,...j->...i",
                      mueller_rotation(a_out) @ M_evt, S_rot)
    # renormalize to unit I (importance weight folded into the branch)
    I_new = xp.maximum(S_new[..., 0:1], 1e-300)
    S_new = S_new / I_new
    return keep, d_new, S_new


def _mixed_scatter(d, S, z, col: WaterColumn, rng, xp):
    """One volume-scattering event off the mixed phase function:
    component (molecular / particulate / bubble) chosen in proportion
    to the local scattering coefficients at depth z."""
    K = d.shape[0]
    b_w = col.rayleigh_scattering * xp.ones(K)
    b_p = col.particulate_scattering * xp.ones(K)
    b_b = col.bubble_scattering * xp.exp(z / col.bubble_efold_m) \
        if col.bubble_scattering > 0 else xp.zeros(K)
    b_tot = xp.maximum(b_w + b_p + b_b, 1e-300)
    u = rng.random(K) * b_tot
    is_w = u < b_w
    is_p = (~is_w) & (u < b_w + b_p)
    is_b = ~(is_w | is_p)

    mu = xp.empty(K)
    P = xp.empty((K, 4, 4))
    n_w = int(xp.sum(is_w))
    n_p = int(xp.sum(is_p))
    n_b = int(xp.sum(is_b))
    if n_w:
        mu_w = sample_rayleigh_scattering(n_w, rng, col.depolarization,
                                          xp=xp)
        mu[is_w] = mu_w
        P[is_w] = rayleigh_phase_mueller(mu_w, col.depolarization)
    if n_p:
        if col.particulate_phase == "ff":
            mu_p = sample_ff_scattering(n_p, rng, col.ff_n, col.ff_mu_junge,
                                        xp=xp)
            P[is_p] = ff_phase_mueller(mu_p, col.ff_n, col.ff_mu_junge,
                                       col.particulate_depol)
        else:
            mu_p = sample_hg_scattering(n_p, rng, col.particulate_g, xp=xp)
            P[is_p] = hg_phase_mueller(mu_p, col.particulate_g,
                                       col.particulate_depol)
        mu[is_p] = mu_p
    if n_b:
        mu_b = sample_hg_scattering(n_b, rng, col.bubble_g, xp=xp)
        mu[is_b] = mu_b
        P[is_b] = hg_phase_mueller(mu_b, col.bubble_g, depol=0.9)
    return polarized_scatter_event(d, S, mu, P[:, 0, 0], P[:, 0, 1], P, rng)


# ---------------------------------------------------------------------------
# Table builder
# ---------------------------------------------------------------------------

def build_upwelling_table(sky_fn, water: WaterColumn,
                          sun: tuple | None = None,
                          n_water: float | None = None,
                          n_photons: int = 200_000,
                          n_mu: int = 16, n_phi: int = 24,
                          n_zen_src: int = 24, n_az_src: int = 48,
                          max_events: int = 400,
                          rng=None,
                          backend: str | None = None,
                          device=None, dtype=None
                          ) -> UpwellingRadianceTable:
    """Monte Carlo sub-surface upwelling radiance table for one band.

    sky_fn   : sky model (dirs -> Stokes), same object handed to the
               renderers (radiance units set the table units)
    water    : WaterColumn with SCALAR (single-band) IOPs
    sun      : optional (zenith_deg, azimuth_deg, E_sun) direct beam --
               the transmitted-beam path; keep the renderer's sun_glint
               on for the reflected path
    n_water  : interface index (defaults to water.n_water)

    Returns an UpwellingRadianceTable; info carries the energy budget
    (E_in entering the water, E_u upwelling below the surface,
    absorbed/terminated/unresolved fractions).
    """
    if water.n_bands != 1:
        raise ValueError("build_upwelling_table takes a single-band "
                         "WaterColumn; use WaterColumn.at_band or loop "
                         "bands (see seapol.spectral)")
    xp = get_xp(backend, device, dtype)
    rng = adapt_rng(rng, xp)
    n_w = float(water.n_water if n_water is None else n_water)

    # --- source spectrum: sky quadrature cells + optional sun beam
    d_cells, S_cells, w_cells = _sky_source_cells(sky_fn, n_zen_src,
                                                  n_az_src)
    E_sky = float(w_cells.sum())
    if sun is not None:
        zen_deg, az_deg, E_sun = sun
        mu_s = max(np.cos(np.deg2rad(zen_deg)), 0.0)
        d_sun = -to_numpy(direction_from_angles(np.deg2rad(zen_deg),
                                                np.deg2rad(az_deg)))
        E_beam = float(E_sun) * mu_s
    else:
        E_beam = 0.0
    E_down = E_sky + E_beam
    if E_down <= 0:
        raise ValueError("source has no downwelling irradiance")
    w_photon = E_down / n_photons

    # sample source: sun beam with probability E_beam / E_down, else a
    # sky cell proportional to its irradiance weight
    cdf = np.cumsum(w_cells) / max(E_sky, 1e-300)
    u = to_numpy(rng.random(n_photons))
    from_sun = u < (E_beam / E_down)
    cell = np.searchsorted(cdf, to_numpy(rng.random(n_photons)),
                           side="left").clip(0, len(w_cells) - 1)
    d0 = d_cells[cell]
    S0 = S_cells[cell]
    if sun is not None and from_sun.any():
        d0[from_sun] = d_sun
        S0[from_sun] = np.array([1.0, 0.0, 0.0, 0.0])
    d = xp.asarray(d0, dtype=float)
    S = xp.asarray(S0, dtype=float)

    # --- transmit through the flat interface (analog branch)
    keep, d, S = _flat_interface_event(d, S, going_up=False,
                                       n_water=n_w, rng=rng, xp=xp)
    alive = keep
    z = xp.zeros(n_photons)

    # --- table accumulators
    mu_edges = xp.linspace(0.0, 1.0, n_mu + 1)
    phi_edges = xp.linspace(-np.pi, np.pi, n_phi + 1)
    S_acc = xp.zeros((n_mu * n_phi, 4))
    absorbed = 0.0
    exited = float(n_photons - int(xp.sum(alive)))  # interface rejections

    c_bulk = float(water.absorption + water.rayleigh_scattering
                   + water.particulate_scattering)
    c_max = c_bulk + float(water.bubble_scattering)

    for _ in range(max_events):
        if not bool(xp.any(alive)):
            break
        idx = xp.flatnonzero(alive)
        dz = d[idx, 2]
        zi = z[idx]

        # Woodcock tentative free path vs distance to the surface
        l_col = rng.exponential(1.0 / c_max, idx.shape[0])
        t_surf = xp.where(dz > 1e-12, -zi / xp.maximum(dz, 1e-12),
                          xp.inf * xp.ones_like(zi))
        hits_surface = l_col >= t_surf

        # ---- surface crossings: deposit, then reflect or terminate
        sidx = idx[hits_surface]
        if sidx.shape[0]:
            ds = d[sidx]
            mu_up = xp.clip(ds[:, 2], 0.0, 1.0)
            phi_up = xp.arctan2(ds[:, 1], ds[:, 0])
            bi = xp.clip(xp.astype(mu_up * n_mu, int), 0, n_mu - 1)
            bj = xp.clip(xp.astype((phi_up + np.pi) / (2 * np.pi) * n_phi,
                                   int), 0, n_phi - 1)
            xp.index_add(S_acc, bi * n_phi + bj, S[sidx] * w_photon)
            keep, d_new, S_new = _flat_interface_event(
                ds, S[sidx], going_up=True, n_water=n_w, rng=rng, xp=xp)
            d[sidx] = d_new
            S[sidx] = S_new
            z[sidx] = -1e-9 * xp.ones(sidx.shape[0])
            exited += float(xp.sum(~keep))
            alive[sidx[~keep]] = False

        # ---- volume: advance to the tentative collision, null-test
        vidx = idx[~hits_surface]
        if vidx.shape[0] == 0:
            continue
        z_new = z[vidx] + l_col[~hits_surface] * d[vidx, 2]
        z[vidx] = z_new
        if water.bubble_scattering > 0:
            b_bub = water.bubble_scattering \
                * xp.exp(z_new / water.bubble_efold_m)
        else:
            b_bub = xp.zeros(vidx.shape[0])
        c_loc = c_bulk + b_bub
        real = rng.random(vidx.shape[0]) < (c_loc / c_max)
        ridx = vidx[real]
        if ridx.shape[0] == 0:
            continue

        # absorption roulette at real collisions
        c_r = c_loc[real]
        die = rng.random(ridx.shape[0]) < (water.absorption
                                           / xp.maximum(c_r, 1e-300))
        absorbed += float(xp.sum(die)) * w_photon
        alive[ridx[die]] = False
        scat = ridx[~die]
        if scat.shape[0]:
            d_new, M_step = _mixed_scatter(d[scat], S[scat], z[scat],
                                           water, rng, xp)
            S_new = xp.einsum("...ij,...j->...i", M_step, S[scat])
            I_new = xp.maximum(S_new[..., 0:1], 1e-300)
            S[scat] = S_new / I_new
            d[scat] = d_new

    unresolved = float(xp.sum(alive)) * w_photon

    # crossing estimator: L(mu, phi) = sum w S / (mu_bin dmu dphi)
    mu_c = 0.5 * (mu_edges[:-1] + mu_edges[1:])
    d_mu = mu_edges[1:] - mu_edges[:-1]
    d_phi = phi_edges[1:] - phi_edges[:-1]
    geom = (xp.maximum(mu_c, 1e-6) * d_mu)[:, None] * d_phi[None, :]
    S_tab = S_acc.reshape(n_mu, n_phi, 4) / geom[..., None]

    info = dict(E_down=E_down, E_sky=E_sky, E_beam=E_beam,
                E_u=float(xp.sum(S_acc[:, 0])),
                absorbed=absorbed, exited=exited * w_photon,
                unresolved=unresolved,
                n_photons=n_photons, max_events=max_events,
                water=dict(a=float(water.absorption),
                           b_w=float(water.rayleigh_scattering),
                           b_p=float(water.particulate_scattering),
                           b_bub=float(water.bubble_scattering)),
                sun=sun)
    return UpwellingRadianceTable(S=S_tab, mu_edges=mu_edges,
                                  phi_edges=phi_edges, n_water=n_w,
                                  info=info)


# ---------------------------------------------------------------------------
# Renderer binding
# ---------------------------------------------------------------------------

def _table_lookup(table: UpwellingRadianceTable, mu, phi, xp):
    """Bilinear interpolation of the table Stokes at (mu, phi), periodic
    in phi, clamped in mu."""
    S = xp.asarray(table.S)
    n_mu, n_phi = S.shape[0], S.shape[1]
    mu_c0 = float(to_numpy(table.mu_edges[0]))
    mu_c1 = float(to_numpy(table.mu_edges[-1]))
    d_mu = (mu_c1 - mu_c0) / n_mu
    d_phi = 2.0 * np.pi / n_phi

    fi = (xp.clip(mu, mu_c0, mu_c1) - mu_c0) / d_mu - 0.5
    fj = (phi + np.pi) / d_phi - 0.5
    i0 = xp.astype(xp.floor(fi), int)
    j0 = xp.astype(xp.floor(fj), int)
    ti = (fi - i0)[..., None]
    tj = (fj - j0)[..., None]
    i0c = xp.clip(i0, 0, n_mu - 1)
    i1c = xp.clip(i0 + 1, 0, n_mu - 1)
    j0c = j0 % n_phi
    j1c = (j0 + 1) % n_phi
    return ((1 - ti) * (1 - tj) * S[i0c, j0c]
            + (1 - ti) * tj * S[i0c, j1c]
            + ti * (1 - tj) * S[i1c, j0c]
            + ti * tj * S[i1c, j1c])


def water_leaving_from_table(d_out, n_hat, table: UpwellingRadianceTable,
                             n_water: float | None = None):
    """Water-leaving Stokes contribution (..., 4) along d_out through
    facets with upward normals n_hat: the in-water view direction is
    found by inverse refraction at each facet, the table is sampled
    there, and the polarized transmission chain (with the n^2 radiance
    law) maps the sub-surface Stokes to the camera frame."""
    xp = xp_of(d_out, n_hat)
    n_w = float(table.n_water if n_water is None else n_water)
    M, d_water, valid = transmission_chain(d_out, n_hat, n_w)
    mu_w = xp.clip(d_water[..., 2], 0.0, 1.0)
    phi_w = xp.arctan2(d_water[..., 1], d_water[..., 0])
    S_u = _table_lookup(table, mu_w, phi_w, xp)
    S = xp.einsum("...ij,...j->...i", M, S_u)
    return xp.where(valid[..., None], S, xp.zeros_like(S))


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def save_table(path, table: UpwellingRadianceTable):
    """Save a table to .npz (numpy arrays; info as a pickle-free dict of
    scalars where possible)."""
    np.savez(path, S=to_numpy(table.S),
             mu_edges=to_numpy(table.mu_edges),
             phi_edges=to_numpy(table.phi_edges),
             n_water=table.n_water,
             E_down=table.info.get("E_down", np.nan),
             E_u=table.info.get("E_u", np.nan))


def load_table(path) -> UpwellingRadianceTable:
    z = np.load(path)
    info = dict(E_down=float(z["E_down"]), E_u=float(z["E_u"]))
    return UpwellingRadianceTable(S=z["S"], mu_edges=z["mu_edges"],
                                  phi_edges=z["phi_edges"],
                                  n_water=float(z["n_water"]), info=info)
