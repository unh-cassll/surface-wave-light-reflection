"""
Forward polarized Monte Carlo ray tracer for a wind-roughened air-water
interface, after Mobley (2015), Appl. Opt. 54, 4828-4849.

Design:
    * Batch-parallel stepping: all live rays advance one bounce at a time
      through vectorized kernels (numpy or torch; with torch tensors the
      whole trace runs on the tensors' device).
    * Heightmap DDA (Amanatides-Woo) grid march for ray/surface
      intersection: only the two triangles of the current (x, y) cell are
      tested, O(sqrt(M)) per ray instead of O(M) triangles.
    * True periodic boundaries: the mesh tiles the full FFT period (the
      last cell row/column wraps to the first), and ray positions shift by
      +/- L when their cell index wraps, so intersection tests remain
      consistent after boundary crossings.
    * Mueller-matrix path accumulation: each ray carries the 4x4 product of
      per-bounce event matrices (meridian frame -> scattering plane ->
      Fresnel -> meridian frame of the new direction).  One trace yields
      the full effective Mueller matrix; no multi-Stokes-input solve.
    * Exact polarization frames at every bounce (no z-meridian shortcut).
    * Branch probabilities and scattering azimuths are importance-sampled
      with the path's unpolarized-launch Stokes column (column 0 of the
      accumulated Mueller), so the I-path weight of every ray stays
      exactly 1 through arbitrarily many polarizing events and energy
      closes to machine precision.  For an unpolarized launch this is
      the analog polarized MC of Mobley (2015); the other Mueller columns
      remain unbiased importance-sampled estimates (up to measure-zero
      degeneracies, e.g. exactly Brewster hits on fully polarized paths).
    * Per-ray medium tracking (air/water) instead of inferring the side
      from the sign of d_z.
    * Optional in-water medium (seapol.water.WaterOptics): exponential
      free paths in c = a + b against the surface-hit distance,
      absorption roulette with survival omega_0 = b / c, and polarized
      Rayleigh-with-depolarization volume scattering.
"""

from __future__ import annotations

import numpy as np

from .backend import adapt_rng, xp_of
from .polarization import (frame_rotation_angle, fresnel_mueller,
                           meridian_frame, mueller_rotation, normalize)
from .water import (WaterOptics, polarized_scatter_event,
                    rayleigh_phase_mueller, sample_rayleigh_scattering)


# ---------------------------------------------------------------------------
# Periodic triangle mesh
# ---------------------------------------------------------------------------

def build_cell_tables(eta, dx: float) -> dict:
    """Per-cell corner vertices of the periodic facet mesh.

    The N x N elevation grid defines N x N cells covering the full FFT
    period [0, L)^2 with L = N dx; the last row/column of cells wraps to
    the first grid line.  Corner ordering [v00, v01, v11, v10] gives
    triangle 1 = (v00, v01, v11), triangle 2 = (v00, v11, v10).
    """
    xp = xp_of(eta)
    eta = xp.asarray(eta, dtype=float)
    N = eta.shape[0]
    if eta.shape[0] != eta.shape[1]:
        raise ValueError("eta must be square")
    idx = xp.arange(N)
    nxt = (idx + 1) % N
    z00 = eta
    z01 = eta[:, nxt]
    z10 = eta[nxt, :]
    z11 = eta[nxt[:, None], nxt[None, :]]
    xf = xp.astype(idx, float) * dx
    ones_col = xp.ones((N, 1))
    ones_row = xp.ones((1, N))
    x0 = xf[None, :] * ones_col
    x1 = (xf + dx)[None, :] * ones_col
    y0 = xf[:, None] * ones_row
    y1 = (xf + dx)[:, None] * ones_row
    v00 = xp.stack([x0, y0, z00], axis=-1)
    v01 = xp.stack([x1, y0, z01], axis=-1)
    v11 = xp.stack([x1, y1, z11], axis=-1)
    v10 = xp.stack([x0, y1, z10], axis=-1)
    V = xp.stack([v00, v01, v11, v10], axis=-2)
    return dict(V=V, N=N, dx=dx, L=N * dx,
                z_min=float(xp.min(eta)), z_max=float(xp.max(eta)))


def _ray_tri(orig, d, v0, v1, v2, eps: float = 1e-12):
    """Moller-Trumbore on B rays against B triangles (elementwise).
    Inputs (B, 3); returns (t, hit) of shape (B,).  Loose barycentric
    tolerances keep edge/diagonal hits from slipping between triangles."""
    xp = xp_of(orig, d)
    e1 = v1 - v0
    e2 = v2 - v0
    h = xp.cross(d, e2)
    a = xp.sum(e1 * h, axis=-1)
    par = xp.abs(a) < eps
    inv = 1.0 / xp.where(par, xp.ones_like(a), a)
    s = orig - v0
    u = inv * xp.sum(s * h, axis=-1)
    q = xp.cross(s, e1)
    v = inv * xp.sum(d * q, axis=-1)
    t = inv * xp.sum(e2 * q, axis=-1)
    hit = (~par & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1.0 + 1e-9))
    return t, hit


def heightmap_intersect(origins, dirs, cell_tbl: dict,
                        max_steps: int | None = None):
    """First intersection of B rays with the periodic facet mesh.

    Returns (t_hit, hit_pos, n_hat):
        t_hit   : (B,) distance along the ray, +inf where no hit
        hit_pos : (B, 3) hit point in wrapped [0, L) coordinates
        n_hat   : (B, 3) upward unit normal of the struck triangle
    """
    xp = xp_of(origins, dirs)
    V = cell_tbl["V"]
    N = cell_tbl["N"]
    dx = cell_tbl["dx"]
    L = cell_tbl["L"]
    z_min = cell_tbl["z_min"]
    z_max = cell_tbl["z_max"]
    if max_steps is None:
        max_steps = 4 * N

    B = origins.shape[0]
    pos = xp.copy(xp.asarray(origins, dtype=float))
    pos[:, 0] %= L
    pos[:, 1] %= L
    dirs = xp.asarray(dirs, dtype=float)
    dxr, dyr, dzr = dirs[:, 0], dirs[:, 1], dirs[:, 2]

    cj = xp.clip(xp.astype(xp.floor(pos[:, 0] / dx), int), 0, N - 1)
    ci = xp.clip(xp.astype(xp.floor(pos[:, 1] / dx), int), 0, N - 1)
    step_j = xp.where(dxr > 0, 1, -1)
    step_i = xp.where(dyr > 0, 1, -1)

    t_base = xp.zeros(B)
    t_hit = xp.full(B, xp.inf)
    hit_pos = xp.zeros((B, 3))
    n_hat = xp.zeros((B, 3))
    n_hat[:, 2] = 1.0
    alive = xp.ones(B, dtype=bool)
    tiny = 1e-30

    for _ in range(max_steps):
        if not bool(xp.any(alive)):
            break

        Vc = V[ci, cj]
        v00, v01, v11, v10 = Vc[:, 0], Vc[:, 1], Vc[:, 2], Vc[:, 3]
        t1, h1 = _ray_tri(pos, dirs, v00, v01, v11)
        t2, h2 = _ray_tri(pos, dirs, v00, v11, v10)

        cjf = xp.astype(cj, float) + xp.astype(dxr > 0, float)
        cif = xp.astype(ci, float) + xp.astype(dyr > 0, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            tx = (cjf * dx - pos[:, 0]) / dxr
            ty = (cif * dx - pos[:, 1]) / dyr
        tx = xp.where(xp.abs(dxr) < tiny, xp.inf, xp.maximum(tx, 0.0))
        ty = xp.where(xp.abs(dyr) < tiny, xp.inf, xp.maximum(ty, 0.0))
        t_exit = xp.minimum(tx, ty)

        ok1 = h1 & (t1 >= -1e-12) & (t1 <= t_exit + 1e-9)
        ok2 = h2 & (t2 >= -1e-12) & (t2 <= t_exit + 1e-9)
        t1v = xp.where(ok1, t1, xp.inf)
        t2v = xp.where(ok2, t2, xp.inf)
        cand = xp.minimum(t1v, t2v)
        newly = alive & xp.isfinite(cand)

        if bool(xp.any(newly)):
            use1 = t1v <= t2v
            n1 = xp.cross(v01 - v00, v11 - v00)
            n2 = xp.cross(v11 - v00, v10 - v00)
            n_tri = xp.where(use1[:, None], n1, n2)
            n_tri = normalize(n_tri)
            n_tri = xp.where(n_tri[:, 2:3] < 0, -n_tri, n_tri)
            cand_safe = xp.where(xp.isfinite(cand), cand,
                                 xp.zeros_like(cand))
            t_hit = xp.where(newly, t_base + cand_safe, t_hit)
            hit_pos = xp.where(newly[:, None],
                               pos + cand_safe[:, None] * dirs, hit_pos)
            n_hat = xp.where(newly[:, None], n_tri, n_hat)
            alive = alive & ~newly
            if not bool(xp.any(alive)):
                break

        # Advance live rays to the next cell boundary
        adv = xp.where(alive & xp.isfinite(t_exit), t_exit,
                       xp.zeros_like(t_exit))
        pos = pos + adv[:, None] * dirs
        t_base = t_base + adv
        go_x = tx <= ty
        zero_i = xp.zeros_like(step_j)
        cj_new = cj + xp.where(alive & go_x, step_j, zero_i)
        ci_new = ci + xp.where(alive & ~go_x, step_i, zero_i)
        # Periodic wrap: shift cell index and physical position together
        lo = cj_new < 0
        hi = cj_new >= N
        pos[:, 0] += xp.where(lo, L, 0.0) - xp.where(hi, L, 0.0)
        cj = xp.where(lo, cj_new + N, xp.where(hi, cj_new - N, cj_new))
        lo = ci_new < 0
        hi = ci_new >= N
        pos[:, 1] += xp.where(lo, L, 0.0) - xp.where(hi, L, 0.0)
        ci = xp.where(lo, ci_new + N, xp.where(hi, ci_new - N, ci_new))

        # Vertical cull: no further hit possible
        gone = (((dzr > 0) & (pos[:, 2] > z_max + 1e-9))
                | ((dzr < 0) & (pos[:, 2] < z_min - 1e-9)))
        alive = alive & ~gone

    return t_hit, hit_pos, n_hat


# ---------------------------------------------------------------------------
# Polarized volume scattering event
# ---------------------------------------------------------------------------

def _scatter_event(d, S_path, depol: float, rng):
    """Rayleigh volume scattering of K in-water rays: sample the
    scattering cosine from the unpolarized marginal (exact), then the
    shared polarized azimuth/event machinery
    (water.polarized_scatter_event)."""
    xp = xp_of(d, S_path)
    K = d.shape[0]
    Delta = (1.0 - depol) / (1.0 + 0.5 * depol)
    mu = sample_rayleigh_scattering(K, rng, depol, xp=xp)
    p00 = Delta * 0.75 * (1.0 + mu**2) + (1.0 - Delta)
    p01 = -Delta * 0.75 * (1.0 - mu**2)
    P = rayleigh_phase_mueller(mu, depol)
    return polarized_scatter_event(d, S_path, mu, p00, p01, P, rng)


# ---------------------------------------------------------------------------
# Forward tracer with Mueller-matrix path accumulation
# ---------------------------------------------------------------------------

def trace_forward(eta, dx: float, origins, dirs,
                  n_water: float = 1.34,
                  max_bounces: int = 12,
                  rng=None,
                  cell_tbl: dict | None = None,
                  water: WaterOptics | None = None) -> dict:
    """Trace B rays until they escape upward (reflectance), escape
    downward (transmittance), or are absorbed in the water body.  Each
    ray accumulates the Mueller matrix mapping its launch Stokes vector
    (meridian frame of the launch direction) to its escape Stokes vector
    (meridian frame of the escape direction).

    With `water` given, in-water rays propagate through a semi-infinite
    homogeneous medium by analog event MC: free paths are exponential in
    the beam attenuation c = a + b, collisions survive with probability
    omega_0 = b / c (else the ray terminates absorbed), and survivors
    scatter through the Rayleigh-with-depolarization phase Mueller matrix
    with explicit meridian <-> scattering-plane frame rotations.

    Surface branch probabilities and scattering azimuths are
    importance-sampled with the path's unpolarized-launch Stokes column,
    so the I-path weight of every ray stays exactly 1 through arbitrarily
    many polarizing events; first bounces reduce to the unpolarized
    reflectance M_R[0, 0].

    max_bounces caps the total event count (surface + volume); rays
    still in flight after the cap are classified by their current
    direction, except in-water rays under a `water` medium, which stay
    unresolved (escaped_up = escaped_down = absorbed = False).

    Returns dict(escaped_up, escaped_down, absorbed, dir, mueller,
    bounces, scatters, in_water).
    """
    xp = xp_of(eta, origins, dirs)
    rng = adapt_rng(rng, xp)
    if cell_tbl is None:
        cell_tbl = build_cell_tables(eta, dx)
    dx = cell_tbl["dx"]
    c_att = water.attenuation if water is not None else 0.0
    omega0 = water.albedo if water is not None else 0.0

    B = origins.shape[0]
    pos = xp.copy(xp.asarray(origins, dtype=float))
    dirn = normalize(xp.asarray(dirs, dtype=float))
    M = xp.copy(xp.broadcast_to(xp.eye(4), (B, 4, 4)))
    alive = xp.ones(B, dtype=bool)
    in_water = xp.zeros(B, dtype=bool)
    bounces = xp.zeros(B, dtype=int)
    scatters = xp.zeros(B, dtype=int)
    esc_up = xp.zeros(B, dtype=bool)
    esc_dn = xp.zeros(B, dtype=bool)
    absorbed = xp.zeros(B, dtype=bool)

    for _ in range(max_bounces):
        if not bool(xp.any(alive)):
            break
        idx = xp.flatnonzero(alive)
        t_hit, hpos, nrm = heightmap_intersect(pos[idx], dirn[idx], cell_tbl)

        # Volume collisions: exponential free path vs surface distance
        l_col = xp.full(idx.shape[0], xp.inf)
        if c_att > 0:
            wsel = in_water[idx]
            l_col[wsel] = rng.exponential(1.0 / c_att,
                                          int(xp.sum(wsel)))
        collide = l_col < t_hit
        nohit = ~xp.isfinite(t_hit) & ~collide

        esc = idx[nohit]
        if esc.shape[0]:
            up = dirn[esc, 2] > 0
            esc_up[esc[up]] = True
            esc_dn[esc[~up]] = True
            alive[esc] = False

        # Collision events: absorption roulette, then polarized scatter
        cidx = idx[collide]
        if cidx.shape[0]:
            pos[cidx] += l_col[collide][:, None] * dirn[cidx]
            die = rng.random(cidx.shape[0]) >= omega0
            absorbed[cidx[die]] = True
            alive[cidx[die]] = False
            sidx = cidx[~die]
            if sidx.shape[0]:
                d_new, M_step = _scatter_event(dirn[sidx], M[sidx][:, :, 0],
                                               water.depolarization, rng)
                M[sidx] = M_step @ M[sidx]
                dirn[sidx] = d_new
                scatters[sidx] += 1

        # Surface events
        hmask = xp.isfinite(t_hit) & ~collide
        hidx = idx[hmask]
        if hidx.shape[0] == 0:
            continue
        d = dirn[hidx]
        n = nrm[hmask]
        hp = hpos[hmask]

        # Orient normal against the incoming ray
        n = xp.where(xp.sum(n * d, axis=-1)[:, None] > 0, -n, n)
        cos_i = xp.clip(-xp.sum(d * n, axis=-1), 1e-9, 1.0)
        n_rel = xp.where(in_water[hidx], 1.0 / n_water, n_water)
        M_R, M_T, tir = fresnel_mueller(cos_i, n_rel)

        # Incident frame rotation: meridian -> scattering plane
        v_in, h_in = meridian_frame(d)
        s_axis = xp.cross(d, n)
        s_n = xp.linalg.norm(s_axis, axis=-1, keepdims=True)
        s_axis = xp.where(s_n > 1e-9, s_axis / xp.maximum(s_n, 1e-300),
                          h_in)
        p_in = xp.cross(s_axis, d)
        a_in = frame_rotation_angle(d, v_in, p_in)

        # Branch on the polarized reflectance of the path Stokes
        # (column 0 of M, rotated into the scattering basis): the event
        # weight M_R / p (or M_T / (1 - p)) then preserves the I-column
        # path weight exactly; first bounces reduce to M_R[0, 0].
        S_p = xp.einsum("...ij,...j->...i",
                        mueller_rotation(a_in), M[hidx][:, :, 0])
        with np.errstate(invalid="ignore", divide="ignore"):
            p_pol = ((M_R[:, 0, 0] * S_p[:, 0] + M_R[:, 0, 1] * S_p[:, 1])
                     / S_p[:, 0])
        p_refl = xp.where(S_p[:, 0] > 1e-12, p_pol, M_R[:, 0, 0])
        p_refl = xp.clip(p_refl, 0.0, 1.0)
        p_refl = xp.where(tir, xp.ones_like(p_refl), p_refl)
        do_refl = tir | (rng.random(hidx.shape[0]) < p_refl)

        # Specular and refracted directions
        d_refl = normalize(d + 2.0 * cos_i[:, None] * n)
        sin_t2 = xp.clip((1.0 - cos_i**2) / n_rel**2, 0.0, 1.0)
        cos_t = xp.sqrt(1.0 - sin_t2)
        d_trans = normalize(d / n_rel[:, None]
                            + (cos_i / n_rel - cos_t)[:, None] * n)
        d_new = xp.where(do_refl[:, None], d_refl, d_trans)

        # Outgoing frame rotation: scattering plane -> new meridian
        v_out, _ = meridian_frame(d_new)
        p_out = xp.cross(s_axis, d_new)
        a_out = frame_rotation_angle(d_new, p_out, v_out)

        w_r = xp.maximum(p_refl, 1e-12)[:, None, None]
        w_t = xp.maximum(1.0 - p_refl, 1e-12)[:, None, None]
        M_evt = xp.where(do_refl[:, None, None], M_R / w_r, M_T / w_t)
        M_step = mueller_rotation(a_out) @ M_evt @ mueller_rotation(a_in)
        M[hidx] = M_step @ M[hidx]
        in_water[hidx] ^= ~do_refl

        # Restart just off the surface on the departure side
        side = xp.sign(xp.sum(n * d_new, axis=-1))[:, None]
        pos[hidx] = hp + (1e-6 * dx) * side * n
        dirn[hidx] = d_new
        bounces[hidx] += 1

    rest = xp.flatnonzero(alive)
    if rest.shape[0]:
        if water is not None:
            rest = rest[~in_water[rest]]  # in-water leftovers unresolved
        up = dirn[rest, 2] > 0
        esc_up[rest[up]] = True
        esc_dn[rest[~up]] = True

    return dict(escaped_up=esc_up, escaped_down=esc_dn, absorbed=absorbed,
                dir=dirn, mueller=M, bounces=bounces, scatters=scatters,
                in_water=in_water)


# ---------------------------------------------------------------------------
# Effective Mueller matrix for one incident direction
# ---------------------------------------------------------------------------

def effective_mueller_for_incident(eta, dx: float,
                                   theta_i_deg: float,
                                   phi_i_deg: float = 0.0,
                                   n_rays: int = 100_000,
                                   n_theta_bins: int = 9,
                                   n_phi_bins: int = 24,
                                   n_water: float = 1.34,
                                   max_bounces: int = 12,
                                   rng=None,
                                   water: WaterOptics | None = None
                                   ) -> dict:
    """Effective reflection Mueller matrix M(theta_r, phi_r) for one sky
    direction, from a single Mueller-accumulating trace.

    With `water` given, transmitted rays propagate through the scattering
    water body and up-escaping rays split into a surface (glint) part
    that never scattered in water and a water-leaving part that did.

    Returns dict with:
        M_eff       : (nT, nP, 4, 4) summed Mueller per launched ray;
                      M_eff[..., 0, 0] sums to the hemispherical
                      reflectance for unpolarized input
        M_eff_glint, M_eff_water : same split by in-water scatter count
        M_bin_mean  : (nT, nP, 4, 4) mean Mueller of rays in each bin
        counts, counts_glint, counts_water
        theta_edges, phi_edges, bin_solid_angle
        R_total, R_glint, R_water, T_total, A_total, W_total
        n_up, n_down, mean_bounces, mean_scatters, launched

    A_total is absorbed energy; W_total is the unresolved in-water
    remainder at the event cap, so R + T + A + W = 1 exactly.
    """
    xp = xp_of(eta)
    rng = adapt_rng(rng, xp)
    cell_tbl = build_cell_tables(eta, dx)
    L = cell_tbl["L"]

    th = np.deg2rad(theta_i_deg)
    ph = np.deg2rad(phi_i_deg)
    d = xp.asarray([np.sin(th) * np.cos(ph),
                    np.sin(th) * np.sin(ph),
                    -np.cos(th)], dtype=float)
    z0 = cell_tbl["z_max"] + 0.5
    origins = xp.stack([rng.uniform(0.0, L, n_rays),
                        rng.uniform(0.0, L, n_rays),
                        xp.full(n_rays, z0)], axis=1)
    dirs = xp.broadcast_to(d, (n_rays, 3))

    res = trace_forward(eta, dx, origins, dirs, n_water=n_water,
                        max_bounces=max_bounces, rng=rng, cell_tbl=cell_tbl,
                        water=water)

    up = res["escaped_up"]
    dn = res["escaped_down"]
    ab = res["absorbed"]
    unresolved = ~(up | dn | ab)
    d_up = res["dir"][up]
    M_up = res["mueller"][up]
    scat_up = res["scatters"][up]

    theta_edges = xp.linspace(0.0, np.pi / 2, n_theta_bins + 1)
    phi_edges = xp.linspace(-np.pi, np.pi, n_phi_bins + 1)
    th_up = xp.arccos(xp.clip(d_up[:, 2], -1.0, 1.0))
    ph_up = xp.arctan2(d_up[:, 1], d_up[:, 0])
    ti = xp.clip(xp.searchsorted(theta_edges, th_up, side="right") - 1,
                 0, n_theta_bins - 1)
    pj = xp.clip(xp.searchsorted(phi_edges, ph_up, side="right") - 1,
                 0, n_phi_bins - 1)
    flat = ti * n_phi_bins + pj

    def _bin(sel):
        M_sum = xp.zeros((n_theta_bins * n_phi_bins, 4, 4))
        xp.index_add(M_sum, flat[sel], M_up[sel])
        cnt = xp.bincount(flat[sel],
                          minlength=n_theta_bins * n_phi_bins)
        return (M_sum.reshape(n_theta_bins, n_phi_bins, 4, 4),
                cnt.reshape(n_theta_bins, n_phi_bins))

    all_up = xp.ones(d_up.shape[0], dtype=bool)
    M_sum, counts = _bin(all_up)
    M_sum_g, counts_g = _bin(scat_up == 0)
    M_sum_w, counts_w = _bin(scat_up > 0)

    with np.errstate(invalid="ignore", divide="ignore"):
        M_bin_mean = M_sum / xp.where(counts > 0, xp.astype(counts, float),
                                      xp.nan)[..., None, None]
    d_omega = ((xp.cos(theta_edges[:-1]) - xp.cos(theta_edges[1:]))[:, None]
               * xp.diff(phi_edges)[None, :])

    return dict(M_eff=M_sum / n_rays,
                M_eff_glint=M_sum_g / n_rays,
                M_eff_water=M_sum_w / n_rays,
                M_bin_mean=M_bin_mean,
                counts=counts, counts_glint=counts_g, counts_water=counts_w,
                theta_edges=theta_edges, phi_edges=phi_edges,
                bin_solid_angle=d_omega,
                R_total=float(xp.sum(M_up[:, 0, 0]) / n_rays),
                R_glint=float(xp.sum(M_up[scat_up == 0, 0, 0]) / n_rays),
                R_water=float(xp.sum(M_up[scat_up > 0, 0, 0]) / n_rays),
                T_total=float(xp.sum(res["mueller"][dn][:, 0, 0]) / n_rays),
                A_total=float(xp.sum(res["mueller"][ab][:, 0, 0]) / n_rays),
                W_total=float(xp.sum(res["mueller"][unresolved][:, 0, 0])
                              / n_rays),
                n_up=int(xp.sum(up)), n_down=int(xp.sum(dn)),
                mean_bounces=float(xp.mean(xp.astype(res["bounces"],
                                                     float))),
                mean_scatters=float(xp.mean(xp.astype(res["scatters"],
                                                      float))),
                launched=n_rays)
