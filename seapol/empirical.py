"""
Bridges from measured (ASIT-style) spectra into the synthesis.

* reduce_kf_cube      : stream a raw Skw(f, kx, ky) cube into S(|k|, f)
                        with forbidden-region noise subtraction and an
                        off-shell (bound) fraction beta_obs(k)
* bound_fraction_from_kf_reduction : fitted beta(k) callable for
                        generate_sea_surface(bound_fraction=...)
* load_asit_run / psi_from_asit / generate_asit_surface : synthesize
                        surfaces directly from a measured directional
                        slope spectrum S_k_theta (validated convention:
                        Cartesian slope density per dkx dky, theta
                        wind-relative; mss = int k S cos^2/sin^2,
                        elevation density Psi_eta = S / k^2)
"""

from __future__ import annotations

import datetime as _dt
import re
import warnings
from pathlib import Path

import numpy as np

from .spectrum import angular_frequency, directional_spectrum

__all__ = ["bound_fraction_from_kf_reduction", "bound_fraction_for_wind",
           "bound_fraction_for_conditions", "run_conditions",
           "smooth_beta_curve", "reduce_kf_cube", "cube_timestamp",
           "match_env_run", "load_asit_run", "psi_from_asit",
           "generate_asit_surface"]


# ---------------------------------------------------------------------------
# Raw-cube reduction (k-f spectrum, noise floor, bound fraction)
# ---------------------------------------------------------------------------

def cube_timestamp(path) -> float:
    """POSIX timestamp parsed from an ASIT cube filename."""
    m = re.search(r"(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})",
                  Path(path).name)
    if not m:
        raise ValueError(f"no timestamp in {path}")
    return _dt.datetime(*map(int, m.groups()),
                        tzinfo=_dt.timezone.utc).timestamp()


def _env_winds(env_path) -> tuple[np.ndarray, np.ndarray]:
    """(t, U10) per run; eddy-covariance wind with COARE fallback."""
    import netCDF4 as nc
    d = nc.Dataset(env_path)
    t = np.array(d["t_seconds_since_January_1_1970"][:])
    U = np.array(d["EC_U_m_s"][:])
    if "COARE_U10" in d.variables:
        U_c = np.array(d["COARE_U10"][:])
        U = np.where(np.isfinite(U), U, U_c)
    d.close()
    return t, U


def match_env_run(env_path, timestamp: float,
                  max_dt: float = 1800.0) -> tuple[int, float]:
    """(run index, U10) of the nearest finite-wind run in the
    environmental file; raises if none within max_dt seconds."""
    t, U = _env_winds(env_path)
    order = np.argsort(np.abs(t - timestamp))
    for i in order:
        if np.isfinite(U[i]):
            if abs(t[i] - timestamp) > max_dt:
                break
            return int(i), float(U[i])
    raise ValueError(f"no finite-wind run within {max_dt} s")


def reduce_kf_cube(path, n_kbin: int = 60, f_stride: int = 3,
                   shell_tol: float = 0.20,
                   k_edges: np.ndarray | None = None,
                   verbose: bool = False) -> dict:
    """Stream one raw Skw(f, kx, ky) cube into a ring-integrated S(k, f)
    with the instrument noise floor (estimated per frequency from the
    forbidden region above the free shell and any physical phase speed)
    subtracted, plus the off-shell bound fraction beta_obs(k)."""
    import netCDF4 as nc
    d = nc.Dataset(path)
    f = np.array(d["f"][:])
    kx = np.array(d["kx"][:])
    ky = np.array(d["ky"][:])
    K = np.hypot(kx, ky).ravel()
    if k_edges is None:
        k_edges = np.geomspace(3.0, 1400.0, n_kbin + 1)
    n_kbin = k_edges.size - 1
    idx = np.digitize(K, k_edges) - 1
    valid = (idx >= 0) & (idx < n_kbin)
    idx_v = idx[valid]
    counts = np.bincount(idx_v, minlength=n_kbin)

    f_sel = np.arange(0, f.size, f_stride)
    S_kf = np.zeros((n_kbin, f_sel.size))
    for j, i_f in enumerate(f_sel):
        plane = np.asarray(d["Skw"][i_f, :, :]).ravel()[valid]
        S_kf[:, j] = np.bincount(idx_v, weights=plane, minlength=n_kbin)
        if verbose and j % 50 == 0:
            print(f"  {j}/{f_sel.size}")
    d.close()
    f_u = f[f_sel]
    k_c = np.sqrt(k_edges[:-1] * k_edges[1:])
    f_disp = angular_frequency(k_c) / (2 * np.pi)

    # noise floor from the forbidden region
    S_mean = S_kf / np.maximum(counts, 1)[:, None]
    KK, FF = np.meshgrid(k_c, f_u, indexing="ij")
    forb = ((FF > 1.5 * f_disp[:, None])
            & (FF > 3.5 * KK / (2 * np.pi)) & (FF > 1.5))
    noise_f = np.array([np.median(S_mean[forb[:, j], j])
                        if forb[:, j].sum() > 5 else np.nan
                        for j in range(f_u.size)])
    ok_f = np.isfinite(noise_f)
    if ok_f.any():
        noise_f = np.interp(np.arange(f_u.size), np.flatnonzero(ok_f),
                            noise_f[ok_f])
        S_kf = np.clip(S_kf - noise_f[None, :] * counts[:, None], 0.0, None)
    else:
        warnings.warn("no forbidden-region (k, f) bins available for "
                      "noise-floor estimation; skipping noise subtraction",
                      stacklevel=2)

    # off-shell fraction where the free shell is inside the band
    beta_obs = np.full(n_kbin, np.nan)
    good_f = f_u > 0.25
    for i in range(n_kbin):
        tot = S_kf[i, good_f].sum()
        if tot <= 0 or counts[i] == 0 or f_disp[i] > 0.85 * f_u.max():
            continue
        on = good_f & (np.abs(f_u - f_disp[i])
                       < np.maximum(shell_tol * f_disp[i],
                                    2 * np.diff(f_u).mean()))
        beta_obs[i] = 1.0 - S_kf[i, on].sum() / tot

    return dict(k=k_c, f=f_u, S_kf=S_kf, beta_obs=beta_obs,
                counts=counts, noise_f=noise_f, shell_tol=shell_tol,
                timestamp=cube_timestamp(path))


# ---------------------------------------------------------------------------
# beta(k) fitting
# ---------------------------------------------------------------------------

def smooth_beta_curve(k: np.ndarray, beta: np.ndarray,
                      n_pass: int = 2) -> np.ndarray:
    """Light 1-2-1 smoothing of a beta(k) curve with NaN gaps filled by
    interpolation in log k."""
    k = np.asarray(k, dtype=float)
    beta = np.asarray(beta, dtype=float).copy()
    ok = np.isfinite(beta)
    if ok.sum() < 4:
        raise ValueError("too few finite beta points")
    beta = np.interp(np.log(k), np.log(k[ok]), beta[ok])
    for _ in range(n_pass):
        beta[1:-1] = 0.25 * beta[:-2] + 0.5 * beta[1:-1] + 0.25 * beta[2:]
    return np.clip(beta, 0.0, 1.0)


def bound_fraction_from_kf_reduction(npz_path,
                                     k_taper: float = 8.0,
                                     beta_cap: float = 0.99,
                                     monotone: bool = True):
    """Callable beta(K) from a saved k-f reduction (asit_kf_reduced.npz).

    The measured off-shell fraction is smoothed, capped at beta_cap
    (a perfectly bound sea is unphysical: some free-wave background
    always survives), tapered to zero below k_taper (the dominant waves
    themselves are free), and held constant beyond the last valid
    measurement.  monotone enforces the non-decreasing-in-k envelope
    seen in every clean reduction, suppressing band-edge
    noise-subtraction dips.
    """
    d = np.load(npz_path)
    k = np.asarray(d["k"], dtype=float)
    beta = smooth_beta_curve(k, np.asarray(d["beta_obs"], dtype=float))
    if monotone:
        beta = np.maximum.accumulate(beta)
    beta = np.minimum(beta, beta_cap)

    log_k = np.log(k)

    def beta_of_k(K):
        K = np.asarray(K, dtype=float)
        with np.errstate(divide="ignore"):
            b = np.interp(np.log(np.maximum(K, 1e-12)), log_k, beta,
                          left=beta[0], right=beta[-1])
        s = np.clip(K / max(k_taper, 1e-9), 0.0, 1.0)
        return b * s * s * (3.0 - 2.0 * s)

    return beta_of_k


def bound_fraction_for_wind(library_dir, U10: float,
                            n_neighbors: int = 3,
                            k_taper: float = 8.0,
                            beta_cap: float = 0.99,
                            monotone: bool = True):
    """Callable beta(K) for a given wind speed, from a library of cube
    reductions (one npz per cube with fields U10, k, beta_obs, as
    produced by reduce_kf_cube).

    The n_neighbors reductions nearest in wind speed are averaged with
    inverse-distance weights, then smoothed/capped/tapered like
    bound_fraction_from_kf_reduction."""
    paths = sorted(Path(library_dir).glob("*.npz"))
    if not paths:
        raise FileNotFoundError(f"no reductions in {library_dir}")
    entries = []
    for p in paths:
        d = np.load(p)
        entries.append((float(d["U10"]), np.asarray(d["k"], dtype=float),
                        np.asarray(d["beta_obs"], dtype=float)))
    entries.sort(key=lambda e: abs(e[0] - U10))
    sel = entries[:max(1, n_neighbors)]
    k = sel[0][1]
    w_sum = 0.0
    beta_acc = np.zeros_like(k)
    for U_j, k_j, b_j in sel:
        b = smooth_beta_curve(k_j, b_j)
        if not np.array_equal(k_j, k):
            b = np.interp(np.log(k), np.log(k_j), b)
        w = 1.0 / (abs(U_j - U10) + 0.5)
        beta_acc += w * b
        w_sum += w
    beta = beta_acc / w_sum
    if monotone:
        beta = np.maximum.accumulate(beta)
    beta = np.minimum(beta, beta_cap)
    log_k = np.log(k)

    def beta_of_k(K):
        K = np.asarray(K, dtype=float)
        b = np.interp(np.log(np.maximum(K, 1e-12)), log_k, beta,
                      left=beta[0], right=beta[-1])
        s = np.clip(K / max(k_taper, 1e-9), 0.0, 1.0)
        return b * s * s * (3.0 - 2.0 * s)

    return beta_of_k


def run_conditions(stats_path, env_path) -> dict:
    """Per-run sea-state metadata: U10, dominant frequency fp, peak phase
    speed cp (deep water), and inverse wave age Omega = U10 / cp, from
    the theta-integrated S_f_theta and the merged winds."""
    import netCDF4 as nc
    g = 9.81
    d = nc.Dataset(stats_path)
    f = np.array(d["f_Hz"][:])
    th = np.array(d["theta_rad"][:])
    Sf = np.array(d["S_f_theta"][:])          # (runs, theta, f)
    d.close()
    dth = np.diff(th).mean()
    spec = Sf.sum(axis=1) * dth               # (runs, f)
    # dominant wind-wave frequency: restrict to f in [0.08, 1.5] Hz;
    # runs with no valid spectrum get NaN
    band = (f >= 0.08) & (f <= 1.5)
    sub = spec[:, band]
    fp = np.full(sub.shape[0], np.nan)
    has = np.isfinite(sub).any(axis=1)
    fp[has] = f[band][np.nanargmax(np.where(np.isfinite(sub), sub,
                                            -np.inf)[has], axis=1)]
    cp = g / (2.0 * np.pi * fp)
    t, U = _env_winds(env_path)
    with np.errstate(invalid="ignore", divide="ignore"):
        Omega = U / cp
    return dict(t=t, U10=U, fp=fp, cp=cp, inverse_wave_age=Omega)


def bound_fraction_for_conditions(library_dir, stats_path, env_path,
                                  U10: float,
                                  inverse_wave_age: float | None = None,
                                  n_neighbors: int = 3,
                                  k_taper: float = 8.0,
                                  beta_cap: float = 0.99,
                                  u_scale: float = 2.0,
                                  age_scale: float = 0.5):
    """Callable beta(K) selected by sea state, not wind alone: library
    reductions are ranked by distance in normalized (U10, Omega) space,
    Omega = U10 / cp the inverse wave age.  Falls back to wind-only
    ranking when inverse_wave_age is None."""
    conds = run_conditions(stats_path, env_path)
    paths = sorted(Path(library_dir).glob("*.npz"))
    if not paths:
        raise FileNotFoundError(f"no reductions in {library_dir}")
    entries = []
    for p in paths:
        d = np.load(p)
        run = int(d["run"])
        U_j = float(d["U10"])
        Om_j = float(conds["inverse_wave_age"][run])
        dist = abs(U_j - U10) / u_scale
        if inverse_wave_age is not None and np.isfinite(Om_j):
            dist = np.hypot(dist,
                            (Om_j - inverse_wave_age) / age_scale)
        entries.append((dist, np.asarray(d["k"], dtype=float),
                        np.asarray(d["beta_obs"], dtype=float)))
    entries.sort(key=lambda e: e[0])
    sel = entries[:max(1, n_neighbors)]
    k = sel[0][1]
    w_sum = 0.0
    beta_acc = np.zeros_like(k)
    for dist, k_j, b_j in sel:
        b = smooth_beta_curve(k_j, b_j)
        if not np.array_equal(k_j, k):
            b = np.interp(np.log(k), np.log(k_j), b)
        w = 1.0 / (dist + 0.25)
        beta_acc += w * b
        w_sum += w
    beta = np.maximum.accumulate(beta_acc / w_sum)
    beta = np.minimum(beta, beta_cap)
    log_k = np.log(k)

    def beta_of_k(K):
        K = np.asarray(K, dtype=float)
        b = np.interp(np.log(np.maximum(K, 1e-12)), log_k, beta,
                      left=beta[0], right=beta[-1])
        s = np.clip(K / max(k_taper, 1e-9), 0.0, 1.0)
        return b * s * s * (3.0 - 2.0 * s)

    return beta_of_k


# ---------------------------------------------------------------------------
# Synthesis from a measured directional slope spectrum
# ---------------------------------------------------------------------------

def load_asit_run(stats_path, run: int,
                  env_path=None) -> dict:
    """One run's measured directional slope spectrum and metadata.

    S_k_theta convention (validated against the stored mss_upwind /
    mss_crosswind to ~20%): Cartesian slope spectral density per
    dkx dky as a function of (theta, k), theta wind-relative, so
        mss_up = int k S cos^2(theta) dk dtheta,
        Psi_eta(kx, ky) = S / k^2.
    """
    import netCDF4 as nc
    d = nc.Dataset(stats_path)
    out = dict(k=np.array(d["k_rad_m"][:]),
               theta=np.array(d["theta_rad"][:]),
               S=np.array(d["S_k_theta"][run]),
               mss_upwind=float(d["mss_upwind"][run]),
               mss_crosswind=float(d["mss_crosswind"][run]),
               run=run)
    d.close()
    if env_path is not None:
        t, U = _env_winds(env_path)
        out["U10"] = float(U[run])
        out["t"] = float(t[run])
    return out


def psi_from_asit(run_data: dict, U10: float | None = None,
                  tail_blend: tuple[float, float] = (300.0, 550.0),
                  k_max: float | None = None):
    """Callable Psi_eta(KX, KY) [m^4/rad^2] from a measured run.

    Inside the measured band, bilinear interpolation of the slope
    density over (theta, log k) divided by k^2.  Beyond the last
    measured wavenumber the Elfouhaily form continues the tail, scaled
    to match the measured azimuthal-mean level over tail_blend.  Below
    the first measured wavenumber the density is zero (use tiles of
    side <~ 2 pi / k_min).

    k_max [rad/m] band-limits the spectrum with a cosine taper to zero
    over [0.85 k_max, k_max].  Use it to cap the synthesis at an
    instrument's reliable resolution: a camera sampling at FS frames/s
    cannot measure wave dynamics above FS/2, which by the gravity-
    capillary dispersion is roughly k ~ 400 (FS/2 ~ 15 Hz) -- content
    above that is guesswork and, when time-evolved, aliases."""
    k_m = run_data["k"]
    th_m = run_data["theta"]
    S = np.asarray(run_data["S"], dtype=float)  # (n_theta, n_k)
    if U10 is None:
        U10 = run_data.get("U10")
    if U10 is None or not np.isfinite(U10):
        raise ValueError("U10 required for the spectral tail")

    log_k = np.log(k_m)
    dth = np.diff(th_m).mean()
    n_th = th_m.size

    # scalar seam factor: measured vs Elfouhaily azimuthal-mean slope
    # density over the blend band
    sel = (k_m >= tail_blend[0]) & (k_m <= tail_blend[1])
    meas_level = float(np.mean(S[:, sel]))
    phi_ring = th_m[:, None]
    k_ring = k_m[None, sel]
    elf = directional_spectrum(k_ring * np.cos(phi_ring),
                               k_ring * np.sin(phi_ring), U10)
    elf_level = float(np.mean(elf * k_ring**2))
    seam = meas_level / max(elf_level, 1e-300)

    def psi(KX, KY):
        K = np.hypot(KX, KY)
        PHI = np.mod(np.arctan2(KY, KX), 2.0 * np.pi)
        out = np.zeros_like(K)

        inside = (K >= k_m[0]) & (K <= k_m[-1])
        if np.any(inside):
            lk = np.log(K[inside])
            ik = np.clip(np.searchsorted(log_k, lk) - 1, 0, k_m.size - 2)
            tk = (lk - log_k[ik]) / (log_k[ik + 1] - log_k[ik])
            # anchor to the measured theta grid origin (grids may start
            # at -pi rather than 0)
            ft = np.mod(PHI[inside] - th_m[0], 2.0 * np.pi) / dth
            it = np.floor(ft).astype(int) % n_th
            it1 = (it + 1) % n_th
            tt = ft - np.floor(ft)
            s_val = ((1 - tt) * ((1 - tk) * S[it, ik]
                                 + tk * S[it, ik + 1])
                     + tt * ((1 - tk) * S[it1, ik]
                             + tk * S[it1, ik + 1]))
            out[inside] = s_val / K[inside] ** 2

        beyond = K > k_m[-1]
        if np.any(beyond):
            out[beyond] = seam * directional_spectrum(
                KX[beyond], KY[beyond], U10)
        if k_max is not None:
            s = np.clip((k_max - K) / (0.15 * k_max), 0.0, 1.0)
            out = out * (s * s * (3.0 - 2.0 * s))
        return out

    return psi


def generate_asit_surface(stats_path, run: int, L: float, N: int,
                          env_path=None, U10: float | None = None,
                          k_max: float | None = None,
                          **kwargs):
    """Sea surface synthesized from a measured ASIT directional spectrum.

    k_max [rad/m] band-limits the spectrum (see psi_from_asit) -- e.g.
    to a camera's reliable resolution for instrument-matched synthesis.
    Remaining keyword arguments (times, bound_fraction, bound_speed,
    rng, backend, ...) pass through to generate_sea_surface."""
    from .surface import generate_sea_surface
    run_data = load_asit_run(stats_path, run, env_path=env_path)
    if U10 is None:
        U10 = run_data.get("U10", np.nan)
    psi = psi_from_asit(run_data, U10=U10, k_max=k_max)
    surf = generate_sea_surface(L, N, U10, psi_override=psi, **kwargs)
    surf.info["asit_run"] = run
    surf.info["mss_measured"] = (run_data["mss_upwind"],
                                 run_data["mss_crosswind"])
    return surf
