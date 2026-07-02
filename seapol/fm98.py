"""
Fedorov & Melville (1998) steady, wind-forced, viscously damped nonlinear
gravity-capillary waves, plus a gauge-fixed Stokes-coefficient table for
the hybrid surface generator.

Reference:
    Fedorov, A. V., Melville, W. K. (1998). "Nonlinear gravity-capillary
    waves with forcing and dissipation." J. Fluid Mech. 354, 1-42.

Formulation (paper Sec. 2, recast):
    Z^i(zeta) = (1/k) [zeta + i S(zeta)],  S = sum_m a_m exp(-i m zeta)
    U^i = -(c/k) |Z^i_zeta|^{-1}
    1/R^i = -Im(Z^i_zetazeta / Z^i_zeta) / |Z^i_zeta|
    Z = Z^i - (2 i nu / (c^2 k)) integral [U^i/R^i - <U^i/R^i>] dzeta
    Bernoulli: U^i^2/2 + g Y + (sigma/rho)/R + P0(X)/rho
               + (nu k / c) d(U^i^2)/dzeta = E
    with forcing P0(X) = p rho c0^2 cos(k X) and amplitude constraint
    (max Y - min Y)/2 = a.  Class 1 (wind-forced, parasitic-capillary-
    bearing) has arg(a_1) ~ pi.

Solved by spectral collocation on N = 4M zeta points and Newton iteration
(MINPACK hybr with Levenberg-Marquardt fallback), with geometric
continuation in (ak, p) from the small-amplitude linear guess.

FM98 uses clean-water constants (rho = 1000, sigma = 0.073, nu = 1e-6);
these are kept local to this module.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.fft import fft, ifft
from scipy.optimize import root
from threadpoolctl import threadpool_limits

G_EARTH = 9.81
RHO_WATER = 1000.0
SIGMA_WATER = 0.073
NU_WATER = 1.0e-6
T_OVER_RHO = SIGMA_WATER / RHO_WATER


def linear_phase_speed(k: float, T: float = T_OVER_RHO,
                       g: float = G_EARTH) -> float:
    """c0 = sqrt(g/k + T k)."""
    return float(np.sqrt(g / k + T * k))


def default_forcing(k: float, ak: float, delta: float = 0.02,
                    T: float = T_OVER_RHO, nu: float = NU_WATER,
                    g: float = G_EARTH) -> float:
    """Wind-forcing amplitude p giving a forcing/profile phase offset
    ~delta rad (linear balance delta = 4 nu k ak / (c0 p)).

    FM98 admit a RANGE of p at fixed (k, ak), parameterized by that
    phase shift.  Marginal forcing (delta ~ 1) lands on weakly forced
    branches whose ripples sit in the trough and creep up the rear
    face; strong forcing selects the branch whose parasitic train
    rides the forward face below the crest and decays through the
    trough.  delta = 0.02 is calibrated to the published example of
    Fedorov, Melville & Rozenberg (1998, Phys. Fluids 10, fig. 1):
    lam = 7 cm, ak = 0.30, p = 0.015 (this rule gives 0.0158)."""
    c0 = linear_phase_speed(k, T=T, g=g)
    return float(4.0 * nu * k * ak / (c0 * delta))


# ---------------------------------------------------------------------------
# Spectral utilities on the 2pi-periodic zeta grid
# ---------------------------------------------------------------------------

def _zeta_grid(N: int) -> np.ndarray:
    return 2.0 * np.pi * np.arange(N) / N


def _spectral_deriv(f: np.ndarray) -> np.ndarray:
    N = f.size
    F = fft(f)
    kk = np.fft.fftfreq(N, d=1.0 / N).astype(float)
    Fp = 1j * kk * F
    if N % 2 == 0:
        Fp[N // 2] = 0.0
    out = ifft(Fp)
    return out.real if np.isrealobj(f) else out


def _spectral_antideriv(f: np.ndarray) -> np.ndarray:
    """Antiderivative of a mean-zero periodic f, shifted so I(0) = 0."""
    N = f.size
    F = fft(f)
    kk = np.fft.fftfreq(N, d=1.0 / N).astype(float)
    Fi = np.zeros_like(F)
    nz = kk != 0
    Fi[nz] = F[nz] / (1j * kk[nz])
    I = ifft(Fi)
    return I - I[0]


def _stokes_series_derivs(a: np.ndarray, zeta: np.ndarray):
    m = np.arange(1, a.size + 1)
    phase = np.exp(-1j * np.outer(zeta, m))
    S = phase @ a
    Sp = phase @ (a * (-1j * m))
    Spp = phase @ (a * (-m * m))
    return S, Sp, Spp


def compute_surface(a: np.ndarray, c: float, k: float, zeta: np.ndarray,
                    T: float = T_OVER_RHO, nu: float = NU_WATER) -> dict:
    """Surface geometry, interface speed, curvature, and damping terms."""
    S, Sp, Spp = _stokes_series_derivs(a, zeta)

    Zi = (zeta + 1j * S) / k
    Zi_zeta = (1.0 + 1j * Sp) / k
    Zi_zetazeta = (1j * Spp) / k

    Zi_zeta_abs2 = (Zi_zeta * np.conj(Zi_zeta)).real
    Ui = -(c / k) / np.sqrt(Zi_zeta_abs2)
    inv_Ri = -np.imag(Zi_zetazeta / Zi_zeta) / np.sqrt(Zi_zeta_abs2)
    Ui2_zeta = _spectral_deriv(Ui * Ui)

    # Viscous boundary-layer correction to the surface position
    integrand = Ui * inv_Ri
    I = _spectral_antideriv(integrand - integrand.mean())
    Z = Zi - (2.0j * nu) / (c * c * k) * I
    X, Y = Z.real, Z.imag

    # Curvature of the corrected surface
    Z_periodic = Z - zeta / k
    Z_zeta = 1.0 / k + _spectral_deriv(Z_periodic)
    Z_zetazeta = _spectral_deriv(Z_zeta)
    inv_R = -np.imag(Z_zetazeta / Z_zeta) \
        / np.sqrt((Z_zeta * np.conj(Z_zeta)).real)

    return dict(Ui=Ui, inv_Ri=inv_Ri, Ui2_zeta=Ui2_zeta,
                Z=Z, X=X, Y=Y, inv_R=inv_R)


# ---------------------------------------------------------------------------
# Residual system
# ---------------------------------------------------------------------------

def _pack(a: np.ndarray, c: float, E: float) -> np.ndarray:
    x = np.empty(2 * a.size + 2)
    x[0:2 * a.size:2] = a.real
    x[1:2 * a.size:2] = a.imag
    x[-2] = c
    x[-1] = E
    return x


def _unpack(x: np.ndarray):
    M = (x.size - 2) // 2
    a = x[0:2 * M:2] + 1j * x[1:2 * M:2]
    return a, float(x[-2]), float(x[-1])


def _smooth_extreme(y: np.ndarray, a_scale: float, sign: float) -> float:
    """Differentiable max (sign=+1) / min (sign=-1) via log-sum-exp.
    beta = 5000/a keeps the bias ~0.1% of the amplitude: parasitic
    ripple generation is exponentially sensitive to crest steepness, so
    a soft constraint that under-delivers ak by even a few percent
    visibly suppresses the train."""
    beta = 5000.0 / max(a_scale, 1e-12)
    z = sign * y
    s = z.max()
    return sign * (s + np.log(np.exp(beta * (z - s)).sum()) / beta)


def residual(x: np.ndarray, k: float, a_amp: float, p: float, N: int,
             T: float = T_OVER_RHO, nu: float = NU_WATER, g: float = G_EARTH,
             c0: float | None = None) -> np.ndarray:
    """2M + 2 residuals: Bernoulli Fourier projections + amplitude.

    The amplitude constraint uses smooth extrema: with parasitic ripples
    in the trough the argmin of Y jumps between ripple troughs, and a
    hard min would make the system non-differentiable there."""
    a, c, E = _unpack(x)
    M = a.size
    zeta = _zeta_grid(N)
    surf = compute_surface(a, c, k, zeta, T=T, nu=nu)

    if c0 is None:
        c0 = linear_phase_speed(k, T=T, g=g)
    P0_over_rho = (p * c0 * c0) * np.cos(k * surf["X"])

    bern = (0.5 * surf["Ui"] ** 2 + g * surf["Y"] + T * surf["inv_R"]
            + P0_over_rho + (nu * k / c) * surf["Ui2_zeta"] - E)

    B = fft(bern) / N
    eqs = np.empty(2 * M + 1)
    eqs[0] = B[0].real
    eqs[1:2 * M + 1:2] = B[1:M + 1].real
    eqs[2:2 * M + 1:2] = B[1:M + 1].imag
    # Relative amplitude misfit scaled to the Bernoulli rows (m^2/s^2):
    # an unscaled (meters) constraint lets Newton park a several-percent
    # amplitude shortfall inside an apparently small residual norm, and
    # parasitic ripple generation is exponentially sensitive to ak
    amp_eq = (0.5 * (_smooth_extreme(surf["Y"], a_amp, 1.0)
                     - _smooth_extreme(surf["Y"], a_amp, -1.0)) / a_amp
              - 1.0) * (c0 * c0)
    return np.concatenate([eqs, [amp_eq]])


def linear_initial_guess(k: float, a_amp: float, p: float, M: int,
                         which_class: int = 1,
                         T: float = T_OVER_RHO, nu: float = NU_WATER,
                         g: float = G_EARTH) -> np.ndarray:
    """Small-amplitude guess from FM98 Sec. 8."""
    c0 = linear_phase_speed(k, T=T, g=g)
    delta_lin = 4.0 * nu * k / c0 * (a_amp * k) / max(p, 1e-12)
    delta_lin = float(np.clip(delta_lin, 0.0, 1.0))
    Theta = (-np.pi + delta_lin) if which_class == 1 else (-delta_lin)
    ak = a_amp * k
    a = np.zeros(M, dtype=complex)
    a[0] = ak * np.exp(-1j * Theta)
    c2 = c0 * c0 * (1.0 + (p / max(ak, 1e-12)) * np.cos(Theta))
    c = np.sqrt(max(c2, 0.25 * c0 * c0))
    return _pack(a, c, 0.5 * c * c)


# ---------------------------------------------------------------------------
# Solver and continuation
# ---------------------------------------------------------------------------

@dataclass
class FMSolution:
    wavelength: float
    k: float
    ak: float
    p: float
    which_class: int
    c0: float
    c: float
    E: float
    a: np.ndarray
    zeta: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    residual_norm: float
    converged: bool


def solve_fm98(wavelength: float, ak: float, p: float,
               M: int = 32, N: int | None = None,
               which_class: int = 1,
               x0: np.ndarray | None = None,
               T: float = T_OVER_RHO, nu: float = NU_WATER,
               g: float = G_EARTH,
               tol: float = 1e-11, maxiter: int = 300) -> FMSolution:
    """Solve the FM98 system at a single (wavelength, ak, p)."""
    k = 2.0 * np.pi / wavelength
    a_amp = ak / k
    if N is None:
        N = max(4 * M, 128)
    if N < 2 * M + 2:
        raise ValueError("need N >= 2M + 2")
    c0 = linear_phase_speed(k, T=T, g=g)
    if x0 is None:
        x0 = linear_initial_guess(k, a_amp, p, M, which_class=which_class,
                                  T=T, nu=nu, g=g)

    def F(x):
        return residual(x, k=k, a_amp=a_amp, p=p, N=N, T=T, nu=nu, g=g,
                        c0=c0)

    # Single-threaded BLAS: the Newton iterations are dominated by tiny
    # (N x M) products where multi-threaded BLAS dispatch overhead
    # swamps the arithmetic (26x slowdown measured on a 32-core host)
    best_x = x0.copy()
    with threadpool_limits(limits=1):
        best_res = float(np.linalg.norm(F(best_x)))
        try:
            sol = root(F, x0, method="hybr",
                       options={"xtol": tol,
                                "maxfev": max(2000, 40 * (x0.size + 1)),
                                "factor": 1.0})
            r = float(np.linalg.norm(F(sol.x)))
            if r < best_res:
                best_res, best_x = r, sol.x.copy()
        except Exception:
            pass
        if best_res > 1e-7:
            try:
                sol_lm = root(F, best_x, method="lm",
                              options={"xtol": 1e-13, "ftol": 1e-13,
                                       "maxiter": maxiter})
                r = float(np.linalg.norm(F(sol_lm.x)))
                if r < best_res:
                    best_res, best_x = r, sol_lm.x.copy()
            except Exception:
                pass

    a, c, E = _unpack(best_x)
    zeta = _zeta_grid(N)
    surf = compute_surface(a, c, k, zeta, T=T, nu=nu)
    return FMSolution(wavelength=wavelength, k=k, ak=ak, p=p,
                      which_class=which_class, c0=c0, c=c, E=E, a=a,
                      zeta=zeta, X=surf["X"], Y=surf["Y"],
                      residual_norm=best_res, converged=best_res < 1e-6)


def solve_fm98_continuation(wavelength: float, ak_target: float,
                            p_target: float,
                            M: int = 32, N: int | None = None,
                            which_class: int = 1,
                            n_steps: int = 10,
                            ak_start: float = 0.02,
                            p_start_ratio: float = 0.2,
                            p_ramp: str = "geometric",
                            res_ok: float = 2e-4,
                            res_cap: float = 2e-3,
                            max_bisect: int = 6,
                            **kwargs) -> FMSolution:
    """Continuation in (ak, p) from the linear guess, with adaptive
    step-size marching and a tolerance-relaxation ladder.

    Branch selection is PATH dependent: the default geometric p ramp
    (p from p_start_ratio * p_target, i.e. over-forced early relative
    to p ~ ak) reliably lands on the well-forced forward-face-train
    branch; p_ramp="proportional" (constant phase offset delta along
    the chain) is kept as an alternative.

    The achievable Newton residual grows with ak (the system passes
    near branch folds), so marching accepts steps at res_ok and, when
    the step size bottoms out, relaxes the acceptance threshold by 5x
    at a time up to res_cap before giving up."""
    # target below the requested start: solve directly, no marching
    direct = ak_target <= ak_start
    ak_start = min(ak_start, ak_target * 0.5)

    if p_ramp == "proportional":
        def p_of(ak1):
            return p_target * ak1 / ak_target
    elif p_ramp == "geometric":
        p_start = max(p_target * p_start_ratio, 1e-6)
        if p_target > p_start and ak_target > ak_start:
            slope = (np.log(p_target / p_start)
                     / np.log(ak_target / ak_start))

            def p_of(ak1):
                return p_start * (ak1 / ak_start) ** slope
        else:
            def p_of(ak1):
                return p_target
    else:
        raise ValueError(f"unknown p_ramp mode: {p_ramp!r}")

    def _solve(ak1, x_init):
        return solve_fm98(wavelength=wavelength, ak=float(ak1),
                          p=float(p_of(ak1)), M=M, N=N,
                          which_class=which_class, x0=x_init, **kwargs)

    if direct:
        return _solve(ak_target, None)

    # Adaptive step-size marching in log ak: grow the step after
    # successes, shrink it on failure; the steep end (where the
    # parasitic train develops fastest) gets the small steps it needs.
    step_full = np.log(ak_target / ak_start) / max(n_steps - 1, 1)
    step_min = step_full / 2 ** max_bisect
    step = step_full
    ak_cur = ak_start
    sol = _solve(ak_start, None)
    if sol.residual_norm >= res_ok:
        return sol
    x0 = _pack(sol.a, sol.c, sol.E)
    best = sol
    res_lim = res_ok

    # Acceptance = converged AND single-crested: 'lumpy' multi-bump
    # steady states coexist with the carrier + parasitic-train branch
    # (residual and c/c0 do not discriminate them), and one accepted
    # lumpy step poisons every warm start after it.
    budget = 8 * n_steps
    while ak_cur < ak_target * (1.0 - 1e-12) and budget > 0:
        budget -= 1
        ak_try = min(ak_cur * np.exp(step), ak_target)
        trial = _solve(ak_try, x0)
        if trial.residual_norm < res_lim and _single_crested(trial.Y):
            ak_cur = ak_try
            x0 = _pack(trial.a, trial.c, trial.E)
            best = trial
            step = min(step * 1.4, step_full)
        else:
            step *= 0.5
            if step < step_min:
                if res_lim < res_cap:
                    res_lim = min(res_lim * 5.0, res_cap)
                    step = 0.25 * step_full
                else:
                    break
    if ak_cur >= ak_target * (1.0 - 1e-12):
        return best
    # stalled below target: a last attempt at the full target
    final = _solve(ak_target, x0)
    if final.residual_norm < res_lim and _single_crested(final.Y):
        return final
    return best


# ---------------------------------------------------------------------------
# Gauge-fixed Stokes-coefficient table
# ---------------------------------------------------------------------------

def gauge_fix(a: np.ndarray, target_arg_a1: float = np.pi) -> np.ndarray:
    """Rotate coefficients so arg(a_1) = target_arg_a1.

    The FM98 system is invariant under zeta -> zeta + dz, which maps
    a_m -> a_m exp(-i m dz); Newton solves land in arbitrary gauges.
    Fixing arg(a_1) (Class 1 analytic limit = pi) makes bilinear
    interpolation of complex a_m across table grid points well-defined
    while preserving the gauge invariants |a_m| and m arg(a_1) - arg(a_m).
    """
    if abs(a[0]) < 1e-14:
        return a.copy()
    dz = np.angle(a[0]) - target_arg_a1
    m = np.arange(1, a.size + 1)
    return a * np.exp(-1j * m * dz)


def _single_crested(Y: np.ndarray, frac: float = 0.6) -> bool:
    """One dominant crest per period: coexisting 'lumpy' steady states
    carry a secondary mid-trough bump at >~ 60% of the crest height
    (measured from the minimum), while genuine parasitic-train ripple
    peaks stay below ~25-40%."""
    Y = np.asarray(Y, dtype=float)
    Yr = np.roll(Y, Y.size // 2 - int(np.argmax(Y)))
    z = Yr - Yr.min()
    peaks = (z[1:-1] > z[:-2]) & (z[1:-1] > z[2:]) \
        & (z[1:-1] > frac * z.max())
    return int(peaks.sum()) <= 1


def is_physically_plausible(a: np.ndarray, ak: float,
                            tol_stokes: float = 1.0) -> bool:
    """Stokes-like plausibility: |a_1| ~ ak, |a_m| decaying through the
    first three harmonics, bounded harmonic ratios."""
    if a.size < 3 or abs(a[0]) < 1e-10:
        return False
    if abs(a[0]) < 0.5 * ak or abs(a[0]) > 2.0 * ak:
        return False
    mags = np.abs(a)
    if not (mags[0] > mags[1] > mags[2]):
        return False
    if mags[1] / mags[0] > tol_stokes:
        return False
    if mags[2] / mags[0] > 0.6 * tol_stokes:
        return False
    return True


@dataclass
class FM98Table:
    """Tabulated gauge-fixed complex Stokes coefficients a_m(k, ak)."""
    k_grid: np.ndarray
    ak_grid: np.ndarray
    a_coeffs: np.ndarray          # (n_k, n_ak, M_keep) complex
    c_over_c0: np.ndarray
    converged: np.ndarray
    M_keep: int
    p: float

    def save(self, path) -> None:
        np.savez(path, k_grid=self.k_grid, ak_grid=self.ak_grid,
                 a_coeffs=self.a_coeffs, c_over_c0=self.c_over_c0,
                 converged=self.converged,
                 meta=np.array([self.M_keep, self.p]))

    @classmethod
    def load(cls, path) -> "FM98Table":
        d = np.load(path)
        return cls(k_grid=d["k_grid"], ak_grid=d["ak_grid"],
                   a_coeffs=d["a_coeffs"], c_over_c0=d["c_over_c0"],
                   converged=d["converged"],
                   M_keep=int(d["meta"][0]), p=float(d["meta"][1]))

    def interp(self, k: float, ak: float) -> tuple[np.ndarray, float]:
        """Bilinear (log k, linear ak) interpolation, rescaled so
        |a_1| = ak.  Below half the smallest tabulated ak the linear
        limit (a_1 = ak, a_m>=2 = 0) is returned."""
        ak = float(ak)
        k = float(k)
        kg, ag = self.k_grid, self.ak_grid
        if ak <= ag[0] * 0.5:
            # linear limit in the table gauge arg(a_1) = pi
            a = np.zeros(self.M_keep, dtype=complex)
            a[0] = -ak + 0j
            return a, 1.0
        ak = float(np.clip(ak, ag[0], ag[-1]))
        k = float(np.clip(k, kg[0], kg[-1]))

        log_k, log_kg = np.log(k), np.log(kg)
        i = int(np.clip(np.searchsorted(log_kg, log_k) - 1, 0, len(kg) - 2))
        j = int(np.clip(np.searchsorted(ag, ak) - 1, 0, len(ag) - 2))
        tk = (log_k - log_kg[i]) / (log_kg[i + 1] - log_kg[i])
        ta = (ak - ag[j]) / (ag[j + 1] - ag[j])

        A = self.a_coeffs
        a = ((1 - tk) * (1 - ta) * A[i, j] + tk * (1 - ta) * A[i + 1, j]
             + (1 - tk) * ta * A[i, j + 1] + tk * ta * A[i + 1, j + 1])
        C = self.c_over_c0
        c_ratio = float((1 - tk) * (1 - ta) * C[i, j]
                        + tk * (1 - ta) * C[i + 1, j]
                        + (1 - tk) * ta * C[i, j + 1]
                        + tk * ta * C[i + 1, j + 1])
        if abs(a[0]) > 1e-12:
            a = a * (ak / abs(a[0]))
        else:
            a[0] = ak + 0j
        return a, c_ratio


def _solve_with_validation(wavelength, ak_target, p_target, M, n_steps,
                           which_class):
    """Continuation with progressively softer restarts; only physically
    plausible solutions are accepted.  Beyond the branch fold the
    marching can land on collapsed-wave-speed states (c/c0 dropping by
    tens of percent, e.g. lam = 5 cm pushed to ak = 0.38), so clean
    propagating solutions also require c/c0 near 1."""
    best = None
    attempts = [
        dict(ak_start=0.02, p_start_ratio=0.2, n_steps=n_steps),
        dict(ak_start=0.01, p_start_ratio=0.1, n_steps=n_steps + 2),
        dict(ak_start=0.005, p_start_ratio=0.05, n_steps=n_steps + 4),
    ]
    for att in attempts:
        try:
            sol = solve_fm98_continuation(
                wavelength=wavelength, ak_target=float(ak_target),
                p_target=p_target, M=M, N=4 * M, which_class=which_class,
                **att)
            if (is_physically_plausible(sol.a, ak_target)
                    and abs(sol.c / sol.c0 - 1.0) < 0.15
                    and _single_crested(sol.Y)):
                if best is None or sol.residual_norm < best.residual_norm:
                    best = sol
                if sol.residual_norm < 1e-5:
                    break
        except Exception:
            continue
    return best


def build_fm98_table(k_grid: np.ndarray, ak_grid: np.ndarray,
                     M_keep: int = 6, M_solve: int = 32,
                     p: float | None = None, which_class: int = 1,
                     n_steps: int = 10, verbose: bool = False) -> FM98Table:
    """Solve FM98 at every (k, ak) grid point.  Failed/implausible points
    fall back to linear Stokes (a_1 = -ak, correct small-amplitude gauge).

    p=None (default) applies the well-forced branch rule per node,
    p = default_forcing(k, ak); a scalar p forces that value
    everywhere."""
    n_k, n_ak = len(k_grid), len(ak_grid)
    a_coeffs = np.zeros((n_k, n_ak, M_keep), dtype=complex)
    c_over_c0 = np.ones((n_k, n_ak))
    converged = np.zeros((n_k, n_ak), dtype=bool)

    for i, k_val in enumerate(k_grid):
        lam = 2.0 * np.pi / k_val
        for j, ak in enumerate(ak_grid):
            p_node = default_forcing(k_val, float(ak)) if p is None else p
            sol = _solve_with_validation(lam, float(ak), p_node, M_solve,
                                         n_steps, which_class)
            if sol is not None:
                a_coeffs[i, j, :] = gauge_fix(sol.a[:M_keep])
                c_over_c0[i, j] = sol.c / sol.c0
                converged[i, j] = sol.residual_norm < 1e-3
            else:
                a_coeffs[i, j, 0] = -ak + 0j
            if verbose:
                tag = "ok" if converged[i, j] else "LIN"
                print(f"  [{tag}] lam={lam * 100:7.3f} cm ak={ak:.3f}  "
                      f"|a2/a1|={abs(a_coeffs[i, j, 1]) / max(abs(a_coeffs[i, j, 0]), 1e-12):.3f}")
    return FM98Table(k_grid=np.asarray(k_grid, dtype=float),
                     ak_grid=np.asarray(ak_grid, dtype=float),
                     a_coeffs=a_coeffs, c_over_c0=c_over_c0,
                     converged=converged, M_keep=M_keep,
                     p=float("nan") if p is None else p)
