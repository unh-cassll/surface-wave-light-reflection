"""
Measured wavenumber-frequency slope spectrum from one raw ASIT 2019 cube:
reduce Skw(f, kx, ky) to S(|k|, f) by streaming frequency planes, show the
observed dispersion diagram, and split each wavenumber's energy into
on-shell (free) and off-shell (bound/advected) parts to obtain the
empirical bound fraction beta_obs(k) and the wavenumber where bound waves
begin to dominate.  Results are cached for use by demo_kw_spectrum.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

from seapol.spectrum import angular_frequency

RAW = Path("/mnt/DATA/raw_ASIT2019_spectra")
CUBE = RAW / "ASIT_day009_2019_10_15_20_00_00_full_mean_dirspect.nc"
ENV = Path("/mnt/DATA/Dropbox/Professional/Github/E-PSS_paper/_data/"
           "ASIT2019_supporting_environmental_observations.nc")
OUT = Path(__file__).parent / "output"
F_STRIDE = 3
N_KBIN = 60
SHELL_TOL = 0.20          # fractional frequency tolerance around the shell


def run_wind_speed():
    """U10 for the cube's timestamp (2019-10-15 20:00 UTC)."""
    import datetime as dt
    d = nc.Dataset(ENV)
    t = np.array(d["t_seconds_since_January_1_1970"][:])
    U = np.array(d["EC_U_m_s"][:])
    d.close()
    target = dt.datetime(2019, 10, 15, 20, 0,
                         tzinfo=dt.timezone.utc).timestamp()
    order = np.argsort(np.abs(t - target))
    for i in order:               # nearest run with finite wind
        if np.isfinite(U[i]):
            return float(U[i]), abs(t[i] - target)
    return np.nan, np.inf


def main():
    OUT.mkdir(exist_ok=True)
    U10, dt_match = run_wind_speed()
    print(f"run wind speed U10 = {U10:.2f} m/s (timestamp match "
          f"{dt_match:.0f} s)")

    d = nc.Dataset(CUBE)
    f = np.array(d["f"][:])
    kx = np.array(d["kx"][:])
    ky = np.array(d["ky"][:])
    K = np.hypot(kx, ky)
    print(f"cube: {d['Skw'].shape}, k to {K.max():.0f} rad/m, "
          f"f to {f.max():.1f} Hz")

    k_edges = np.geomspace(3.0, 1400.0, N_KBIN + 1)
    idx = np.digitize(K.ravel(), k_edges) - 1
    valid = (idx >= 0) & (idx < N_KBIN)
    idx_v = idx[valid]
    counts = np.bincount(idx_v, minlength=N_KBIN)

    f_sel = np.arange(0, f.size, F_STRIDE)
    S_kf = np.zeros((N_KBIN, f_sel.size))
    print(f"streaming {f_sel.size} frequency planes ...")
    for j, i_f in enumerate(f_sel):
        plane = np.asarray(d["Skw"][i_f, :, :]).ravel()[valid]
        S_kf[:, j] = np.bincount(idx_v, weights=plane, minlength=N_KBIN)
        if j % 50 == 0:
            print(f"  {j}/{f_sel.size}")
    d.close()
    f_u = f[f_sel]
    k_c = np.sqrt(k_edges[:-1] * k_edges[1:])
    f_disp = angular_frequency(k_c) / (2 * np.pi)

    # Instrument noise floor: estimated per frequency from the
    # "forbidden" region (well above the free shell and faster than any
    # physical phase speed), then subtracted.  Without this, the floor
    # integrated over the large off-shell area inflates the bound
    # fraction and the slow-phase-speed energy several-fold.
    S_mean = S_kf / np.maximum(counts, 1)[:, None]
    KK, FF = np.meshgrid(k_c, f_u, indexing="ij")
    forb = ((FF > 1.5 * f_disp[:, None])
            & (FF > 3.5 * KK / (2 * np.pi)) & (FF > 1.5))
    noise_f = np.array([np.median(S_mean[forb[:, j], j])
                        if forb[:, j].sum() > 5 else np.nan
                        for j in range(f_u.size)])
    ok_f = np.isfinite(noise_f)
    noise_f = np.interp(np.arange(f_u.size), np.flatnonzero(ok_f),
                        noise_f[ok_f])
    S_kf = np.clip(S_kf - noise_f[None, :] * counts[:, None], 0.0, None)

    # on-shell (free) vs off-shell split per wavenumber; only where the
    # free shell lies inside the measured frequency band
    good_f = f_u > 0.25
    beta_obs = np.full(N_KBIN, np.nan)
    for i in range(N_KBIN):
        tot = S_kf[i, good_f].sum()
        if (tot <= 0 or counts[i] == 0
                or f_disp[i] > 0.85 * f_u.max()):
            continue
        on = good_f & (np.abs(f_u - f_disp[i])
                       < np.maximum(SHELL_TOL * f_disp[i],
                                    2 * np.diff(f_u).mean()))
        beta_obs[i] = 1.0 - S_kf[i, on].sum() / tot
    ok = np.isfinite(beta_obs)
    cross = k_c[ok][np.argmax(beta_obs[ok] > 0.5)] if np.any(
        beta_obs[ok] > 0.5) else np.nan
    print(f"bound fraction: beta(50 rad/m) = "
          f"{np.interp(50, k_c, beta_obs):.2f}, beta(300) = "
          f"{np.interp(300, k_c, beta_obs):.2f}; "
          f"crossover beta > 0.5 at k ~ {cross:.0f} rad/m")

    np.savez(OUT / "asit_kf_reduced.npz", k=k_c, f=f_u, S_kf=S_kf,
             beta_obs=beta_obs, counts=counts, noise_f=noise_f,
             U10=U10, shell_tol=SHELL_TOL)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), dpi=110)
    ax = axes[0]
    with np.errstate(divide="ignore"):
        lp = np.log10(np.maximum(S_kf, S_kf[S_kf > 0].min()))
    im = ax.pcolormesh(k_c, f_u, lp.T, cmap="magma",
                       vmin=lp.max() - 5, vmax=lp.max())
    ax.plot(k_c, f_disp, "c--", lw=1.2, label="linear dispersion")
    for c_ref, ls in [(2.0, ":"), (0.8, "-.")]:
        ax.plot(k_c, c_ref * k_c / (2 * np.pi), "w", ls=ls, lw=0.9,
                label=f"c = {c_ref} m/s")
    ax.set_xscale("log")
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("f [Hz]")
    ax.set_ylim(0, f_u.max())
    ax.set_title(f"ASIT measured slope k-f spectrum (U10 = {U10:.1f} m/s)")
    ax.legend(fontsize=7, loc="upper left")
    plt.colorbar(im, ax=ax, fraction=0.046, label="log10 S")

    ax = axes[1]
    ax.semilogx(k_c, beta_obs, "ko-", ms=3)
    if np.isfinite(cross):
        ax.axvline(cross, color="r", ls="--", lw=0.8,
                   label=f"bound dominance k ~ {cross:.0f} rad/m")
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("off-shell (bound) fraction")
    ax.set_ylim(0, 1)
    ax.set_title("empirical bound fraction beta_obs(k)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = OUT / "demo_asit_dispersion.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
