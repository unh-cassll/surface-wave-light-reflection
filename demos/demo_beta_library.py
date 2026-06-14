"""
The empirical bound-fraction library: beta(k) curves from the batch cube
reductions, colored by wind speed, plus the wind-interpolated curves used
by the synthesis.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (bound_fraction_for_wind, run_conditions,
                    smooth_beta_curve)

LIB = Path(__file__).parent / "output" / "asit_beta_library"
OUT = Path(__file__).parent / "output"
DATA = Path("/mnt/DATA/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = DATA / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = DATA / "ASIT2019_supporting_environmental_observations.nc"


def main():
    paths = sorted(LIB.glob("*.npz"))
    print(f"{len(paths)} reductions")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4), dpi=110)

    entries = []
    for p in paths:
        d = np.load(p)
        entries.append((float(d["U10"]), np.asarray(d["k"]),
                        np.asarray(d["beta_obs"])))
    entries.sort(key=lambda e: e[0])
    U_lo, U_hi = entries[0][0], entries[-1][0]

    ax = axes[0]
    cmap = plt.cm.viridis
    for U10, k, b in entries:
        try:
            bs = smooth_beta_curve(k, b)
        except ValueError:
            continue
        ax.semilogx(k, bs, color=cmap((U10 - U_lo) / (U_hi - U_lo)),
                    lw=1.0, alpha=0.8)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(U_lo, U_hi))
    plt.colorbar(sm, ax=ax, label="U10 [m/s]")
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("off-shell (bound) fraction")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"measured beta(k), {len(entries)} runs")
    ax.grid(alpha=0.3)

    ax = axes[1]
    K = np.geomspace(3, 1300, 200)
    for U10 in [3.0, 6.0, 9.0, 13.0, 17.0]:
        beta_fn = bound_fraction_for_wind(LIB, U10)
        ax.semilogx(K, beta_fn(K),
                    color=cmap((U10 - U_lo) / (U_hi - U_lo)), lw=1.6,
                    label=f"U10 = {U10:.0f}")
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("beta(k) for synthesis")
    ax.set_ylim(0, 1.05)
    ax.set_title("wind-interpolated bound fraction "
                 "(smoothed, tapered, capped)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # binding onset vs sea state: beta at a short-gravity wavenumber
    if STATS.exists():
        conds = run_conditions(STATS, ENV)
        Om_all = conds["inverse_wave_age"]
        ax = axes[2]
        pts = []
        for p in paths:
            d = np.load(p)
            run = int(d["run"])
            try:
                b30 = float(np.interp(np.log(30.0),
                                      np.log(np.asarray(d["k"])),
                                      smooth_beta_curve(
                                          np.asarray(d["k"]),
                                          np.asarray(d["beta_obs"]))))
            except ValueError:
                continue
            if np.isfinite(Om_all[run]):
                pts.append((Om_all[run], b30, float(d["U10"])))
        pts = np.array(pts)
        sc = ax.scatter(pts[:, 0], pts[:, 1], c=pts[:, 2], cmap=cmap,
                        s=18)
        plt.colorbar(sc, ax=ax, label="U10 [m/s]")
        ax.set_xlabel("inverse wave age  U10 / cp")
        ax.set_ylabel("beta at k = 30 rad/m")
        ax.set_ylim(0, 1.05)
        ax.set_title("binding strength vs sea state "
                     f"({pts.shape[0]} runs)")
        ax.grid(alpha=0.3)

    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_beta_library.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
