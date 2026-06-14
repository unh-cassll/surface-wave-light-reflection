"""
Slope statistics of linear vs hybrid surfaces against Cox & Munk (1954):
mean-square slope components vs wind, slope PDFs with their Gaussian
references, and the skewness measures that quantify what the synthesis
does and does not yet reproduce.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import generate_hybrid_surface, generate_sea_surface
from seapol.fm98 import FM98Table

OUT = Path(__file__).parent / "output"
TABLE_CACHE = OUT / "fm98_table_deep.npz"

L, N = 2.0, 2048                       # dx ~ 1 mm
WINDS = [5.0, 7.0, 10.0]
K_SPLIT = 2 * np.pi / 0.5
MTF = 6.5


def cox_munk(U10):
    """Clean-surface Cox-Munk fits: (sigma_u^2, sigma_c^2, c21, c03)."""
    return (3.16e-3 * U10,
            0.003 + 1.92e-3 * U10,
            0.01 - 0.0086 * U10,
            0.04 - 0.033 * U10)


def main():
    OUT.mkdir(exist_ok=True)
    table = FM98Table.load(TABLE_CACHE)

    rows = []
    pdfs = {}
    for U10 in WINDS:
        lin = generate_sea_surface(L, N, U10,
                                   rng=np.random.default_rng(11))
        hyb = generate_hybrid_surface(L=L, N=N, U10=U10, table=table,
                                      long_wave_mtf=MTF, k_split=K_SPLIT,
                                      rng=np.random.default_rng(11))
        sa2, sc2 = lin.info["sigma_a2_cut"], lin.info["sigma_c2_cut"]
        for tag, s in [("lin", lin), ("hyb", hyb)]:
            sx, sy = s.slope_x, s.slope_y
            skew_x = float(np.mean((sx - sx.mean()) ** 3) / sx.std() ** 3)
            skew_y = float(np.mean((sy - sy.mean()) ** 3) / sy.std() ** 3)
            kurt_x = float(np.mean((sx - sx.mean()) ** 4) / sx.var() ** 2
                           - 3.0)
            rows.append((U10, tag, sx.var() + sa2, sy.var() + sc2,
                         skew_x, skew_y, kurt_x))
        if U10 == 7.0:
            pdfs["lin"] = lin.slope_x
            pdfs["hyb"] = hyb.slope_x

    su_cm, sc_cm, c21, c03 = np.array([cox_munk(u) for u in WINDS]).T

    print(f"{'U10':>4} {'model':>6} | {'mss_up':>8} {'CM':>8} | "
          f"{'mss_cr':>8} {'CM':>8} | {'skew_up':>8} {'skew_cr':>8} "
          f"{'xkurt_up':>9}")
    for (U10, tag, mu, mc, kx, ky, qx) in rows:
        i = WINDS.index(U10)
        print(f"{U10:4.0f} {tag:>6} | {mu:8.4f} {su_cm[i]:8.4f} | "
              f"{mc:8.4f} {sc_cm[i]:8.4f} | {kx:+8.3f} {ky:+8.3f} "
              f"{qx:+9.3f}")
    print(f"\nCox-Munk (classical field reference; non-Gaussian shape "
          f"emerges from Elfouhaily + FM98, not fitted): "
          f"c21 = {dict(zip(WINDS, c21.round(3)))}, "
          f"c03 = {dict(zip(WINDS, c03.round(3)))}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), dpi=110)

    ax = axes[0]
    u_fine = np.linspace(2, 14, 50)
    ax.plot(u_fine, 3.16e-3 * u_fine, "k-", lw=1, label="CM upwind")
    ax.plot(u_fine, 0.003 + 1.92e-3 * u_fine, "k--", lw=1,
            label="CM crosswind")
    for tag, marker in [("lin", "o"), ("hyb", "s")]:
        sel = [r for r in rows if r[1] == tag]
        ax.plot([r[0] for r in sel], [r[2] for r in sel], marker + "-",
                color="C0" if tag == "lin" else "C3",
                label=f"{tag} upwind")
        ax.plot([r[0] for r in sel], [r[3] for r in sel], marker + "--",
                color="C0" if tag == "lin" else "C3", alpha=0.6,
                label=f"{tag} crosswind")
    ax.set_xlabel("U10 [m/s]")
    ax.set_ylabel("mean-square slope")
    ax.set_title("MSS (resolved + sub-grid) vs Cox-Munk")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1]
    edges = np.linspace(-5, 5, 121)
    centers = 0.5 * (edges[:-1] + edges[1:])
    for tag, color in [("lin", "C0"), ("hyb", "C3")]:
        s = pdfs[tag]
        z = (s - s.mean()) / s.std()
        hist, _ = np.histogram(z.ravel(), bins=edges, density=True)
        ax.semilogy(centers, hist, color=color, label=f"{tag} along-wind")
    ax.semilogy(centers, np.exp(-centers**2 / 2) / np.sqrt(2 * np.pi),
                "k:", label="Gaussian")
    ax.set_ylim(1e-6, 1)
    ax.set_xlabel("normalized along-wind slope")
    ax.set_ylabel("PDF")
    ax.set_title(f"slope PDF, U10 = 7 m/s")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = OUT / "demo_slope_statistics.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
