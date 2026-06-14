"""
Validate the forward Monte Carlo tracer: hemispherical reflectance vs
incidence angle for several wind speeds, compared with the flat-surface
Fresnel curve.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import effective_mueller_for_incident, generate_sea_surface
from seapol.polarization import fresnel_mueller

OUT = Path(__file__).parent / "output"


def main():
    U10_values = [2.0, 5.0, 10.0]
    theta_deg = np.array([10, 20, 30, 40, 50, 60, 70, 80])
    n_rays = 20_000

    fig, ax = plt.subplots(figsize=(8, 5), dpi=110)
    th_fine = np.linspace(0, 85, 200)
    R_flat = fresnel_mueller(np.cos(np.deg2rad(th_fine)), 1.34)[0][..., 0, 0]
    ax.plot(th_fine, R_flat, "k-", lw=2, label="Flat Fresnel (analytic)")

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(U10_values)))
    for U10, col in zip(U10_values, colors):
        surf = generate_sea_surface(L=32.0, N=64, U10=U10,
                                    rng=np.random.default_rng(0))
        R = []
        for th in theta_deg:
            res = effective_mueller_for_incident(
                surf.eta, surf.info["dx"], float(th), n_rays=n_rays,
                rng=np.random.default_rng(1))
            R.append(res["R_total"])
        ax.plot(theta_deg, R, "o-", color=col, label=f"MC, U10 = {U10:.0f} m/s")
        print(f"U10={U10:4.1f}: R = {np.array(R).round(4)}")

    ax.set_xlabel("Incidence angle [deg]")
    ax.set_ylabel("Hemispherical reflectance R")
    ax.set_title("Forward-MC reflectance vs flat Fresnel")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_mc_reflectance.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
