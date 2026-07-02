"""
In-water polarized scattering in the Monte Carlo tracer: energy budget
vs single-scattering albedo, the angular shape of the emergent
water-leaving radiance against the first-order isotropic model of
seapol.water, and the damping of above-water DoLP by the (weakly
polarized) water-leaving component.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (WaterBody, WaterOptics, effective_mueller_for_incident,
                    generate_sea_surface, water_leaving_stokes)
from seapol.polarization import fresnel_mueller

OUT = Path(__file__).parent / "output"

THETA_I = 40.0
N_WATER = 1.34


def run(eta, dx, water, n_rays, seed, **kw):
    return effective_mueller_for_incident(
        eta, dx, THETA_I, n_rays=n_rays, water=water, max_bounces=80,
        n_theta_bins=18, rng=np.random.default_rng(seed), **kw)


def main():
    rng = np.random.default_rng(0)
    surf = generate_sea_surface(L=32.0, N=64, U10=5.0, rng=rng)
    eta, dx = surf.eta, surf.info["dx"]
    flat = np.zeros((16, 16))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), dpi=110)

    # --- energy budget vs single-scattering albedo (rough surface)
    c_att = 0.15
    omegas = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.95])
    budget = []
    for om in omegas:
        w = WaterOptics(absorption=c_att * (1 - om), scattering=c_att * om)
        r = run(eta, dx, w, n_rays=20_000, seed=1)
        budget.append([r["R_glint"], r["R_water"], r["A_total"],
                       r["W_total"]])
        print(f"omega0={om:.2f}: glint={r['R_glint']:.4f} "
              f"water={r['R_water']:.4f} absorbed={r['A_total']:.4f} "
              f"unresolved={r['W_total']:.5f}")
    budget = np.array(budget).T
    ax = axes[0]
    labels = ["surface glint", "water-leaving", "absorbed", "unresolved"]
    ax.stackplot(omegas, budget, labels=labels, alpha=0.85,
                 colors=["#c8b273", "#3d85c6", "#1b3a57", "#999999"])
    ax.set_xlabel(r"single-scattering albedo $\omega_0$")
    ax.set_ylabel("fraction of incident energy")
    ax.set_title(f"Energy budget, U10 = 5 m/s, "
                 f"$\\theta_i$ = {THETA_I:.0f}$^\\circ$")
    ax.set_ylim(0, 1)
    ax.legend(loc="center left", fontsize=8)

    # --- water-leaving radiance shape vs first-order model (flat surface)
    water = WaterOptics(absorption=0.10, scattering=0.05)
    r = run(flat, 0.5, water, n_rays=200_000, seed=2)
    th_c = 0.5 * (r["theta_edges"][:-1] + r["theta_edges"][1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        # energy per bin / (solid angle * cos theta): radiance, comparable
        # to the first-order model (without cos theta it is only the
        # projected flux and fakes a limb darkening)
        L_mc = (r["M_eff_water"][..., 0, 0].sum(axis=1)
                / (r["bin_solid_angle"].sum(axis=1) * np.cos(th_c)))
    S_w = water_leaving_stokes(
        np.stack([np.sin(th_c), np.zeros_like(th_c), np.cos(th_c)], axis=-1),
        np.broadcast_to([0.0, 0.0, 1.0], (th_c.size, 3)),
        WaterBody(case=1), n_water=N_WATER)
    ok = L_mc > 0
    ax = axes[1]
    ax.plot(np.rad2deg(th_c[ok]), L_mc[ok] / L_mc[ok][0], "o-",
            label="MC water-leaving")
    ax.plot(np.rad2deg(th_c), S_w[:, 0] / S_w[0, 0], "k--",
            label="first-order isotropic model")
    ax.set_xlabel("view zenith angle [deg]")
    ax.set_ylabel("radiance / nadir radiance")
    ax.set_title("Water-leaving radiance shape, flat surface")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- DoLP of up-escaping light: glint only vs with water-leaving
    r = run(eta, dx, water, n_rays=200_000, seed=3)
    phi_c = 0.5 * (r["phi_edges"][:-1] + r["phi_edges"][1:])
    fwd = np.abs(phi_c) < np.pi / 3          # forward (specular) azimuths
    S_g = r["M_eff_glint"][:, fwd][..., :, 0].sum(axis=1)
    S_t = (r["M_eff_glint"] + r["M_eff_water"])[:, fwd][..., :, 0].sum(axis=1)

    def dolp(S):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.hypot(S[:, 1], S[:, 2]) / np.where(S[:, 0] > 0,
                                                         S[:, 0], np.nan)
    ax = axes[2]
    ax.plot(np.rad2deg(th_c), dolp(S_g), "o-", label="glint only")
    ax.plot(np.rad2deg(th_c), dolp(S_t), "s-",
            label="glint + water-leaving")
    ax.set_xlabel("view zenith angle [deg]")
    ax.set_ylabel("DoLP")
    ax.set_title("DoLP damping by water-leaving light")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    R_F = fresnel_mueller(np.array(np.cos(np.deg2rad(THETA_I))),
                          N_WATER)[0][0, 0]
    print(f"flat Fresnel R({THETA_I:.0f} deg) = {R_F:.4f}; "
          f"MC glint = {r['R_glint']:.4f}, water-leaving = "
          f"{r['R_water']:.4f}, mean scatters = {r['mean_scatters']:.2f}")

    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_mc_water_body.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
