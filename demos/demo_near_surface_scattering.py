"""
Near-surface scattering: what the in-water polarized Monte Carlo adds
to the rendered scene.

Top row: S0 with (a) specular reflection only -- the shiny-mirror look,
(b) the first-order isotropic WaterBody term, (c) the directional
upwelling-radiance table (sky + transmitted sun beam scattered by the
water column).  Middle row: the corresponding DoLP images (in-water
light damps the reflected polarization, more where the water is
bright).  Bottom row: the sub-surface light field itself, L_u(mu_w,
phi), and its azimuth structure (anti-solar brightening of the
scattered beam).

The table is cached in output/ -- delete it to rebuild (e.g. after
changing the water type).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (PinholeCamera, SubpixelSlopes, WATER_TYPES, WaterBody,
                    build_upwelling_table, generate_sea_surface, load_table,
                    make_clear_sky, render_camera_image, save_table)
from seapol.polarization import stokes_dolp

OUT = Path(__file__).parent / "output"
TABLE_CACHE = OUT / "upwelling_coastal_550.npz"

L, N, U10 = 32.0, 512, 6.5
SUN_ZEN, SUN_AZ = 45.0, 195.0
E_SUN = 30.0
WATER_TYPE = "coastal_case2"


def get_table(sky):
    if TABLE_CACHE.exists():
        print(f"using cached table {TABLE_CACHE.name}")
        return load_table(TABLE_CACHE)
    print("building upwelling table (5e5 photons) ...")
    col = WATER_TYPES[WATER_TYPE].column()
    tab = build_upwelling_table(sky, col, sun=(SUN_ZEN, SUN_AZ, E_SUN),
                                n_photons=500_000,
                                rng=np.random.default_rng(0))
    save_table(TABLE_CACHE, tab)
    return tab


def main():
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(2)
    surf = generate_sea_surface(L=L, N=N, U10=U10, rng=rng)
    dx = surf.info["dx"]
    sub = SubpixelSlopes.from_cox_munk(U10, surf.info["k_cutoff"])
    cam = PinholeCamera(altitude_m=200.0, zenith_deg=42.0, azimuth_deg=0.0,
                        hfov_deg=5.0, img_shape=(360, 360))
    sky = make_clear_sky(SUN_ZEN, SUN_AZ, I_sky=1.0, turbidity=0.15)
    glint = (SUN_ZEN, SUN_AZ, E_SUN)
    tab = get_table(sky)
    print(f"table: E_down={tab.info.get('E_down', float('nan')):.2f}, "
          f"E_u={tab.info.get('E_u', float('nan')):.3f}")

    waters = [("specular only", None),
              ("first-order WaterBody", WaterBody(case=2)),
              ("scattering table", tab)]
    renders = []
    for name, water in waters:
        print(f"rendering: {name} ...")
        S = render_camera_image(surf.eta, dx, camera=cam, sky=sky,
                                slope_x=surf.slope_x, slope_y=surf.slope_y,
                                subpixel=sub, n_subpixel=12, water=water,
                                sun_glint=glint,
                                rng=np.random.default_rng(5))
        renders.append((name, S))

    fig = plt.figure(figsize=(13.5, 13), dpi=110)
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.75])
    # scale to the non-glint background so the water term is visible;
    # the sun glint saturates (its dynamic range is the point of
    # demo_glint_foam_current, not this one)
    s0_max = 1.3 * np.nanquantile(renders[2][1][..., 0], 0.85)
    for col_i, (name, S) in enumerate(renders):
        ax = fig.add_subplot(gs[0, col_i])
        im = ax.imshow(S[..., 0], cmap="gray", vmin=0, vmax=s0_max)
        ax.set_title(f"{name}\nS0 (shared scale, glint saturated)")
        ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046)

        ax = fig.add_subplot(gs[1, col_i])
        im = ax.imshow(stokes_dolp(S), cmap="viridis", vmin=0, vmax=1)
        ax.set_title("DoLP")
        ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046)

    # sub-surface light field
    ax = fig.add_subplot(gs[2, 0:2])
    S_tab = np.asarray(tab.S)
    mu_c = 0.5 * (np.asarray(tab.mu_edges[:-1])
                  + np.asarray(tab.mu_edges[1:]))
    phi_c = 0.5 * (np.asarray(tab.phi_edges[:-1])
                   + np.asarray(tab.phi_edges[1:]))
    im = ax.pcolormesh(np.rad2deg(phi_c), mu_c, S_tab[..., 0],
                       shading="nearest", cmap="magma")
    ax.axvline(SUN_AZ - 360.0, color="c", ls="--", lw=1,
               label="sun azimuth")
    ax.set_xlabel("azimuth of in-water propagation [deg]")
    ax.set_ylabel(r"$\mu_w$ (upward)")
    ax.set_title(r"sub-surface upwelling $L_u(\mu_w, \phi)$"
                 f"  [{WATER_TYPE}]")
    ax.legend(loc="upper right", fontsize=8)
    plt.colorbar(im, ax=ax, fraction=0.03)

    ax = fig.add_subplot(gs[2, 2])
    ax.plot(np.rad2deg(phi_c), S_tab[..., 0].mean(axis=0), "k-")
    ax.axvline(SUN_AZ - 360.0, color="c", ls="--", lw=1)
    ax.set_xlabel("azimuth [deg]")
    ax.set_ylabel(r"$\langle L_u\rangle_{\mu}$")
    ax.set_title("anti-solar brightening")
    ax.grid(alpha=0.3)

    fig.suptitle("Near-surface scattering: from mirror to water "
                 f"(U10 = {U10} m/s, {WATER_TYPE})", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "demo_near_surface_scattering.png",
                bbox_inches="tight")
    print(f"wrote {OUT / 'demo_near_surface_scattering.png'}")

    for name, S in renders:
        s0 = S[..., 0]
        print(f"{name:24s} S0 mean {np.nanmean(s0):.4f}  "
              f"cv {np.nanstd(s0) / np.nanmean(s0):.3f}  "
              f"DoLP {np.nanmean(stokes_dolp(S)):.3f}")


if __name__ == "__main__":
    main()
