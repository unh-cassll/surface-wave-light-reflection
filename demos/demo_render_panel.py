"""
Render a polarization image of a wind-driven sea surface under a Rayleigh
sky and save Stokes I/Q/U, DoLP, AoP, and the input elevation as a panel.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (PinholeCamera, SubpixelSlopes, generate_sea_surface,
                    make_rayleigh_sky, render_camera_image)
from seapol.polarization import stokes_aop, stokes_dolp

OUT = Path(__file__).parent / "output"


def main():
    L, N, U10 = 32.0, 512, 5.0
    rng = np.random.default_rng(1)
    surf = generate_sea_surface(L=L, N=N, U10=U10, rng=rng)
    dx = surf.info["dx"]
    print(f"surface: Hs = {surf.info['Hs_realized']:.3f} m, "
          f"subpixel slope var = "
          f"{surf.info['sigma_a2_cut'] + surf.info['sigma_c2_cut']:.4f}")

    cam = PinholeCamera(altitude_m=300.0, zenith_deg=45.0, azimuth_deg=0.0,
                        hfov_deg=3.0, img_shape=(400, 400))
    sky = make_rayleigh_sky(sun_zenith_deg=45.0, sun_azimuth_deg=90.0)
    sub = SubpixelSlopes(sigma_a2=surf.info["sigma_a2_cut"],
                         sigma_c2=surf.info["sigma_c2_cut"])

    print("rendering...")
    S = render_camera_image(surf.eta, dx, camera=cam, sky=sky,
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=sub, n_subpixel=16, shadowing=True,
                            rng=rng)
    I, Q, U = S[..., 0], S[..., 1], S[..., 2]
    DoLP = stokes_dolp(S)
    AoP = np.rad2deg(stokes_aop(S))

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), dpi=110)
    panels = [
        (I, "Stokes I", "gray", None),
        (Q, "Stokes Q", "RdBu_r", "sym"),
        (U, "Stokes U", "RdBu_r", "sym"),
        (DoLP, "DoLP", "viridis", (0, 1)),
        (AoP, "AoP [deg]", "twilight", (-90, 90)),
        (surf.eta, "Input sea surface eta [m]", "RdBu_r", "sym"),
    ]
    for ax, (img, title, cmap, scale) in zip(axes.ravel(), panels):
        if scale == "sym":
            vmax = np.nanpercentile(np.abs(img), 99)
            im = ax.imshow(img, cmap=cmap, vmin=-vmax, vmax=vmax)
        elif scale is None:
            im = ax.imshow(img, cmap=cmap)
        else:
            im = ax.imshow(img, cmap=cmap, vmin=scale[0], vmax=scale[1])
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.045)

    fig.suptitle(f"Sea-surface polarization image: U10 = {U10} m/s, "
                 f"SZA = 45 deg, view zenith = {cam.zenith_deg} deg",
                 fontsize=12)
    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_panel.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
