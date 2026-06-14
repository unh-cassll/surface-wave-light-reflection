"""
Rendered Stokes imagery across sky types and water bodies: clear sky with
sun glint, partly cloudy, and overcast, each with Case 1 water-leaving
radiance.  A qualitative gap-finder for the radiometric realism of the
pipeline (glint dynamic range, cloud texture, DoLP suppression by
upwelling light).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (PinholeCamera, SubpixelSlopes, WaterBody,
                    generate_sea_surface, make_clear_sky,
                    make_overcast_sky, make_partly_cloudy_sky,
                    render_camera_image)
from seapol.polarization import stokes_dolp

OUT = Path(__file__).parent / "output"

L, N, U10 = 32.0, 512, 6.0
# camera sits at azimuth 0 looking back along -x; the specular sky
# direction therefore has azimuth ~180, so put the sun there for glint
SUN_ZEN, SUN_AZ = 40.0, 195.0


def main():
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(2)
    surf = generate_sea_surface(L=L, N=N, U10=U10, rng=rng)
    dx = surf.info["dx"]
    sub = SubpixelSlopes(surf.info["sigma_a2_cut"],
                         surf.info["sigma_c2_cut"])
    cam = PinholeCamera(altitude_m=200.0, zenith_deg=42.0, azimuth_deg=0.0,
                        hfov_deg=5.0, img_shape=(360, 360))
    water = WaterBody(case=1)

    # the direct beam enters via the analytic Cox-Munk glint term, not a
    # sun disk in the sky model (avoids Monte Carlo glitter speckle)
    skies = [
        ("clear + sun", make_clear_sky(SUN_ZEN, SUN_AZ, I_sky=1.0,
                                       turbidity=0.15),
         (SUN_ZEN, SUN_AZ, 50.0)),
        ("partly cloudy", make_partly_cloudy_sky(SUN_ZEN, SUN_AZ,
                                                 cloud_fraction=0.45,
                                                 cloud_brightness=4.0,
                                                 rng=np.random.default_rng(8)),
         None),
        ("overcast", make_overcast_sky(I_zenith=1.0), None),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(9.5, 13), dpi=110)
    for row, (name, sky, glint) in enumerate(skies):
        print(f"rendering {name} ...")
        S = render_camera_image(surf.eta, dx, camera=cam, sky=sky,
                                slope_x=surf.slope_x, slope_y=surf.slope_y,
                                subpixel=sub, n_subpixel=12,
                                shadowing=True, water=water,
                                sun_glint=glint,
                                rng=np.random.default_rng(3))
        I = S[..., 0]
        DoLP = stokes_dolp(S)

        ax = axes[row, 0]
        vmax = np.nanpercentile(I, 99.5)
        im = ax.imshow(I, cmap="gray", vmin=0, vmax=vmax)
        ax.set_title(f"{name}: I (clipped at 99.5%)")
        plt.colorbar(im, ax=ax, fraction=0.046)

        ax = axes[row, 1]
        im = ax.imshow(DoLP, cmap="viridis", vmin=0, vmax=1)
        ax.set_title(f"{name}: DoLP")
        plt.colorbar(im, ax=ax, fraction=0.046)

        med = np.nanmedian(DoLP)
        print(f"  median DoLP = {med:.3f}, I p50/p99.5 = "
              f"{np.nanpercentile(I, 50):.3f}/{vmax:.3f}")

    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"Sky/water gallery: U10 = {U10} m/s, Case 1 water, "
                 f"view zenith {cam.zenith_deg} deg", fontsize=12)
    plt.tight_layout()
    out = OUT / "demo_sky_water_gallery.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
