"""
Time-evolving sea surface with per-frame polarized reflection (the Python
counterpart of produce_simulated_sea_surface_modeled_reflection.m).
Saves an elevation/Stokes-I strip across frames and an animated GIF.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (CameraGeometry, SubpixelSlopes, generate_sea_surface,
                    make_rayleigh_sky, render_facet_stokes_stack)

OUT = Path(__file__).parent / "output"


def main():
    L, N, U10 = 64.0, 128, 7.0
    fps = 4.0
    n_frames = 8
    times = np.arange(n_frames) / fps

    rng = np.random.default_rng(0)
    surf = generate_sea_surface(L=L, N=N, U10=U10, times=times, rng=rng)
    dx = surf.info["dx"]
    print(f"Hs = {surf.info['Hs_realized']:.2f} m, dx = {dx} m, "
          f"{n_frames} frames at {fps} Hz")

    cam = CameraGeometry(incidence_deg=30.0, azimuth_deg=0.0, height_m=100.0)
    sky = make_rayleigh_sky(sun_zenith_deg=40.0, sun_azimuth_deg=90.0)
    sub = SubpixelSlopes(sigma_a2=surf.info["sigma_a2_cut"],
                         sigma_c2=surf.info["sigma_c2_cut"])

    print("rendering stack...")
    S = render_facet_stokes_stack(surf.eta, dx, slope_x=surf.slope_x,
                                  slope_y=surf.slope_y, camera=cam, sky=sky,
                                  subpixel=sub, n_subpixel=16,
                                  rng=np.random.default_rng(1))

    # Strip of frames: eta on top, Stokes I below
    show = [0, n_frames // 2, n_frames - 1]
    fig, axes = plt.subplots(2, len(show), figsize=(11, 7), dpi=110)
    e_max = 3 * surf.eta.std()
    I_lim = np.nanpercentile(S[:, :, 0, :], [2, 98])
    for col, it in enumerate(show):
        axes[0, col].imshow(surf.eta[:, :, it], cmap="RdBu_r",
                            vmin=-e_max, vmax=e_max)
        axes[0, col].set_title(f"eta, t = {times[it]:.2f} s")
        axes[1, col].imshow(S[:, :, 0, it], cmap="gray",
                            vmin=I_lim[0], vmax=I_lim[1])
        axes[1, col].set_title(f"Stokes I, t = {times[it]:.2f} s")
        for ax in axes[:, col]:
            ax.set_xticks([])
            ax.set_yticks([])
    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_time_evolution.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")

    # Animated GIF of Stokes I
    try:
        from matplotlib import animation
        fig2, ax2 = plt.subplots(figsize=(5, 5), dpi=100)
        im = ax2.imshow(S[:, :, 0, 0], cmap="gray",
                        vmin=I_lim[0], vmax=I_lim[1])
        ax2.set_xticks([])
        ax2.set_yticks([])

        def update(i):
            im.set_data(S[:, :, 0, i])
            ax2.set_title(f"Stokes I, t = {times[i]:.2f} s")
            return [im]

        ani = animation.FuncAnimation(fig2, update, frames=n_frames)
        gif = OUT / "demo_time_evolution.gif"
        ani.save(gif, writer=animation.PillowWriter(fps=int(fps)))
        print(f"saved -> {gif}")
    except Exception as exc:
        print(f"GIF skipped: {exc}")


if __name__ == "__main__":
    main()
