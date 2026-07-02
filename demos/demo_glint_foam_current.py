"""
Rendered light inheriting the surface physics: Gram-Charlier glint
asymmetry (the Cox-Munk upwind/downwind glitter difference), whitecap
foam riding the steepest facets, and the current Doppler in the k-f
dispersion ridge.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (Foam, PinholeCamera, SubpixelSlopes,
                    generate_sea_surface, make_clear_sky, monahan_coverage,
                    render_camera_image, sun_glint_stokes)
from seapol.diagnostics import kf_slope_spectrum
from seapol.polarization import stokes_dolp

OUT = Path(__file__).parent / "output"

U10 = 12.0


def main():
    fig = plt.figure(figsize=(15, 8.5), dpi=110)

    # --- glint transect: Gaussian vs Gram-Charlier slope PDF
    k_cut = 2 * np.pi / 0.5
    gauss = SubpixelSlopes.from_spectrum(U10, k_cut)
    gc = SubpixelSlopes.from_cox_munk(U10, k_cut)
    th_v = np.deg2rad(np.linspace(-60, 60, 241))   # +: down-sun (+x)
    d_out = np.stack([np.sin(th_v), np.zeros_like(th_v), np.cos(th_v)],
                     axis=-1)
    zeros = np.zeros(th_v.shape)
    ax = fig.add_subplot(2, 3, 1)
    for sub, ls, lbl in ((gauss, "--", "Gaussian"),
                         (gc, "-", "Gram-Charlier (Cox-Munk)")):
        S = sun_glint_stokes(d_out, zeros, zeros, sub, 30.0, 0.0, 1.0)
        ax.semilogy(np.rad2deg(th_v), S[:, 0], ls, label=lbl)
    # specular direction: sun at +30 deg zenith mirrors to view -30 deg
    ax.axvline(-30.0, color="k", lw=0.5, alpha=0.5)
    ax.set_xlabel("view zenith toward +x (sun azimuth) [deg]")
    ax.set_ylabel("glint radiance / E_sun")
    ax.set_title(f"Principal-plane glint, U10 = {U10:.0f} m/s\n"
                 "(skewed slopes shift and tilt the glitter lobe)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # asymmetry summary
    S_gc = sun_glint_stokes(d_out, zeros, zeros, gc, 30.0, 0.0, 1.0)[:, 0]
    i_pk = int(np.argmax(S_gc))
    print(f"GC glint peak at {np.rad2deg(th_v[i_pk]):+.1f} deg "
          f"(Gaussian specular = -30.0); c03 = {gc.c03:+.3f}")

    # --- rendered images with and without foam
    rng = np.random.default_rng(0)
    surf = generate_sea_surface(L=64.0, N=512, U10=U10, rng=rng)
    sub = SubpixelSlopes.from_cox_munk(U10, surf.info["k_cutoff"])
    cam = PinholeCamera(altitude_m=60.0, zenith_deg=50.0, azimuth_deg=15.0,
                        hfov_deg=14.0, img_shape=(300, 300))
    sky = make_clear_sky(40.0, 200.0, turbidity=0.15)
    cov = monahan_coverage(U10)
    print(f"Monahan whitecap coverage at U10 = {U10:.0f}: {cov*100:.2f}%")
    common = dict(camera=cam, sky=sky, slope_x=surf.slope_x,
                  slope_y=surf.slope_y, subpixel=sub, n_subpixel=8,
                  shadowing=True, sun_glint=(40.0, 200.0, 30.0))
    S0 = render_camera_image(surf.eta, surf.info["dx"],
                             rng=np.random.default_rng(1), **common)
    S1 = render_camera_image(surf.eta, surf.info["dx"],
                             foam=Foam(coverage=cov),
                             rng=np.random.default_rng(1), **common)

    vmax = np.nanpercentile(S1[..., 0], 99.5)
    for i, (S, tag) in enumerate(((S0, "no foam"),
                                  (S1, f"foam W = {cov*100:.1f}%"))):
        ax = fig.add_subplot(2, 3, 2 + i)
        ax.imshow(S[..., 0], cmap="gray", vmin=0, vmax=vmax)
        ax.set_title(f"Stokes I, {tag}")
        ax.set_xticks([])
        ax.set_yticks([])

    ax = fig.add_subplot(2, 3, 5)
    ax.imshow(stokes_dolp(S1), cmap="viridis", vmin=0, vmax=1)
    ax.set_title("DoLP with foam (whitecaps depolarize)")
    ax.set_xticks([])
    ax.set_yticks([])
    print(f"mean DoLP: {np.nanmean(stokes_dolp(S0)):.3f} -> "
          f"{np.nanmean(stokes_dolp(S1)):.3f} with foam")

    # --- current Doppler in the k-f ridge
    N, L = 128, 32.0
    FS, NT = 12.0, 256
    ax = fig.add_subplot(2, 3, 6)
    for cur, color, lbl in ((None, "c", "U = 0"),
                            ((1.0, 0.0), "r", "U = 1 m/s downwind")):
        buf = np.empty((32, N, NT))

        def cb(it, t, e, gx, gy):
            buf[:, :, it] = e[:32]
        generate_sea_surface(L, N, 8.0, times=np.arange(NT) / FS,
                             current=cur, compute_slopes=False,
                             frame_callback=cb,
                             rng=np.random.default_rng(2))
        kx, f, P = kf_slope_spectrum(buf, L / N, FS, n_seg=2)
        pos = kx > 0
        ridge = f[np.argmax(P[pos], axis=1)]
        ax.plot(kx[pos], ridge, color + ".", ms=3, label=lbl)
    from seapol.spectrum import angular_frequency
    kk = np.linspace(0.4, kx.max(), 100)
    ax.plot(kk, angular_frequency(kk) / (2 * np.pi), "k-", lw=1,
            label="free dispersion")
    ax.plot(kk, (angular_frequency(kk) + kk * 1.0) / (2 * np.pi), "k--",
            lw=1, label="+ k U / 2 pi")
    ax.set_xlabel("kx [rad/m]")
    ax.set_ylabel("ridge f [Hz]")
    ax.set_ylim(0, FS / 2)
    ax.set_title("Dispersion ridge with uniform current")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_glint_foam_current.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
