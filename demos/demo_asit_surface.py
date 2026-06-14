"""
End-to-end empirical synthesis: a sea surface generated from a measured
ASIT 2019 directional slope spectrum, validated against the same run's
omnidirectional spectrum, slope histogram, and MSS, then rendered as
polarimetric imagery.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

from seapol import (PinholeCamera, SubpixelSlopes, WaterBody,
                    generate_asit_surface, load_asit_run, make_clear_sky,
                    match_env_run, render_camera_image)

DATA = Path("/mnt/DATA/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = DATA / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = DATA / "ASIT2019_supporting_environmental_observations.nc"
OUT = Path(__file__).parent / "output"

RUN = 84                  # moderate wind with finite EC measurement
L, N = 2.9, 2048          # dx = 1.4 mm; k = 2.2 .. 2218 rad/m


def main():
    OUT.mkdir(exist_ok=True)
    run = RUN
    run_data = load_asit_run(STATS, run, env_path=ENV)
    U10 = run_data["U10"]
    print(f"run {run}, U10 = {U10:.2f} m/s")

    surf = generate_asit_surface(STATS, run, L=L, N=N, env_path=ENV,
                                 rng=np.random.default_rng(0))
    dx = surf.info["dx"]
    mss_real = surf.slope_x.var() + surf.slope_y.var()
    print(f"MSS: realized {mss_real:.4f} | grid target "
          f"{surf.info['mss_resolved']:.4f} | instrument "
          f"{sum(surf.info['mss_measured']):.4f}")

    # realized ring-integrated slope variance density vs measured
    Z = np.fft.fft2(surf.eta)
    kxg = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    K = np.hypot(kxg[None, :], kxg[:, None])
    dk = 2 * np.pi / L
    P_slope = (K**2) * (np.abs(Z) / N**2) ** 2   # slope power per bin
    k_edges = np.geomspace(2.2, 2000.0, 50)
    k_c = np.sqrt(k_edges[:-1] * k_edges[1:])
    D_real = np.empty(k_c.size)
    for i in range(k_c.size):
        ring = (K >= k_edges[i]) & (K < k_edges[i + 1])
        D_real[i] = P_slope[ring].sum() / (k_edges[i + 1] - k_edges[i])
    k_m, th_m, S_m = run_data["k"], run_data["theta"], run_data["S"]
    dth = np.diff(th_m).mean()
    D_meas = k_m * S_m.sum(0) * dth              # slope density per dk

    # slope PDFs vs the instrument histogram
    d = nc.Dataset(STATS)
    centers = np.array(d["slope_centers"][:])
    hist2d = np.array(d["slope_histogram_crosswind_upwind"][run])
    d.close()
    dc = np.diff(centers).mean()
    pdf_up_meas = hist2d.sum(axis=0)
    pdf_up_meas = pdf_up_meas / (pdf_up_meas.sum() * dc)
    pdf_cr_meas = hist2d.sum(axis=1)
    pdf_cr_meas = pdf_cr_meas / (pdf_cr_meas.sum() * dc)

    def pdf(s):
        h, edges = np.histogram(s.ravel(), bins=120,
                                range=(centers[0], centers[-1]),
                                density=True)
        return 0.5 * (edges[:-1] + edges[1:]), h

    print("rendering ...")
    sub = SubpixelSlopes(surf.info["sigma_a2_cut"],
                         surf.info["sigma_c2_cut"])
    cam = PinholeCamera(altitude_m=30.0, zenith_deg=40.0, azimuth_deg=0.0,
                        hfov_deg=2.2, img_shape=(360, 360))
    sky = make_clear_sky(45.0, 195.0, turbidity=0.15, I_sun=200.0)
    S_img = render_camera_image(surf.eta, dx, camera=cam, sky=sky,
                                slope_x=surf.slope_x,
                                slope_y=surf.slope_y,
                                subpixel=sub, n_subpixel=12,
                                shadowing=True, water=WaterBody(case=2),
                                rng=np.random.default_rng(2))

    fig = plt.figure(figsize=(13, 9.5), dpi=110)
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    ax.loglog(k_m, D_meas, "k-", lw=1.4, label="measured")
    ax.loglog(k_c, D_real, "r--", lw=1.1, label="realized")
    ax.axvline(k_m[-1], color="b", ls=":", lw=0.8)
    ax.text(k_m[-1] * 1.1, D_meas.max() * 0.3, "tail blend", rotation=90,
            fontsize=7)
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("slope variance density [1/(rad/m)]")
    ax.set_title("omnidirectional slope spectrum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = fig.add_subplot(gs[0, 1])
    for s, lab, color in [(surf.slope_x, "realized up", "r"),
                          (surf.slope_y, "realized cross", "m")]:
        c, h = pdf(s)
        ax.semilogy(c, h, color=color, lw=1.0, label=lab)
    ax.semilogy(centers, pdf_up_meas, "k-", lw=1.4, label="measured up")
    ax.semilogy(centers, pdf_cr_meas, "k--", lw=1.4,
                label="measured cross")
    ax.set_ylim(1e-4, 20)
    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel("slope")
    ax.set_ylabel("PDF")
    ax.set_title("slope PDFs (instrument histogram for context only:\n"
                 "scale-limited, not a fitting target)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 2])
    im = ax.imshow(surf.eta * 100, cmap="RdBu_r",
                   extent=[0, L * 100, 0, L * 100])
    ax.set_title(f"synthesized eta [cm], run {run}")
    ax.set_xlabel("x [cm]")
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = fig.add_subplot(gs[1, 0])
    vmax = np.nanpercentile(S_img[..., 0], 99.5)
    im = ax.imshow(S_img[..., 0], cmap="gray", vmin=0, vmax=vmax)
    ax.set_title("rendered Stokes I (clear + sun)")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = fig.add_subplot(gs[1, 1])
    dolp = np.hypot(S_img[..., 1], S_img[..., 2]) / S_img[..., 0]
    im = ax.imshow(dolp, cmap="viridis", vmin=0, vmax=1)
    ax.set_title("rendered DoLP")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = fig.add_subplot(gs[1, 2])
    zoom = slice(0, N // 8)
    im = ax.imshow(surf.slope_x[zoom, zoom], cmap="gray",
                   extent=[0, L / 8 * 100, 0, L / 8 * 100])
    ax.set_title("along-wind slope, 36 cm zoom")
    ax.set_xlabel("x [cm]")
    plt.colorbar(im, ax=ax, fraction=0.046)

    fig.suptitle(f"Surface synthesized from measured ASIT spectrum: "
                 f"run {run}, U10 = {U10:.1f} m/s, L = {L} m, "
                 f"dx = {dx * 1e3:.2f} mm", fontsize=12)
    plt.tight_layout()
    out = OUT / "demo_asit_surface.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
