"""
Capstone: the full seapol pipeline driven end-to-end by ASIT 2019
measurements for a single run.

    measured S(k, theta)  ->  Psi(kx, ky)            (empirical spectrum)
    measured beta(k; U10) ->  free/bound partition   (cube-reduction
                              library when present; scalar ramp fallback)
    resolved Psi          ->  carrier-speed weights  (bound_speed="spectrum")
    + FM98 bound harmonics, long-wave MTF binding, orbital advection
    -> time-evolving surface -> polarimetric imagery (analytic glint)
    -> slope Q(nu) validated against the same run's measured Q(nu, theta).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

from seapol import (PinholeCamera, SubpixelSlopes, WaterBody,
                    bound_fraction_for_conditions, generate_hybrid_surface,
                    load_asit_run, make_clear_sky, psi_from_asit,
                    render_camera_image, run_conditions)
from seapol.fm98 import FM98Table

DATA = Path("/home/nathanlaxague/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = DATA / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = DATA / "ASIT2019_supporting_environmental_observations.nc"
OUT = Path(__file__).parent / "output"
LIB = OUT / "asit_beta_library"

RUN = 84
L, N = 2.9, 512                # dx = 5.7 mm, k_nyq = 555 rad/m
FS, NT = 36.0, 720             # 20 s at 36 Hz (streamed transects)
N_ROWS = 48
K_SPLIT = 2 * np.pi / 0.5


from seapol.diagnostics import kf_slope_spectrum, q_nu, slow_nu_fraction


def main():
    OUT.mkdir(exist_ok=True)
    run_data = load_asit_run(STATS, RUN, env_path=ENV)
    U10 = run_data["U10"]
    print(f"run {RUN}, U10 = {U10:.2f} m/s")

    psi = psi_from_asit(run_data)
    # bound fraction keyed by sea state (wind AND inverse wave age) when
    # a measured beta library exists; scalar ramp fallback otherwise
    if LIB.is_dir() and any(LIB.glob("*.npz")):
        Om = float(run_conditions(STATS, ENV)["inverse_wave_age"][RUN])
        print(f"inverse wave age Omega = {Om:.2f}")
        beta = bound_fraction_for_conditions(LIB, STATS, ENV, U10,
                                             inverse_wave_age=Om)
        beta_tag = "measured beta(k; U10, Omega)"
    else:
        beta = 0.9
        beta_tag = f"default ramp beta_max = {beta} (no beta library)"
    print(f"bound fraction: {beta_tag}")
    table = FM98Table.load(OUT / "fm98_table_deep.npz")

    print(f"hybrid record (streamed): L = {L} m, N = {N}, {NT} frames "
          f"at {FS} Hz ...")
    buf = np.empty((N_ROWS, N, NT))

    def cb(it, t, e, gx, gy):
        buf[:, :, it] = e[:N_ROWS]

    surf = generate_hybrid_surface(
        L=L, N=N, U10=U10, times=np.arange(NT) / FS, table=table,
        long_wave_mtf=6.5, k_split=K_SPLIT, orbital_advection=True,
        psi_override=psi, bound_fraction=beta, bound_speed="spectrum",
        frame_callback=cb, rng=np.random.default_rng(11), verbose=True)
    dx = surf.info["dx"]
    mss_real = surf.slope_x.var() + surf.slope_y.var()
    print(f"MSS realized {mss_real:.4f} (frame 0) | stored run mss "
          f"{run_data['mss_upwind'] + run_data['mss_crosswind']:.4f} "
          f"(context only)")

    # Q(nu) against the same run's measured inverse-phase-speed spectrum
    kx, f, P = kf_slope_spectrum(buf, dx, FS, n_seg=3)
    nu_s, Q_s = q_nu(kx, f, P)
    d = nc.Dataset(STATS)
    nu_m = np.array(d["nu_s_m"][:])
    th = np.array(d["theta_rad"][:])
    Q_m = np.array(d["Qs_nu_theta"][RUN]).sum(axis=0) * np.diff(th).mean()
    d.close()
    for tag, nn, qq in [("synthetic", nu_s, Q_s), ("measured", nu_m, Q_m)]:
        dq = np.gradient(nn)
        frac = np.nansum((qq * dq)[nn > 2.0]) / np.nansum(qq * dq)
        print(f"{tag}: slope-Q fraction at nu > 2 s/m = {frac:.3f}")

    # rendered movie frames: a short stack regenerated with the same
    # seed (deterministic draw) so the streamed record stays stack-free
    print("rendering movie frames ...")
    movie_times = np.arange(0, 144, 6) / FS
    surf_m = generate_hybrid_surface(
        L=L, N=N, U10=U10, times=movie_times, table=table,
        long_wave_mtf=6.5, k_split=K_SPLIT, orbital_advection=True,
        psi_override=psi, bound_fraction=beta, bound_speed="spectrum",
        rng=np.random.default_rng(11))
    sub = SubpixelSlopes(surf.info["sigma_a2_cut"],
                         surf.info["sigma_c2_cut"])
    cam = PinholeCamera(altitude_m=30.0, zenith_deg=42.0, azimuth_deg=0.0,
                        hfov_deg=2.4, img_shape=(280, 280))
    sky = make_clear_sky(45.0, 195.0, turbidity=0.15)
    frames = []
    f_idx = range(movie_times.size)
    for it in f_idx:
        S = np.asarray(
            render_camera_image(
                surf_m.eta[:, :, it], dx, camera=cam, sky=sky,
                slope_x=surf_m.slope_x[:, :, it],
                slope_y=surf_m.slope_y[:, :, it],
                subpixel=sub, n_subpixel=8, shadowing=True,
                water=WaterBody(case=2), sun_glint=(45.0, 195.0, 50.0),
                rng=np.random.default_rng(3)))
        frames.append(S[..., 0])
    frames = np.stack(frames, axis=-1)
    vmax = np.nanpercentile(frames, 99.0)

    # GIF
    from matplotlib import animation
    fig_g, ax_g = plt.subplots(figsize=(4.6, 4.6), dpi=100)
    im = ax_g.imshow(frames[..., 0], cmap="gray", vmin=0, vmax=vmax)
    ax_g.set_xticks([])
    ax_g.set_yticks([])

    def update(i):
        im.set_data(frames[..., i])
        ax_g.set_title(f"Stokes I, t = {movie_times[i]:.2f} s", fontsize=9)
        return [im]

    ani = animation.FuncAnimation(fig_g, update, frames=frames.shape[-1])
    gif = OUT / "demo_full_pipeline.gif"
    ani.save(gif, writer=animation.PillowWriter(fps=6))
    print(f"saved -> {gif}")
    plt.close(fig_g)

    # summary panel
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), dpi=110)
    ax = axes[0]
    # streamed surf carries frame 0 only (2-D eta)
    im = ax.imshow(surf.eta * 100, cmap="RdBu_r",
                   extent=[0, L * 100, 0, L * 100])
    ax.set_title(f"eta [cm], measured-spectrum run {RUN}")
    ax.set_xlabel("x [cm]")
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1]
    im = ax.imshow(frames[..., 0], cmap="gray", vmin=0, vmax=vmax)
    ax.set_title("rendered I (glint + clear sky + water)")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[2]
    ax.semilogy(nu_s, Q_s / np.nanmax(Q_s), "r-", label="synthetic")
    ax.semilogy(nu_m, Q_m / np.nanmax(Q_m), "k--", label="measured (same run)")
    ax.set_xlabel("nu = 1/c [s/m]")
    ax.set_ylabel("slope Q(nu), peak-normalized")
    ax.set_title("inverse-phase-speed validation")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(f"Full pipeline from ASIT run {RUN} "
                 f"(U10 = {U10:.1f} m/s): measured spectrum + measured "
                 f"bound fraction + FM98/MTF/advection", fontsize=11)
    plt.tight_layout()
    out = OUT / "demo_full_pipeline.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
