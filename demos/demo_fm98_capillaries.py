"""
Phase-locked parasitic capillaries, two ways:

  1. Direct FM98 steady profiles at a 5 cm carrier for increasing
     steepness.  The ripple train (harmonics m ~ 10, lambda ~ 5 mm, set
     by the c(k_ripple) = c(carrier) resonance) is ~0.1 mm in elevation
     but ~half the carrier amplitude in slope, so the profiles are shown
     in both elevation and slope.
  2. A capillary-resolving hybrid surface (dx ~ 0.5 mm) with the deep
     (M_keep = 28) gauge-fixed FM98 table: bound-harmonic curvature
     reveals ripple packets phase-locked to the steep wave-group crests.

The deep table build (one-time, ~30 min) is cached to demos/output/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import generate_hybrid_surface
from seapol.fm98 import FM98Table, build_fm98_table, solve_fm98_continuation
from seapol.hybrid import default_carrier_band

OUT = Path(__file__).parent / "output"
TABLE_CACHE = OUT / "fm98_table_deep.npz"

L, N, U10 = 2.0, 4096, 9.0
WAVELENGTH = 0.05
MTF, MTF_PHASE_DEG = 6.5, 30.0      # long-wave binding (hydrodynamic MTF)
K_SPLIT = 2 * np.pi / 0.5           # waves longer than 50 cm are "long"


def deep_table() -> FM98Table:
    """Deep gauge-fixed FM98 table on the well-forced branch
    (default_forcing rule): ak to 0.40 (MF15 figure-2 steepness),
    M_keep = 28 so carriers down to ~70 rad/m reach their capillary
    resonance m* = (k_m / k)^2."""
    if TABLE_CACHE.exists():
        print(f"loading cached table {TABLE_CACHE}")
        return FM98Table.load(TABLE_CACHE)
    print("building deep FM98 table (one-time, ~30 min)...")
    k_grid = np.geomspace(60.0, 380.0, 8)
    ak_grid = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                        0.34, 0.375, 0.40])
    tbl = build_fm98_table(k_grid, ak_grid, M_keep=28, M_solve=64,
                           n_steps=14, verbose=True)
    OUT.mkdir(exist_ok=True)
    # tmp name must end in .npz or np.savez appends it and os.replace
    # then targets a nonexistent path
    tmp = TABLE_CACHE.with_suffix(".tmp.npz")
    tbl.save(tmp)
    import os
    os.replace(tmp, TABLE_CACHE)
    return tbl


def fm98_profile_panels(ax_eta, ax_slope):
    """Direct FM98 solves in elevation and slope; slope makes the
    parasitic train obvious.  Well-forced branch (default_forcing):
    the train rides the forward face and decays into the trough, as in
    Melville & Fedorov (2015) figure 2; steepness up to their range."""
    from seapol.fm98 import default_forcing
    k_c = 2 * np.pi / WAVELENGTH
    # targets stay below the lam = 5 cm branch fold (~0.32): pushing
    # past it returns lumpy multi-crested states instead of the clean
    # carrier + train
    for ak, color in zip([0.22, 0.28, 0.32],
                         plt.cm.viridis([0.15, 0.5, 0.85])):
        sol = solve_fm98_continuation(WAVELENGTH, ak,
                                      p_target=default_forcing(k_c, ak),
                                      M=64, n_steps=28)
        order = np.argsort(sol.X)
        X1, Y1 = sol.X[order], sol.Y[order]
        # phase-align: solution phase origins are arbitrary, so put
        # every crest at lambda/4 (overlapping unaligned curves read
        # as spurious double crests)
        X1 = (X1 - X1[np.argmax(Y1)] + 0.25 * WAVELENGTH) % WAVELENGTH
        order = np.argsort(X1)
        X1, Y1 = X1[order], Y1[order]
        X = np.concatenate([X1, X1 + WAVELENGTH])
        Y = np.concatenate([Y1, Y1])
        slope = np.gradient(Y, X)
        ak_real = 0.5 * (sol.Y.max() - sol.Y.min()) * k_c
        ax_eta.plot(X * 100, Y * 1000, color=color, lw=1.2,
                    label=f"ak = {ak_real:.2f} (target {ak:.2f})")
        ax_slope.plot(X * 100, slope, color=color, lw=1.0)
    ax_eta.set_ylabel("elevation [mm]")
    ax_eta.set_title(f"Direct FM98 Class-1 profiles, lambda = "
                     f"{WAVELENGTH * 100:.0f} cm (wave moves +x)")
    ax_eta.legend(fontsize=8, loc="upper right")
    ax_eta.grid(alpha=0.3)
    ax_slope.set_xlabel("x [cm]")
    ax_slope.set_ylabel("slope")
    ax_slope.set_title("surface slope: parasitic ripple train on the "
                       "forward face")
    ax_slope.grid(alpha=0.3)


def _bandpass_field(eta, dx, k_lo, k_hi):
    Nn = eta.shape[0]
    kx = 2.0 * np.pi * np.fft.fftfreq(Nn, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    Z = np.fft.fft2(eta)
    return np.fft.ifft2(np.where((K >= k_lo) & (K <= k_hi), Z, 0.0)).real


def main():
    OUT.mkdir(exist_ok=True)
    table = deep_table()
    print(f"table: {int(table.converged.sum())}/{table.converged.size} "
          f"points converged")

    print(f"hybrid surface: L = {L} m, N = {N}, dx = {L / N * 1e3:.3f} mm, "
          f"U10 = {U10} m/s, long-wave MTF = {MTF} ...")
    surf = generate_hybrid_surface(L=L, N=N, U10=U10, table=table,
                                   long_wave_mtf=MTF,
                                   mtf_phase_deg=MTF_PHASE_DEG,
                                   k_split=K_SPLIT,
                                   rng=np.random.default_rng(3),
                                   verbose=True)
    dx = surf.info["dx"]
    k_lo, k_hi = surf.info["carrier_band"]

    eh = surf.eta_high - surf.eta_high.mean()
    ex = np.gradient(surf.eta_high, dx, axis=1)
    print(f"skew(eta_high) = {np.mean(eh**3) / eh.std()**3:+.3f} "
          f"(FM98: > 0), skew(d eta_high/dx) = "
          f"{np.mean((ex - ex.mean())**3) / ex.std()**3:+.3f} (FM98: < 0)")

    # Carrier-band field (wave groups) and bound-harmonic curvature
    eta_band = _bandpass_field(surf.eta, dx, k_lo, k_hi)
    curv_high = (np.gradient(np.gradient(surf.eta_high, dx, axis=1), dx,
                             axis=1)
                 + np.gradient(np.gradient(surf.eta_high, dx, axis=0), dx,
                               axis=0))

    # Window centered on the strongest ripple packet
    from scipy.ndimage import uniform_filter
    win = int(0.30 / dx)
    energy = uniform_filter(np.abs(surf.eta_high), size=win // 4)
    pad = win // 2 + 1
    e_in = energy[pad:-pad, pad:-pad]
    iy, ix = np.unravel_index(np.argmax(e_in), e_in.shape)
    iy += pad
    ix += pad
    sl_y = slice(iy - win // 2, iy + win // 2)
    sl_x = slice(ix - win // 2, ix + win // 2)
    extent = [0, win * dx * 100, 0, win * dx * 100]

    fig = plt.figure(figsize=(13, 12), dpi=110)
    gs = fig.add_gridspec(4, 2, height_ratios=[0.9, 0.9, 1.5, 0.9])

    ax_eta = fig.add_subplot(gs[0, :])
    ax_slope = fig.add_subplot(gs[1, :], sharex=ax_eta)
    fm98_profile_panels(ax_eta, ax_slope)

    ax1 = fig.add_subplot(gs[2, 0])
    vband = 3 * eta_band[sl_y, sl_x].std()
    im1 = ax1.imshow(eta_band[sl_y, sl_x] * 1000, cmap="RdBu_r",
                     origin="lower", extent=extent, vmin=-vband * 1000,
                     vmax=vband * 1000)
    ax1.set_title("carrier-band elevation (2-8 cm waves) [mm]")
    ax1.set_xlabel("x [cm] (wind ->)")
    ax1.set_ylabel("y [cm]")
    plt.colorbar(im1, ax=ax1, fraction=0.046)

    # Tighter sub-window so the ~5 mm ripple stripes survive rasterization
    win2 = 2 * (int(0.12 / dx) // 2)
    sl2_y = slice(iy - win2 // 2, iy + win2 // 2)
    sl2_x = slice(ix - win2 // 2, ix + win2 // 2)
    extent2 = [0, win2 * dx * 100, 0, win2 * dx * 100]

    ax2 = fig.add_subplot(gs[2, 1])
    cz = curv_high[sl2_y, sl2_x]
    vmax = np.percentile(np.abs(cz), 99.5)
    im2 = ax2.imshow(cz, cmap="gray", origin="lower", extent=extent2,
                     vmin=-vmax, vmax=vmax)
    ax2.contour(np.linspace(*extent2[:2], win2),
                np.linspace(*extent2[2:], win2),
                eta_band[sl2_y, sl2_x] * 1000, levels=4,
                colors="r", linewidths=0.6, alpha=0.7)
    ax2.set_title("bound-harmonic curvature, 12 cm zoom:\n"
                  "ripple packets locked to carrier crests (red)")
    ax2.set_xlabel("x [cm] (wind ->)")
    plt.colorbar(im2, ax=ax2, fraction=0.046)
    # mark the curvature sub-window on the carrier-band panel
    x0 = (sl2_x.start - sl_x.start) * dx * 100
    y0 = (sl2_y.start - sl_y.start) * dx * 100
    ax1.add_patch(plt.Rectangle((x0, y0), win2 * dx * 100, win2 * dx * 100,
                                fill=False, edgecolor="k", lw=1.0))

    # Detrended slice through the packet: carrier band + ripples
    ax3 = fig.add_subplot(gs[3, 0])
    row = iy
    xs = (np.arange(sl_x.start, sl_x.stop) - sl_x.start) * dx * 100
    ax3.plot(xs, eta_band[row, sl_x] * 1000, "c-", lw=1.0,
             label="carrier band (2-8 cm)")
    ax3.plot(xs, (eta_band[row, sl_x] + surf.eta_high[row, sl_x]) * 1000,
             "k-", lw=0.7, alpha=0.8, label="carrier + bound harmonics")
    ax3.plot(xs, surf.eta_high[row, sl_x] * 1000 * 5, "r-", lw=0.7,
             label="bound harmonics x 5")
    ax3.set_xlabel("x [cm] (wind ->)")
    ax3.set_ylabel("elevation [mm]")
    ax3.set_title("slice through the strongest ripple packet (long waves "
                  "removed)")
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)

    # Phase-resolved binding: bound-harmonic energy vs long-wave phase
    from seapol.hybrid import _analytic_signal_along_wind
    kxf = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    Kf = np.hypot(kxf[None, :], kxf[:, None])
    Zf = np.fft.fft2(surf.eta_lin)
    eta_L = np.fft.ifft2(np.where((Kf > 0) & (Kf <= K_SPLIT), Zf, 0)).real
    phi_L = np.angle(_analytic_signal_along_wind(eta_L, kxf, kxf, 0.0))
    energy = uniform_filter(surf.eta_high**2, size=int(0.02 / dx))
    bins = np.linspace(-np.pi, np.pi, 25)
    idx = np.digitize(phi_L.ravel(), bins) - 1
    prof = np.array([energy.ravel()[idx == i].mean() for i in range(24)])
    centers = 0.5 * (bins[:-1] + bins[1:])

    ax4 = fig.add_subplot(gs[3, 1])
    ax4.plot(np.rad2deg(centers), prof / energy.mean(), "ko-", ms=3)
    ax4.axvline(0, color="b", lw=0.8, alpha=0.6)
    ax4.axvline(MTF_PHASE_DEG, color="r", ls="--", lw=0.8,
                label=f"MTF phase ({MTF_PHASE_DEG:.0f} deg)")
    ax4.text(2, ax4.get_ylim()[0] + 0.05, "crest", color="b", fontsize=7,
             rotation=90)
    ax4.set_xlabel("long-wave phase [deg]  (+ = forward face)")
    ax4.set_ylabel("capillary energy / mean")
    ax4.set_title("bound capillaries locked to the long wave")
    ax4.legend(fontsize=8)
    ax4.grid(alpha=0.3)

    fig.suptitle(
        f"Fedorov-Melville parasitic capillaries: direct solutions and "
        f"phase-locked hybrid surface (U10 = {U10} m/s, dx = "
        f"{dx * 1e3:.2f} mm)", fontsize=12)
    plt.tight_layout()
    out = OUT / "demo_fm98_capillaries.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
