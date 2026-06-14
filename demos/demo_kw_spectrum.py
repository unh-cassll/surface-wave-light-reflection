"""
Wavenumber-frequency diagnostics of the evolving hybrid surface against
the ASIT 2019 observations.

A time-evolving surface stack is reduced to the along-wind slope (kx, f)
spectrum and the inverse-phase-speed spectrum Q(nu), nu = kx/(2 pi f).
ASIT measurements show short-wave slope variance concentrated at the
dominant-wave phase speed (only a few percent rides the free dispersion
shell), so the synthesis is run with the free/bound spectral partition
(bound_fraction): the bound share of high-k power advects rigidly at the
dominant speed while bin powers (Hs, MSS) are unchanged.  The measured
Q(nu) from the ASIT stats file is overlaid when present in
demos/output/.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import generate_hybrid_surface, generate_sea_surface
from seapol.fm98 import FM98Table
from seapol.spectrum import angular_frequency

OUT = Path(__file__).parent / "output"
TABLE_CACHE = OUT / "fm98_table_deep.npz"

# Small fine tile sampled fast: bound harmonics of the 3-8 cm carriers
# (m = 2, 3 -> 12-26 Hz) stay below the frequency Nyquist instead of
# aliasing, which is what hides them at coarser dx / slower fs.
L, N = 2.0, 512                    # dx = 3.9 mm, k_nyq = 804 rad/m
# L = 2 m fits ~6 long-wave modes below k_split: orbital-advection
# Doppler sidebands then merge into the smooth broadened ridge seen in
# the measured spectra (a 1 m tile has only ~3 modes, and the same
# sidebands appear as a sharp discrete comb fanning across (k, f))
U10 = 7.0
FS, NT = 60.0, 2048                # 34 s record at 60 Hz (streamed)
N_ROWS = 48                        # transect block retained per frame
K_SPLIT = 2 * np.pi / 0.30
MTF = 6.5
BETA = 0.9                         # bound fraction at high k (ASIT: most
C_BOUND = 2.0                      # short-wave slope rides near ~2 m/s)
# bound waves ride a spectrum of carriers: dominant + intermediate scales
BOUND_SPEEDS = ((2.0, 1.0, 0.55), (0.55, 0.30, 0.15))


from seapol.diagnostics import kf_slope_spectrum, q_nu, slow_nu_fraction


def main():
    OUT.mkdir(exist_ok=True)
    table = FM98Table.load(TABLE_CACHE)
    times = np.arange(NT) / FS
    rng_seed = 7

    def transect_buffer():
        buf = np.empty((N_ROWS, N, NT))

        def cb(it, t, e, gx, gy):
            buf[:, :, it] = e[:N_ROWS]
        return buf, cb

    print(f"free-only record (streamed): L={L} m, N={N}, {NT} frames "
          f"at {FS} Hz ...")
    buf_lin, cb = transect_buffer()
    lin = generate_sea_surface(L, N, U10, times=times, compute_slopes=False,
                               frame_callback=cb,
                               rng=np.random.default_rng(rng_seed))
    dx = lin.info["dx"]

    # bound fraction: fitted from the measured ASIT k-f reduction when
    # available, else the scalar default
    asit_npz = OUT / "asit_kf_reduced.npz"
    if asit_npz.exists():
        from seapol import bound_fraction_from_kf_reduction
        beta = bound_fraction_from_kf_reduction(asit_npz)
        beta_tag = "measured beta(k)"
    else:
        beta = BETA
        beta_tag = f"beta = {BETA}"
    print(f"hybrid stack ({beta_tag} at c = {C_BOUND} m/s + orbital "
          f"advection + FM98 + binding) ...")
    buf_hyb, cb = transect_buffer()
    # f_nyq: harmonics with m f(k_loc) above the record Nyquist would
    # alias into sharp folded ridges absent from the measured spectra
    hyb = generate_hybrid_surface(L=L, N=N, U10=U10, times=times,
                                  table=table, long_wave_mtf=MTF,
                                  k_split=K_SPLIT, orbital_advection=True,
                                  bound_fraction=beta,
                                  bound_speed=BOUND_SPEEDS,
                                  f_nyq=FS / 2.0,
                                  frame_callback=cb,
                                  rng=np.random.default_rng(rng_seed))
    print(f"  m_max = {hyb.info['m_max']}, carrier band = "
          f"{hyb.info['carrier_band'][0]:.0f}-"
          f"{hyb.info['carrier_band'][1]:.0f} rad/m")

    kx, f, P_lin = kf_slope_spectrum(buf_lin, dx, FS, n_seg=3)
    _, _, P_hyb = kf_slope_spectrum(buf_hyb, dx, FS, n_seg=3)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), dpi=110)
    k_fine = np.linspace(kx[0], kx[-1], 400)
    f_disp = angular_frequency(k_fine) / (2 * np.pi)

    for ax, P, title in [(axes[0], P_lin, "free-only synthesis"),
                         (axes[1], P_hyb,
                          f"hybrid ({beta_tag} at {C_BOUND} m/s "
                          f"+ advection)")]:
        lp = np.log10(np.maximum(P, P[P > 0].min()))
        im = ax.pcolormesh(kx, f, lp.T, cmap="magma",
                           vmin=lp.max() - 6, vmax=lp.max())
        ax.plot(k_fine, f_disp, "c--", lw=1.0, label="linear dispersion")
        for c_b, ls in [(C_BOUND, ":"), (0.29, "-.")]:
            ax.plot(k_fine, c_b * k_fine / (2 * np.pi), "w", ls=ls, lw=0.8,
                    label=f"c = {c_b} m/s")
        ax.set_xlim(0, kx[-1])
        ax.set_ylim(0, FS / 2)
        ax.set_xlabel("kx [rad/m]")
        ax.set_title(title)
        if ax is axes[0]:
            ax.set_ylabel("f [Hz]")
            ax.legend(fontsize=7, loc="upper left")
        plt.colorbar(im, ax=ax, fraction=0.046, label="log10 P")

    nu, Q_lin = q_nu(kx, f, P_lin)
    _, Q_hyb = q_nu(kx, f, P_hyb)
    ax = axes[2]
    ax.semilogy(nu, Q_lin / Q_lin.max(), "b-", label="free-only")
    ax.semilogy(nu, Q_hyb / Q_hyb.max(), "r-", label="hybrid (bound)")
    asit = OUT / "asit_Q_nu_U7.npy"
    if asit.exists():
        nu_o, Q_o = np.load(asit)
        ax.semilogy(nu_o, Q_o / np.nanmax(Q_o), "k--", lw=1.2,
                    label="ASIT measured (U10 ~ 7)")
    ax.axvline(1.0 / C_BOUND, color="r", ls=":", lw=0.8)
    ax.text(1.0 / C_BOUND + 0.05, 3e-3, "bound c", rotation=90, fontsize=7)
    ax.axvline(1.0 / 0.23, color="b", ls=":", lw=0.8)
    ax.text(1.0 / 0.23 - 0.35, 3e-3, "slowest free", rotation=90,
            fontsize=7)
    ax.set_xlabel("nu = 1/c [s/m]")
    ax.set_ylabel("slope Q(nu) (peak-normalized)")
    ax.set_title("inverse-phase-speed spectrum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(f"slope k-f diagnostics, U10 = {U10} m/s: bound waves "
                 f"ride c = {C_BOUND} m/s, off the dispersion shell")
    plt.tight_layout()
    out = OUT / "demo_kw_spectrum.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")
    # free-fraction metric: slope-Q at slow free-wave phase speeds.
    # References: 0.078 = noise-subtracted cube restricted to this tile's
    # k support (6.3-800 rad/m); 0.032 = the processed ASIT stats over
    # the full support including the dominant band.
    for tag, Q in [("free-only", Q_lin), ("hybrid", Q_hyb)]:
        frac = slow_nu_fraction(nu, Q)
        print(f"{tag}: fraction of slope Q at nu > 2 s/m = {frac:.3f} "
              f"(measured: 0.078 tile-support / 0.032 full)")


if __name__ == "__main__":
    main()
