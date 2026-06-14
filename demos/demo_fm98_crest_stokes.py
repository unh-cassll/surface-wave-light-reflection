"""
Polarimetric signature of parasitic capillaries on a steep wave crest.

The exact FM98 steady solution (lambda = 5 cm) is laid out as a
long-crested surface at ~0.1 mm resolution and run through the polarized
facet renderer under a Rayleigh sky: Stokes I, Q/I, U/I, DoLP, and AoP
across the crest, for a steep carrier (ak = 0.28, pronounced parasitic
train) against a gentle one (ak = 0.15, none).  The capillary train
imprints a distinctive polarization banding on the forward face --
the observable the polarimetric slope-sensing technique exploits.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import CameraGeometry, make_clear_sky, render_facet_stokes
from seapol.fm98 import solve_fm98_continuation
from seapol.polarization import stokes_aop, stokes_dolp

OUT = Path(__file__).parent / "output"

LAM = 0.05
N_LAM = 4                  # wavelengths along x
NX, NY = 2048, 128         # dx ~ 0.098 mm
CAM = CameraGeometry(incidence_deg=53.0, azimuth_deg=0.0, height_m=50.0)
SKY = make_clear_sky(sun_zenith_deg=45.0, sun_azimuth_deg=120.0,
                     turbidity=0.1)


def fm98_surface(ak: float):
    """Long-crested FM98 surface eta(y, x) over N_LAM wavelengths
    (well-forced branch: published forward-face train structure)."""
    from seapol.fm98 import default_forcing
    sol = solve_fm98_continuation(LAM, ak,
                                  p_target=default_forcing(2 * np.pi / LAM,
                                                           ak),
                                  M=64, n_steps=24)
    x_p = sol.X - sol.X.min()
    order = np.argsort(x_p)
    x_s, y_s = x_p[order], sol.Y[order]
    # periodic closure for interpolation
    x_s = np.concatenate([x_s, [x_s[0] + LAM]])
    y_s = np.concatenate([y_s, [y_s[0]]])

    L = N_LAM * LAM
    dx = L / NX
    xg = np.arange(NX) * dx
    eta_1d = np.interp(np.mod(xg, LAM), x_s, y_s)
    kx = 2 * np.pi * np.fft.fftfreq(NX, d=dx)
    sx_1d = np.fft.ifft(1j * kx * np.fft.fft(eta_1d)).real
    eta = np.tile(eta_1d, (NY, 1))
    sx = np.tile(sx_1d, (NY, 1))
    sy = np.zeros_like(sx)
    return eta, sx, sy, dx, xg, sol


def render(ak: float):
    eta, sx, sy, dx, xg, sol = fm98_surface(ak)
    S = render_facet_stokes(eta, dx, camera=CAM, sky=SKY,
                            slope_x=sx, slope_y=sy)
    row = NY // 2
    return dict(x=xg, eta=eta[row], sx=sx[row], S=S, S_row=S[row],
                dx=dx, res=sol.residual_norm)


def main():
    OUT.mkdir(exist_ok=True)
    steep = render(0.28)
    gentle = render(0.15)
    print(f"FM98 residuals: steep {steep['res']:.1e}, "
          f"gentle {gentle['res']:.1e}")

    x_cm = steep["x"] * 100
    win = (x_cm >= 5.0) & (x_cm <= 10.5)   # one wavelength + crest margin

    fig = plt.figure(figsize=(12.5, 12), dpi=110)
    gs = fig.add_gridspec(6, 1, height_ratios=[1.1, 0.9, 0.9, 0.9, 0.9,
                                               1.0])

    ax = fig.add_subplot(gs[0])
    ax.plot(x_cm[win], steep["eta"][win] * 1000, "k-", lw=1.2,
            label="ak = 0.28 (parasitic train)")
    ax.plot(x_cm[win], gentle["eta"][win] * 1000, "c--", lw=1.0,
            label="ak = 0.15")
    ax.set_ylabel("elevation [mm]")
    ax.set_title("FM98 steady profiles, lambda = 5 cm (wave moves +x; "
                 "camera at 53 deg incidence)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    panels = [
        ("I", lambda r: r["S_row"][:, 0], "radiance"),
        ("Q/I", lambda r: r["S_row"][:, 1] / r["S_row"][:, 0], None),
        ("U/I", lambda r: r["S_row"][:, 2] / r["S_row"][:, 0], None),
        ("DoLP", lambda r: stokes_dolp(r["S_row"]), None),
    ]
    for i, (name, fn, ylab) in enumerate(panels):
        ax = fig.add_subplot(gs[1 + i])
        ax.plot(x_cm[win], fn(steep)[win], "k-", lw=1.0)
        ax.plot(x_cm[win], fn(gentle)[win], "c--", lw=0.9)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.set_title("Stokes transects across the crest "
                         "(black: steep, dashed cyan: gentle)")

    ax = fig.add_subplot(gs[5])
    # AoP is defined mod 180 deg; map to [0, 180) to avoid +/-90 wraps
    aop = np.mod(np.rad2deg(stokes_aop(steep["S_row"])), 180.0)
    aop_g = np.mod(np.rad2deg(stokes_aop(gentle["S_row"])), 180.0)
    ax.plot(x_cm[win], aop[win], "k-", lw=1.0)
    ax.plot(x_cm[win], aop_g[win], "c--", lw=0.9)
    ax.set_ylabel("AoP [deg, mod 180]")
    ax.set_xlabel("x [cm]")
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = OUT / "demo_fm98_crest_stokes.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")

    # 2D strips: DoLP and AoP around one steep crest
    S = steep["S"]
    dolp2 = stokes_dolp(S)
    aop2 = np.rad2deg(stokes_aop(S))
    j = win.nonzero()[0]
    fig2, axes = plt.subplots(2, 1, figsize=(11, 4.4), dpi=110)
    ext = [x_cm[j[0]], x_cm[j[-1]], 0, NY * steep["dx"] * 100]
    im = axes[0].imshow(dolp2[:, j], cmap="viridis", origin="lower",
                        extent=ext, aspect="auto")
    axes[0].set_title("DoLP strip, ak = 0.28: capillary banding on the "
                      "forward faces")
    plt.colorbar(im, ax=axes[0], fraction=0.03)
    im = axes[1].imshow(aop2[:, j], cmap="twilight", origin="lower",
                        extent=ext, aspect="auto", vmin=-90, vmax=90)
    axes[1].set_title("AoP strip [deg]")
    axes[1].set_xlabel("x [cm]")
    plt.colorbar(im, ax=axes[1], fraction=0.03)
    plt.tight_layout()
    out2 = OUT / "demo_fm98_crest_stokes_strips.png"
    plt.savefig(out2, bbox_inches="tight")
    print(f"saved -> {out2}")


if __name__ == "__main__":
    main()
