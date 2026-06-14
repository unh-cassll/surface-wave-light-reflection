"""
Closing the loop: polarimetric reconstruction of wave slope and height
statistics from rendered imagery.

A synthetic surface is rendered under an overcast (unpolarized) sky;
per facet, DoLP gives the Fresnel incidence angle and AoP the plane of
incidence, from which the facet normal -- and hence the slope field --
is reconstructed (seapol.inversion).  Spectral integration of the
recovered slopes returns the height field.  Reconstructed slope PDFs,
MSS, height PDF, Hs, and the omnidirectional spectrum are compared
against the known truth: the framework demonstrates the full
measure-from-light workflow used by polarimetric slope-sensing
instruments, with the ground truth that field data never has.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (CameraGeometry, generate_sea_surface,
                    height_from_slopes, make_overcast_sky,
                    render_facet_stokes, slopes_from_stokes)
from seapol.polarization import normalize, stokes_aop, stokes_dolp

OUT = Path(__file__).parent / "output"

L, N, U10 = 32.0, 256, 6.0
INC_DEG = 30.0          # camera incidence: keeps facets sub-Brewster


def omni_spectrum(eta, dx, n_bins=48):
    n = eta.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    P = np.abs(np.fft.fft2(eta - eta.mean()))**2
    edges = np.geomspace(2 * np.pi / (n * dx), kx.max(), n_bins + 1)
    idx = np.clip(np.digitize(K.ravel(), edges) - 1, 0, n_bins - 1)
    pow_k = np.bincount(idx, weights=P.ravel(), minlength=n_bins)
    cnt = np.maximum(np.bincount(idx, minlength=n_bins), 1)
    kc = np.sqrt(edges[:-1] * edges[1:])
    return kc, pow_k / cnt


def main():
    rng = np.random.default_rng(0)
    surf = generate_sea_surface(L=L, N=N, U10=U10, rng=rng)
    dx = surf.info["dx"]

    print("rendering overcast-sky polarized imagery ...")
    cam = CameraGeometry(incidence_deg=INC_DEG, azimuth_deg=0.0,
                         height_m=100.0)
    S = render_facet_stokes(surf.eta, dx, camera=cam,
                            sky=make_overcast_sky(1.0),
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=None, shadowing=False)

    # inversion: view directions from the flat-surface approximation
    xg = np.arange(N) * dx
    X, Y = np.meshgrid(xg, xg)
    P3 = np.stack([X, Y, np.zeros_like(X)], axis=-1)
    cpos = cam.position(((N - 1) * dx / 2.0, (N - 1) * dx / 2.0))
    d_out = normalize(cpos[None, None, :] - P3)
    sx_r, sy_r, valid = slopes_from_stokes(S, d_out)
    eta_r = height_from_slopes(sx_r, sy_r, dx)
    eta_t = surf.eta - surf.eta.mean()
    print(f"valid pixels: {valid.mean()*100:.1f}%")

    mss_t = surf.slope_x.var() + surf.slope_y.var()
    mss_r = np.nanvar(sx_r) + np.nanvar(sy_r)
    c_eta = np.corrcoef(eta_r.ravel(), eta_t.ravel())[0, 1]
    print(f"MSS: truth {mss_t:.5f} -> reconstructed {mss_r:.5f} "
          f"({mss_r/mss_t*100:.1f}%)")
    print(f"Hs:  truth {4*eta_t.std():.3f} m -> reconstructed "
          f"{4*eta_r.std():.3f} m (corr {c_eta:.4f})")

    fig = plt.figure(figsize=(15, 9), dpi=110)

    # measurements
    ax = fig.add_subplot(2, 4, 1)
    im = ax.imshow(stokes_dolp(S), cmap="viridis")
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("measured DoLP")
    ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(2, 4, 2)
    im = ax.imshow(np.rad2deg(stokes_aop(S)), cmap="twilight",
                   vmin=-90, vmax=90)
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("measured AoP [deg]")
    ax.set_xticks([]); ax.set_yticks([])

    # slope field: truth vs reconstruction
    v = np.nanpercentile(np.abs(surf.slope_x), 99)
    for i, (f, t) in enumerate(((surf.slope_x, "slope x (truth)"),
                                (sx_r, "slope x (reconstructed)"))):
        ax = fig.add_subplot(2, 4, 3 + i)
        ax.imshow(f, cmap="RdBu_r", vmin=-v, vmax=v)
        ax.set_title(t)
        ax.set_xticks([]); ax.set_yticks([])

    # slope PDFs
    ax = fig.add_subplot(2, 4, 5)
    bins = np.linspace(-3.5 * np.sqrt(mss_t / 2), 3.5 * np.sqrt(mss_t / 2),
                       80)
    for tru, rec, lbl in ((surf.slope_x, sx_r, "along"),
                          (surf.slope_y, sy_r, "cross")):
        h_t, _ = np.histogram(tru[valid], bins=bins, density=True)
        h_r, _ = np.histogram(rec[valid & np.isfinite(rec)], bins=bins,
                              density=True)
        cb = 0.5 * (bins[:-1] + bins[1:])
        ax.semilogy(cb, h_t, "-", label=f"{lbl} truth")
        ax.semilogy(cb, h_r, ".", ms=3, label=f"{lbl} recon")
    ax.set_xlabel("slope")
    ax.set_ylabel("PDF")
    ax.set_title("slope PDFs")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # height field + PDF
    vh = np.nanpercentile(np.abs(eta_t), 99.5)
    ax = fig.add_subplot(2, 4, 6)
    ax.imshow(eta_r, cmap="RdBu_r", vmin=-vh, vmax=vh)
    ax.set_title(f"reconstructed eta (corr {c_eta:.3f})")
    ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(2, 4, 7)
    hb = np.linspace(-3.5 * eta_t.std(), 3.5 * eta_t.std(), 60)
    for f, st, lbl in ((eta_t, "-", "truth"), (eta_r, "--", "recon")):
        h, _ = np.histogram(f.ravel(), bins=hb, density=True)
        ax.plot(0.5 * (hb[:-1] + hb[1:]), h, st, label=lbl)
    ax.set_xlabel("eta [m]")
    ax.set_title(f"height PDF: Hs {4*eta_t.std():.3f} -> "
                 f"{4*eta_r.std():.3f} m")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # omnidirectional elevation spectrum
    ax = fig.add_subplot(2, 4, 8)
    for f, st, lbl in ((eta_t, "-", "truth"), (eta_r, "--", "recon")):
        kc, pw = omni_spectrum(f, dx)
        ax.loglog(kc, pw, st, label=lbl)
    ax.set_xlabel("k [rad/m]")
    ax.set_ylabel("ring power")
    ax.set_title("omnidirectional spectrum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_polarimetric_reconstruction.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
