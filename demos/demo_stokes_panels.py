"""
Stokes panels of the full Elfouhaily + FM98 capillary blend across wind
speed.

For each U10 in 3 .. 15 m/s (step 2) a 1024 x 1024 surface over a
1 m x 1 m patch is synthesized with generate_hybrid_surface (linear
Elfouhaily field plus the phase-locked FM98 parasitic-capillary
augmentation), rendered to a per-facet Stokes image under a clear
polarized sky, and shown as a 1 x 3 grayscale panel:

    S0        reflected radiance (no colorbar; color limits set to a
              robust percentile range for visible dynamic range)
    S1 / S0   = Q/I, normalized polarization, color -1 .. 1
    S2 / S0   = U/I, normalized polarization, color -1 .. 1

Runs on the GPU when torch + CUDA are available (backend dispatch);
falls back to numpy otherwise.  The FM98 coefficient table is loaded
from the cached deep table when present (its k 60-380 rad/m grid covers
the 2-8 cm carrier band of this grid) and built + cached otherwise.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import backend
from seapol.fm98 import FM98Table
from seapol.hybrid import default_fm98_table, generate_hybrid_surface
from seapol.render import (CameraGeometry, make_clear_sky,
                           render_facet_stokes)

OUT = Path(__file__).parent / "output"
PANELS = OUT / "stokes_panels"
DEEP_CACHE = OUT / "fm98_table_deep.npz"     # reused read-only if it covers the band
OWN_CACHE = OUT / "fm98_table_panels.npz"    # built here if the deep table is absent

L, N = 1.0, 1024                       # 1 m x 1 m, 1024 px (dx ~ 0.98 mm)
WINDS = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
SEED = 20240                            # same realization across winds
S0_PCTL = (1.0, 99.0)                   # robust S0 color limits

# fixed illumination / view, so wind is the only variable
SUN_ZEN, SUN_AZ = 50.0, 90.0
CAMERA = CameraGeometry(incidence_deg=40.0, azimuth_deg=0.0, height_m=50.0)


def get_table(k_nyq):
    band = (78.5, 314.2)                      # 2-8 cm carrier band of this grid
    for cache in (DEEP_CACHE, OWN_CACHE):
        if cache.exists():
            t = FM98Table.load(cache)
            if t.k_grid[0] <= band[0] and t.k_grid[-1] >= band[1]:
                print(f"using cached FM98 table {cache.name} "
                      f"(k {t.k_grid[0]:.0f}-{t.k_grid[-1]:.0f}, "
                      f"M_keep={t.M_keep})")
                return t
    print("building FM98 table (one-time, cached) ...")
    t = default_fm98_table(k_nyq, m_keep=20, verbose=True)
    t.save(OWN_CACHE)                          # never clobber the deep table
    return t


def render_stokes(U10, table, xp, bk, dev):
    surf = generate_hybrid_surface(L, N, U10, table=table,
                                   rng=np.random.default_rng(SEED),
                                   backend=bk, device=dev)
    sky = make_clear_sky(SUN_ZEN, SUN_AZ, I_sky=1.0, turbidity=0.1)
    S = render_facet_stokes(surf.eta, surf.info["dx"], CAMERA, sky=sky,
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=None,
                            rng=backend.default_rng(0, xp))
    return backend.to_numpy(S), surf.info


def panel_figure(U10, S):
    """1 x 3 grayscale panel; returns (fig, s0_lims)."""
    I = S[..., 0]
    eps = 1e-6 * np.nanmax(np.abs(I))
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(np.abs(I) > eps, S[..., 1] / I, np.nan)
        u = np.where(np.abs(I) > eps, S[..., 2] / I, np.nan)
    q = np.clip(q, -1.0, 1.0)
    u = np.clip(u, -1.0, 1.0)
    vmin, vmax = np.nanpercentile(I, S0_PCTL)

    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.0), dpi=220)
    ext = [0.0, 100.0, 0.0, 100.0]            # cm
    titles = ["$S_0$ (reflected radiance)", "$S_1/S_0$", "$S_2/S_0$"]
    imgs = [(I, vmin, vmax), (q, -1.0, 1.0), (u, -1.0, 1.0)]
    norm_im = None
    for ax, (data, lo, hi), title in zip(axes, imgs, titles):
        im = ax.imshow(data, cmap="gray", vmin=lo, vmax=hi, origin="lower",
                       extent=ext, interpolation="nearest")
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("x [cm]")
        if title.startswith("$S_1"):
            norm_im = im
    axes[0].set_ylabel("y [cm]  (wind along +x)")
    # one shared colorbar for the normalized panels (-1..1); S0 has none
    cb = fig.colorbar(norm_im, ax=[axes[1], axes[2]], fraction=0.025,
                      pad=0.02)
    cb.set_label("normalized Stokes")
    fig.suptitle(f"Elfouhaily + FM98 capillary blend  |  "
                 f"U10 = {U10:.0f} m/s  |  1 m x 1 m, {N}$^2$ px",
                 fontsize=13, y=1.02)
    return fig, (vmin, vmax)


def main():
    PANELS.mkdir(parents=True, exist_ok=True)
    if backend.has_torch():
        import torch
        if torch.cuda.is_available():
            bk, dev = "torch", "cuda"
            print(f"backend: torch / {torch.cuda.get_device_name(0)}")
        else:
            bk, dev = "torch", "cpu"
            print("backend: torch / cpu")
    else:
        bk, dev = None, None
        print("backend: numpy")
    xp = backend.get_xp(bk, dev)

    table = get_table(np.pi / (L / N))
    for U10 in WINDS:
        import time
        t0 = time.perf_counter()
        S, info = render_stokes(U10, table, xp, bk, dev)
        fig, (vmin, vmax) = panel_figure(U10, S)
        path = PANELS / f"stokes_panel_U{int(U10):02d}.png"
        fig.savefig(path, bbox_inches="tight")
        if U10 == 11.0:        # representative figure for the gallery
            fig.savefig(OUT / "demo_stokes_panels.png", bbox_inches="tight")
        plt.close(fig)
        finite = np.isfinite(S[..., 0])
        dolp = (np.hypot(S[..., 1], S[..., 2])
                / np.where(finite, S[..., 0], np.nan))
        print(f"U10={U10:>4.0f}  Hs={info['Hs_realized']*1e3:5.1f} mm  "
              f"S0 lims=[{vmin:.4f}, {vmax:.4f}]  "
              f"mean DoLP={np.nanmean(dolp):.3f}  "
              f"({time.perf_counter()-t0:.2f} s)  -> {path.name}")
    print(f"\nwrote {len(WINDS)} panels to {PANELS}")


if __name__ == "__main__":
    main()
