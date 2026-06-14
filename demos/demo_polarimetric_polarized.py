"""
Polarimetric slope reconstruction under a CLEAR (polarized) sky at an
off-Brewster viewing angle -- the geometry a division-of-focal-plane
(DoFP) camera actually sees in the field (e.g. 30 deg incidence).

The standard DoLP-based inversion (slopes_from_stokes) assumes
unpolarized illumination: it reads the facet incidence angle off the
Fresnel DoLP, which is only valid under an overcast sky.  Under a clear
sky the incident skylight is already polarized, so reflected Q/U mix the
sky polarization with the Fresnel rotation and the DoLP no longer maps
to facet orientation.  slopes_from_stokes_polarized inverts the full
Mueller chain against the known sky model instead, recovering the slopes
correctly.

This demo renders a synthetic surface under a clear sky at 30 deg, then
compares both inversions against the known truth and integrates the
polarized-inversion slopes to height.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (CameraGeometry, generate_sea_surface,
                    height_from_slopes, make_clear_sky, render_facet_stokes,
                    slopes_from_stokes, slopes_from_stokes_polarized)
from seapol.polarization import normalize, stokes_dolp

OUT = Path(__file__).parent / "output"
L, N, U10 = 12.0, 384, 6.0
INCIDENCE, SUN_ZEN, SUN_AZ = 30.0, 50.0, 90.0


def _view_dirs(eta, dx, cam):
    H, W = eta.shape
    xg = np.arange(W) * dx
    X, Y = np.meshgrid(xg, np.arange(H) * dx)
    P = np.stack([X, Y, eta], axis=-1)
    cpos = cam.position(((W - 1) * dx / 2.0, (H - 1) * dx / 2.0))
    return normalize(cpos[None, None, :] - P)


def main():
    OUT.mkdir(exist_ok=True)
    surf = generate_sea_surface(L=L, N=N, U10=U10,
                                rng=np.random.default_rng(0))
    dx = surf.info["dx"]
    cam = CameraGeometry(incidence_deg=INCIDENCE, azimuth_deg=0.0,
                         height_m=120.0)
    sky = make_clear_sky(SUN_ZEN, SUN_AZ, 1.0, turbidity=0.05)
    S = render_facet_stokes(surf.eta, dx, camera=cam, sky=sky,
                            slope_x=surf.slope_x, slope_y=surf.slope_y,
                            subpixel=None)
    d_out = _view_dirs(surf.eta, dx, cam)

    sxu, syu, vu = slopes_from_stokes(S, d_out)             # unpolarized model
    sxp, syp, vp = slopes_from_stokes_polarized(S, d_out, sky)  # full Mueller
    eta_p = height_from_slopes(sxp, syp, dx)
    eta_t = surf.eta - surf.eta.mean()

    def corr(a, b, m):
        return np.corrcoef(a[m], b[m])[0, 1]

    mu = vu & np.isfinite(sxu)
    mp = vp & np.isfinite(sxp)
    cpx = corr(sxp, surf.slope_x, mp)
    cpy = corr(syp, surf.slope_y, mp)
    cux = corr(sxu, surf.slope_x, mu)
    ceta = np.corrcoef(eta_p.ravel(), eta_t.ravel())[0, 1]
    print(f"clear sky, {INCIDENCE:.0f} deg incidence (DoFP geometry):")
    print(f"  unpolarized model  : slope_x corr {cux:.3f}")
    print(f"  polarized model    : slope_x corr {cpx:.3f}, "
          f"slope_y corr {cpy:.3f}")
    print(f"  height corr {ceta:.4f}, Hs {4 * eta_p.std():.3f} -> "
          f"{4 * eta_t.std():.3f} m")

    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.6), dpi=120)
    sl = dict(cmap="RdBu_r", vmin=-0.4, vmax=0.4, origin="lower")
    ax[0, 0].imshow(surf.slope_x, **sl); ax[0, 0].set_title("truth slope_x")
    ax[0, 1].imshow(np.where(mp, sxp, np.nan), **sl)
    ax[0, 1].set_title(f"polarized inversion\nslope_x (corr {cpx:.3f})")
    ax[0, 2].imshow(np.where(mu, sxu, np.nan), **sl)
    ax[0, 2].set_title(f"unpolarized inversion\nslope_x (corr {cux:.3f})")
    for a in ax[0]:
        a.axis("off")

    ax[1, 0].imshow(stokes_dolp(S), cmap="viridis", vmin=0, vmax=1,
                    origin="lower")
    ax[1, 0].set_title("observed DoLP (clear sky)")
    ax[1, 0].axis("off")
    a = ax[1, 1]
    a.scatter(surf.slope_x[mp][::13], sxp[mp][::13], s=2, alpha=0.3,
              label=f"polarized (corr {cpx:.3f})")
    a.scatter(surf.slope_x[mu][::13], sxu[mu][::13], s=2, alpha=0.3,
              color="C3", label=f"unpolarized (corr {cux:.3f})")
    a.plot([-0.5, 0.5], [-0.5, 0.5], "k--", lw=0.8)
    a.set_xlabel("true slope_x"); a.set_ylabel("recovered slope_x")
    a.set_xlim(-0.5, 0.5); a.set_ylim(-0.5, 0.5)
    a.legend(fontsize=8, markerscale=4); a.set_title("slope recovery")
    a.grid(alpha=0.3)
    im = ax[1, 2].imshow(eta_p, cmap="coolwarm", origin="lower")
    ax[1, 2].set_title(f"reconstructed height\n(corr {ceta:.4f})")
    ax[1, 2].axis("off")
    plt.colorbar(im, ax=ax[1, 2], fraction=0.046, label="eta [m]")

    fig.suptitle("Polarimetric reconstruction under a clear (polarized) "
                 f"sky at {INCIDENCE:.0f} deg -- DoFP field geometry", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "demo_polarimetric_polarized.png", bbox_inches="tight")
    print(f"wrote {OUT / 'demo_polarimetric_polarized.png'}")


if __name__ == "__main__":
    main()
