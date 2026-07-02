"""Compare seapol against available external renders and emit metrics + figures.

Loads the seapol baseline (always) and whatever external renders are present:
    validation/output/seapol.npz    (required)
    validation/output/blender.exr   (optional, intensity)
    validation/output/mitsuba.npz   (optional, Stokes)

Degrades gracefully: with neither external render present it still reports the
seapol baseline statistics, so the harness runs end to end on a bare install.

Usage:
    uv run validation/compare.py [--out validation/output]
"""

from __future__ import annotations

import argparse
import os

import numpy as np


def _valid(*arrs):
    m = np.ones(arrs[0].shape, dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    return m


def _norm_unit_mean(x, mask):
    mu = np.nanmean(x[mask])
    return x / mu if mu else x


def _circ_aop_err(a, b):
    """Angle-of-polarization error folded to [0, pi/2] (pi-periodic axis)."""
    d = np.abs(a - b) % np.pi
    return np.minimum(d, np.pi - d)


def _aop_to_meridian(out, aop_sp, aop_mi, mask):
    """Rotate Mitsuba's camera-frame AoP into seapol's meridian frame using the
    per-pixel meridian angle from the scene camera geometry, with a 2-DOF
    (sign, global constant) convention fit.  Returns the rotated AoP, or None if
    scene.json is unavailable."""
    sj = os.path.join(out, "scene.json")
    if not os.path.exists(sj):
        return None
    from aop_align import camera_dirs_and_meridian
    from scene import SceneSpec
    gamma = camera_dirs_and_meridian(SceneSpec.from_json(sj))
    cs = np.deg2rad(np.arange(0, 180))
    best = None
    for s in (+1.0, -1.0):
        base = aop_mi + s * gamma
        for c in cs:
            e = np.nanmean(_circ_aop_err(aop_sp[mask], base[mask] + c))
            if best is None or e < best[0]:
                best = (e, base + c)
    return best[1]


def load_blender_intensity(out):
    """Scalar Blender radiance image, or None.  Prefers the top-down .npy that
    render_blender.py writes (no external EXR reader needed); falls back to the
    EXR via imageio/cv2 if only the EXR is present."""
    npy = os.path.join(out, "blender.npy")
    if os.path.exists(npy):
        rgb = np.load(npy).astype(np.float64)
        return rgb[..., :3].mean(axis=-1) if rgb.ndim == 3 else rgb
    exr = os.path.join(out, "blender.exr")
    if not os.path.exists(exr):
        return None
    try:
        import imageio.v3 as iio
        rgb = np.asarray(iio.imread(exr), dtype=np.float64)
    except Exception:
        try:
            import cv2
            rgb = cv2.imread(exr, cv2.IMREAD_UNCHANGED | cv2.IMREAD_ANYDEPTH)
            rgb = np.asarray(rgb, dtype=np.float64)[..., ::-1]
        except Exception as e:
            print(f"  [blender] cannot read EXR ({e}); install blender.npy "
                  f"route or imageio[freeimage]; skipping")
            return None
    return rgb[..., :3].mean(axis=-1) if rgb.ndim == 3 else rgb


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/output")
    args = ap.parse_args()
    out = args.out

    sp_path = os.path.join(out, "seapol.npz")
    if not os.path.exists(sp_path):
        raise SystemExit(f"missing {sp_path}; run render_seapol.py first")
    sp = np.load(sp_path)
    I0, dolp0, aop0 = sp["I"], sp["dolp"], sp["aop"]
    base_valid = _valid(I0)
    print(f"seapol baseline: {I0.shape}, {base_valid.mean()*100:.1f}% valid; "
          f"I mean={np.nanmean(I0):.4g}, DoLP mean={np.nanmean(dolp0):.3f}")

    panels = [("seapol I", I0), ("seapol DoLP", dolp0)]

    # --- Blender (intensity) ----------------------------------------------
    Ib = load_blender_intensity(out)
    if Ib is None:
        print("[blender] no render present; skipping intensity comparison")
    elif Ib.shape != I0.shape:
        print(f"[blender] shape {Ib.shape} != seapol {I0.shape}; skipping")
    else:
        m = _valid(I0, Ib) & (Ib > 0)
        a = _norm_unit_mean(I0, m)
        b = _norm_unit_mean(Ib, m)
        rmse = np.sqrt(np.nanmean((a[m] - b[m]) ** 2))
        corr = np.corrcoef(a[m], b[m])[0, 1]
        print(f"[blender] normalized-radiance RMSE={rmse:.3f}, corr={corr:.3f} "
              f"over {m.sum()} px (INDICATIVE: Blender is unpolarized and uses a "
              f"different dielectric/illumination model -- qualitative I check)")
        panels.append(("blender I", Ib))

    # --- Mitsuba (Stokes) -------------------------------------------------
    mi_path = os.path.join(out, "mitsuba.npz")
    if not os.path.exists(mi_path):
        print("[mitsuba] no render present; skipping polarized comparison")
    else:
        ms = np.load(mi_path)
        Im, dolpm, aopm = ms["I"], ms["dolp"], ms["aop"]
        if Im.shape != I0.shape:
            print(f"[mitsuba] shape {Im.shape} != seapol {I0.shape}; skipping")
        else:
            # drop near-zero-radiance pixels (DoLP/AoP are ill-defined there)
            thr = 1e-3 * np.nanmean(Im[np.isfinite(Im)])
            m = _valid(I0, Im, dolp0, dolpm) & (Im > thr) & (I0 > 0)
            dd = np.abs(dolp0 - dolpm)
            corrI = np.corrcoef(_norm_unit_mean(I0, m)[m],
                                _norm_unit_mean(Im, m)[m])[0, 1]
            print(f"[mitsuba] DoLP |err| mean={np.nanmean(dd[m]):.4f} "
                  f"median={np.nanmedian(dd[m]):.4f} max={np.nanmax(dd[m]):.3f} "
                  f"over {m.sum()} px")
            # the DoLP comparison validates Fresnel only when both codes
            # see the same (unpolarized) sky; Mitsuba's environment is
            # always unpolarized, seapol's default is a Rayleigh sky
            sj = os.path.join(out, "scene.json")
            unpol = False
            if os.path.exists(sj):
                from scene import SceneSpec
                unpol = SceneSpec.from_json(sj).unpolarized_sky
            tag = ("frame-invariant -> validates Fresnel" if unpol else
                   "skies differ: seapol Rayleigh vs Mitsuba unpolarized; "
                   "set unpolarized_sky=true for a Fresnel validation")
            print(f"  seapol DoLP mean={np.nanmean(dolp0[m]):.3f}  "
                  f"mitsuba DoLP mean={np.nanmean(dolpm[m]):.3f}  ({tag})")
            # AoP: rotate Mitsuba's camera-frame AoP into seapol's meridian
            # frame using the per-pixel geometric meridian angle (a fixed global
            # offset is insufficient -- the rotation varies per pixel).
            aop_panel = aopm
            aopm_merid = _aop_to_meridian(out, aop0, aopm, m)
            if aopm_merid is not None:
                da = _circ_aop_err(aop0, aopm_merid)
                print(f"  AoP |err| raw={np.rad2deg(np.nanmean(_circ_aop_err(aop0, aopm)[m])):.1f} deg"
                      f" -> {np.rad2deg(np.nanmean(da[m])):.1f} deg after meridian rotation"
                      f" (validates AoP)")
                aop_panel = (aopm_merid + np.pi / 2) % np.pi - np.pi / 2
            else:
                print("  AoP: scene.json absent; skipping meridian rotation")
            print(f"  I corr={corrI:.3f} (low when incidence varies little -> "
                  f"near-flat I field, not a disagreement)")
            panels += [("mitsuba DoLP", dolpm),
                       ("mitsuba AoP->meridian", np.where(m, aop_panel, np.nan))]

    _save_panels(panels, os.path.join(out, "comparison.png"))


def _save_panels(panels, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping figure")
        return
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.2))
    axes = np.atleast_1d(axes)
    for ax, (title, img) in zip(axes, panels):
        # AoP is a pi-periodic axis -> cyclic colormap so +/-90 deg read alike
        cmap = "twilight" if "AoP" in title else "viridis"
        im = ax.imshow(np.asarray(img), cmap=cmap)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
