"""Rotate Mitsuba's AoP into seapol's meridian frame and compare (out_flat).

seapol reports the angle of polarization in its meridian Stokes frame
(e_perp = z x d, e_par = e_perp x d, +Q along e_par); Mitsuba reports it in its
sensor/camera frame.  The two differ by a PER-PIXEL rotation -- the angle of the
meridian e_par measured in the camera (right, up) frame -- not a single global
offset, which is why a constant shift leaves a large residual.

This computes that per-pixel meridian angle gamma from the camera geometry
(replicating seapol.render._camera_rays), then tests the standard frame
relationship  AoP_seapol = AoP_mitsuba + s*gamma + c  with a sign s in {+1,-1}
and ONE global convention constant c (grid-searched).  A small residual after
this 2-DOF fit confirms the AoP agrees once expressed in a common frame.

Usage:
    uv run validation/aop_align.py --out validation/out_flat
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from scene import SceneSpec


def camera_dirs_and_meridian(spec):
    """Per-pixel observed direction d (toward camera) and meridian basis,
    replicating seapol.render._camera_rays (d_out = -ray)."""
    cb = spec.camera_basis()
    look, right, up = cb["look"], cb["right"], cb["up"]
    H, W, half = cb["H"], cb["W"], cb["half"]
    step = 2.0 * half / max(H, W)   # edge-aligned pixel centers, as in seapol
    ys = (np.arange(H) - (H - 1) / 2.0) * step
    xs = (np.arange(W) - (W - 1) / 2.0) * step
    XS, YS = np.meshgrid(xs, -ys)
    dirs = (look[None, None, :] + right[None, None, :] * XS[..., None]
            + up[None, None, :] * YS[..., None])
    dirs /= np.linalg.norm(dirs, axis=-1, keepdims=True)
    d = -dirs                                    # observed/outgoing direction
    z = np.array([0.0, 0.0, 1.0])
    e_perp = np.cross(z, d)
    e_perp /= np.linalg.norm(e_perp, axis=-1, keepdims=True)
    e_par = np.cross(e_perp, d)
    # meridian angle measured in the camera (right, up) frame
    gamma = np.arctan2(e_par @ up, e_par @ right)
    return gamma


def _fold(a):
    """Fold an AoP difference (pi-periodic) to [0, pi/2]."""
    d = np.abs(a) % np.pi
    return np.minimum(d, np.pi - d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/out_flat")
    args = ap.parse_args()
    out = args.out

    spec = SceneSpec.from_json(os.path.join(out, "scene.json"))
    sp = np.load(os.path.join(out, "seapol.npz"))
    ms = np.load(os.path.join(out, "mitsuba.npz"))
    aop_sp, aop_mi = sp["aop"], ms["aop"]
    m = np.isfinite(aop_sp) & np.isfinite(aop_mi) & np.isfinite(sp["I"]) \
        & (ms["I"] > 1e-3 * np.nanmean(ms["I"]))

    gamma = camera_dirs_and_meridian(spec)

    # 2-DOF fit: sign s in {+1,-1}, global constant c (1-deg grid)
    cs = np.deg2rad(np.arange(0, 180))
    best = None
    for s in (+1.0, -1.0):
        base = aop_mi + s * gamma
        for c in cs:
            err = np.rad2deg(np.nanmean(_fold(aop_sp[m] - (base[m] + c))))
            if best is None or err < best[0]:
                best = (err, s, c, base)
    res, s, c, base = best
    aop_mi_merid = base + c
    raw = np.rad2deg(np.nanmean(_fold(aop_sp[m] - aop_mi[m])))
    print(f"AoP raw |err| mean = {raw:.1f} deg")
    print(f"after per-pixel meridian rotation (s={s:+.0f}, c={np.rad2deg(c):.1f} deg): "
          f"|err| mean = {res:.1f} deg over {int(m.sum())} px")

    _figure(out, aop_sp, aop_mi, aop_mi_merid, m, res)


def _figure(out, aop_sp, aop_mi, aop_mi_merid, m, res):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping figure")
        return
    # wrap rotated AoP back to [-pi/2, pi/2] for display
    disp = (aop_mi_merid + np.pi / 2) % np.pi - np.pi / 2
    disp = np.where(m, disp, np.nan)
    sp_d = np.where(m, aop_sp, np.nan)
    mi_d = np.where(m, aop_mi, np.nan)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, img, t in zip(
            ax, [sp_d, mi_d, disp],
            ["seapol AoP (meridian)", "Mitsuba AoP (camera frame)",
             f"Mitsuba AoP -> meridian  (|err| {res:.1f} deg)"]):
        im = a.imshow(np.rad2deg(img), cmap="twilight", vmin=-90, vmax=90)
        a.set_title(t, fontsize=10)
        a.axis("off")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="AoP [deg]")
    fig.suptitle("Angle of polarization: seapol vs Mitsuba in a common frame "
                 "(out_flat)", fontsize=11)
    path = os.path.join(out, "aop_align.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
