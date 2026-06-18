"""Color-aware S0 animation: seapol vs Mitsuba over an evolving wave field.

Renders a time-evolving seapol surface's reflected-sky intensity (S0) in color,
both with seapol (spectral pipeline) and Mitsuba (reflection-only dielectric),
side by side, and encodes an mp4.  Default 10 s at 30 fps (300 frames).

Temporal stability: the absolute-radiance scale (alpha) and the display exposure
are fit ONCE on the first frame and held fixed, so brightness does not flicker.

Usage:
    uv run validation/s0_color_anim.py [--out validation/anim] [--seconds 10]
        [--fps 30] [--spp 2048] [--n 256] [--res 320] [--frames N(dry-run)]
"""

from __future__ import annotations

import argparse
import os
import subprocess

import numpy as np

from scene import SceneSpec, export_obj

LUM = np.array([0.2126, 0.7152, 0.0722])


def _srgb(x):
    x = np.clip(x, 0.0, 1.0)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * x ** (1 / 2.4) - 0.055)


def seapol_frame(bands, skies, eta, dx, info, cam, n_water, sub):
    from seapol import render_camera_image_spectral
    S = render_camera_image_spectral(eta, dx, bands, skies, camera=cam,
                                     n_water=n_water, slope_x=info["sx"],
                                     slope_y=info["sy"], subpixel=sub,
                                     n_subpixel=8, seed=5)
    S0 = np.asarray(S)[..., 0]                          # (H, W, 3) per band
    return np.stack([S0[..., 2], S0[..., 1], S0[..., 0]], axis=-1)  # R,G,B


def mitsuba_frame(mi, spec, obj, spp):
    cb = spec.camera_basis()
    H, W = cb["H"], cb["W"]
    to_world = mi.ScalarTransform4f().look_at(
        origin=mi.ScalarPoint3f(*cb["origin"]),
        target=mi.ScalarPoint3f(*cb["center"]),
        up=mi.ScalarVector3f(*cb["up"]))
    scene = mi.load_dict({
        "type": "scene",
        "integrator": {"type": "path", "max_depth": 4},
        "sensor": {"type": "perspective", "fov_axis": "x",
                   "fov": float(2.0 * spec.hfov_deg), "to_world": to_world,
                   "film": {"type": "hdrfilm", "width": W, "height": H,
                            "pixel_format": "rgb", "rfilter": {"type": "box"}},
                   "sampler": {"type": "independent", "sample_count": spp}},
        "surface": {"type": "obj", "filename": obj,
                    "bsdf": {"type": "dielectric",
                             "int_ior": spec.refractive_index(), "ext_ior": 1.0,
                             "specular_transmittance": 0.0}},
        "sky": {"type": "constant",
                "radiance": {"type": "rgb", "value": list(spec.sky_rgb)}},
    })
    return np.asarray(mi.render(scene, spp=spp), dtype=np.float64)[..., :3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/anim")
    ap.add_argument("--seconds", type=float, default=10.0)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--spp", type=int, default=2048)
    ap.add_argument("--n", type=int, default=256)       # surface grid
    ap.add_argument("--res", type=int, default=320)      # image size
    ap.add_argument("--frames", type=int, default=0)     # >0: dry-run subset
    args = ap.parse_args()
    out = args.out
    fdir = os.path.join(out, "frames")
    os.makedirs(fdir, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from seapol import (PinholeCamera, SpectralBands, SubpixelSlopes,
                        generate_sea_surface, spectral_sky_factories)

    nframes = int(round(args.seconds * args.fps))
    if args.frames:
        nframes = args.frames
    spec = SceneSpec(L=8.0, N=args.n, U10=7.0, seed=0,
                     altitude_m=12.0, zenith_deg=40.0, azimuth_deg=0.0,
                     hfov_deg=18.0, img_h=args.res, img_w=args.res,
                     sun_zenith_deg=45.0, sun_azimuth_deg=90.0)
    n_water = spec.refractive_index()
    cam = PinholeCamera(altitude_m=spec.altitude_m, zenith_deg=spec.zenith_deg,
                        azimuth_deg=spec.azimuth_deg, hfov_deg=spec.hfov_deg,
                        img_shape=spec.img_shape)
    bands = SpectralBands.rgb()
    skies = spectral_sky_factories(bands, "clear",
                                   sun_zenith_deg=spec.sun_zenith_deg,
                                   sun_azimuth_deg=spec.sun_azimuth_deg,
                                   turbidity=0.1)

    print(f"generating evolving surface N={args.n}, {nframes} frames ...")
    times = np.arange(nframes) / args.fps
    surf = generate_sea_surface(L=spec.L, N=spec.N, U10=spec.U10,
                                times=times, rng=np.random.default_rng(spec.seed))
    dx = surf.info["dx"]
    sub = SubpixelSlopes(surf.info.get("sigma_a2_cut", 0.0),
                         surf.info.get("sigma_c2_cut", 0.0))

    import mitsuba as mi
    for v in ("scalar_rgb", "llvm_ad_rgb"):
        if v in mi.variants():
            mi.set_variant(v)
            break

    obj = os.path.join(out, "_frame.obj")
    alpha = exposure = None
    for i in range(nframes):
        eta = surf.eta[:, :, i]
        info = {"sx": surf.slope_x[:, :, i], "sy": surf.slope_y[:, :, i]}
        sp = seapol_frame(bands, skies, eta, dx, info, cam, n_water, sub)
        export_obj(eta, dx, obj)
        miL = mitsuba_frame(mi, spec, obj, args.spp)

        if alpha is None:                       # fix scale/exposure on frame 0
            valid = np.isfinite(sp[..., 0]) & (miL @ LUM > 0)
            Lsp, Lmi = sp @ LUM, miL @ LUM
            alpha = float(np.sum(Lsp[valid] * Lmi[valid])
                          / np.sum(Lmi[valid] ** 2))
            exposure = float(np.quantile(Lsp[valid], 0.99))
            print(f"alpha={alpha:.3g}, exposure={exposure:.3g}")

        sp_d = _srgb(np.nan_to_num(sp / max(exposure, 1e-9)))
        mi_d = _srgb(np.nan_to_num(miL * alpha / max(exposure, 1e-9)))
        fig, ax = plt.subplots(1, 2, figsize=(8.2, 4.2))
        ax[0].imshow(np.clip(sp_d, 0, 1)); ax[0].set_title("seapol S0", fontsize=10)
        ax[1].imshow(np.clip(mi_d, 0, 1))
        ax[1].set_title(f"Mitsuba S0 (x{alpha:.2g})", fontsize=10)
        for a in ax:
            a.axis("off")
        fig.suptitle(f"Color S0 off an evolving wave field   t = {times[i]:5.2f} s",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(fdir, f"f{i:04d}.png"), dpi=110)
        plt.close(fig)
        if i % 10 == 0 or i == nframes - 1:
            print(f"  frame {i+1}/{nframes}")

    mp4 = os.path.join(out, "s0_color.mp4")
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-i", os.path.join(fdir, "f%04d.png"),
           # libx264/yuv420p require even dimensions
           "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", mp4]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    print(f"wrote {mp4} ({nframes} frames @ {args.fps} fps)")


if __name__ == "__main__":
    main()
