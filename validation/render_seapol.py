"""Render the seapol reference image for a SceneSpec.

Produces the baseline every external renderer is compared against:
    <out>/seapol.npz     I, Stokes (H,W,4), DoLP, AoP
    <out>/surface.obj    the exported mesh (for Blender / Mitsuba)
    <out>/scene.json     the spec (copied so the run is self-contained)

Usage:
    uv run validation/render_seapol.py [scene.json] [--out validation/output]
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from scene import SceneSpec, export_obj, generate_surface


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scene", nargs="?", default=None,
                    help="SceneSpec JSON (default: built-in defaults)")
    ap.add_argument("--out", default="validation/output")
    args = ap.parse_args()

    spec = SceneSpec.from_json(args.scene) if args.scene else SceneSpec()
    os.makedirs(args.out, exist_ok=True)

    from seapol import (PinholeCamera, SubpixelSlopes, make_rayleigh_sky,
                        make_unpolarized_sky, render_camera_image, stokes_aop,
                        stokes_dolp)

    eta, dx, info = generate_surface(spec)
    n_water = spec.refractive_index()

    cam = PinholeCamera(altitude_m=spec.altitude_m, zenith_deg=spec.zenith_deg,
                        azimuth_deg=spec.azimuth_deg, hfov_deg=spec.hfov_deg,
                        img_shape=spec.img_shape)
    if spec.unpolarized_sky:
        sky = make_unpolarized_sky(spec.I_sky)
    else:
        sky = make_rayleigh_sky(sun_zenith_deg=spec.sun_zenith_deg,
                                sun_azimuth_deg=spec.sun_azimuth_deg,
                                I_sky=spec.I_sky)
    sub = SubpixelSlopes(info.get("sigma_a2_cut", 0.0),
                         info.get("sigma_c2_cut", 0.0))

    # subpixel ensemble seeded (seed + 1; the surface draw uses seed) so
    # repeated runs of the same scene.json are bit-reproducible
    S = render_camera_image(eta, dx, camera=cam, sky=sky,
                            slope_x=info.get("slope_x"),
                            slope_y=info.get("slope_y"),
                            n_water=n_water, subpixel=sub, n_subpixel=16,
                            shadowing=True,
                            rng=np.random.default_rng(spec.seed + 1))
    S = np.asarray(S)
    I = S[..., 0]
    dolp = np.asarray(stokes_dolp(S))
    aop = np.asarray(stokes_aop(S))

    np.savez(os.path.join(args.out, "seapol.npz"),
             I=I, S=S, dolp=dolp, aop=aop, n_water=n_water)
    export_obj(eta, dx, os.path.join(args.out, "surface.obj"))
    # persist the resolved n so Blender/Mitsuba (which cannot import
    # seapol's Quan & Fry) use the identical refractive index
    spec.n_water = n_water
    spec.to_json(os.path.join(args.out, "scene.json"))

    valid = np.isfinite(I)
    print(f"seapol render: {S.shape[0]}x{S.shape[1]}, "
          f"{valid.mean()*100:.1f}% valid pixels, n_water={n_water:.4f}")
    print(f"  I (valid):    mean={np.nanmean(I):.4g}  max={np.nanmax(I):.4g}")
    print(f"  DoLP (valid): mean={np.nanmean(dolp):.3f}  max={np.nanmax(dolp):.3f}")
    print(f"wrote seapol.npz, surface.obj, scene.json to {args.out}/")


if __name__ == "__main__":
    main()
