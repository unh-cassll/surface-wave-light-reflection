"""Render the exported surface in Mitsuba 3 (polarized) for a STOKES
cross-check of seapol's polarized reflection.

Mitsuba 3 with a *_polarized variant traces the full Stokes vector through
polarized Fresnel at the dielectric interface, so the reflected DoLP/AoP it
produces should match seapol's Mueller chain.

Usage:
    uv run validation/render_mitsuba.py validation/output/scene.json \
        validation/output/surface.obj --out validation/output

Sky handling:
    The environment here is UNPOLARIZED (constant radiance I_sky).  Under an
    unpolarized sky the reflected polarization is generated entirely by the
    surface Fresnel reflection -- which is exactly the quantity to validate,
    and for a flat surface (SceneSpec.flat=True) reduces to analytic Fresnel
    in both codes.  Matching seapol's *polarized* Coulson sky would require a
    custom polarized environment emitter (precompute rayleigh_sky_stokes onto
    an envmap with per-texel Stokes); that is a documented extension, not done
    here, because it does not affect the flat-surface Fresnel sanity check.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from scene import SceneSpec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scene")
    ap.add_argument("obj")
    ap.add_argument("--out", default="validation/output")
    # High spp by default: the reflection-only dielectric suppresses the
    # transmission lobe, so only ~R (a few percent) of BSDF samples carry the
    # reflected signal -- low spp speckles badly, worst in DoLP (a ratio over a
    # noisy small S0).  spp also averages sub-pixel facets, like seapol.
    ap.add_argument("--spp", type=int, default=2048)
    args = ap.parse_args()

    try:
        import mitsuba as mi
    except ImportError:
        raise SystemExit("mitsuba not installed: pip install mitsuba "
                         "(or `uv sync --extra validation`)")
    # spectral + polarized; fall back to llvm if scalar is unavailable
    for v in ("scalar_spectral_polarized", "llvm_ad_spectral_polarized"):
        if v in mi.variants():
            mi.set_variant(v)
            break
    else:
        raise SystemExit(f"no polarized variant available; have {mi.variants()}")

    spec = SceneSpec.from_json(args.scene)
    cb = spec.camera_basis()
    n_water = spec.n_water if spec.n_water > 0 else 1.34
    H, W = cb["H"], cb["W"]

    origin = cb["origin"]
    target = cb["center"]
    up = cb["up"]
    to_world = mi.ScalarTransform4f().look_at(
        origin=mi.ScalarPoint3f(*origin),
        target=mi.ScalarPoint3f(*target),
        up=mi.ScalarVector3f(*up))

    scene = mi.load_dict({
        "type": "scene",
        "integrator": {"type": "stokes", "nested": {"type": "path",
                                                    "max_depth": 8}},
        "sensor": {
            "type": "perspective",
            "fov_axis": "larger",
            "fov": float(2.0 * spec.hfov_deg),
            "to_world": to_world,
            "film": {"type": "hdrfilm", "width": W, "height": H,
                     "pixel_format": "rgb", "rfilter": {"type": "box"}},
            "sampler": {"type": "independent", "sample_count": args.spp},
        },
        "surface": {
            "type": "obj",
            "filename": args.obj,
            # Reflection-only dielectric: polarized Fresnel reflection Mueller
            # matrix, transmitted radiance suppressed (specular_transmittance=0)
            # to match seapol's render_camera_image with water=None (dark-water
            # sky-reflection model).  The reflection Mueller is unaffected.
            "bsdf": {"type": "dielectric", "int_ior": n_water, "ext_ior": 1.0,
                     "specular_transmittance": 0.0},
        },
        "sky": {"type": "constant",
                "radiance": {"type": "uniform", "value": float(spec.I_sky)}},
    })

    img = mi.render(scene, spp=args.spp)
    # The stokes integrator film has 15 channels: the nested RGB output (S0)
    # followed by S0, S1, S2, S3 each as an RGB AOV.  Reshape to 5 groups of
    # RGB and average each group to a scalar; groups 1..4 are S0..S3.
    arr = np.asarray(img, dtype=np.float64).reshape(H, W, 5, 3).mean(axis=-1)
    S0, S1, S2, S3 = arr[..., 1], arr[..., 2], arr[..., 3], arr[..., 4]
    S = np.stack([S0, S1, S2, S3], axis=-1)

    with np.errstate(divide="ignore", invalid="ignore"):
        dolp = np.sqrt(S1**2 + S2**2) / S0
        aop = 0.5 * np.arctan2(S2, S1)

    os.makedirs(args.out, exist_ok=True)
    np.savez(os.path.join(args.out, "mitsuba.npz"),
             I=S0, S=S, dolp=dolp, aop=aop, n_water=n_water)
    print(f"mitsuba render: {H}x{W}, I mean={np.nanmean(S0):.4g}, "
          f"DoLP mean={np.nanmean(dolp):.3f}")
    print(f"wrote mitsuba.npz to {args.out}/")


if __name__ == "__main__":
    main()
