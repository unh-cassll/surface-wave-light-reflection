"""Color-aware S0 cross-comparison (seapol vs Mitsuba) off one synthetic wave field.

Renders one seapol surface's reflected-sky intensity (S0) in color two ways and
compares them on a common ABSOLUTE radiance scale:
    seapol   -- spectral pipeline (SpectralBands.rgb), reflection-only
    Mitsuba  -- spectral_polarized, reflection-only dielectric, colored sky

Both compute the same physics (dielectric Fresnel reflection of the sky), so
after a single absolute-radiance scale factor they should agree closely; the
figure shows seapol | Mitsuba(scaled) | luminance difference.

Notes:
  * High spp is used because the reflection-only dielectric suppresses the
    transmission lobe, so only ~R (a few percent) of BSDF samples carry the
    reflected signal -- low spp therefore speckles badly.  High spp also
    averages sub-pixel facets, matching seapol's sub-pixel slope averaging.
  * Absolute radiance is matched by a single least-squares luminance scale
    (alpha) fit over the overlapping valid pixels, so the comparison is about
    spatial/chromatic structure on a shared scale, not arbitrary exposure.
  * Blender is intentionally not in this figure (unpolarized, different sky/
    dielectric); use render_blender.py for a separate qualitative intensity view.
  * DoLP is NOT compared here -- polarization validation stays restricted to the
    Mitsuba run in render_mitsuba.py / compare.py.

Usage:
    uv run validation/s0_color.py [--out validation/output_color] [--spp 4096]
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from scene import SceneSpec, export_obj, generate_surface

LUM = np.array([0.2126, 0.7152, 0.0722])


def _srgb(x):
    x = np.clip(x, 0.0, 1.0)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * x ** (1 / 2.4) - 0.055)


def render_seapol_linear(spec, eta, dx, info, out):
    """seapol reflection-only spectral render -> linear RGB radiance (H,W,3),
    NaN where the camera ray misses the patch."""
    from seapol import (PinholeCamera, SpectralBands, SubpixelSlopes,
                        render_camera_image_spectral, spectral_sky_factories)
    bands = SpectralBands.rgb()                       # [450, 550, 650] nm
    cam = PinholeCamera(altitude_m=spec.altitude_m, zenith_deg=spec.zenith_deg,
                        azimuth_deg=spec.azimuth_deg, hfov_deg=spec.hfov_deg,
                        img_shape=spec.img_shape)
    skies = spectral_sky_factories(bands, "clear",
                                   sun_zenith_deg=spec.sun_zenith_deg,
                                   sun_azimuth_deg=spec.sun_azimuth_deg,
                                   turbidity=0.1)
    sub = SubpixelSlopes(info.get("sigma_a2_cut", 0.0),
                         info.get("sigma_c2_cut", 0.0))
    S = render_camera_image_spectral(eta, dx, bands, skies, camera=cam,
                                     n_water=spec.refractive_index(),
                                     slope_x=info.get("slope_x"),
                                     slope_y=info.get("slope_y"),
                                     subpixel=sub, n_subpixel=10, seed=5)
    S0 = np.asarray(S)[..., 0]                         # (H, W, 3) per band
    lin = np.stack([S0[..., 2], S0[..., 1], S0[..., 0]], axis=-1)  # R,G,B
    np.save(os.path.join(out, "seapol_lin.npy"), lin)
    print(f"seapol color render done (valid {np.isfinite(lin[...,0]).mean()*100:.0f}%)")
    return lin


def render_mitsuba_linear(spec, obj, out, spp):
    """Mitsuba reflection-only render -> linear RGB radiance (H,W,3)."""
    try:
        import mitsuba as mi
    except ImportError:
        print("[mitsuba] not installed; skipping")
        return None
    for v in ("scalar_spectral_polarized", "scalar_rgb", "llvm_ad_rgb"):
        if v in mi.variants():
            mi.set_variant(v)
            break
    cb = spec.camera_basis()
    H, W = cb["H"], cb["W"]
    to_world = mi.ScalarTransform4f().look_at(
        origin=mi.ScalarPoint3f(*cb["origin"]),
        target=mi.ScalarPoint3f(*cb["center"]),
        up=mi.ScalarVector3f(*cb["up"]))
    integ = ({"type": "stokes", "nested": {"type": "path", "max_depth": 4}}
             if "polarized" in mi.variant() else {"type": "path", "max_depth": 4})
    scene = mi.load_dict({
        "type": "scene",
        "integrator": integ,
        "sensor": {"type": "perspective", "fov_axis": "x",
                   "fov": float(2.0 * spec.hfov_deg), "to_world": to_world,
                   "film": {"type": "hdrfilm", "width": W, "height": H,
                            "pixel_format": "rgb", "rfilter": {"type": "box"}},
                   "sampler": {"type": "independent", "sample_count": spp}},
        "surface": {"type": "obj", "filename": obj,
                    "bsdf": {"type": "dielectric", "int_ior": spec.refractive_index(),
                             "ext_ior": 1.0, "specular_transmittance": 0.0}},
        "sky": {"type": "constant",
                "radiance": {"type": "rgb", "value": list(spec.sky_rgb)}},
    })
    arr = np.asarray(mi.render(scene, spp=spp), dtype=np.float64)
    # polarized stokes film: 15 channels -> S0 is group 1 (RGB)
    lin = arr.reshape(H, W, 5, 3)[:, :, 1, :] if arr.shape[-1] == 15 else arr[..., :3]
    np.save(os.path.join(out, "mitsuba_lin.npy"), lin)
    print(f"mitsuba color render done (spp={spp})")
    return lin


def make_figure(sp, mi, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping figure")
        return
    # absolute-radiance scale: single LS luminance factor over valid overlap
    valid = np.isfinite(sp[..., 0]) & np.isfinite(mi[..., 0]) & (mi[..., 0] >= 0)
    Lsp = sp @ LUM
    Lmi = mi @ LUM
    m = valid & (Lmi > 0)
    alpha = float(np.sum(Lsp[m] * Lmi[m]) / np.sum(Lmi[m] ** 2))
    mi_s = mi * alpha
    Lmi_s = Lmi * alpha
    rrmse = float(np.sqrt(np.mean((Lsp[m] - Lmi_s[m]) ** 2)) / np.mean(Lsp[m]))
    corr = float(np.corrcoef(Lsp[m], Lmi_s[m])[0, 1])
    print(f"absolute-radiance scale alpha={alpha:.3g}; after scaling: "
          f"luminance rel-RMSE={rrmse:.3f}, corr={corr:.3f} over {int(m.sum())} px")

    # shared exposure from seapol's 99th-pct luminance, then sRGB
    scale = np.quantile(Lsp[m], 0.99)
    sp_d = _srgb(np.nan_to_num(sp / max(scale, 1e-9)))
    mi_d = _srgb(np.nan_to_num(mi_s / max(scale, 1e-9)))
    diff = np.abs(Lsp - Lmi_s) / max(scale, 1e-9)
    diff = np.where(valid, diff, np.nan)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.4))
    ax[0].imshow(np.clip(sp_d, 0, 1)); ax[0].set_title("seapol S0 (color)")
    ax[1].imshow(np.clip(mi_d, 0, 1))
    ax[1].set_title(f"Mitsuba S0 (color, x{alpha:.2g})")
    im = ax[2].imshow(diff, cmap="magma", vmin=0, vmax=np.nanquantile(diff, 0.99))
    ax[2].set_title(f"|luminance diff|  (rel-RMSE {rrmse:.2f})")
    fig.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04)
    for a in ax:
        a.axis("off")
    fig.suptitle("Color-aware S0 off one synthetic wave field: seapol vs Mitsuba "
                 "(shared absolute scale)", fontsize=11)
    fig.tight_layout()
    path = os.path.join(out, "s0_color.png")
    fig.savefig(path, dpi=130)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/output_color")
    ap.add_argument("--spp", type=int, default=4096)
    args = ap.parse_args()
    out = args.out
    os.makedirs(out, exist_ok=True)

    spec = SceneSpec(L=8.0, N=512, U10=7.0, seed=0,
                     altitude_m=12.0, zenith_deg=40.0, azimuth_deg=0.0,
                     hfov_deg=18.0, img_h=384, img_w=384,
                     sun_zenith_deg=45.0, sun_azimuth_deg=90.0)
    eta, dx, info = generate_surface(spec)
    obj = os.path.join(out, "surface.obj")
    export_obj(eta, dx, obj)
    spec.to_json(os.path.join(out, "scene.json"))
    print(f"synthetic surface L={spec.L} N={spec.N} U10={spec.U10}; exported OBJ")

    sp = render_seapol_linear(spec, eta, dx, info, out)
    mi = render_mitsuba_linear(spec, obj, out, args.spp)
    if mi is not None:
        make_figure(sp, mi, out)


if __name__ == "__main__":
    main()
