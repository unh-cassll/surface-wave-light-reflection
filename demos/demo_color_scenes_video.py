"""
Video form of demo_color_scenes: the five color-aware sky/water scenes
evolving over a SHARED, time-stepping wave field, so the only thing
differing between panels is the illumination and water optics.

A 2 x 3 grid of 30 fps RGB panels (450/550/650 nm bands -> sRGB):
    clear sky / clear ocean        clear sky / productive Case 1
    overcast / coastal Case 2      partly cloudy / Case 1 + bubbles
    low sun / clear ocean          (band reflectance-spectra legend)

One surface is synthesized and evolved under free dispersion; each frame
is rendered through every scene's per-band upwelling tables (cached from
demo_color_scenes; built here if missing).  Default 10 s.

Wind speed is set with --u10 (sets the surface roughness only; the sky/
water tables are wind-independent and shared), so a low- and a moderate-
high-wind pair can reuse the same cached tables, e.g.:
    demo_color_scenes_video.py --u10 3      # calm
    demo_color_scenes_video.py --u10 11     # rough

Frames are piped to ffmpeg (libx264).  Runs on the GPU when available.
Output: output/videos/color_scenes_video_U<u10>.mp4 (atomic; resumable).
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import _video
from seapol import backend
from seapol.render import PinholeCamera, SubpixelSlopes
from seapol.spectral import (SpectralBands, build_spectral_tables,
                             render_camera_image_spectral,
                             spectral_sky_factories)
from seapol.surface import generate_sea_surface
from seapol.color import stokes_bands_to_rgb
from seapol.water import WATER_TYPES, WaterType

OUT = Path(__file__).parent / "output"
VIDEOS = OUT / "videos"
CACHE = OUT / "color_tables"

L, N = 32.0, 512
DEFAULT_U10 = 6.5
PX = 400                            # per-panel render size (higher: see features)
FPS = 30
SEED = 2
BANDS = SpectralBands.rgb()
N_PHOT = 200_000
CAM = PinholeCamera(altitude_m=200.0, zenith_deg=42.0, azimuth_deg=0.0,
                    hfov_deg=5.0, img_shape=(PX, PX))

_BUBBLY = WaterType(chlorophyll_mg_m3=0.3, bubble_scattering=1.5)
# titles kept short so they fit one panel width without colliding
SCENES = [
    ("clear / clear ocean",
     dict(tag="clear_clear", sky_kind="clear",
          water=WATER_TYPES["clear"], sun_zen=40.0, E_sun=30.0)),
    ("clear / productive C1",
     dict(tag="clear_productive", sky_kind="clear",
          water=WATER_TYPES["productive_case1"], sun_zen=40.0, E_sun=30.0)),
    ("overcast / coastal C2",
     dict(tag="overcast_coastal", sky_kind="overcast",
          water=WATER_TYPES["coastal_case2"], sun_zen=40.0, E_sun=0.0)),
    ("cloudy / C1 + bubbles",
     dict(tag="cloudy_bubbly", sky_kind="partly_cloudy", water=_BUBBLY,
          sun_zen=40.0, E_sun=30.0,
          cloud=dict(cloud_fraction=0.45, cloud_brightness=4.0,
                     cloud_seed=8))),
    ("low sun / clear ocean",
     dict(tag="lowsun_clear", sky_kind="clear",
          water=WATER_TYPES["clear"], sun_zen=78.0, E_sun=60.0)),
]
SUN_AZ, TURB = 195.0, 0.15

SEP, HEADER = 2, 22                 # tight grid: minimal gutter and title bar
_FS = _video.font(16)


def _scene_setup(kw):
    skies = spectral_sky_factories(
        BANDS, kw["sky_kind"], sun_zenith_deg=kw["sun_zen"],
        sun_azimuth_deg=SUN_AZ, turbidity=TURB, **(kw.get("cloud") or {}))
    sun = ((kw["sun_zen"], SUN_AZ, kw["E_sun"]) if kw["E_sun"] > 0
           else None)
    paths = [CACHE / f"{kw['tag']}_{int(wl)}.npz"
             for wl in BANDS.wavelengths_nm]
    if all(p.exists() for p in paths):
        from seapol.scattering import load_table
        tabs = [load_table(p) for p in paths]
    else:
        print(f"  building tables: {kw['tag']} ...")
        tabs = build_spectral_tables(skies, kw["water"], BANDS, sun=sun,
                                     turbidity=TURB, n_photons=N_PHOT)
        CACHE.mkdir(parents=True, exist_ok=True)
        from seapol.scattering import save_table
        for p, t in zip(paths, tabs):
            save_table(p, t)
    return skies, tabs, sun


def _label(img_arr, titles):
    img = Image.fromarray(img_arr)
    draw = ImageDraw.Draw(img)
    cols = 3
    for i, t in enumerate(titles):
        r, c = divmod(i, cols)
        x = c * (PX + SEP)
        y = r * (HEADER + PX)
        draw.text((x + 6, y + 8), t, fill=(255, 255, 255), font=_FS)
    return np.asarray(img)




def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--duration", type=float, default=10.0)
    ap.add_argument("--u10", type=float, default=DEFAULT_U10,
                    help="wind speed (m/s); sets surface roughness")
    args = ap.parse_args()
    U10 = args.u10
    VIDEOS.mkdir(parents=True, exist_ok=True)
    final = VIDEOS / f"color_scenes_video_U{U10:g}.mp4"
    if final.exists():
        print(f"exists, skipping -> {final.name}")
        return

    if backend.has_torch():
        import torch
        bk, dev = ("torch", "cuda") if torch.cuda.is_available() \
            else ("torch", "cpu")
        print(f"backend: torch / "
              f"{torch.cuda.get_device_name(0) if dev == 'cuda' else 'cpu'}")
    else:
        bk, dev = None, None
        print("backend: numpy")
    xp = backend.get_xp(bk, dev)

    print("setting up scenes (tables cached from demo_color_scenes) ...")
    setups = [( _scene_setup(kw), kw) for _, kw in SCENES]

    # grid geometry: 2 rows x 3 cols, last cell blank
    cols, rows = 3, 2
    cw, ch = PX + SEP, HEADER + PX
    W, H = cols * cw - SEP, rows * ch
    if W % 2:
        W += 1
    if H % 2:
        H += 1
    titles = [t for t, _ in SCENES]

    log = open(VIDEOS / f"color_scenes_video_U{U10:g}.log", "w")
    tmp = final.with_suffix(".mp4.tmp")
    proc = _video.open_ffmpeg(tmp, W, H, FPS, log, crf=20)
    n_frames = int(round(args.duration * FPS))
    times = np.arange(n_frames) / FPS
    sub = SubpixelSlopes.from_cox_munk(U10, np.pi / (L / N))
    t0 = time.perf_counter()

    def cb(it, t, eta, sx, sy):
        canvas = np.zeros((H, W, 3), np.uint8)
        for i, ((skies, tabs, sun), kw) in enumerate(setups):
            S = render_camera_image_spectral(
                eta, L / N, BANDS, skies, camera=CAM, water=tabs,
                sun_glint=sun, turbidity=TURB, slope_x=sx, slope_y=sy,
                subpixel=sub, n_subpixel=1, seed=SEED)
            rgb = (np.clip(stokes_bands_to_rgb(
                S, BANDS.wavelengths_nm, mode="direct",
                expose_quantile=0.995), 0, 1) * 255).astype(np.uint8)
            r, c = divmod(i, cols)
            y0 = r * ch + HEADER
            x0 = c * cw
            canvas[y0:y0 + PX, x0:x0 + PX, :] = rgb
        frame = _label(canvas, titles)
        proc.stdin.write(frame.tobytes())

    generate_sea_surface(L, N, U10, times=times, compute_slopes=True,
                         frame_callback=cb, rng=np.random.default_rng(SEED),
                         backend=bk, device=dev)
    proc.stdin.close()
    rc = proc.wait()
    log.close()
    dt = time.perf_counter() - t0
    if rc != 0:
        print(f"ffmpeg FAILED (rc={rc}); see log")
        return
    os.replace(tmp, final)
    print(f"{n_frames} frames in {dt:.0f}s ({n_frames / dt:.2f} fps)  "
          f"-> {final.name}")


if __name__ == "__main__":
    main()
