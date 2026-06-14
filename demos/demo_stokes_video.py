"""
Stokes videos of the time-evolving Elfouhaily + FM98 capillary blend
across wind speed.

For each U10 in 3 .. 15 m/s (step 2) a 1024 x 1024 surface over a
1 m x 1 m patch is evolved at 30 fps for 30 s (900 frames) with the
streaming hybrid synthesis (linear Elfouhaily field plus the
phase-locked FM98 parasitic-capillary augmentation, re-derived every
frame so the bound capillaries ride their carriers), rendered to a
per-facet Stokes image under a clear polarized sky, and written as a
1 x 3 grayscale panel video:

    S0        reflected radiance (color limits fixed per video to a
              robust percentile range for visible dynamic range)
    S1 / S0   = Q/I, normalized polarization, color -1 .. 1
    S2 / S0   = U/I, normalized polarization, color -1 .. 1

Frames are composed as raw RGB and piped to ffmpeg (libx264), which is
far cheaper than per-frame matplotlib.  Runs on the GPU when torch +
CUDA are available.

NOTE on temporal sampling: at 30 fps the fastest FM98 capillaries are
undersampled (they oscillate well above the 15 Hz frame Nyquist), so
the finest ripple texture shimmers frame to frame -- each frame is a
valid snapshot, but the fast-wave motion aliases.  Pass --fnyq to mask
harmonics above the frame Nyquist for temporally clean (but
capillary-reduced) motion.

Output: output/videos/stokes_video_U{NN}.mp4 (written atomically; an
existing file is skipped, so the run is resumable).
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import _video
from seapol import backend
from seapol.fm98 import FM98Table
from seapol.hybrid import default_fm98_table, generate_hybrid_surface
from seapol.render import (CameraGeometry, make_clear_sky,
                           render_facet_stokes)

OUT = Path(__file__).parent / "output"
VIDEOS = OUT / "videos"
DEEP_CACHE = OUT / "fm98_table_deep.npz"
OWN_CACHE = OUT / "fm98_table_panels.npz"

L, N = 1.0, 1024
SIM_FPS = 30                # simulation sampling rate (wave time step = 1/SIM_FPS)
PLAY_FPS = 30               # video playback rate; SIM_FPS/PLAY_FPS = slow-mo factor
ALL_WINDS = [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
SEED = 20240
S0_PCTL = (1.0, 99.0)
SUN_ZEN, SUN_AZ = 50.0, 90.0
CAMERA = CameraGeometry(incidence_deg=40.0, azimuth_deg=0.0, height_m=50.0)

SEP = 8                     # px gap between panels
HEADER = 64                 # px label strip
_FONT = _video.font(30)
_FONT_S = _video.font(22)


def get_table(k_nyq):
    band = (78.5, 314.2)
    for cache in (DEEP_CACHE, OWN_CACHE):
        if cache.exists():
            t = FM98Table.load(cache)
            if t.k_grid[0] <= band[0] and t.k_grid[-1] >= band[1]:
                print(f"using cached FM98 table {cache.name} "
                      f"(M_keep={t.M_keep})")
                return t
    print("building FM98 table (one-time, cached) ...")
    t = default_fm98_table(k_nyq, m_keep=20, verbose=True)
    t.save(OWN_CACHE)
    return t


def _stokes_panels(S, s0_lo, s0_hi):
    I = S[..., 0]
    eps = 1e-6 * np.nanmax(np.abs(I))
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(np.abs(I) > eps, S[..., 1] / I, np.nan)
        u = np.where(np.abs(I) > eps, S[..., 2] / I, np.nan)
    return (_video.to_u8(I, s0_lo, s0_hi, 0.0),
            _video.to_u8(q, -1.0, 1.0, 0.5),
            _video.to_u8(u, -1.0, 1.0, 0.5))


def _compose(panels, label):
    h, w = N, 3 * N + 2 * SEP
    canvas = np.zeros((HEADER + h, w, 3), np.uint8)
    x = 0
    for p in panels:
        canvas[HEADER:HEADER + h, x:x + N, :] = p[..., None]
        x += N + SEP
    img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(img)
    for i, c in enumerate(["S0", "S1 / S0", "S2 / S0"]):
        _video.draw_centered(draw, i * (N + SEP) + N / 2.0, 18, c, _FONT)
    draw.text((12, 6), label, fill=(190, 190, 190), font=_FONT_S)
    return np.asarray(img)


def _tag():
    """Filename suffix encoding non-default domain / rates."""
    s = ""
    if abs(L - 1.0) > 1e-9:
        s += f"_L{int(round(L * 100)):02d}cm"
    if SIM_FPS != PLAY_FPS:
        s += f"_sim{SIM_FPS}_slowmo{SIM_FPS // PLAY_FPS}x"
    elif SIM_FPS != 30:
        s += f"_F{SIM_FPS}"
    return s


def render_video(U10, table, xp, bk, dev, n_frames, f_nyq):
    name = f"stokes_video_U{int(U10):02d}{_tag()}"
    final = VIDEOS / f"{name}.mp4"
    if final.exists():
        print(f"U10={U10:>4.0f}  exists, skipping  -> {final.name}")
        return
    sky = make_clear_sky(SUN_ZEN, SUN_AZ, I_sky=1.0, turbidity=0.1)

    # fixed S0 color limits from the t=0 frame
    s0 = generate_hybrid_surface(L, N, U10, table=table,
                                 rng=np.random.default_rng(SEED),
                                 backend=bk, device=dev)
    S0 = backend.to_numpy(render_facet_stokes(
        s0.eta, s0.info["dx"], CAMERA, sky=sky, slope_x=s0.slope_x,
        slope_y=s0.slope_y, subpixel=None, rng=backend.default_rng(0, xp)))
    s0_lo, s0_hi = np.nanpercentile(S0[..., 0], S0_PCTL)
    del s0, S0

    tmp = final.with_suffix(".mp4.tmp")
    log = open(VIDEOS / f"{name}.log", "w")
    h, w = HEADER + N, 3 * N + 2 * SEP
    proc = _video.open_ffmpeg(tmp, w, h, PLAY_FPS, log, crf=18)
    times = np.arange(n_frames) / SIM_FPS
    t0 = time.perf_counter()

    def cb(it, t, eta, sx, sy):
        S = backend.to_numpy(render_facet_stokes(
            eta, L / N, CAMERA, sky=sky, slope_x=sx, slope_y=sy,
            subpixel=None, rng=backend.default_rng(it + 1, xp)))
        panels = _stokes_panels(S, s0_lo, s0_hi)
        frame = _compose(panels, f"U10 = {U10:.0f} m/s    t = {t:5.2f} s")
        proc.stdin.write(frame.tobytes())

    generate_hybrid_surface(L, N, U10, table=table, times=times,
                            f_nyq=f_nyq, frame_callback=cb,
                            rng=np.random.default_rng(SEED),
                            backend=bk, device=dev)
    proc.stdin.close()
    rc = proc.wait()
    log.close()
    dt = time.perf_counter() - t0
    if rc != 0:
        print(f"U10={U10:>4.0f}  ffmpeg FAILED (rc={rc}); see log")
        return
    os.replace(tmp, final)
    slowmo = SIM_FPS / PLAY_FPS
    print(f"U10={U10:>4.0f}  {n_frames} frames in {dt:.0f}s "
          f"({n_frames / dt:.1f} fps render)  "
          f"{n_frames / SIM_FPS:.1f}s sim -> {n_frames / PLAY_FPS:.0f}s video "
          f"({slowmo:.0f}x slow-mo)  -> {final.name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--duration", type=float, default=30.0,
                    help="simulated seconds per video (default 30)")
    ap.add_argument("--winds", type=float, nargs="+", default=ALL_WINDS)
    ap.add_argument("--domain", type=float, default=1.0,
                    help="patch side length L [m] (default 1.0; "
                         "0.5 doubles spatial resolution)")
    ap.add_argument("--sim-fps", type=int, default=30,
                    help="simulation sampling rate [Hz] (default 30)")
    ap.add_argument("--play-fps", type=int, default=30,
                    help="video playback rate [fps] (default 30); "
                         "sim-fps/play-fps is the slow-motion factor")
    ap.add_argument("--fnyq", action="store_true",
                    help="mask FM98 harmonics above the simulation Nyquist "
                         "(temporally clean, capillary-reduced)")
    args = ap.parse_args()

    global L, SIM_FPS, PLAY_FPS
    L = args.domain
    SIM_FPS = args.sim_fps
    PLAY_FPS = args.play_fps

    VIDEOS.mkdir(parents=True, exist_ok=True)
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
    n_frames = int(round(args.duration * SIM_FPS))
    f_nyq = (SIM_FPS / 2.0) if args.fnyq else None
    print(f"{len(args.winds)} video(s), L={L} m (dx={L / N * 1e3:.3f} mm), "
          f"{n_frames} frames, sim {SIM_FPS} Hz -> play {PLAY_FPS} fps "
          f"({SIM_FPS / PLAY_FPS:.0f}x slow-mo), "
          f"capillaries: {'Nyquist-masked' if args.fnyq else 'full'}")
    for U10 in args.winds:
        render_video(U10, table, xp, bk, dev, n_frames, f_nyq)
    print(f"\ndone -> {VIDEOS}")


if __name__ == "__main__":
    main()
