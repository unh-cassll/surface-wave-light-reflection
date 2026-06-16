"""
ASIT-matched Stokes videos: time-evolving sea surfaces driven by the
MEASURED ASIT 2019 directional slope spectra, rendered at the
instrument's geometry so the synthetic imagery is directly comparable
to the measured S0/S1/S2 frames.

For each ASIT Stokes run the measured directional slope spectrum
S(k, theta) (E-PSS stats file) is turned into Psi(kx, ky)
(empirical.psi_from_asit, Elfouhaily tail beyond the measured band),
a surface is synthesized and time-evolved under free dispersion at
30 fps, and rendered to a per-facet Stokes image under a clear
polarized sky at the instrument incidence (30 deg above nadir, which is
what gives the measured S1 its nonzero bias).

Geometry matches the camera: 2.903 m FOV (5.67 mm/px on the real 512px
sensor) synthesized at 2048 px (dx = 1.42 mm, 4x the sensor detail),
sub-pixel slope roughness from the Cox-Munk tail.  Panels: S0 (fixed
robust limits), S1/S0, S2/S0 (robust symmetric limits).

Frames are piped to ffmpeg (libx264).  Runs on the GPU when available.
Output: output/videos/stokes_video_asit_<run>.mp4 (atomic, resumable).
"""

import argparse
import os
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

import _video
from seapol import backend
from seapol.empirical import load_asit_run, psi_from_asit
from seapol.render import (CameraGeometry, SubpixelSlopes, make_clear_sky,
                           render_facet_stokes)
from seapol.surface import generate_sea_surface

OUT = Path(__file__).parent / "output"
VIDEOS = OUT / "videos"
EPSS = Path("/home/nathanlaxague/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = EPSS / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = EPSS / "ASIT2019_supporting_environmental_observations.nc"

N = 2048
FOV = 0.00567 * 512                 # 2.903 m, the instrument field of view
# Slow-motion: the simulation is sampled at SIM_FPS and the video plays
# at PLAY_FPS, so SIM_FPS/PLAY_FPS is the slow-mo factor.  Bound short
# waves advect at the carrier speed (~2 m/s), so at 30 fps real-time a
# few-cm wave would move ~1.5 wavelength/frame and alias; sampling at
# 120 Hz and playing at 30 fps (4x slow-mo) shows that motion coherently.
SIM_FPS = 120
PLAY_FPS = 30
SEED = 20240
S0_PCTL = (1.0, 99.0)
N_SUBPIXEL = 0
# Band-limit to the camera's reliable measurement: at 30 fps the data
# above 15 Hz -- ~100-200 rad/m by experience -- is guesswork, so the
# spectrum is tapered to ~150 rad/m (lambda ~ 4 cm) rather than trusting
# the higher-k tail.  Set to None for the full 2048px detail.
BAND_LIMIT_K_RADM = 150.0
# Bound short waves: above the dominant scales the short-wave variance
# rides the longer waves rather than free-dispersing (the seapol bound-
# wave / Q(nu,theta) picture).  A simple monotone ramp (no measured
# beta(k), which is unreliable past 15 Hz) makes the fraction BOUND_BETA
# advect with the carrier instead of on the free shell.
BOUND_BETA = 0.9
# camera 30 deg above nadir; azimuth offsets the view from the wind axis
CAMERA = CameraGeometry(incidence_deg=30.0, azimuth_deg=20.0, height_m=50.0)
SUN_ZEN, SUN_AZ, TURBIDITY = 45.0, 120.0, 0.12

# the six ASIT Stokes runs -> (label, stats/env run index)
RUNS = [("day025_U10.7", 121), ("day032_U10.1", 129), ("day040_U8.3", 143),
        ("day052_U4.6", 146), ("day059_U5.6", 174), ("day061_U8.6", 181)]

SEP, HEADER = 8, 64
_F = _video.font(30)
_FS = _video.font(22)


def _panels(S, s0_lo, s0_hi, q_lim, u_lim):
    I = S[..., 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(I > 1e-6, S[..., 1] / I, np.nan)
        u = np.where(I > 1e-6, S[..., 2] / I, np.nan)
    return (_video.to_u8(I, s0_lo, s0_hi, 0.0),
            _video.to_u8(q, -q_lim, q_lim, 0.5),
            _video.to_u8(u, -u_lim, u_lim, 0.5))


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
        _video.draw_centered(draw, i * (N + SEP) + N / 2.0, 18, c, _F)
    draw.text((12, 6), label, fill=(190, 190, 190), font=_FS)
    return np.asarray(img)


def _psi_array(run, U10, dx):
    """Measured-spectrum Psi(kx, ky) on the synthesis grid (numpy; the
    psi callable is numpy-only, so evaluate once and hand the array to
    the GPU synthesis).  Band-limited to the instrument resolution
    (k_max) so sub-cm waves that alias under 30 fps evolution are never
    synthesized."""
    rd = load_asit_run(str(STATS), run, env_path=str(ENV))
    psi = psi_from_asit(rd, U10=U10, k_max=BAND_LIMIT_K_RADM)
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, kx)
    return psi(KX, KY), rd


def render_run(label, run, xp, bk, dev, n_frames):
    final = VIDEOS / f"stokes_video_asit_{label}.mp4"
    if final.exists():
        print(f"{label:14s} exists, skipping")
        return
    dx = FOV / N
    rd = load_asit_run(str(STATS), run, env_path=str(ENV))
    U10 = float(rd["U10"])
    Psi, _ = _psi_array(run, U10, dx)
    sky = make_clear_sky(SUN_ZEN, SUN_AZ, 1.0, turbidity=TURBIDITY)
    sub = (SubpixelSlopes.from_cox_munk(U10, np.pi / dx)
           if N_SUBPIXEL > 0 else None)

    # bound short waves ride the carrier; the rest free-disperse
    bound_kw = dict(bound_fraction=BOUND_BETA) if BOUND_BETA else {}

    # fixed color limits from the t=0 frame
    s0 = generate_sea_surface(FOV, N, U10, psi_override=Psi,
                              rng=np.random.default_rng(SEED),
                              backend=bk, device=dev, **bound_kw)
    S = backend.to_numpy(render_facet_stokes(
        s0.eta, dx, CAMERA, sky=sky, slope_x=s0.slope_x,
        slope_y=s0.slope_y, subpixel=sub, n_subpixel=N_SUBPIXEL,
        rng=backend.default_rng(0, xp)))
    I0 = S[..., 0]
    s0_lo, s0_hi = np.nanpercentile(I0, S0_PCTL)
    q_lim = float(np.nanpercentile(np.abs(S[..., 1] / np.maximum(I0, 1e-6)),
                                   98))
    u_lim = float(np.nanpercentile(np.abs(S[..., 2] / np.maximum(I0, 1e-6)),
                                   98))
    del s0, S, I0

    tmp = final.with_suffix(".mp4.tmp")
    log = open(VIDEOS / f"stokes_video_asit_{label}.log", "w")
    proc = _video.open_ffmpeg(tmp, 3 * N + 2 * SEP, HEADER + N, PLAY_FPS,
                              log, crf=18)
    times = np.arange(n_frames) / SIM_FPS
    t0 = time.perf_counter()

    def cb(it, t, eta, sx, sy):
        S = backend.to_numpy(render_facet_stokes(
            eta, dx, CAMERA, sky=sky, slope_x=sx, slope_y=sy,
            subpixel=sub, n_subpixel=N_SUBPIXEL,
            rng=backend.default_rng(it + 1, xp)))
        slowmo = SIM_FPS // PLAY_FPS
        note = (f"measured spectrum, bound short waves, {slowmo}x slow-mo"
                if BOUND_BETA else "measured spectrum")
        frame = _compose(_panels(S, s0_lo, s0_hi, q_lim, u_lim),
                         f"ASIT {label}   U10 = {U10:.1f} m/s   "
                         f"t = {t:5.2f} s   ({note})")
        proc.stdin.write(frame.tobytes())

    generate_sea_surface(FOV, N, U10, psi_override=Psi, times=times,
                         compute_slopes=True, frame_callback=cb,
                         rng=np.random.default_rng(SEED),
                         backend=bk, device=dev, **bound_kw)
    proc.stdin.close()
    rc = proc.wait()
    log.close()
    dt = time.perf_counter() - t0
    if rc != 0:
        print(f"{label:14s} ffmpeg FAILED (rc={rc})")
        return
    os.replace(tmp, final)
    print(f"{label:14s} {n_frames} frames in {dt:.0f}s "
          f"({n_frames / dt:.2f} fps)  -> {final.name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--duration", type=float, default=30.0)
    ap.add_argument("--runs", nargs="+", default=None,
                    help="subset of run labels (default all six)")
    args = ap.parse_args()

    if not STATS.exists():
        print(f"missing ASIT stats file {STATS}")
        return
    VIDEOS.mkdir(parents=True, exist_ok=True)
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

    n_frames = int(round(args.duration * SIM_FPS))
    runs = RUNS if args.runs is None else [r for r in RUNS
                                           if r[0] in args.runs]
    print(f"{len(runs)} ASIT video(s), {N}px / {FOV:.3f} m "
          f"(dx={FOV / N * 1e3:.2f} mm), {n_frames} frames @ {SIM_FPS}fps")
    for label, run in runs:
        render_run(label, run, xp, bk, dev, n_frames)
    print(f"\ndone -> {VIDEOS}")


if __name__ == "__main__":
    main()
