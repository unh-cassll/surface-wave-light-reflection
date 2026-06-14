"""
Shared video helpers for the demo renderers (demo_stokes_video,
demo_asit_stokes_video, demo_color_scenes_video).

These centralize the two things that are easy to get wrong when piping
raw frames to ffmpeg:

  * the muxer flag -- the atomic temp path ends in `.tmp`, which ffmpeg
    cannot map to a container by extension, so `-f mp4` is mandatory
    (otherwise the writer dies and the frame pipe breaks); and
  * grayscale -> uint8 mapping with an explicit NaN fill (off-scene
    pixels are NaN and must not propagate into the encoder).

Frame composition (panel layout, labels) stays in each demo, since the
layouts differ.
"""

from __future__ import annotations

import subprocess

import matplotlib.font_manager as fm
import numpy as np
from PIL import ImageFont

_FONT_PATH = fm.findfont("DejaVu Sans")


def font(size: int) -> ImageFont.FreeTypeFont:
    """DejaVu Sans at the given pixel size."""
    return ImageFont.truetype(_FONT_PATH, size)


def to_u8(data, lo, hi, nan_fill: float = 0.0) -> np.ndarray:
    """Map data in [lo, hi] to uint8 [0, 255]; non-finite -> nan_fill
    (in normalized [0, 1] units) so off-scene NaNs render as a chosen
    gray instead of poisoning the encoder."""
    x = (np.asarray(data, dtype=float) - lo) / (hi - lo)
    x = np.where(np.isfinite(x), x, nan_fill)
    return (np.clip(x, 0.0, 1.0) * 255.0).astype(np.uint8)


def draw_centered(draw, cx: float, y: float, text: str, fnt, fill=(255,) * 3):
    """Draw text horizontally centered on cx (Pillow-version agnostic --
    uses textlength rather than the anchor kwarg)."""
    draw.text((cx - draw.textlength(text, font=fnt) / 2.0, y), text,
              fill=fill, font=fnt)


def open_ffmpeg(path, width: int, height: int, fps: int, log,
                crf: int = 18, preset: str = "medium"):
    """Open an ffmpeg subprocess consuming raw rgb24 frames on stdin and
    writing H.264/mp4 to `path`.

    `-f mp4` forces the muxer because the atomic temp path ends in
    `.tmp`; `+faststart` front-loads the moov atom for streaming.  Width
    and height must be even (yuv420p).  Returns the Popen; write frames
    with `proc.stdin.write(frame.tobytes())`, then `proc.stdin.close()`
    and `proc.wait()`."""
    if width % 2 or height % 2:
        raise ValueError(f"ffmpeg yuv420p needs even dims, got "
                         f"{width}x{height}")
    cmd = ["ffmpeg", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
           "-s", f"{width}x{height}", "-r", str(fps), "-i", "-", "-an",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(crf),
           "-preset", preset, "-movflags", "+faststart", "-f", "mp4",
           str(path)]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=log,
                            stderr=log)
