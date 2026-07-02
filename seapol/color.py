"""
Spectral-to-color conversion for the color-aware rendering mode.

Spectral Stokes images (H, W, B, 4) carry per-band radiance; this module
maps the intensity bands to displayable sRGB:

    * "cie"    : trapezoid integration of the band radiances against the
                 CIE 1931 2-degree color matching functions ->
                 XYZ -> linear sRGB (D65 primaries) -> gamma encoding.
                 Quantitative when the band grid resolves the visible
                 (>= ~5 bands).
    * "direct" : three bands taken directly as (B, G, R) channels --
                 the nominal fast mode for wavelengths ~ (450, 550, 650).

Exposure: radiances are relative (sky models are normalized at 550 nm),
so a display scale must be chosen; auto-exposure maps a high quantile
of luminance to white.
"""

from __future__ import annotations

import numpy as np

from .backend import to_numpy

__all__ = ["cie_cmf", "stokes_bands_to_rgb", "srgb_encode",
           "xyz_to_linear_srgb"]

# CIE 1931 2-degree color matching functions, 10 nm, 380-730 nm
_CMF_WL = np.arange(380.0, 731.0, 10.0)
_CMF = np.array([
    [0.0014, 0.0000, 0.0065], [0.0042, 0.0001, 0.0201],
    [0.0143, 0.0004, 0.0679], [0.0435, 0.0012, 0.2074],
    [0.1344, 0.0040, 0.6456], [0.2839, 0.0116, 1.3856],
    [0.3483, 0.0230, 1.7471], [0.3362, 0.0380, 1.7721],
    [0.2908, 0.0600, 1.6692], [0.1954, 0.0910, 1.2876],
    [0.0956, 0.1390, 0.8130], [0.0320, 0.2080, 0.4652],
    [0.0049, 0.3230, 0.2720], [0.0093, 0.5030, 0.1582],
    [0.0633, 0.7100, 0.0782], [0.1655, 0.8620, 0.0422],
    [0.2904, 0.9540, 0.0203], [0.4334, 0.9950, 0.0087],
    [0.5945, 0.9950, 0.0039], [0.7621, 0.9520, 0.0021],
    [0.9163, 0.8700, 0.0017], [1.0263, 0.7570, 0.0011],
    [1.0622, 0.6310, 0.0008], [1.0026, 0.5030, 0.0003],
    [0.8544, 0.3810, 0.0002], [0.6424, 0.2650, 0.0000],
    [0.4479, 0.1750, 0.0000], [0.2835, 0.1070, 0.0000],
    [0.1649, 0.0610, 0.0000], [0.0874, 0.0320, 0.0000],
    [0.0468, 0.0170, 0.0000], [0.0227, 0.0082, 0.0000],
    [0.0114, 0.0041, 0.0000], [0.0058, 0.0021, 0.0000],
    [0.0029, 0.0010, 0.0000], [0.0014, 0.0005, 0.0000]])

# linear sRGB (D65) from XYZ
_XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986],
                         [-0.9689, 1.8758, 0.0415],
                         [0.0557, -0.2040, 1.0570]])


def cie_cmf(wavelength_nm):
    """CIE 1931 2-deg color matching functions (..., 3) at the given
    wavelengths [nm]."""
    wl = np.asarray(wavelength_nm, dtype=float)
    return np.stack([np.interp(wl, _CMF_WL, _CMF[:, i]) for i in range(3)],
                    axis=-1)


def xyz_to_linear_srgb(xyz):
    """XYZ (..., 3) -> linear sRGB (..., 3), unclipped."""
    return np.einsum("ij,...j->...i", _XYZ_TO_SRGB, np.asarray(xyz))


def srgb_encode(linear):
    """Linear -> gamma-encoded sRGB in [0, 1]."""
    c = np.clip(np.asarray(linear, dtype=float), 0.0, 1.0)
    return np.where(c <= 0.0031308, 12.92 * c,
                    1.055 * c ** (1.0 / 2.4) - 0.055)


def stokes_bands_to_rgb(S_bands, wavelengths_nm, stokes_index: int = 0,
                        mode: str | None = None,
                        exposure: float | None = None,
                        expose_quantile: float = 0.99,
                        return_exposure: bool = False):
    """Displayable sRGB image (H, W, 3) from a spectral Stokes image
    (H, W, B, 4) (or any leading shape + (B, 4)).

    stokes_index selects the channel (0 = intensity; DoLP images should
    be colormapped instead).  mode "cie" (>= 5 bands default) or
    "direct" (3 bands ~ 450/550/650 as B/G/R).  exposure scales linear
    RGB before gamma; None auto-exposes so the expose_quantile of
    luminance maps to 1.  return_exposure=True returns (image, exposure)
    so video loops can freeze the frame-0 gain.  NaN pixels (off-scene)
    come back black."""
    S = to_numpy(S_bands)
    wl = np.asarray(wavelengths_nm, dtype=float)
    band = S[..., stokes_index]
    if mode is None:
        mode = "cie" if wl.size >= 5 else "direct"

    if mode == "cie":
        cmf = cie_cmf(wl)                              # (B, 3)
        xyz = np.trapezoid(band[..., None] * cmf, wl, axis=-2)
        rgb = xyz_to_linear_srgb(xyz / max(np.trapezoid(cmf[:, 1], wl),
                                           1e-12))
    elif mode == "direct":
        if wl.size != 3:
            raise ValueError("direct mode needs exactly 3 bands")
        order = np.argsort(wl)[::-1]                   # R, G, B
        rgb = band[..., order]
    else:
        raise ValueError(f"unknown mode: {mode!r}")

    finite = np.isfinite(rgb).all(axis=-1)
    rgb = np.where(finite[..., None], rgb, 0.0)
    rgb = np.maximum(rgb, 0.0)
    if exposure is None:
        lum = rgb @ np.array([0.2126, 0.7152, 0.0722])
        ref = np.quantile(lum[finite], expose_quantile) if finite.any() \
            else 1.0
        exposure = 1.0 / max(ref, 1e-12)
    img = srgb_encode(rgb * exposure)
    return (img, exposure) if return_exposure else img
