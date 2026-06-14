"""
Color-aware rendering: nominal sky colors and water-type absorption
spectra through the spectral pipeline.

Five RGB scenes (450/550/650 nm bands, direct display mode):
    1. clear sky, high sun, clear ocean        -- deep blue water
    2. clear sky, high sun, productive Case 1  -- green-shifted water
    3. overcast, coastal Case 2                -- gray sky, green-brown
    4. partly cloudy, Case 1 with a bubble layer
    5. clear sky, low sun, clear ocean         -- warm glint, blue shadow
plus the quantitative panel: sub-surface reflectance spectra
E_u/E_down(lambda) per water type from 9-band scattering tables.

Per-band upwelling tables are cached under output/color_tables/ --
delete the directory to rebuild.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from seapol import (PinholeCamera, SpectralBands, SubpixelSlopes,
                    WATER_TYPES, build_spectral_tables, generate_sea_surface,
                    load_table, render_camera_image_spectral, save_table,
                    spectral_sky_factories, stokes_bands_to_rgb)
from seapol.water import WaterType

OUT = Path(__file__).parent / "output"
CACHE = OUT / "color_tables"

L, N, U10 = 32.0, 384, 6.5
BANDS = SpectralBands.rgb()
N_PHOT = 200_000


def cached_tables(tag, skies, water, bands, sun, turbidity):
    paths = [CACHE / f"{tag}_{int(wl)}.npz" for wl in bands.wavelengths_nm]
    if all(p.exists() for p in paths):
        print(f"  using cached tables: {tag}")
        return [load_table(p) for p in paths]
    print(f"  building tables: {tag} ({bands.n_bands} bands x "
          f"{N_PHOT} photons) ...")
    tabs = build_spectral_tables(skies, water, bands, sun=sun,
                                 turbidity=turbidity, n_photons=N_PHOT)
    CACHE.mkdir(parents=True, exist_ok=True)
    for p, t in zip(paths, tabs):
        save_table(p, t)
    return tabs


def render_scene(surf, sub, cam, *, tag, sky_kind, water_type,
                 sun_zen=40.0, sun_az=195.0, E_sun=30.0, turbidity=0.15,
                 cloud_kwargs=None):
    skies = spectral_sky_factories(
        BANDS, sky_kind, sun_zenith_deg=sun_zen, sun_azimuth_deg=sun_az,
        turbidity=turbidity, I_sun_550=0.0,
        **(cloud_kwargs or {}))
    sun = (sun_zen, sun_az, E_sun) if sky_kind != "overcast" else None
    tabs = cached_tables(tag, skies, water_type, BANDS, sun, turbidity)
    S = render_camera_image_spectral(
        surf.eta, surf.info["dx"], BANDS, skies, camera=cam, water=tabs,
        sun_glint=sun, turbidity=turbidity,
        slope_x=surf.slope_x, slope_y=surf.slope_y,
        subpixel=sub, n_subpixel=10, seed=5)
    return stokes_bands_to_rgb(S, BANDS.wavelengths_nm, mode="direct",
                               expose_quantile=0.995)


def reflectance_spectra():
    """E_u/E_down(lambda) per water type, 9 bands, isotropic sky."""
    from seapol import build_upwelling_table, make_unpolarized_sky
    bands9 = SpectralBands.visible(9)
    sky = make_unpolarized_sky(1.0)
    curves = {}
    for name, wt in WATER_TYPES.items():
        path = CACHE / f"spectrum_{name}.npz"
        if path.exists():
            z = np.load(path)
            curves[name] = z["R"]
            continue
        print(f"  reflectance spectrum: {name} ...")
        col = wt.column(bands9.wavelengths_nm)
        R = []
        for b in range(bands9.n_bands):
            t = build_upwelling_table(sky, col.at_band(b),
                                      n_photons=60_000, max_events=500,
                                      rng=np.random.default_rng(b))
            R.append(t.info["E_u"] / t.info["E_down"])
        CACHE.mkdir(parents=True, exist_ok=True)
        np.savez(path, R=np.array(R), wl=bands9.wavelengths_nm)
        curves[name] = np.array(R)
    return bands9.wavelengths_nm, curves


def main():
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(2)
    surf = generate_sea_surface(L=L, N=N, U10=U10, rng=rng)
    sub = SubpixelSlopes.from_cox_munk(U10, surf.info["k_cutoff"])
    cam = PinholeCamera(altitude_m=200.0, zenith_deg=42.0, azimuth_deg=0.0,
                        hfov_deg=5.0, img_shape=(300, 300))

    bubbly = WaterType(chlorophyll_mg_m3=0.3, bubble_scattering=1.5)
    scenes = [
        ("clear sky / clear ocean", dict(
            tag="clear_clear", sky_kind="clear",
            water_type=WATER_TYPES["clear"])),
        ("clear sky / productive Case 1", dict(
            tag="clear_productive", sky_kind="clear",
            water_type=WATER_TYPES["productive_case1"])),
        ("overcast / coastal Case 2", dict(
            tag="overcast_coastal", sky_kind="overcast",
            water_type=WATER_TYPES["coastal_case2"])),
        ("partly cloudy / Case 1 + bubbles", dict(
            tag="cloudy_bubbly", sky_kind="partly_cloudy",
            water_type=bubbly,
            cloud_kwargs=dict(cloud_fraction=0.45, cloud_brightness=4.0,
                              cloud_seed=8))),
        ("low sun / clear ocean", dict(
            tag="lowsun_clear", sky_kind="clear", sun_zen=78.0,
            E_sun=60.0, water_type=WATER_TYPES["clear"])),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 9.5), dpi=110)
    for ax, (title, kw) in zip(axes.ravel()[:5], scenes):
        print(f"scene: {title}")
        rgb = render_scene(surf, sub, cam, **kw)
        ax.imshow(rgb)
        ax.set_title(title, fontsize=10)
        ax.axis("off")

    ax = axes.ravel()[5]
    wl, curves = reflectance_spectra()
    for name, R in curves.items():
        ax.plot(wl, 100.0 * R, marker="o", ms=3, label=name)
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel(r"$E_u / E_d$ below surface [%]")
    ax.set_title("water-type reflectance spectra (MC)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle("Color-aware rendering: sky color spectra x water-type "
                 "absorption spectra (450/550/650 nm)", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "demo_color_scenes.png", bbox_inches="tight")
    print(f"wrote {OUT / 'demo_color_scenes.png'}")


if __name__ == "__main__":
    main()
