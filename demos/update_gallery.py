"""Sync the curated gallery figures from output/ (gitignored) into
gallery/ (tracked), so demos/README.md can embed them on GitHub.
Run after regenerating demos: uv run demos/update_gallery.py"""

import shutil
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "output"
GALLERY = HERE / "gallery"

FIGURES = [
    "demo_panel.png",
    "demo_time_evolution.png",
    "demo_polarimetric_reconstruction.png",
    "demo_polarimetric_polarized.png",
    "demo_glint_foam_current.png",
    "demo_sky_water_gallery.png",
    "demo_near_surface_scattering.png",
    "demo_color_scenes.png",
    "demo_stokes_panels.png",
    "demo_color_video_still.png",
    "demo_mc_reflectance.png",
    "demo_mc_water_body.png",
    "demo_fm98_capillaries.png",
    "demo_fm98_crest_stokes.png",
    "demo_fm98_crest_stokes_strips.png",
    "demo_fm98_3d_placement.png",
    "demo_slope_statistics.png",
    "demo_kw_spectrum.png",
    "demo_full_pipeline.png",
]


def main():
    GALLERY.mkdir(exist_ok=True)
    missing = []
    for name in FIGURES:
        src = OUT / name
        if src.exists():
            shutil.copy2(src, GALLERY / name)
            print(f"  synced {name}")
        else:
            missing.append(name)
    if missing:
        print("missing (run the demos first):", ", ".join(missing))


if __name__ == "__main__":
    main()
