# Cross-validation harness

Compare `seapol` renders against third-party renderers, to confirm the
polarized reflection and glint geometry agree with independent implementations.

| Renderer | Validates | Notes |
|---|---|---|
| **seapol** | reference | `render_seapol.py` — the baseline everything is compared to |
| **Blender (Cycles)** | radiance `I` only | no polarization; sky model differs (uniform background used) |
| **Mitsuba 3 (polarized)** | full Stokes / DoLP / AoP | the rigorous polarimetric cross-check |

Unity is intentionally excluded: as a real-time rasterizer it produces neither
spectral nor polarized radiance and is not a quantitative ground truth.

## One source of truth

All renderers read a single [`SceneSpec`](scene.py) JSON (surface, camera, sun,
water refractive index). The camera basis reproduces `seapol.render._camera_rays`
and the exported OBJ reuses the `seapol.montecarlo.build_cell_tables`
triangulation, so an external render of the exported mesh is geometrically
registered to the seapol render of the same surface.

## Run

```bash
# 1. seapol baseline (writes seapol.npz, surface.obj, scene.json)
uv run validation/scene.py validation/output/scene.json     # optional: edit it
uv run validation/render_seapol.py validation/output/scene.json

# 2a. Blender intensity render (needs a `blender` binary on PATH)
blender -b --python validation/render_blender.py -- \
    validation/output/scene.json validation/output/surface.obj \
    validation/output/blender.exr

# 2b. Mitsuba polarized render (needs `pip install mitsuba`)
uv run validation/render_mitsuba.py \
    validation/output/scene.json validation/output/surface.obj

# 3. metrics + side-by-side figure (skips any absent external render)
uv run validation/compare.py
```

`compare.py` runs with neither Blender nor Mitsuba installed — it reports the
seapol baseline and skips the missing stages with a message.

## Color-aware S0: seapol vs Mitsuba

`s0_color.py` renders one synthetic wave field's reflected-sky intensity (S0) in
**color** with seapol and Mitsuba and compares them on a **shared absolute
radiance scale**:

```bash
uv run validation/s0_color.py --out validation/output_color --spp 4096
```

- **seapol** — spectral pipeline (`SpectralBands.rgb`), reflection-only.
- **Mitsuba** — `*_polarized`, reflection-only dielectric, colored constant sky.

Both compute the same physics (dielectric Fresnel reflection of the sky), so
after a single least-squares luminance scale factor `alpha` they agree closely
(latest run: luminance correlation ≈0.7, rel-RMSE ≈0.27). The figure is
seapol | Mitsuba(scaled) | luminance-difference. Blender is intentionally not in
this figure — it is unpolarized and uses a different sky/dielectric model; use
`render_blender.py` for a separate qualitative intensity-only view.

**Why high spp.** The reflection-only dielectric suppresses the transmission
lobe, so only ~R (a few percent at these angles) of the BSDF samples carry the
reflected signal; the rest sample the killed transmission lobe and contribute
zero. Low spp therefore speckles badly (relative noise ≈ √((1−R)/(R·N))), worst
in DoLP, which divides by a noisy small S0. High spp removes the speckle and, as
a bonus, averages sub-pixel facets the way seapol's sub-pixel slope sampling
does — which is why the seapol↔Mitsuba agreement improves with spp.

**Absolute radiance scaling.** Each renderer carries its own radiometric units
(seapol's spectral sky, Mitsuba's emitter radiance, RGB conversion). A single
multiplicative `alpha`, fit by least squares on luminance over the overlapping
valid pixels, puts them on a common absolute scale, so the comparison is about
spatial/chromatic structure rather than arbitrary exposure.

Two registration details the Blender path gets right (easy to miss): the OBJ is
imported with `up_axis="Z", forward_axis="Y"` (Blender's default Y-up conversion
would rotate the patch into a sliver), and the mesh normals are forced to +z (an
inverted normal makes Blender's Fresnel/Glass treat the surface as seen from
inside the water → spurious total reflection).

## Demonstrated result

With Mitsuba 3.8 (polarized) on a flat, unpolarized-sky scene, the reflected
**degree of linear polarization matches seapol to a mean error of 0.001**
(seapol 0.577 vs Mitsuba 0.577) — an independent confirmation of seapol's
polarized Fresnel reflection.

The **angle of polarization also matches** once both are in a common frame.
seapol reports AoP in its meridian Stokes frame; Mitsuba in its camera frame, so
they differ by a *per-pixel* rotation (the meridian angle in the camera frame) —
not a single global offset, which is why a constant shift left ~21° residual.
`aop_align.py` computes that per-pixel meridian angle from the camera geometry
and rotates Mitsuba's AoP into seapol's frame: the error drops from **81.4° to
0.5°** (sign −1, zero constant — i.e. exactly the geometric meridian rotation):

```bash
uv run validation/aop_align.py --out validation/out_flat
```

`compare.py` applies the same rotation, so its AoP panel and metric are reported
in the meridian frame. The Blender intensity render is a qualitative
glint/radiance visualization only (Blender is unpolarized).

## Flat-surface sanity check

Set both `"flat": true` and `"unpolarized_sky": true` in the scene JSON. The
unpolarized sky matches Mitsuba's constant environment, so the reflected
DoLP/AoP are generated purely by surface Fresnel in both codes and should agree
to within a few percent — the cleanest cross-code check. (With the default
polarized Rayleigh sky the incident field differs from Mitsuba's, so DoLP would
not match; keep them aligned for the sanity test.)

## Install

```bash
uv sync --extra validation     # adds mitsuba + imageio (EXR reading)
```

Blender is a standalone application, not a Python package; install it separately
and ensure `blender` is on `PATH`.

## Caveats

- **Sky model.** seapol uses a Coulson Rayleigh polarized sky; Blender uses a
  uniform background here and Mitsuba an unpolarized constant. Comparisons are
  about the *surface* (Fresnel reflection, glint geometry), not the sky. A
  matching polarized environment emitter for Mitsuba (precomputed from
  `rayleigh_sky_stokes`) is a documented extension, not yet implemented.
- **Mesh wrap.** The OBJ is the open (non-wrapped) patch; seapol's periodic
  wrap is an internal Monte Carlo detail and is omitted so external renderers
  see no domain-spanning facet.
