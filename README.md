# surface-wave-light-reflection / seapol

Simulation of wind-driven sea surfaces and the polarized light reflected from the air-sea interface: the Elfouhaily et al. (1997) directional spectrum with FFT synthesis and time evolution, Mobley (2015) polarized reflection, a Rayleigh polarized sky, single-bounce renderers, near-surface scattering (an in-water polarized Monte Carlo that supplies the directional sub-surface light field, so scenes are water rather than bare mirror), a color-aware spectral mode (sky color spectra, water-type absorption spectra, CIE sRGB output), a forward Monte Carlo tracer with an in-water scattering medium, Fedorov & Melville (1998) phase-locked parasitic capillaries bound to their carriers and to the longer waves, polarimetric inversion back to slope and height, and an optional torch backend that runs the whole array pipeline on a GPU.

The Python package (`seapol`) is the primary implementation.  The original MATLAB codes that seeded the project live in [_original_MATLAB_codes/](_original_MATLAB_codes/).

## Install

```bash
uv sync                   # .venv: seapol (editable) + everything for tests/demos
uv run pytest             # test suite
uv run demos/run_all.py   # demo gallery (see demos/README.md)
```

The default `dev` dependency group carries pytest/matplotlib/netCDF4, so a plain `uv sync` (or any `uv run`) is fully provisioned; `uv.lock` pins the resolved environment.  Library consumers get the lean core (numpy + scipy + threadpoolctl) with optional extras `demos` (matplotlib), `asit` (netCDF4), and `gpu` (torch — a multi-GB download with bundled CUDA libraries); without uv: `pip install -e ".[demos,asit]" pytest`.

```bash
uv sync --extra gpu       # adds torch for the GPU backend
```

## Quick start

```python
import numpy as np
from seapol import (generate_sea_surface, make_rayleigh_sky, PinholeCamera,
                    SubpixelSlopes, render_camera_image)

# Time-evolving surface: eta(x, y, t) plus spectral slope fields
surf = generate_sea_surface(L=64.0, N=512, U10=7.0,
                            times=np.arange(32) / 4.0,
                            rng=np.random.default_rng(0))

# Polarization image of frame 0 under a Rayleigh sky
sky = make_rayleigh_sky(sun_zenith_deg=45.0, sun_azimuth_deg=90.0)
cam = PinholeCamera(altitude_m=300.0, zenith_deg=45.0, hfov_deg=3.0,
                    img_shape=(400, 400))
sub = SubpixelSlopes(surf.info["sigma_a2_cut"], surf.info["sigma_c2_cut"])
S = render_camera_image(surf.eta[:, :, 0], surf.info["dx"], camera=cam,
                        sky=sky, slope_x=surf.slope_x[:, :, 0],
                        slope_y=surf.slope_y[:, :, 0],
                        subpixel=sub, n_subpixel=16, shadowing=True)
I, Q, U, V = np.moveaxis(S, -1, 0)
```

Monte Carlo validation of the effective surface Mueller matrix:

```python
from seapol import effective_mueller_for_incident
res = effective_mueller_for_incident(surf.eta[:, :, 0], surf.info["dx"],
                                     theta_i_deg=40.0, n_rays=100_000)
print(res["R_total"], res["M_bin_mean"].shape)   # (9, 24, 4, 4)
```

## Modules

| Module | Purpose |
|---|---|
| `spectrum.py` | ECKV omnidirectional spectrum and spreading, fetch law, two drag laws (logistic Edson/Curcic-Haus fit and the ECKV form), sub-grid slope variances. |
| `surface.py` | FFT surface synthesis with exact variance closure, time evolution under gravity-capillary dispersion, free/bound spectral partition, uniform-current Doppler, spectral (exact) slope fields. |
| `polarization.py` | 4x4 Fresnel reflection/transmission Mueller matrices (incl. TIR phase), Stokes frame rotations, the full meridian-to-scattering-plane reflection chain. |
| `skylight.py` | Rayleigh single-scattering polarized sky (Coulson DoLP max), vector-based so U keeps its sign across the sun meridian. |
| `render.py` | Facet and pinhole-camera Stokes renderers; sub-pixel slope ensembles with optional Gram-Charlier (Cox-Munk) non-Gaussian statistics; analytic sun glint; whitecap foam; Smith/Saunders bistatic shadowing. |
| `montecarlo.py` | Forward polarized MC tracer: batch DDA grid march with true periodic boundaries, per-ray Mueller-matrix accumulation, and an optional in-water scattering medium. |
| `fm98.py` | Fedorov & Melville (1998) steady nonlinear gravity-capillary solver (spectral collocation + Newton continuation) and a gauge-fixed Stokes-coefficient table. |
| `hybrid.py` | Phase-locked hybrid surfaces: envelope demodulation injects FM98 bound harmonics (parasitic capillaries) onto the linear synthesis at the local carrier wavevector, with long-wave MTF binding, orbital advection, per-harmonic Nyquist truncation, and bin-wise variance redistribution. |
| `water.py` | First-order Case 1 / Case 2 water-leaving radiance; spectral IOPs from nominal bio-optics (Pope & Fry pure-water absorption, Bricaud-type chlorophyll, CDOM, sediment, near-surface bubble layers) via `WaterType`/`WaterColumn`; Quan-Fry n(lambda); Rayleigh, Henyey-Greenstein, and Fournier-Forand particulate phase Mueller matrices (with inverse-CDF samplers) and the shared polarized volume-scattering event. |
| `scattering.py` | Near-surface scattering: plane-parallel polarized MC for the water column under the actual sky (+ transmitted sun beam), Woodcock tracking through bubble-layer IOPs, producing the cached directional sub-surface Stokes radiance table the renderers consume as `water=`. |
| `spectral.py` | Color-aware orchestration: per-band skies/water/n(lambda) over the scalar renderers, (H, W, B, 4) spectral Stokes images, per-band scattering tables. |
| `color.py` | CIE 1931 spectral-to-sRGB conversion (and the direct 3-band mode) for spectral Stokes images. |
| `backend.py` | numpy/torch array dispatch: numpy by default, identical code paths on torch tensors (CPU or CUDA/ROCm), with cross-backend reproducible RNG. |
| `diagnostics.py` | Wavenumber-frequency diagnostics for evolving stacks: k-f slope spectra (Blackman-Harris, Welch segment averaging), inverse-phase-speed spectra Q(nu), slow-fraction helpers. |
| `empirical.py` | ASIT 2019 bridges: raw k-f cube reduction with noise-floor subtraction, measured bound-fraction curves beta(k) and the wind-interpolated library, and synthesis directly from measured directional slope spectra (`generate_asit_surface`). |
| `inversion.py` | Polarimetric slope sensing: reconstruct facet slopes from reflected Stokes imagery and wave height by spectral integration. `slopes_from_stokes` (DoLP -> Fresnel incidence angle, AoP -> plane of incidence) for unpolarized/overcast skies; `slopes_from_stokes_polarized` inverts the full Mueller chain against a known sky model (Gauss-Newton), valid under polarized skies at off-Brewster geometry (e.g. a DoFP camera at 30 deg). |

Demos in `demos/` write figures to `demos/output/`; see the **[illustrated gallery](demos/README.md)** for every figure, including the polarimetric slope/height reconstruction round trip. Tests: `uv run pytest`.

## Conventions

* Coordinates: x east, y north, z up; `eta[i, j]` has row `i` = y. Wind direction measured from +x, counter-clockwise.
* Stokes vectors are (I, Q, U, V) with trailing axis 4; propagation direction `d` carries the meridian frame `e_perp = unit(z x d)`, `e_par = e_perp x d`. **+Q is parallel to the meridian plane (vertical)**; reflected glint therefore has Q < 0.
* Mueller matrices act in the (p, s) scattering-plane basis and are rotated into/out of meridian frames explicitly.
* The upward normal of z = eta(x, y) is `(-deta/dx, -deta/dy, 1)`.

## Validation (seed 0, L = 1024 m, N = 1024)

| U10 [m/s] | Hs realized [m] | Hs target [m] | Hs PM [m] | MSS total | MSS Cox-Munk |
|---|---|---|---|---|---|
| 3 | 0.222 | 0.222 | 0.193 | 0.023 | 0.018 |
| 5 | 0.639 | 0.637 | 0.535 | 0.032 | 0.029 |
| 7 | 1.227 | 1.232 | 1.049 | 0.040 | 0.039 |
| 10 | 2.354 | 2.408 | 2.141 | 0.060 | 0.054 |
| 13 | 3.701 | 3.785 | 3.618 | 0.075 | 0.070 |

* Realized variance closes on the spectral target to <2%; Hs sits 5-15% above Pierson-Moskowitz (the ECKV long-wave branch is slightly more energetic); total MSS tracks Cox-Munk within ~10-20%.
* MC hemispherical reflectance matches flat Fresnel within ~3% at moderate angles and shows the correct rough-surface departure at grazing incidence (`demos/gallery/demo_mc_reflectance.png`).
* Flat-surface MC bin-mean Mueller equals the analytic Fresnel matrix to 1e-9 (exercises every frame rotation); R + T = 1 to machine precision.
* Polarimetric inversion round trip: slope fields, MSS, Hs, and the omnidirectional spectrum are recovered from rendered overcast-sky imagery to better than 1% (height correlation 0.9999).

## Hybrid FM98 surfaces (parasitic capillaries)

`generate_hybrid_surface` adds the nonlinear phase coherence that random-phase synthesis cannot supply: bound harmonics of short-gravity carriers, phase-locked to the forward face per the Fedorov-Melville Class-1 solution.  Design points:

1. **Physical envelope amplitudes.**  The slope map ak(r) = A(r) k(r) comes from a synthesis that closes variance on the spectral target; one-sided spreading keeps the whole carrier band in the downwind half-plane captured by the analytic signal.  Parasitic capillaries ignite only for ak >~ 0.2, so envelope amplitude fidelity controls whether trains appear at all.
2. **Deep harmonic tables.**  The ripple train of a 5 cm carrier lives at harmonic numbers m ~ 10 (lambda ~ 5 mm, the c(k_ripple) = c(carrier) resonance; the FM98 solver places the |a_m| bump at m ~ 11 for ak = 0.28).  Tables default to M_keep = 12 (28 for capillary-resolving grids) with per-harmonic Nyquist truncation.
3. **Bin-wise redistribution (not addition).**  The Elfouhaily spectrum is an empirical total — free plus bound variance, tuned to Cox-Munk — so the augmentation must never add power.  The bound-harmonic power is capped bin-by-bin at the available empirical power (bin phases kept, preserving the phase-locking) and exactly that much is removed from the random-phase field: total spectrum and MSS are preserved by construction (verified to 1-2% at U10 = 5-10 m/s, with emergent along-wind slope skewness -0.27 to -0.36 and excess kurtosis +0.5 to +1.5 tracking the Cox-Munk c03 wind trend -0.13 to -0.29 with no fitting).
4. **Local carrier placement (`carrier_k="local"`, default).**  The local carrier wavevector k_loc(r) = |grad phi| comes from the analytic-signal phase gradient, and the FM98 coefficients are bilinearly interpolated at the local (k_loc, ak) per pixel.  Each wave group then carries the capillary train of its own scale and direction: ripple wavelength follows the resonance k_rip ~ k_m^2 / k_loc (longer carriers bear *shorter* ripples), and ripple crests stay parallel to their local carrier.  For a monochromatic carrier this reduces exactly to a single-carrier evaluation at the true wavenumber (`demos/demo_fm98_3d_placement.py`).
5. **Well-forced branch selection.**  FM98 admit a range of wind forcings p at fixed (k, ak), parameterized by the forcing/profile phase offset delta = 4 nu k ak / (c0 p).  Weakly forced members put the parasitic train in the trough; `fm98.default_forcing(k, ak)` (delta ~ 0.02, calibrated to the published Fedorov, Melville & Rozenberg 1998 example p = 0.015 at lam = 7 cm, ak = 0.30) selects the branch whose train rides the forward face below the crest and decays through the trough.  At those exact published parameters the solver reproduces their figure 1: realized crest-trough ak 0.3000, elevation extrema +0.46/-0.21 cm (published +0.45/-0.22), and the ~15-ripple train with crest-front slope extrema of -1.0.  The continuation uses a smooth (log-sum-exp) amplitude constraint scaled to the Bernoulli rows, adaptive step-size marching, and two physical-state gates (single-crested profile, wave speed within 15% of linear) that keep the marching on the carrier-plus-train branch. Realized steepness reaches crest-trough ak = 0.40 for lam >= 7 cm; the lam = 5 cm branch folds near ak ~ 0.32.

```python
import numpy as np
from seapol import generate_hybrid_surface
from seapol.hybrid import default_fm98_table

table = default_fm98_table(k_nyq=np.pi / 0.0005, m_keep=20)  # one-time, cacheable
surf = generate_hybrid_surface(L=2.0, N=4096, U10=9.0, table=table,
                               rng=np.random.default_rng(0))
# surf.eta = surf.eta_lin + surf.eta_shrink + surf.eta_high
```

`demos/demo_fm98_capillaries.py` shows direct FM98 profiles (ripple train emerging with steepness at a 5 cm carrier) alongside the hybrid surface curvature field.  Resolving the ripples requires dx <~ 1 mm; the FM98 augmentation degrades gracefully (fewer harmonics) on coarser grids and raises if even m = 2 is unresolvable.

### Binding to the longer waves

`long_wave_modulation` (enabled via `generate_hybrid_surface(..., long_wave_mtf=...)`) adds the first stage of the binding hierarchy long waves -> short-gravity carriers -> parasitic capillaries.  The field is split at `k_split`; the short-wave part is amplitude-modulated with the standard hydrodynamic MTF form

    E'(r)/E = 1 + mtf * eps_L(r) * cos(phi_L(r) - theta),

where (A_L, phi_L) come from the long-wave analytic signal along the wind, eps_L = A_L k_L is the local long-wave steepness, mtf ~ 4-12 (empirical magnitude), and theta places the enhancement (0 = crest, +90 deg = mid forward face).  The FM98 envelope inherits the modulation, so the capillary packets concentrate on the chosen long-wave phase (~20x trough-to-face contrast at mtf = 6.5).  Because the modulation is recomputed from the instantaneous long-wave phase each frame, the bound content travels at the long-wave phase speed — the non-dispersive signature seen in inverse-phase-speed spectra Q(nu, theta). `apply_orbital_advection` adds the geometric coupling: the short-wave field is warped by the long-wave horizontal parcel displacement, broadening the free dispersion ridge as observed.

## Near-surface scattering (optional)

Specular reflection alone makes the scene a shiny mirror. `seapol.scattering` adds the in-water light: a plane-parallel polarized Monte Carlo propagates the actual sky (and optionally the transmitted solar beam — a different light path from the surface sun glint, so both can be on together) through the water column and tabulates the directional sub-surface upwelling Stokes radiance L_u(mu_w, phi).  The renderers accept the table as `water=`: each facet's view ray is inverse-refracted into the water, the table is sampled there, and the polarized transmission chain (with the n^2 radiance law) carries it to the camera.  Compared to the first-order isotropic `WaterBody` term, the scene gains view-dependent water brightness (anti-solar brightening of the scattered beam, grazing falloff toward the critical angle), water-type dependence, and physically placed DoLP damping.

```python
from seapol import (WATER_TYPES, build_upwelling_table, make_clear_sky,
                    render_camera_image, save_table, load_table)

sky = make_clear_sky(45.0, 195.0, I_sky=1.0, turbidity=0.15)
column = WATER_TYPES["coastal_case2"].column()       # 550 nm bulk IOPs
table = build_upwelling_table(sky, column, sun=(45.0, 195.0, 30.0),
                              n_photons=500_000)     # expensive: cache it
save_table("upwelling_coastal.npz", table)
S = render_camera_image(eta, dx, sky=sky, water=table,
                        sun_glint=(45.0, 195.0, 30.0), subpixel=sub)
```

The tracer is validated against the analytic single-scattering albedo of a semi-infinite Rayleigh medium, conservative-medium closure, and an exact energy budget (every photon ends absorbed, exited, or reported unresolved).  `WaterColumn` carries molecular + particulate (Henyey-Greenstein, Petzold-like g = 0.924, Voss-Fry polarization) + near-surface bubble-layer IOPs (exponential profile, Woodcock tracking).  Tables are horizontally homogeneous (flat mean interface for the column; facet tilts act at the exit refraction) and should be rebuilt when the sky, sun, or water type changes.

## Color-aware rendering (optional)

`seapol.spectral` runs the renderers per wavelength band with consistent physics: nominal sky color spectra (Rayleigh lambda^-4 blue whitened by an Angstrom-law aerosol haze, neutral overcast, Beer-reddened low sun), water-type absorption/scattering spectra (Pope & Fry pure water, Bricaud-type chlorophyll with the 675 nm peak, CDOM exponential, sediment), and Quan-Fry n(lambda) dispersion in both the Fresnel chain and the water-leaving transmission.  Band loops re-seed the sub-pixel ensemble and reuse one cloud field, so spectral noise stays luminance noise.  `seapol.color` maps spectral Stokes images to sRGB through the CIE 1931 color matching functions (or a direct three-band mode at 450/550/650 nm).

```python
from seapol import (SpectralBands, WATER_TYPES, build_spectral_tables,
                    render_camera_image_spectral, spectral_sky_factories,
                    stokes_bands_to_rgb)

bands = SpectralBands.rgb()                  # or SpectralBands.visible(13)
skies = spectral_sky_factories(bands, "clear", sun_zenith_deg=45.0,
                               sun_azimuth_deg=195.0, turbidity=0.15)
tables = build_spectral_tables(skies, WATER_TYPES["productive_case1"],
                               bands, sun=(45.0, 195.0, 30.0))  # cacheable
S = render_camera_image_spectral(eta, dx, bands, skies, water=tables,
                                 sun_glint=(45.0, 195.0, 30.0),
                                 subpixel=sub)        # (H, W, B, 4)
rgb = stokes_bands_to_rgb(S, bands.wavelengths_nm)    # (H, W, 3) sRGB
```

Cost scales linearly with bands; the per-band scattering tables are the expensive part (`demos/demo_color_scenes.py` caches them under `demos/output/color_tables/`).

## GPU backend (optional)

Every array-consuming stage dispatches on its inputs through `seapol.backend`: numpy arrays run plain numpy (the default path is unchanged), torch tensors run the same code on the tensors' device. Entry points that create arrays take `backend`/`device`/`dtype`, and everything downstream — hybrid FM98 augmentation, renderers, both Monte Carlos, diagnostics, inversion — stays on-device.

```python
import numpy as np
from seapol import backend, generate_hybrid_surface, render_camera_image

surf = generate_hybrid_surface(L=2.0, N=2048, U10=8.0, table=table,
                               rng=np.random.default_rng(0),
                               backend="torch", device="cuda")
S = render_camera_image(surf.eta, surf.info["dx"], camera=cam, sky=sky,
                        subpixel=sub, n_subpixel=64,
                        rng=backend.default_rng(1, backend.xp_of(surf.eta)))
S_np = backend.to_numpy(S)
```

Passing a numpy `Generator` into a torch-backed call draws on the CPU and copies — slower, but bitwise-reproducible across backends (the parity tests assert synthesis/render/MC equality at 1e-10). Device-native draws use `backend.default_rng(seed, xp)` / `backend.TorchRNG`.  On CUDA the default dtype is float32 (consumer GPUs run fp64 at a small fraction of fp32 throughput); CPU torch defaults to float64.  Not dispatched: the FM98 Newton solver and netCDF ingest (CPU setup work whose products convert on entry).

## Feature notes

* **Free/bound spectral partition.**  `generate_sea_surface(..., bound_fraction=beta, bound_speed=c)` splits each bin's power into a free part on the dispersion shell and a bound part advected rigidly downwind (bin powers, Hs, and MSS unchanged).  beta may be a scalar, a callable, or a measured curve via `empirical.bound_fraction_from_kf_reduction`; bound_speed may be a single speed, (speeds, weights), or `"spectrum"` (carrier speeds weighted by ring-integrated k^2 Psi).  With the measured beta(k), the synthetic slope k-f spectrum reproduces the measured bound ridges and Q(nu) decay (`demos/demo_kw_spectrum.py`: slow-nu slope fraction 0.49 free-only vs 0.082 hybrid, against 0.078 measured on the same support).
* **In-water polarized scattering.**  `WaterOptics(absorption, scattering, depolarization)` propagates transmitted MC rays through a semi-infinite water body: exponential free paths, absorption roulette, polarized Rayleigh-with-depolarization volume scattering with exact frame rotations.  Branch decisions and scattering azimuths are importance-sampled with the path's unpolarized-launch Stokes column, so per-ray I-weights stay exactly 1 and R + T + A closes to 1e-9 through hundreds of scattering events; the up-escaping field splits into glint and water-leaving parts.  Validated against the analytic single-scattering albedo of a semi-infinite Rayleigh medium and the exact single-scatter DoLP pattern (`demos/demo_mc_water_body.py`).
* **Non-Gaussian sun glint.**  `SubpixelSlopes` carries optional Gram-Charlier shape coefficients (skew(along-wind) = c03, cross term c21, kurtosis c40/c22/c04; `from_cox_munk(U10, k_cutoff)` loads the Cox-Munk 1954 clean-surface fits).  The analytic glint PDF and the sub-pixel ensemble (importance-weighted Gaussian draws) both honor them, reproducing the observed upwind/downwind glitter asymmetry. The correction is tapered beyond ~4 sigma, where the asymptotic GC series would otherwise produce spurious bright wings (`demos/demo_glint_foam_current.py`).
* **Whitecap foam.**  `render.Foam(coverage, albedo)` covers the steepest-slope facets (a breaking proxy, so foam rides the crests the synthesis sharpens) with unpolarized Lambertian radiance albedo E_d / pi; `monahan_coverage(U10)` gives the classical W = 3.84e-6 U10^3.41 fraction.
* **Uniform current.**  `generate_sea_surface(..., current=(Ux, Uy))` Doppler-shifts every component (free and bound) to omega + k.U, preserving Hermitian symmetry and per-frame variance exactly; relevant for k-f comparisons at tidal sites.
* **Streaming synthesis.**  `generate_sea_surface` and `generate_hybrid_surface` accept `frame_callback(it, t, eta, sx, sy)` so long records stream through the full FM98/MTF/advection pipeline without materializing stacks.
* **Wave-age keying.**  `empirical.run_conditions` (U10, fp, cp, inverse wave age per run) and `bound_fraction_for_conditions(..., inverse_wave_age=...)` rank bound-fraction library neighbors in (U10, Omega) space; binding strength rises with both wind and inverse wave age across 126 ASIT runs.

## Empirical (ASIT 2019) operation

With the ASIT data on disk, the pipeline runs end-to-end from measurements (`demos/demo_full_pipeline.py`):

* `load_asit_run` / `psi_from_asit` turn a run's measured directional slope spectrum S(k, theta) (Cartesian slope density per dkx dky, theta wind-relative) into Psi(kx, ky) for `generate_sea_surface(psi_override=...)`, with an Elfouhaily tail beyond the measured band.  Scope note: the measured spectrum's *form* and the k-f bound-wave structure are the calibration targets; instrument slope statistics (PDFs, MSS) are scale-limited and are shown for context only — slope statistics are built from the Elfouhaily + FM98 physics.
* **Instrument-matched synthesis.**  `psi_from_asit(...,  k_max=...)` (and `generate_asit_surface(..., k_max=...)`) band-limit the spectrum with a cosine taper to a camera's reliable resolution.  A sensor at FS frames/s cannot measure wave dynamics above FS/2 (~15 Hz for a 30 fps camera, ~100-200 rad/m by the gravity-capillary dispersion), so for *time-evolving* instrument-matched imagery the short waves should be band-limited there (the higher-k content is unresolved and, evolved forward, aliases) and carried as **bound** rather than free (`bound_fraction=...`, a simple monotone ramp), since cm-scale waves ride the longer waves instead of free-dispersing.  See `demos/demo_asit_stokes_video.py`.
* `reduce_kf_cube` streams a raw 6.9 GB Skw(f, kx, ky) cube into S(|k|, f) with the noise floor (estimated in the forbidden region of (k, f) space) subtracted; batch-reducing the ASIT runs builds a library of measured beta(k) curves (130 runs, U10 = 2.6-18 m/s: bound fraction rises monotonically with k toward ~1 at 100-200 rad/m and increases with wind), and `bound_fraction_for_wind` interpolates the library (monotone envelope, smoothed/tapered/capped) for synthesis.
* The renderers take `sun_glint=(zen, az, E_sun)`: the analytic Cox-Munk glint (sub-pixel slope PDF at the specular slope through the full Mueller chain), which avoids the glitter speckle of Monte Carlo sun-disk sampling.

## Notes relative to the original MATLAB codes

* The upward normal of z = eta(x, y) is `(-deta/dx, -deta/dy, 1)` (the MATLAB driver used `(+gx, +gy, 1)`).
* The Mueller retardance element is computed from the amplitude coefficients (rp rs), which is negative near normal incidence and correctly mirror-flips U and V on reflection (the tan/sin form is positive there).
* Sub-pixel slopes perturb the slope field (then renormalize the normal) and are drawn directly from the analytic distribution rather than a discretized PDF pool.
* The sky-weighting placeholder in the MATLAB driver is replaced by a polarized Rayleigh sky evaluated per (possibly perturbed) facet.
* All four Stokes components are carried throughout, including TIR phase shifts.

## Known simplifications and future work

* The water column is semi-infinite, horizontally homogeneous, and layered only through the bubble profile — no bottom, no stratification, no Raman/fluorescence.  The upwelling table assumes a flat mean interface (facet tilts act at the exit refraction) and a static sky: rebuild it when sun, sky, or water type changes. Particulate scattering offers a Henyey-Greenstein intensity or the Fournier-Forand phase function (`WaterColumn.particulate_phase="ff"`, tunable particle index and Junge slope -- realistic, decoupled backscatter that a single HG term cannot reach); the polarization is a depolarized-Rayleigh shape on either, not a fully measured Mueller matrix (Petzold/Voss-Fry tables remain a refinement candidate).
* The sky spectra are nominal single-scattering shapes (smoothed solar curve, lambda^-4 molecular + Angstrom aerosol, neutral clouds) for scene color, not an atmospheric RT solution; the polarization pattern is Rayleigh with an empirical DoLP cap, with no aerosol phase matrix or path radiance.
* Foam is a steepness-keyed Lambertian cover with a flat spectrum, not a breaking simulation; no wake physics.
* Long waves are linear (no Stokes/Creamer crest sharpening at the gravity-wave scale); the MTF is a scalar with fixed phase; dispersion is deep-water.
* The polarimetric inversion is single-reflection (no foam or water-leaving contamination).  `slopes_from_stokes` assumes an unpolarized sky on the sub-Brewster DoLP branch; `slopes_from_stokes_polarized` lifts the unpolarized assumption by inverting the full Mueller chain against a known sky, but still needs that sky model and a warm start to resolve the rare two-facet ambiguity.
* The torch backend leaves the FM98 Newton solver and netCDF ingest on the CPU; orbital advection uses Catmull-Rom (not spline-prefiltered) interpolation on torch, equal to the scipy path to interpolation order.
* Other candidates: division-of-focal-plane (DoFP) sensor model, per-run beta(k) fits binned by wind and wave age, MSS sensitivity to the drag law (the default logistic fit runs ~30% below Cox-Munk MSS; `drag_model="elfouhaily"` matches the published tuning).
