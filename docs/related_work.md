# seapol in context: constituent elements and prior work

`seapol` is not a new physical theory of air-sea light interaction.  It is an
**integration of independently validated models**, each drawn from the published
literature, assembled into a single render-ready (and invertible) pipeline.  This
document decomposes the package into its constituent elements and names the source each
one implements, so that the physics can be checked component by component against the
expectation set by prior work.

Read this way, almost every block in `seapol` is a textbook or peer-reviewed model.  The
genuinely new contributions are narrow and specific:

1. **Phase-locked FM98 parasitic-capillary surfaces** integrated into a render-ready
   synthesis — the Fedorov & Melville (1998) nonlinear gravity-capillary solution injected
   onto a linear random-phase field by envelope demodulation, conserving spectral variance
   ([seapol/fm98.py](../seapol/fm98.py), [seapol/hybrid.py](../seapol/hybrid.py)).
2. **A closed forward → inverse loop** — polarimetric retrieval of facet slope and wave
   height from the same Stokes imagery the forward model renders
   ([seapol/inversion.py](../seapol/inversion.py)).

Everything else is a published model, used as published.

---

## 1. Sea-surface wave field

| Element | Source | Implementation |
|---|---|---|
| Directional wave spectrum (ECKV unified gravity + capillary) | Elfouhaily, Chapron, Katsaros & Vandemark (1997) | [spectrum.py](../seapol/spectrum.py) `elfouhaily_omni`, `directional_spectrum` |
| Wind-stress / drag laws | Edson et al. (2013); Curcic & Haus (2020); ECKV (1997) form | [spectrum.py](../seapol/spectrum.py) `drag_coefficient` |
| FFT random-phase synthesis, exact variance closure, gravity-capillary dispersion | Standard linear wave theory | [surface.py](../seapol/surface.py) `generate_sea_surface` |
| Gaussian slope statistics | Cox & Munk (1954) | [render.py](../seapol/render.py), [spectrum.py](../seapol/spectrum.py) `cutoff_slope_variances` |
| Non-Gaussian (Gram-Charlier) slope skew / kurtosis | Cox & Munk (1954, 1956) | [render.py](../seapol/render.py) `SubpixelSlopes` |
| Nonlinear bound / parasitic-capillary harmonics | Fedorov & Melville (1998); phase-locking after Melville & Fedorov (2015) | [fm98.py](../seapol/fm98.py), [hybrid.py](../seapol/hybrid.py) |
| Empirical directional spectra & bound-fraction library | ASIT field measurements | [empirical.py](../seapol/empirical.py) |

## 2. Polarized surface reflection / refraction

| Element | Source | Implementation |
|---|---|---|
| Stokes-vector / Mueller-matrix conventions | Mobley (2015) | [polarization.py](../seapol/polarization.py) |
| Polarized Fresnel coefficients + TIR retardance | Classical electromagnetics (Born & Wolf) | [polarization.py](../seapol/polarization.py) `fresnel_mueller` |
| Facet (Cox-Munk) Mueller reflection chain | Cox & Munk (1954) + Mobley (2015) frame algebra | [polarization.py](../seapol/polarization.py) `reflection_chain`, [render.py](../seapol/render.py) |

## 3. Sky / incident light field

| Element | Source | Implementation |
|---|---|---|
| Rayleigh single-scattering polarized sky (Coulson DoLP max) | Coulson (1988); Xue et al. (2021) | [skylight.py](../seapol/skylight.py) `rayleigh_sky_stokes` |
| Rayleigh optical depth | Bodhaine et al. (1999) | [skylight.py](../seapol/skylight.py) `rayleigh_optical_depth` |
| Aerosol optical depth (Angstrom turbidity) | Angstrom (1929) | [skylight.py](../seapol/skylight.py) `aerosol_optical_depth` |
| Overcast luminance gradation | Moon & Spencer (1942) | [skylight.py](../seapol/skylight.py) `overcast_sky_stokes` |

## 4. Water-body optics

| Element | Source | Implementation |
|---|---|---|
| Seawater refractive index n(S, T, lambda) | Quan & Fry (1995) | [water.py](../seapol/water.py) `water_refractive_index` |
| Pure-water absorption | Pope & Fry (1997) | [water.py](../seapol/water.py) `pure_water_absorption` |
| Seawater Rayleigh scattering + depolarization | Morel (1974) | [water.py](../seapol/water.py) `pure_water_scattering` |
| Phytoplankton (chlorophyll) absorption | Bricaud et al. (Case 1) | [water.py](../seapol/water.py) `WaterType` |
| Particulate scattering scaling | Loisel & Morel (1998) | [water.py](../seapol/water.py) `WaterType` |
| Rayleigh phase Mueller matrix | Hansen & Travis (1974) | [water.py](../seapol/water.py) `rayleigh_phase_mueller` |
| Particulate phase functions (HG, Fournier-Forand; Petzold g) | Henyey & Greenstein (1941); Fournier & Forand (1994); Petzold (1972) | [water.py](../seapol/water.py) `hg_phase_mueller`, `fournier_forand_phase` |
| Particulate DoLP cap | Voss & Fry (1984) | [water.py](../seapol/water.py) |

## 5. Radiative transfer / Monte Carlo

| Element | Source | Implementation |
|---|---|---|
| Polarized MC with Stokes importance sampling | Mobley (2015) | [montecarlo.py](../seapol/montecarlo.py) `trace_forward`, [scattering.py](../seapol/scattering.py) `build_upwelling_table` |
| Null-collision (Woodcock) tracking | Woodcock et al. (1965) | [scattering.py](../seapol/scattering.py) |
| DDA voxel/grid traversal | Amanatides & Woo (1987) | [montecarlo.py](../seapol/montecarlo.py) `heightmap_intersect` |
| Ray-triangle intersection | Moller & Trumbore (1997) | [montecarlo.py](../seapol/montecarlo.py) `_ray_tri` |

## 6. Scene / rendering

| Element | Source | Implementation |
|---|---|---|
| Whitecap / foam coverage | Monahan & O'Muircheartaigh (1980) | [render.py](../seapol/render.py) `monahan_coverage` |
| Surface shadowing (bistatic) | Smith (1967); Saunders | [render.py](../seapol/render.py) `smith_lambda`, `saunders_bistatic` |
| Spectral-to-color conversion | CIE 1931 CMFs; IEC sRGB | [color.py](../seapol/color.py), [spectral.py](../seapol/spectral.py) |

## 7. Polarimetric inversion (epss)

| Element | Source | Implementation |
|---|---|---|
| DoLP -> incidence, AoP -> plane-of-incidence retrieval | Polarized Fresnel inversion (sub-Brewster branch) | [inversion.py](../seapol/inversion.py) `slopes_from_stokes` |
| Mueller-chain slope fit | Levenberg-Marquardt over the forward model | [inversion.py](../seapol/inversion.py) `slopes_from_stokes_polarized` |
| Slope -> height spectral integration | Fourier slope-to-elevation | [inversion.py](../seapol/inversion.py) `height_from_slopes` |

## 8. Compute

| Element | Source | Implementation |
|---|---|---|
| numpy / torch array dispatch (CPU + GPU) | — | [backend.py](../seapol/backend.py) |

---

## Integration-level precedents

The five works below each combine a subset of the elements above in essentially the way
`seapol` does, establishing that the combination matches accepted practice.  `seapol` is a
superset of any one of them rather than a departure.

| Work | Element subset combined | What it establishes |
|---|---|---|
| Hieronymi (2016) | 1-2 | Cox-Munk facets + Mueller Fresnel -> polarized reflectance/transmittance distribution functions (pBRDF/pBTDF). Validated reference tables. |
| You et al. (2011) | 1, 4-5 | Dynamic surface + polarized in-water RT, compared with underwater measurements. The closest analog to seapol's surface + in-water MC path. |
| D'Alimonte & Kajiyama (2016) | 1-2 | Polarization x (non-)Gaussian slope statistics -> sea-surface reflectance factor for ocean-color glint correction. |
| Xue et al. (2021) | 2-3, 6 | Polarized sky + facet reflection + airborne camera imaging chain (submarine Kelvin wakes). |
| Chu et al. (2020) | 2-3 | Rayleigh sky + Fresnel facets in the skylight-dominated (low-glint) regime. |

Where these are individually deeper than `seapol`: Hieronymi's tables are a validated
reference product (with whitecap and multiple-surface contributions); You et al. couple a
spatially-resolved underwater field to the explicit dynamic surface and validate against
in-water data, whereas `seapol`'s near-surface scattering ([scattering.py](../seapol/scattering.py))
is plane-parallel (facet tilt applied only at exit refraction); D'Alimonte & Kajiyama
report a focused, validated reflectance factor.  `seapol` trades some of that
specialty depth for breadth, the nonlinear surface microstructure, spectral color, and the
inversion loop.

---

## Appendix: capability matrix

| Capability | Hieronymi 2016 | You 2011 | D'Alimonte 2016 | Xue 2021 | Chu 2020 | seapol |
|---|---|---|---|---|---|---|
| 4x4 Mueller / polarized Fresnel | yes | yes | yes | yes | yes | yes |
| Explicit phase-resolved surface | - | yes | - | yes (wake) | yes | yes |
| Dynamic / time-evolving + currents | - | yes | - | - | - | yes |
| Nonlinear bound waves (FM98 parasitic capillaries) | - | - | - | - | - | **yes** |
| Non-Gaussian (Gram-Charlier) slopes | - | - | yes | - | - | yes |
| In-water polarized RT/MC | tables | yes | - | - | - | yes (plane-parallel) |
| Polarized sky | yes | yes | yes | yes | yes | yes |
| Spectral / color (CIE) rendering | - | - | - | partial | - | **yes** |
| Camera imaging chain | - | - | - | yes (+sensor MTF) | - | yes (no sensor MTF) |
| Polarimetric inversion (slope/height) | - | - | - | - | - | **yes** |
| GPU backend | - | - | - | - | - | yes |

---

## References

- Amanatides, J., & Woo, A. (1987). A fast voxel traversal algorithm for ray tracing. *Eurographics '87*, 3-10.
- Angstrom, A. (1929). On the atmospheric transmission of sun radiation and on dust in the air. *Geografiska Annaler*, 11, 156-166.
- Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999). On Rayleigh optical depth calculations. *J. Atmos. Oceanic Technol.*, 16(11), 1854-1861.
- Born, M., & Wolf, E. (1999). *Principles of Optics* (7th ed.). Cambridge University Press.
- Bricaud, A., Babin, M., Morel, A., & Claustre, H. (1995). Variability in the chlorophyll-specific absorption coefficients of natural phytoplankton. *J. Geophys. Res.*, 100(C7), 13321-13332.
- Chu, J., et al. (2020). Simulation of polarization distribution model under wavy water surfaces dominated by skylight. *Acta Optica Sinica* (Sept. 2020). [Volume/pages to confirm.]
- Cox, C., & Munk, W. (1954). Measurement of the roughness of the sea surface from photographs of the sun's glitter. *J. Opt. Soc. Am.*, 44(11), 838-850.
- Cox, C., & Munk, W. (1956). Slopes of the sea surface deduced from photographs of sun glitter. *Bull. Scripps Inst. Oceanogr.*, 6(9), 401-488.
- Curcic, M., & Haus, B. K. (2020). Revised estimates of ocean surface drag in strong winds. *Geophys. Res. Lett.*, 47, e2020GL087647.
- D'Alimonte, D., & Kajiyama, T. (2016). Effects of light polarization and waves slope statistics on the reflectance factor of the sea surface. *Optics Express*, 24(8), 7922-7942.
- Edson, J. B., et al. (2013). On the exchange of momentum over the open ocean. *J. Phys. Oceanogr.*, 43(8), 1589-1610.
- Elfouhaily, T., Chapron, B., Katsaros, K., & Vandemark, D. (1997). A unified directional spectrum for long and short wind-driven waves. *J. Geophys. Res.*, 102(C7), 15781-15796.
- Fedorov, A. V., & Melville, W. K. (1998). Nonlinear gravity-capillary waves with forcing and dissipation. *J. Fluid Mech.*, 354, 1-42.
- Fournier, G. R., & Forand, J. L. (1994). Analytic phase function for ocean water. *Proc. SPIE*, 2258, 194-201.
- Hansen, J. E., & Travis, L. D. (1974). Light scattering in planetary atmospheres. *Space Sci. Rev.*, 16, 527-610.
- Henyey, L. G., & Greenstein, J. L. (1941). Diffuse radiation in the galaxy. *Astrophys. J.*, 93, 70-83.
- Hieronymi, M. (2016). Polarized reflectance and transmittance distribution functions of the ocean surface. *Optics Express*, 24(14), A1045-A1068.
- Loisel, H., & Morel, A. (1998). Light scattering and chlorophyll concentration in case 1 waters. *Limnol. Oceanogr.*, 43(5), 847-858.
- Melville, W. K., & Fedorov, A. V. (2015). The equilibrium dynamics and statistics of gravity-capillary waves. *J. Fluid Mech.*, 767, 449-466.
- Mobley, C. D. (2015). Polarized reflectance and transmittance properties of windblown sea surfaces. *Applied Optics*, 54(15), 4828-4849.
- Moller, T., & Trumbore, B. (1997). Fast, minimum storage ray-triangle intersection. *J. Graphics Tools*, 2(1), 21-28.
- Monahan, E. C., & O'Muircheartaigh, I. (1980). Optimal power-law description of oceanic whitecap coverage dependence on wind speed. *J. Phys. Oceanogr.*, 10(12), 2094-2099.
- Moon, P., & Spencer, D. E. (1942). Illumination from a non-uniform sky. *Illum. Eng.*, 37, 707-726.
- Morel, A. (1974). Optical properties of pure water and pure sea water. In *Optical Aspects of Oceanography* (pp. 1-24). Academic Press.
- Petzold, T. J. (1972). *Volume scattering functions for selected ocean waters*. SIO Ref. 72-78, Scripps Institution of Oceanography.
- Pope, R. M., & Fry, E. S. (1997). Absorption spectrum (380-700 nm) of pure water. II. Integrating cavity measurements. *Applied Optics*, 36(33), 8710-8723.
- Quan, X., & Fry, E. S. (1995). Empirical equation for the index of refraction of seawater. *Applied Optics*, 34(18), 3477-3480.
- Smith, B. G. (1967). Geometrical shadowing of a random rough surface. *IEEE Trans. Antennas Propag.*, 15(5), 668-671.
- Voss, K. J., & Fry, E. S. (1984). Measurement of the Mueller matrix for ocean water. *Applied Optics*, 23(23), 4427-4439.
- Woodcock, E. R., Murphy, T., Hemmings, P. J., & Longworth, T. C. (1965). Techniques used in the GEM code for Monte Carlo neutronics calculations. *Proc. Conf. Applications of Computing Methods to Reactor Problems*, ANL-7050.
- Xue, F., Jin, W., Qiu, S., & Yang, J. (2021). Airborne optical polarization imaging for observation of submarine Kelvin wakes on the sea surface: imaging chain and simulation. *ISPRS J. Photogramm. Remote Sens.*, 178, 136-154. [Cited by some sources as 2022; published in vol. 178, 2021.]
- You, Y., et al. (2011). Polarized light field under dynamic ocean surfaces: numerical modeling compared with measurements. *J. Geophys. Res.*, 116, C00H05.
