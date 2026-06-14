"""
seapol: sea-surface wave simulation and polarized light reflection.

Elfouhaily directional spectrum with FFT synthesis and time evolution,
Fedorov-Melville phase-locked parasitic capillaries, Mobley (2015)
polarized reflection, Rayleigh skylight, single-bounce renderers with
non-Gaussian sun glint and whitecap foam, near-surface scattering (an
in-water polarized Monte Carlo producing the directional sub-surface
light field), a color-aware spectral mode (sky color spectra, water-type
absorption spectra, CIE sRGB output), a forward Monte Carlo tracer with
an in-water scattering medium, k-f spectral diagnostics, polarimetric
inversion back to slope and height, and an optional torch backend that
runs the whole pipeline on a GPU.
"""

from .backend import (TorchRNG, default_rng, get_xp, has_torch, to_numpy,
                      xp_of)
from .spectrum import (angular_frequency, cutoff_slope_variances,
                       directional_spectrum, directional_spread,
                       drag_coefficient, elfouhaily_delta, elfouhaily_omni,
                       inverse_wave_age, phase_speed,
                       total_mean_square_slope)
from .surface import SeaSurface, default_bound_ramp, generate_sea_surface
from .polarization import (apply_mueller, brewster_angle,
                           frame_rotation_angle, fresnel_mueller,
                           meridian_frame, mueller_rotation, normalize,
                           reflection_chain, stokes_aop, stokes_dolp,
                           transmission_chain)
from .skylight import (aerosol_optical_depth, direction_from_angles,
                       dolp_max, overcast_sky_stokes, overcast_spectrum,
                       rayleigh_optical_depth, rayleigh_sky_stokes,
                       rayleigh_sky_stokes_angles, sky_radiance_spectrum,
                       solar_spectrum, sun_beam_spectrum, sun_disk_stokes)
from .render import (CameraGeometry, Foam, PinholeCamera, SubpixelSlopes,
                     compose_skies, downwelling_irradiance, make_clear_sky,
                     make_overcast_sky, make_partly_cloudy_sky,
                     make_rayleigh_sky, make_unpolarized_sky,
                     monahan_coverage, render_camera_image,
                     render_facet_stokes, render_facet_stokes_stack,
                     saunders_bistatic, smith_lambda, sun_glint_stokes)
from .water import (WATER_TYPES, WaterBody, WaterColumn, WaterOptics,
                    WaterType, ff_phase_mueller, fournier_forand_phase,
                    hg_phase_mueller, polarized_scatter_event,
                    pure_water_absorption, pure_water_scattering,
                    rayleigh_phase_mueller, sample_ff_scattering,
                    sample_hg_scattering, sample_rayleigh_scattering,
                    water_leaving_stokes, water_refractive_index)
from .scattering import (UpwellingRadianceTable, build_upwelling_table,
                         load_table, save_table, water_leaving_from_table)
from .spectral import (SpectralBands, build_spectral_tables,
                       render_camera_image_spectral,
                       render_facet_stokes_spectral,
                       spectral_sky_factories, spectral_sun_irradiance)
from .color import cie_cmf, srgb_encode, stokes_bands_to_rgb
from .montecarlo import (build_cell_tables, effective_mueller_for_incident,
                         heightmap_intersect, trace_forward)
from .fm98 import (FM98Table, FMSolution, build_fm98_table, gauge_fix,
                   linear_phase_speed, solve_fm98, solve_fm98_continuation)
from .hybrid import (HybridSeaSurface, apply_orbital_advection,
                     augment_fm98, default_carrier_band,
                     default_fm98_table, generate_hybrid_surface,
                     long_wave_modulation)
from .diagnostics import kf_slope_spectrum, q_nu, slow_nu_fraction
from .inversion import (fresnel_dolp_curve, height_from_slopes,
                        slopes_from_stokes, slopes_from_stokes_polarized)
from .empirical import (bound_fraction_for_conditions,
                        bound_fraction_for_wind,
                        bound_fraction_from_kf_reduction, cube_timestamp,
                        generate_asit_surface, load_asit_run,
                        match_env_run, psi_from_asit, reduce_kf_cube,
                        run_conditions, smooth_beta_curve)

__version__ = "0.2.0"
