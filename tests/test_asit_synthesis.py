"""Synthesis from measured ASIT directional slope spectra.  Skipped when
the data files are not present."""

from pathlib import Path

import numpy as np
import pytest

DATA = Path("/home/nathanlaxague/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = DATA / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = DATA / "ASIT2019_supporting_environmental_observations.nc"

pytestmark = pytest.mark.skipif(not STATS.exists(),
                                reason="ASIT data not available")

RUN = 80


@pytest.fixture(scope="module")
def run_data():
    from seapol import load_asit_run
    return load_asit_run(STATS, RUN, env_path=ENV)


def test_convention_against_stored_mss(run_data):
    """Cartesian slope-density convention reproduces the stored MSS
    components within the gain/processing scatter."""
    k, th, S = run_data["k"], run_data["theta"], run_data["S"]
    dth = np.diff(th).mean()
    m_u = np.trapezoid(k * ((np.cos(th) ** 2)[:, None] * S).sum(0) * dth, k)
    m_c = np.trapezoid(k * ((np.sin(th) ** 2)[:, None] * S).sum(0) * dth, k)
    assert 0.5 < m_u / run_data["mss_upwind"] < 1.6
    assert 0.5 < m_c / run_data["mss_crosswind"] < 1.6


def test_psi_interpolation_conserves_mss(run_data):
    """Psi_eta on a fine grid integrates back to the measured-band MSS."""
    from seapol import psi_from_asit
    psi = psi_from_asit(run_data)
    N, L = 1024, 2.9
    kx = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, kx, indexing="xy")
    P = psi(KX, KY)
    dk = 2 * np.pi / L
    K = np.hypot(KX, KY)
    band = K <= run_data["k"][-1]
    mss_grid = float(np.sum((K**2 * P)[band]) * dk * dk)
    k, th, S = run_data["k"], run_data["theta"], run_data["S"]
    dth = np.diff(th).mean()
    mss_meas = np.trapezoid(k * S.sum(0) * dth, k)
    assert abs(mss_grid / mss_meas - 1.0) < 0.25, (mss_grid, mss_meas)


def test_generate_asit_surface(run_data):
    """End-to-end synthesis from the measured spectrum: realized slope
    variance matches the spectral target on the resolved band."""
    from seapol import generate_asit_surface
    surf = generate_asit_surface(STATS, RUN, L=2.9, N=512, env_path=ENV,
                                 rng=np.random.default_rng(0))
    assert np.isfinite(surf.eta).all()
    mss_real = surf.slope_x.var() + surf.slope_y.var()
    ratio = mss_real / surf.info["mss_resolved"]
    assert abs(ratio - 1.0) < 0.1, ratio
    assert surf.info["asit_run"] == RUN
    # in the ballpark of the instrument-reported MSS
    mss_meas = sum(surf.info["mss_measured"])
    assert 0.4 < mss_real / mss_meas < 2.5


def test_asit_surface_with_bound_partition(run_data):
    from seapol import generate_asit_surface
    surf = generate_asit_surface(STATS, RUN, L=2.9, N=256, env_path=ENV,
                                 bound_fraction=0.9,
                                 bound_speed="spectrum",
                                 times=np.array([0.0, 0.1]),
                                 rng=np.random.default_rng(1))
    assert surf.eta.shape == (256, 256, 2)
    assert np.isfinite(surf.eta).all()


def test_psi_band_limit(run_data):
    """k_max tapers the spectrum to zero above the cutoff (instrument-
    resolution band-limit), removing the high-k slope variance that
    aliases under time evolution while preserving the low-k field."""
    from seapol import psi_from_asit
    N, L = 512, 2.9
    kx = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY = np.meshgrid(kx, kx, indexing="xy")
    K = np.hypot(KX, KY)
    k_cut = 150.0
    full = psi_from_asit(run_data)(KX, KY)
    lim = psi_from_asit(run_data, k_max=k_cut)(KX, KY)
    # nothing survives above the cutoff; the low-k band is untouched
    assert lim[K > k_cut].max() == 0.0
    low = K < 0.8 * k_cut
    np.testing.assert_allclose(lim[low], full[low], rtol=1e-9)
    # it removes real slope variance (a cut, not a no-op)
    dk2 = (2 * np.pi / L) ** 2
    assert float(np.sum(K**2 * lim) * dk2) < float(np.sum(K**2 * full) * dk2)


def test_generate_asit_surface_k_max(run_data):
    from seapol import generate_asit_surface
    surf = generate_asit_surface(STATS, RUN, L=2.9, N=256, env_path=ENV,
                                 k_max=150.0, rng=np.random.default_rng(2))
    assert np.isfinite(surf.eta).all()
    # the realized field carries no resolved content above the cutoff
    N = 256
    kx = 2 * np.pi * np.fft.fftfreq(N, d=2.9 / N)
    K = np.hypot(*np.meshgrid(kx, kx, indexing="xy"))
    P = np.abs(np.fft.fft2(surf.eta)) ** 2
    assert P[K > 170.0].sum() / P[K > 1e-9].sum() < 1e-6
