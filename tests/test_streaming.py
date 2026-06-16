"""Streaming (frame_callback) synthesis and conditions-keyed beta."""

from pathlib import Path

import numpy as np
import pytest

from seapol import generate_hybrid_surface, generate_sea_surface

DATA = Path("/home/nathanlaxague/Dropbox/Professional/Github/E-PSS_paper/_data")
STATS = DATA / "ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc"
ENV = DATA / "ASIT2019_supporting_environmental_observations.nc"
LIB = Path(__file__).parent.parent / "demos/output/asit_beta_library"


def test_surface_streaming_matches_stack():
    """Streaming frames must be identical to the materialized stack."""
    times = np.array([0.0, 0.2, 0.4])
    kw = dict(L=8.0, N=64, U10=6.0)
    stack = generate_sea_surface(times=times,
                                 rng=np.random.default_rng(3), **kw)
    got = {}

    def cb(it, t, e, gx, gy):
        got[it] = (t, e.copy(), gx.copy(), gy.copy())

    out = generate_sea_surface(times=times, frame_callback=cb,
                               rng=np.random.default_rng(3), **kw)
    assert len(got) == 3
    for it in range(3):
        np.testing.assert_allclose(got[it][1], stack.eta[:, :, it],
                                   atol=1e-12)
        np.testing.assert_allclose(got[it][2], stack.slope_x[:, :, it],
                                   atol=1e-12)
    # returned object carries frame 0
    np.testing.assert_allclose(out.eta, stack.eta[:, :, 0], atol=1e-12)


def test_hybrid_streaming(tiny_table_path=None):
    from seapol import build_fm98_table
    table = build_fm98_table(np.geomspace(70.0, 350.0, 2),
                             np.array([0.05, 0.2]), M_keep=3, M_solve=12,
                             n_steps=4)
    times = np.array([0.0, 0.1])
    seen = []

    def cb(it, t, eta, sx, sy):
        seen.append((it, eta.shape, np.isfinite(eta).all()))

    surf = generate_hybrid_surface(L=1.0, N=128, U10=7.0, table=table,
                                   times=times, frame_callback=cb,
                                   rng=np.random.default_rng(4))
    assert len(seen) == 2
    assert all(s[1] == (128, 128) and s[2] for s in seen)
    assert surf.info.get("streamed") is True
    assert surf.eta.shape == (128, 128)


@pytest.mark.skipif(not (STATS.exists() and LIB.exists()),
                    reason="ASIT data/library not available")
def test_conditions_keyed_beta():
    from seapol import bound_fraction_for_conditions, run_conditions
    conds = run_conditions(STATS, ENV)
    assert np.isfinite(conds["inverse_wave_age"]).sum() > 100
    K = np.array([5.0, 30.0, 200.0])
    for om in (None, 1.0, 3.0):
        b = bound_fraction_for_conditions(LIB, STATS, ENV, 8.0,
                                          inverse_wave_age=om)(K)
        assert np.all((b >= 0) & (b <= 0.99))
        assert b[0] < b[1] <= b[2] + 1e-9   # monotone-ish in k
