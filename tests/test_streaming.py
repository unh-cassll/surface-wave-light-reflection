"""Streaming (frame_callback) synthesis and conditions-keyed beta."""

import numpy as np
import pytest

from seapol import generate_hybrid_surface, generate_sea_surface


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


def _conditions_fixture(tmp_path, n_runs=6):
    """Synthetic ASIT-format stats/env netCDFs and beta-reduction library
    covering the fields run_conditions and bound_fraction_for_conditions
    read."""
    nc = pytest.importorskip("netCDF4")
    f = np.linspace(0.05, 2.0, 40)
    th = np.linspace(-np.pi, np.pi, 24, endpoint=False)
    fp_true = np.linspace(0.15, 0.4, n_runs)

    stats = tmp_path / "stats.nc"
    d = nc.Dataset(stats, "w")
    d.createDimension("run", n_runs)
    d.createDimension("f", f.size)
    d.createDimension("theta", th.size)
    d.createVariable("f_Hz", "f8", ("f",))[:] = f
    d.createVariable("theta_rad", "f8", ("theta",))[:] = th
    Sf = d.createVariable("S_f_theta", "f8", ("run", "theta", "f"))
    for r in range(n_runs):
        Sf[r] = np.broadcast_to(
            np.exp(-0.5 * ((f - fp_true[r]) / 0.05) ** 2), (th.size, f.size))
    d.close()

    env = tmp_path / "env.nc"
    d = nc.Dataset(env, "w")
    d.createDimension("run", n_runs)
    d.createVariable("t_seconds_since_January_1_1970", "f8",
                     ("run",))[:] = np.arange(n_runs) * 3600.0
    d.createVariable("EC_U_m_s", "f8",
                     ("run",))[:] = np.linspace(4.0, 12.0, n_runs)
    d.close()

    lib = tmp_path / "beta_lib"
    lib.mkdir()
    k = np.geomspace(3.0, 1400.0, 40)
    for r in range(n_runs):
        U = 4.0 + 8.0 * r / (n_runs - 1)
        beta = (0.9 * U / 12.0) / (1.0 + (35.0 / k) ** 2)
        np.savez(lib / f"run{r:03d}.npz", run=r, U10=U, k=k, beta_obs=beta)
    return stats, env, lib


def test_conditions_keyed_beta(tmp_path):
    from seapol import bound_fraction_for_conditions, run_conditions
    stats, env, lib = _conditions_fixture(tmp_path)
    conds = run_conditions(stats, env)
    assert np.isfinite(conds["inverse_wave_age"]).sum() == 6
    K = np.array([5.0, 30.0, 200.0])
    for om in (None, 1.0, 3.0):
        b = bound_fraction_for_conditions(lib, stats, env, 8.0,
                                          inverse_wave_age=om)(K)
        assert np.all((b >= 0) & (b <= 0.99))
        assert b[0] < b[1] <= b[2] + 1e-9   # monotone-ish in k
