"""Surface synthesis checks: variance closure, slopes, time evolution."""

import numpy as np

from seapol import generate_sea_surface


def test_variance_matches_spectral_target():
    surf = generate_sea_surface(L=256.0, N=512, U10=7.0,
                                rng=np.random.default_rng(0))
    ratio = surf.eta.var() / surf.info["var_target"]
    assert 0.9 < ratio < 1.1, ratio


def test_slope_variance_matches_target():
    surf = generate_sea_surface(L=128.0, N=512, U10=7.0,
                                rng=np.random.default_rng(1))
    mss = surf.slope_x.var() + surf.slope_y.var()
    ratio = mss / surf.info["mss_resolved"]
    assert 0.9 < ratio < 1.1, ratio


def test_one_sided_preserves_variance():
    a = generate_sea_surface(L=128.0, N=256, U10=7.0, one_sided=True,
                             rng=np.random.default_rng(2))
    b = generate_sea_surface(L=128.0, N=256, U10=7.0, one_sided=False,
                             rng=np.random.default_rng(2))
    # unpaired Nyquist modes leave a tiny residual on even grids
    np.testing.assert_allclose(a.info["var_target"], b.info["var_target"],
                               rtol=1e-4)


def test_time_evolution_consistency():
    rng_a = np.random.default_rng(3)
    rng_b = np.random.default_rng(3)
    static = generate_sea_surface(L=64.0, N=128, U10=6.0, rng=rng_a)
    movie = generate_sea_surface(L=64.0, N=128, U10=6.0,
                                 times=np.array([0.0, 0.5, 1.0]), rng=rng_b)
    np.testing.assert_allclose(movie.eta[:, :, 0], static.eta, atol=1e-12)
    # stationarity of frame variance
    v = movie.eta.reshape(-1, 3).var(axis=0)
    assert np.all(np.abs(v / v[0] - 1.0) < 0.05)
    # frames actually evolve
    assert np.abs(movie.eta[:, :, 1] - movie.eta[:, :, 0]).max() > 1e-3


def test_waves_propagate_downwind():
    """With one_sided spreading the dominant waves move toward +x."""
    dt = 0.4
    surf = generate_sea_surface(L=128.0, N=128, U10=7.0, one_sided=True,
                                times=np.array([0.0, dt]),
                                rng=np.random.default_rng(4))
    f0 = np.fft.fft2(surf.eta[:, :, 0])
    f1 = np.fft.fft2(surf.eta[:, :, 1])
    xcorr = np.fft.ifft2(f1 * np.conj(f0)).real
    iy, ix = np.unravel_index(np.argmax(xcorr), xcorr.shape)
    shift_x = ix if ix < 64 else ix - 128
    # peak-wave phase speed ~ 8.4 m/s -> ~3.4 m in 0.4 s, dx = 1 m
    assert 1 <= shift_x <= 8, shift_x


def test_wind_direction_rotates_field():
    surf = generate_sea_surface(L=128.0, N=256, U10=8.0,
                                wind_dir_rad=np.pi / 2,
                                rng=np.random.default_rng(5))
    # along-wind (+y now) slope variance should dominate
    assert surf.slope_y.var() > surf.slope_x.var()


def test_bound_partition_preserves_power():
    """The free/bound split must leave variance and MSS on the spectral
    target (independent populations, power partitioned per bin)."""
    free = generate_sea_surface(L=64.0, N=512, U10=7.0,
                                rng=np.random.default_rng(7))
    part = generate_sea_surface(L=64.0, N=512, U10=7.0, bound_fraction=0.9,
                                rng=np.random.default_rng(7))
    # peak-band variance is realization-limited on a 64 m tile, so compare
    # against the free realization (identical low-k draw); MSS is high-k
    # dominated and must hold tightly
    assert abs(part.eta.var() / free.eta.var() - 1.0) < 0.05
    mss_f = free.slope_x.var() + free.slope_y.var()
    mss_p = part.slope_x.var() + part.slope_y.var()
    assert abs(mss_p / mss_f - 1.0) < 0.05


def test_bound_part_advects_at_bound_speed():
    """With beta = 1 above the ramp, the high-pass field translates
    rigidly downwind at bound_speed."""
    c_b, dt = 2.0, 0.12
    surf = generate_sea_surface(L=4.0, N=256, U10=7.0, bound_fraction=1.0,
                                bound_speed=c_b,
                                times=np.array([0.0, dt]),
                                rng=np.random.default_rng(8))
    dx = surf.info["dx"]
    kx = 2 * np.pi * np.fft.fftfreq(256, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    hp = []
    for it in range(2):
        Z = np.fft.fft2(surf.eta[:, :, it])
        hp.append(np.fft.ifft2(np.where(K > 80.0, Z, 0)).real)
    xc = np.fft.ifft2(np.fft.fft2(hp[1]) * np.conj(np.fft.fft2(hp[0]))).real
    iy, ix = np.unravel_index(np.argmax(xc), xc.shape)
    shift_cells = ix if ix < 128 else ix - 256
    expected = c_b * dt / dx
    assert abs(shift_cells - expected) <= 2, (shift_cells, expected)
    assert iy in (0, 1, 255)  # no crosswind drift
