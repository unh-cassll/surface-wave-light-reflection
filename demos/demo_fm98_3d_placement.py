"""
Local (3D-appropriate) placement of Fedorov-Melville parasitic
capillaries: each wave group carries the ripple train of its own scale
and direction.

The FM98 capillary bump sits at the harmonic number m* where
c(m* k) ~ c(k), i.e. at ripple wavenumber k_rip ~ k_m^2 / k_carrier
(k_m = sqrt(g rho / sigma) ~ 363 rad/m): SHORTER carriers bear LONGER
ripples.  carrier_k="local" reproduces that inversion across the field;
the band-averaged "representative" scheme pins every ripple train to
one scale.  Panels: local carrier wavenumber map, bound-harmonic
curvature zoom, and conditional ripple spectra by k_loc tercile for
both schemes.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demo_fm98_capillaries import deep_table
from seapol import augment_fm98, generate_sea_surface, long_wave_modulation
from seapol.hybrid import default_carrier_band

OUT = Path(__file__).parent / "output"

L, N, U10 = 2.0, 2048, 9.0
MTF, MTF_PHASE_DEG = 6.5, 30.0
K_SPLIT = 2 * np.pi / 0.5
K_M = 363.0                         # clean-water gravity-capillary minimum


def carrier_envelope_fields(eta, dx, k_lo, k_hi):
    """Carrier-band envelope A, phase phi, and local wavenumber k_loc
    (same construction as augment_fm98)."""
    n = eta.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    K2 = kx[None, :] ** 2 + kx[:, None] ** 2
    Z = np.fft.fft2(eta)
    Z_band = np.where((K2 >= k_lo**2) & (K2 <= k_hi**2), Z, 0.0)
    F_a = Z_band * np.where(kx[None, :] > 0, 2.0,
                            np.where(kx[None, :] < 0, 0.0, 1.0))
    Z_a = np.fft.ifft2(F_a)
    A = np.abs(Z_a)
    A2 = np.maximum(A**2, 1e-30)
    dZx = np.fft.ifft2(1j * kx[None, :] * F_a)
    dZy = np.fft.ifft2(1j * kx[:, None] * F_a)
    k_loc = np.hypot((np.conj(Z_a) * dZx).imag / A2,
                     (np.conj(Z_a) * dZy).imag / A2)
    return A, np.angle(Z_a), np.clip(k_loc, k_lo, k_hi)


def radial_spectrum(field, dx, n_bins=160):
    n = field.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    K = np.hypot(kx[None, :], kx[:, None])
    P = np.abs(np.fft.fft2(field)) ** 2
    edges = np.linspace(0, kx.max(), n_bins + 1)
    idx = np.clip(np.digitize(K.ravel(), edges) - 1, 0, n_bins - 1)
    pow_k = np.bincount(idx, weights=P.ravel(), minlength=n_bins)
    return 0.5 * (edges[:-1] + edges[1:]), pow_k


def main():
    table = deep_table()
    dx = L / N
    k_lo, k_hi = default_carrier_band(np.pi / dx)

    print("synthesizing linear surface + long-wave MTF binding...")
    lin = generate_sea_surface(L, N, U10, rng=np.random.default_rng(4))
    eta_mod, _ = long_wave_modulation(lin.eta, dx, k_split=K_SPLIT,
                                      mtf=MTF, mtf_phase_deg=MTF_PHASE_DEG)

    A, _, k_loc = carrier_envelope_fields(eta_mod, dx, k_lo, k_hi)

    results = {}
    for mode in ("local", "representative"):
        hi, _, info = augment_fm98(eta_mod, dx, table, carrier_k=mode)
        results[mode] = hi
        print(f"{mode:>14s}: k_rep={info['k_rep']:.0f} rad/m, "
              f"k_loc p50/p90 = {info['k_loc_p50']:.0f}/"
              f"{info['k_loc_p90']:.0f}, var_high={info['var_high']:.2e}")

    # Carrier-scale groups: explicit k_loc windows whose resonance
    # harmonic m* = (k_m / k_loc)^2 stays inside the table depth
    # (M_keep = 28), restricted to the strongest 30% of the envelope
    groups = {"long carriers": (90.0, 110.0),
              "short carriers": (150.0, 200.0)}
    akm = A * k_loc
    ak_thresh = np.percentile(akm, 70)

    fig = plt.figure(figsize=(15, 4.8), dpi=110)

    # --- local carrier wavenumber map
    ax = fig.add_subplot(1, 3, 1)
    show = np.where(A * k_loc > 0.02, k_loc, np.nan)
    im = ax.imshow(show[:512, :512], origin="lower", cmap="viridis",
                   extent=[0, 512 * dx, 0, 512 * dx],
                   vmin=k_lo, vmax=k_hi)
    fig.colorbar(im, ax=ax, label=r"$k_{loc}$ [rad/m]")
    ax.set_title("Local carrier wavenumber (ak > 0.02)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    # --- bound-harmonic curvature zoom at the steepest group; window
    # clamped inside the domain so the group stays in view near edges
    iy, ix = np.unravel_index(np.argmax(A * k_loc), A.shape)
    half = 192
    y0 = int(np.clip(iy - half, 0, N - 2 * half))
    x0 = int(np.clip(ix - half, 0, N - 2 * half))
    sl = (slice(y0, y0 + 2 * half), slice(x0, x0 + 2 * half))
    hi = results["local"]
    curv = (np.gradient(np.gradient(hi, dx, axis=1), dx, axis=1)
            + np.gradient(np.gradient(hi, dx, axis=0), dx, axis=0))
    ax = fig.add_subplot(1, 3, 2)
    z = curv[sl]
    vmax = np.percentile(np.abs(z), 99)
    ax.imshow(z, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
              extent=[0, z.shape[1] * dx, 0, z.shape[0] * dx])
    ax.contour(np.arange(z.shape[1]) * dx, np.arange(z.shape[0]) * dx,
               eta_mod[sl], levels=6, colors="k", linewidths=0.5,
               alpha=0.5)
    ax.set_title("Bound-harmonic curvature (steepest group)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")

    # --- conditional ripple spectra by local carrier scale
    from scipy.ndimage import gaussian_filter
    ax = fig.add_subplot(1, 3, 3)
    colors = plt.cm.viridis([0.2, 0.8])
    for (name, (ka, kb)), col in zip(groups.items(), colors):
        sel = (k_loc >= ka) & (k_loc < kb) & (akm > ak_thresh)
        k_cen = float(np.average(k_loc[sel], weights=A[sel] ** 2))
        m_s = gaussian_filter(sel.astype(float), sigma=12)
        for mode, ls in (("local", "-"), ("representative", "--")):
            kr, pw = radial_spectrum(results[mode] * m_s, dx)
            band = (kr > 1.5 * k_hi) & (kr < 2500.0)
            pw = pw / np.mean(m_s**2)
            lbl = (f"{name} ($k_{{loc}}$~{k_cen:.0f})"
                   if mode == "local" else None)
            ax.semilogx(kr[band], (pw * kr**3)[band], ls, color=col,
                        lw=1.5, label=lbl)
        ax.axvline(K_M**2 / k_cen, color=col, ls=":", alpha=0.9)
        print(f"{name}: <k_loc> = {k_cen:.0f} rad/m, "
              f"resonance k_rip = {K_M**2 / k_cen:.0f} rad/m")
    ax.plot([], [], "k-", label="local")
    ax.plot([], [], "k--", label="representative")
    ax.set_xlabel(r"$k$ [rad/m]")
    ax.set_ylabel(r"$k^3 \times$ ripple power (per masked area)")
    ax.set_title("Conditional ripple spectra by carrier scale\n"
                 r"(dotted: resonance $k_m^2/k_{loc}$)")
    ax.legend(fontsize=7, loc="upper left")

    plt.tight_layout()
    OUT.mkdir(exist_ok=True)
    out = OUT / "demo_fm98_3d_placement.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
