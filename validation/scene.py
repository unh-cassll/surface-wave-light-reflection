"""Shared scene definition for cross-validation against external renderers.

A SceneSpec is the single source of truth read by every renderer
(render_seapol.py, render_blender.py, render_mitsuba.py): surface, camera,
sun, and water refractive index.  Only numpy is needed for JSON I/O and the
camera basis, so this module imports cleanly inside Blender's bundled Python;
seapol is imported lazily only where the surface is actually synthesized.

Camera and OBJ conventions mirror seapol exactly:
    * camera basis reproduces seapol.render._camera_rays
    * OBJ triangulation reproduces seapol.montecarlo.build_cell_tables
      (corner order [v00, v01, v11, v10] -> tris (v00,v01,v11),(v00,v11,v10))
so an external render of the exported mesh is geometrically registered to the
seapol render of the same surface.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

import numpy as np


@dataclass
class SceneSpec:
    """Full description of a validation scene."""
    # surface synthesis
    L: float = 4.0
    N: int = 256
    U10: float = 7.0
    wind_dir_deg: float = 0.0
    seed: int = 0
    flat: bool = False                 # flat z=0 sheet (sanity-check geometry)
    # camera (mirrors seapol.render.PinholeCamera)
    altitude_m: float = 30.0
    zenith_deg: float = 35.0
    azimuth_deg: float = 0.0
    hfov_deg: float = 6.0
    img_h: int = 256
    img_w: int = 256
    # illumination
    sun_zenith_deg: float = 45.0
    sun_azimuth_deg: float = 90.0
    I_sky: float = 1.0
    unpolarized_sky: bool = False      # match Mitsuba's constant sky (Fresnel-only)
    sky_rgb: tuple = (0.20, 0.45, 1.0)  # blue sky color for Mitsuba/Blender color S0
    # water
    wavelength_nm: float = 550.0
    salinity_psu: float = 35.0
    temperature_c: float = 20.0
    n_water: float = field(default=0.0)   # 0 -> filled from Quan & Fry on demand

    # ---- serialization -------------------------------------------------
    def to_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(asdict(self), f, indent=2)

    @classmethod
    def from_json(cls, path: str) -> "SceneSpec":
        with open(path) as f:
            return cls(**json.load(f))

    # ---- derived quantities -------------------------------------------
    @property
    def img_shape(self) -> tuple[int, int]:
        return (self.img_h, self.img_w)

    def refractive_index(self) -> float:
        """Seawater n; uses the stored value if set, else Quan & Fry (1995)."""
        if self.n_water and self.n_water > 0:
            return float(self.n_water)
        from seapol import water_refractive_index
        return float(water_refractive_index(self.wavelength_nm,
                                             self.salinity_psu,
                                             self.temperature_c))

    def camera_basis(self) -> dict:
        """Pinhole origin/look/right/up and FOV, identical to
        seapol.render._camera_rays.  Pure numpy (Blender-safe)."""
        L = self.N * (self.L / self.N)   # = self.L; written so dx*N is explicit
        th = np.deg2rad(self.zenith_deg)
        ph = np.deg2rad(self.azimuth_deg)
        z0 = self.altitude_m
        center = np.array([L / 2.0, L / 2.0, 0.0])
        origin = center + np.array([np.tan(th) * z0 * np.cos(ph),
                                    np.tan(th) * z0 * np.sin(ph), z0])
        look = (center - origin) / np.linalg.norm(center - origin)
        right = np.cross(look, np.array([0.0, 0.0, 1.0]))
        right = right / np.linalg.norm(right)
        up = np.cross(right, look)
        H, W = self.img_shape
        # seapol places the outer pixel EDGES of the larger axis at
        # +/- tan(hfov) (edge-aligned, matching Mitsuba/Blender)
        half = np.tan(np.deg2rad(self.hfov_deg))
        return dict(origin=origin, look=look, right=right, up=up,
                    center=center, half=half, H=H, W=W, L=L)


def generate_surface(spec: SceneSpec):
    """Synthesize the static surface for a spec.  Returns (eta, dx, info)."""
    import numpy as _np
    from seapol import generate_sea_surface
    dx = spec.L / spec.N
    if spec.flat:
        eta = _np.zeros((spec.N, spec.N))
        return eta, dx, {"dx": dx, "sigma_a2_cut": 0.0, "sigma_c2_cut": 0.0}
    surf = generate_sea_surface(L=spec.L, N=spec.N, U10=spec.U10,
                                wind_dir_rad=_np.deg2rad(spec.wind_dir_deg),
                                rng=_np.random.default_rng(spec.seed))
    eta = surf.eta if surf.eta.ndim == 2 else surf.eta[:, :, 0]
    sx = surf.slope_x if surf.slope_x.ndim == 2 else surf.slope_x[:, :, 0]
    sy = surf.slope_y if surf.slope_y.ndim == 2 else surf.slope_y[:, :, 0]
    info = dict(surf.info)
    info["slope_x"], info["slope_y"] = sx, sy
    return eta, dx, info


def export_obj(eta, dx: float, path: str) -> None:
    """Write the surface as a Wavefront OBJ, reusing seapol's facet
    triangulation convention (open patch; the periodic wrap is a seapol MC
    internal and is omitted so external renderers see no domain-spanning
    facet).  Vertex (i, j) is at (x=j*dx, y=i*dx, z=eta[i, j])."""
    eta = np.asarray(eta, dtype=float)
    N = eta.shape[0]
    if eta.shape[0] != eta.shape[1]:
        raise ValueError("eta must be square")
    jj, ii = np.meshgrid(np.arange(N), np.arange(N))
    xs = (jj * dx).ravel()
    ys = (ii * dx).ravel()
    zs = eta.ravel()

    def vid(i, j):                      # 1-based OBJ vertex id
        return i * N + j + 1

    i0 = np.arange(N - 1)[:, None]
    j0 = np.arange(N - 1)[None, :]
    v00 = (i0 * N + j0 + 1).ravel()
    v01 = (i0 * N + (j0 + 1) + 1).ravel()
    v11 = ((i0 + 1) * N + (j0 + 1) + 1).ravel()
    v10 = ((i0 + 1) * N + j0 + 1).ravel()

    with open(path, "w") as f:
        f.write("# seapol validation surface mesh\n")
        for x, y, z in zip(xs, ys, zs):
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        # tri 1 = (v00, v01, v11), tri 2 = (v00, v11, v10)
        for a, b, c in zip(v00, v01, v11):
            f.write(f"f {a} {b} {c}\n")
        for a, b, c in zip(v00, v11, v10):
            f.write(f"f {a} {b} {c}\n")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Write a default SceneSpec JSON")
    ap.add_argument("path", nargs="?", default="scene.json")
    args = ap.parse_args()
    import os
    d = os.path.dirname(args.path)
    if d:
        os.makedirs(d, exist_ok=True)
    SceneSpec().to_json(args.path)
    print(f"wrote {args.path}")
