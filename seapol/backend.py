"""
Array-backend dispatch: numpy (default) or torch (CPU / CUDA / ROCm).

Every array-consuming seapol function obtains its array namespace from
its inputs via xp_of(*arrays): numpy inputs run through numpy unchanged;
torch tensors run through a numpy-compatible shim bound to the tensors'
device and float dtype, so the whole pipeline (synthesis -> hybrid
augmentation -> rendering -> Monte Carlo -> diagnostics) stays on-device
end to end.  Entry points that create arrays accept backend/device/dtype
and build a creation namespace with get_xp().

Conventions:
    * The numpy shim forwards every attribute to numpy, so the default
      path is numerically identical to plain numpy.
    * The torch shim maps numpy call signatures (axis -> dim,
      keepdims -> keepdim, ddof-0 std/var, endpoint-aware linspace, ...)
      and carries (device, dtype); zeros/ones/empty/full/linspace/
      geomspace default to the shim float dtype.
    * RNG: default_rng(seed, xp) returns a numpy Generator or a TorchRNG
      on the shim device.  Passing a numpy Generator into a torch-backed
      call wraps it so draws happen in numpy and are copied to the
      device -- slower, but bitwise-reproducible across backends, which
      is what the parity tests use.

Not dispatched: the FM98 Newton solver (tiny dense systems, scipy) and
netCDF ingest in seapol.empirical stay on the CPU; their products
(coefficient tables, measured spectra) are converted on entry to the
array pipeline.
"""

from __future__ import annotations

import sys

import numpy as np

__all__ = ["xp_of", "get_xp", "default_rng", "adapt_rng", "to_numpy",
           "has_torch", "TorchRNG", "NUMPY_XP"]

_TORCH = None


def _torch():
    global _TORCH
    if _TORCH is None:
        import torch
        _TORCH = torch
    return _TORCH


def has_torch() -> bool:
    try:
        _torch()
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# numpy shim: forwards everything, adds the few cross-backend extras
# ---------------------------------------------------------------------------

class _NumpyShim:
    is_torch = False
    name = "numpy"
    device = "cpu"

    def __getattr__(self, attr):
        return getattr(np, attr)

    @staticmethod
    def asarray(x, dtype=None):
        if dtype is not None:
            return np.asarray(x, dtype=dtype)
        return np.asarray(x)

    @staticmethod
    def astype(x, dtype):
        return np.asarray(x).astype(dtype)

    @staticmethod
    def copy(x):
        return np.array(x)

    @staticmethod
    def erfc(x):
        from scipy.special import erfc
        return erfc(x)

    @staticmethod
    def index_add(dest, index, src):
        """dest[index] += src along axis 0, with repeated indices."""
        np.add.at(dest, index, src)
        return dest

    @staticmethod
    def to_device(x):
        return np.asarray(x)


NUMPY_XP = _NumpyShim()


# ---------------------------------------------------------------------------
# torch shim
# ---------------------------------------------------------------------------

class _TorchFFT:
    def __init__(self, shim):
        self._s = shim

    def fft(self, x, axis=-1):
        return _torch().fft.fft(self._s.asarray(x), dim=axis)

    def ifft(self, x, axis=-1):
        return _torch().fft.ifft(self._s.asarray(x), dim=axis)

    def fft2(self, x):
        return _torch().fft.fft2(self._s.asarray(x))

    def ifft2(self, x):
        return _torch().fft.ifft2(self._s.asarray(x))

    def fftfreq(self, n, d=1.0):
        return _torch().fft.fftfreq(n, d=d, device=self._s.device,
                                    dtype=self._s.dtype)

    def fftshift(self, x, axes=None):
        return _torch().fft.fftshift(x, dim=axes)


class _TorchLinalg:
    def __init__(self, shim):
        self._s = shim

    def norm(self, x, axis=None, keepdims=False):
        t = _torch()
        if axis is None:
            return t.linalg.vector_norm(x)
        return t.linalg.vector_norm(x, dim=axis, keepdim=keepdims)

    def lstsq(self, a, b, rcond=None):
        res = _torch().linalg.lstsq(a, b)
        return res.solution, res.residuals, None, None


class _TorchShim:
    """numpy-compatible namespace over torch, bound to (device, dtype)."""

    is_torch = True
    name = "torch"

    pi = float(np.pi)
    inf = float("inf")
    nan = float("nan")
    newaxis = None

    def __init__(self, device, dtype):
        t = _torch()
        self.device = t.device(device)
        self.dtype = dtype
        self.cdtype = (t.complex64 if dtype == t.float32 else t.complex128)
        self.fft = _TorchFFT(self)
        self.linalg = _TorchLinalg(self)

    # -- dtype mapping ------------------------------------------------------
    def _map_dtype(self, dtype):
        t = _torch()
        if dtype is None:
            return None
        if isinstance(dtype, t.dtype):
            return dtype
        if dtype in (float, "float", np.float64, np.float32, "float64",
                     "float32"):
            return self.dtype
        if dtype in (complex, "complex", np.complex128, np.complex64):
            return self.cdtype
        if dtype in (int, "int", np.int64, np.int32, "int64"):
            return t.int64
        if dtype in (bool, "bool", np.bool_):
            return t.bool
        raise TypeError(f"unsupported dtype for torch backend: {dtype!r}")

    def _result_float(self, x):
        """Promote integer/bool tensors to the shim float dtype."""
        if x.dtype.is_floating_point or x.dtype.is_complex:
            return x
        return x.to(self.dtype)

    # -- creation -----------------------------------------------------------
    def asarray(self, x, dtype=None):
        t = _torch()
        d = self._map_dtype(dtype)
        if isinstance(x, t.Tensor):
            out = x.to(self.device) if x.device != self.device else x
            return out.to(d) if d is not None and out.dtype != d else out
        if d is None:
            arr = np.asarray(x)
            if arr.dtype.kind == "f":
                d = self.dtype
            elif arr.dtype.kind == "c":
                d = self.cdtype
        return t.as_tensor(np.asarray(x), dtype=d, device=self.device)

    def array(self, x, dtype=None):
        out = self.asarray(x, dtype=dtype)
        return out.clone()

    def copy(self, x):
        return self.asarray(x).clone()

    def astype(self, x, dtype):
        return self.asarray(x).to(self._map_dtype(dtype))

    def to_device(self, x):
        return self.asarray(x)

    def _size(self, shape):
        if isinstance(shape, (int, np.integer)):
            return (int(shape),)
        return tuple(int(s) for s in shape)

    def zeros(self, shape, dtype=float):
        return _torch().zeros(self._size(shape), dtype=self._map_dtype(dtype),
                              device=self.device)

    def ones(self, shape, dtype=float):
        return _torch().ones(self._size(shape), dtype=self._map_dtype(dtype),
                             device=self.device)

    def empty(self, shape, dtype=float):
        return _torch().empty(self._size(shape), dtype=self._map_dtype(dtype),
                              device=self.device)

    def full(self, shape, fill, dtype=None):
        if dtype is None:
            dtype = bool if isinstance(fill, bool) else float
        return _torch().full(self._size(shape), fill,
                             dtype=self._map_dtype(dtype), device=self.device)

    def zeros_like(self, x, dtype=None):
        return _torch().zeros_like(x, dtype=self._map_dtype(dtype))

    def ones_like(self, x, dtype=None):
        return _torch().ones_like(x, dtype=self._map_dtype(dtype))

    def empty_like(self, x, dtype=None):
        return _torch().empty_like(x, dtype=self._map_dtype(dtype))

    def full_like(self, x, fill, dtype=None):
        return _torch().full_like(x, fill, dtype=self._map_dtype(dtype))

    def arange(self, *args, dtype=None):
        t = _torch()
        d = self._map_dtype(dtype)
        if d is None and any(isinstance(a, float) for a in args):
            d = self.dtype
        return t.arange(*args, dtype=d, device=self.device)

    def linspace(self, start, stop, num=50, endpoint=True, dtype=float):
        t = _torch()
        d = self._map_dtype(dtype)
        if endpoint:
            return t.linspace(start, stop, num, dtype=d, device=self.device)
        return t.linspace(start, stop, num + 1, dtype=d,
                          device=self.device)[:-1]

    def geomspace(self, start, stop, num=50):
        t = _torch()
        return t.exp(t.linspace(float(np.log(start)), float(np.log(stop)),
                                num, dtype=self.dtype, device=self.device))

    def eye(self, n, dtype=float):
        return _torch().eye(n, dtype=self._map_dtype(dtype),
                            device=self.device)

    def meshgrid(self, *arrays, indexing="xy"):
        return _torch().meshgrid(*arrays, indexing=indexing)

    def atleast_1d(self, x):
        return _torch().atleast_1d(self.asarray(x))

    def broadcast_to(self, x, shape):
        return _torch().broadcast_to(self.asarray(x), self._size(shape))

    # -- elementwise math ----------------------------------------------------
    def _two(self, a, b):
        """Coerce the second operand of a binary torch op."""
        t = _torch()
        if not isinstance(a, t.Tensor):
            a = self.asarray(a)
        if not isinstance(b, t.Tensor):
            b = t.as_tensor(b, dtype=a.dtype if a.dtype.is_floating_point
                            else self.dtype, device=self.device)
        return a, b

    def sqrt(self, x):
        return _torch().sqrt(self._result_float(self.asarray(x)))

    def exp(self, x):
        return _torch().exp(self._result_float(self.asarray(x)))

    def log(self, x):
        return _torch().log(self._result_float(self.asarray(x)))

    def log10(self, x):
        return _torch().log10(self._result_float(self.asarray(x)))

    def sin(self, x):
        return _torch().sin(self._result_float(self.asarray(x)))

    def cos(self, x):
        return _torch().cos(self._result_float(self.asarray(x)))

    def tan(self, x):
        return _torch().tan(self._result_float(self.asarray(x)))

    def tanh(self, x):
        return _torch().tanh(self._result_float(self.asarray(x)))

    def arccos(self, x):
        return _torch().arccos(self._result_float(self.asarray(x)))

    def arcsin(self, x):
        return _torch().arcsin(self._result_float(self.asarray(x)))

    def arctan(self, x):
        return _torch().arctan(self._result_float(self.asarray(x)))

    def arctan2(self, y, x):
        a, b = self._two(self.asarray(y), x)
        return _torch().arctan2(self._result_float(a), self._result_float(b))

    def hypot(self, a, b):
        x, y = self._two(self.asarray(a), b)
        return _torch().hypot(self._result_float(x), self._result_float(y))

    def abs(self, x):
        return _torch().abs(self.asarray(x))

    def sign(self, x):
        return _torch().sign(x)

    def floor(self, x):
        return _torch().floor(self._result_float(self.asarray(x)))

    def ceil(self, x):
        return _torch().ceil(self._result_float(self.asarray(x)))

    def cbrt(self, x):
        t = _torch()
        x = self._result_float(self.asarray(x))
        return t.sign(x) * t.pow(t.abs(x), 1.0 / 3.0)

    def deg2rad(self, x):
        return _torch().deg2rad(self._result_float(self.asarray(x)))

    def rad2deg(self, x):
        return _torch().rad2deg(self._result_float(self.asarray(x)))

    degrees = rad2deg

    def real(self, x):
        return x.real if x.is_complex() else x

    def imag(self, x):
        return x.imag if x.is_complex() else _torch().zeros_like(x)

    def conj(self, x):
        return _torch().conj(x)

    def angle(self, x):
        return _torch().angle(x)

    def maximum(self, a, b):
        a, b = self._two(a, b)
        return _torch().maximum(a, b)

    def minimum(self, a, b):
        a, b = self._two(a, b)
        return _torch().minimum(a, b)

    def clip(self, x, lo, hi):
        return _torch().clamp(self.asarray(x), min=lo, max=hi)

    def mod(self, a, b):
        a, b = self._two(a, b)
        return _torch().remainder(a, b)

    def where(self, cond, a=None, b=None):
        t = _torch()
        if a is None:
            return t.nonzero(cond, as_tuple=True)
        if not isinstance(a, t.Tensor) and not isinstance(b, t.Tensor):
            a = self.asarray(a)
        return t.where(cond, a, b)

    def isfinite(self, x):
        return _torch().isfinite(x)

    def isnan(self, x):
        return _torch().isnan(x)

    # -- reductions ----------------------------------------------------------
    def _red(self, fn, x, axis=None, keepdims=False):
        if axis is None:
            out = fn(x)
            return out.reshape((1,) * x.dim()) if keepdims else out
        return fn(x, dim=axis, keepdim=keepdims)

    def sum(self, x, axis=None, keepdims=False):
        return self._red(_torch().sum, self.asarray(x), axis, keepdims)

    def mean(self, x, axis=None, keepdims=False):
        return self._red(_torch().mean, self._result_float(self.asarray(x)),
                         axis, keepdims)

    def std(self, x, axis=None, keepdims=False):
        t = _torch()
        x = self._result_float(self.asarray(x))
        if axis is None:
            out = t.std(x, correction=0)
            return out.reshape((1,) * x.dim()) if keepdims else out
        return t.std(x, dim=axis, keepdim=keepdims, correction=0)

    def var(self, x, axis=None, keepdims=False):
        t = _torch()
        x = self._result_float(self.asarray(x))
        if axis is None:
            out = t.var(x, correction=0)
            return out.reshape((1,) * x.dim()) if keepdims else out
        return t.var(x, dim=axis, keepdim=keepdims, correction=0)

    def min(self, x, axis=None, keepdims=False):
        return self._red(_torch().amin, x, axis, keepdims)

    def max(self, x, axis=None, keepdims=False):
        return self._red(_torch().amax, x, axis, keepdims)

    def nansum(self, x, axis=None, keepdims=False):
        return self._red(_torch().nansum, x, axis, keepdims)

    def nanmean(self, x, axis=None, keepdims=False):
        return self._red(_torch().nanmean, x, axis, keepdims)

    def any(self, x, axis=None):
        t = _torch()
        return t.any(x) if axis is None else t.any(x, dim=axis)

    def all(self, x, axis=None):
        t = _torch()
        return t.all(x) if axis is None else t.all(x, dim=axis)

    def argmax(self, x, axis=None):
        return _torch().argmax(x) if axis is None \
            else _torch().argmax(x, dim=axis)

    def argsort(self, x, axis=-1):
        return _torch().argsort(x, dim=axis)

    def cumsum(self, x, axis=None):
        if axis is None:
            return _torch().cumsum(self.asarray(x).reshape(-1), dim=0)
        return _torch().cumsum(x, dim=axis)

    def diff(self, x, axis=-1):
        return _torch().diff(x, dim=axis)

    def median(self, x):
        return _torch().quantile(self._result_float(self.asarray(x)), 0.5)

    def percentile(self, x, q):
        x = self._result_float(self.asarray(x)).reshape(-1)
        return _torch().quantile(x, float(q) / 100.0)

    def quantile(self, x, q):
        x = self._result_float(self.asarray(x)).reshape(-1)
        return _torch().quantile(x, float(q))

    def nanquantile(self, x, q):
        x = self._result_float(self.asarray(x)).reshape(-1)
        return _torch().nanquantile(x, float(q))

    # -- structure -----------------------------------------------------------
    def stack(self, arrays, axis=0):
        return _torch().stack([self.asarray(a) for a in arrays], dim=axis)

    def concatenate(self, arrays, axis=0):
        return _torch().cat([self.asarray(a) for a in arrays], dim=axis)

    def roll(self, x, shift, axis=None):
        return _torch().roll(x, shift, dims=axis)

    def outer(self, a, b):
        return _torch().outer(self.asarray(a), self.asarray(b))

    def cross(self, a, b, axis=-1):
        a, b = self._two(self.asarray(a), b)
        a, b = _torch().broadcast_tensors(a, b)
        return _torch().linalg.cross(a, b, dim=axis)

    def einsum(self, eq, *ops):
        return _torch().einsum(eq, *[self.asarray(o) for o in ops])

    def searchsorted(self, a, v, side="left"):
        return _torch().searchsorted(self.asarray(a), self.asarray(v),
                                     right=(side == "right"))

    def digitize(self, x, bins, right=False):
        # torch.bucketize's `right` is the inverse of np.digitize's
        return _torch().bucketize(self.asarray(x), self.asarray(bins),
                                  right=not right)

    def bincount(self, x, minlength=0):
        return _torch().bincount(x, minlength=minlength)

    def flatnonzero(self, x):
        return _torch().nonzero(self.asarray(x).reshape(-1),
                                as_tuple=False).reshape(-1)

    def interp(self, x, xp_pts, fp):
        t = _torch()
        x = self._result_float(self.asarray(x))
        xp_pts = self._result_float(self.asarray(xp_pts))
        fp = self.asarray(fp)
        i = t.clamp(t.searchsorted(xp_pts, x, right=True), 1,
                    xp_pts.numel() - 1)
        x0, x1 = xp_pts[i - 1], xp_pts[i]
        f0, f1 = fp[i - 1], fp[i]
        w = (x - x0) / (x1 - x0)
        out = f0 + w * (f1 - f0)
        return t.where(x <= xp_pts[0], fp[0],
                       t.where(x >= xp_pts[-1], fp[-1], out))

    def gradient(self, f, *spacing, axis=None):
        t = _torch()
        f = self._result_float(self.asarray(f))
        sp = spacing if spacing else (1.0,)
        out = t.gradient(f, spacing=sp[0] if len(sp) == 1 else list(sp),
                         dim=axis)
        if axis is not None and isinstance(axis, int):
            return out[0]
        return out[0] if f.dim() == 1 else list(out)

    def flip(self, x, axis=None):
        t = _torch()
        if axis is None:
            axis = tuple(range(self.asarray(x).dim()))
        if isinstance(axis, int):
            axis = (axis,)
        return t.flip(self.asarray(x), dims=axis)

    def trapezoid(self, y, x=None, axis=-1):
        if x is None:
            return _torch().trapezoid(y, dim=axis)
        return _torch().trapezoid(y, x=self.asarray(x), dim=axis)

    def index_add(self, dest, index, src):
        return dest.index_add_(0, index, self.asarray(src, dtype=dest.dtype))

    def erfc(self, x):
        return _torch().special.erfc(self._result_float(self.asarray(x)))

    @staticmethod
    def errstate(**kwargs):
        return np.errstate(**kwargs)

    def ndim(self, x):
        return _torch().as_tensor(x).dim() if isinstance(x, _torch().Tensor) \
            else np.ndim(x)


_SHIM_CACHE: dict = {}


def _torch_shim(device, dtype):
    t = _torch()
    if dtype is None:
        dtype = t.float32 if t.device(device).type == "cuda" else t.float64
    elif not isinstance(dtype, t.dtype):
        dtype = {"float32": t.float32, "float64": t.float64,
                 np.float32: t.float32, np.float64: t.float64,
                 float: t.float64}[dtype]
    key = (str(t.device(device)), dtype)
    if key not in _SHIM_CACHE:
        _SHIM_CACHE[key] = _TorchShim(device, dtype)
    return _SHIM_CACHE[key]


# ---------------------------------------------------------------------------
# dispatch
# ---------------------------------------------------------------------------

def xp_of(*arrays):
    """Array namespace for the given arrays: the torch shim of the first
    tensor argument (its device; float dtype of the first floating
    tensor), else the numpy shim.

    Fast path: a torch tensor cannot exist unless `torch` is imported, so
    when it is absent from sys.modules we skip straight to numpy without
    type inspection (this dispatch runs on nearly every array op)."""
    if _TORCH is not None or "torch" in sys.modules:
        try:
            t = _torch()
        except ImportError:
            return NUMPY_XP
        first = None
        dtype = None
        for a in arrays:
            if isinstance(a, t.Tensor):
                if first is None:
                    first = a
                if dtype is None and a.dtype.is_floating_point:
                    dtype = a.dtype
                elif dtype is None and a.dtype.is_complex:
                    dtype = (t.float32 if a.dtype == t.complex64
                             else t.float64)
        if first is not None:
            return _torch_shim(first.device, dtype or t.float64)
    return NUMPY_XP


def get_xp(backend=None, device=None, dtype=None):
    """Creation namespace: backend 'numpy' (default) or 'torch'.
    For torch, dtype defaults to float32 on CUDA devices and float64 on
    the CPU (consumer GPUs run fp64 at a fraction of fp32 throughput)."""
    if backend is None or backend == "numpy":
        return NUMPY_XP
    if backend == "torch":
        return _torch_shim(device if device is not None else "cpu", dtype)
    raise ValueError(f"unknown backend: {backend!r}")


def to_numpy(x):
    """Copy/convert an array (or nested tuple/list) to numpy."""
    if isinstance(x, (tuple, list)):
        return type(x)(to_numpy(v) for v in x)
    if _TORCH is not None and isinstance(x, _TORCH.Tensor):
        return x.detach().cpu().numpy()
    if "torch" in str(type(x)):
        return x.detach().cpu().numpy()
    return np.asarray(x)


# ---------------------------------------------------------------------------
# random number generation
# ---------------------------------------------------------------------------

class TorchRNG:
    """numpy.random.Generator-like wrapper over torch.Generator, drawing
    on the shim's device in the shim's float dtype."""

    def __init__(self, seed=None, xp=None):
        t = _torch()
        self.xp = xp if xp is not None else _torch_shim("cpu", None)
        self.gen = t.Generator(device=self.xp.device)
        if seed is not None:
            self.gen.manual_seed(int(seed))

    def _shape(self, size):
        if size is None:
            return ()
        if isinstance(size, (int, np.integer)):
            return (int(size),)
        return tuple(int(s) for s in size)

    def standard_normal(self, size=None):
        return _torch().randn(self._shape(size), generator=self.gen,
                              dtype=self.xp.dtype, device=self.xp.device)

    def normal(self, loc=0.0, scale=1.0, size=None):
        out = self.standard_normal(size)
        return out * scale + loc

    def random(self, size=None):
        return _torch().rand(self._shape(size), generator=self.gen,
                             dtype=self.xp.dtype, device=self.xp.device)

    def uniform(self, low=0.0, high=1.0, size=None):
        return self.random(size) * (high - low) + low

    def exponential(self, scale=1.0, size=None):
        u = self.random(size)
        return -scale * _torch().log1p(-u)


class _NumpyDrawRNG:
    """Draws with a numpy Generator, returns tensors on the shim device:
    bitwise-reproducible across backends (used by the parity tests)."""

    def __init__(self, np_rng, xp):
        self.np_rng = np_rng
        self.xp = xp

    def standard_normal(self, size=None):
        return self.xp.asarray(self.np_rng.standard_normal(size))

    def normal(self, loc=0.0, scale=1.0, size=None):
        return self.xp.asarray(self.np_rng.normal(loc, scale, size))

    def random(self, size=None):
        return self.xp.asarray(self.np_rng.random(size))

    def uniform(self, low=0.0, high=1.0, size=None):
        return self.xp.asarray(self.np_rng.uniform(low, high, size))

    def exponential(self, scale=1.0, size=None):
        return self.xp.asarray(self.np_rng.exponential(scale, size))


def default_rng(seed=None, xp=None):
    """RNG matching the backend of xp (numpy Generator or TorchRNG)."""
    if xp is None or not getattr(xp, "is_torch", False):
        return np.random.default_rng(seed)
    return TorchRNG(seed, xp)


def adapt_rng(rng, xp):
    """Coerce rng to the backend of xp.  None -> fresh default;
    numpy Generator under torch -> numpy-drawing device wrapper."""
    if rng is None:
        return default_rng(None, xp)
    if not getattr(xp, "is_torch", False):
        if isinstance(rng, np.random.Generator):
            return rng
        raise TypeError("numpy backend requires a numpy Generator; got "
                        f"{type(rng).__name__}")
    if isinstance(rng, TorchRNG):
        return rng
    if isinstance(rng, np.random.Generator):
        return _NumpyDrawRNG(rng, xp)
    return rng
