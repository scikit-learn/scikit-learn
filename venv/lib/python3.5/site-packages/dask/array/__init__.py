from __future__ import absolute_import, division, print_function

from ..utils import ignoring
from .core import (Array, block, concatenate, stack, from_array, store,
                   map_blocks, atop, to_hdf5, to_npy_stack, from_npy_stack,
                   from_delayed, asarray, asanyarray,
                   broadcast_arrays, broadcast_to)
from .routines import (take, choose, argwhere, where, coarsen, insert,
                       ravel, roll, unique, squeeze, ptp, diff, ediff1d,
                       gradient, bincount, digitize, histogram, cov, array,
                       dstack, vstack, hstack, compress, extract, round,
                       count_nonzero, flatnonzero, nonzero, around, isin,
                       isnull, notnull, isclose, allclose, corrcoef, swapaxes,
                       tensordot, transpose, dot, vdot, matmul,
                       apply_along_axis, apply_over_axes, result_type,
                       atleast_1d, atleast_2d, atleast_3d, piecewise, flip,
                       flipud, fliplr, einsum)
from .reshape import reshape
from .ufunc import (add, subtract, multiply, divide, logaddexp, logaddexp2,
        true_divide, floor_divide, negative, power, remainder, mod, conj, exp,
        exp2, log, log2, log10, log1p, expm1, sqrt, square, cbrt, reciprocal,
        sin, cos, tan, arcsin, arccos, arctan, arctan2, hypot, sinh, cosh,
        tanh, arcsinh, arccosh, arctanh, deg2rad, rad2deg, greater,
        greater_equal, less, less_equal, not_equal, equal, maximum,
        bitwise_and, bitwise_or, bitwise_xor, bitwise_not, minimum,
        logical_and, logical_or, logical_xor, logical_not, fmax, fmin,
        isreal, iscomplex, isfinite, isinf, isnan, signbit, copysign,
        nextafter, spacing, ldexp, fmod, floor, ceil, trunc, degrees,
        radians, rint, fix, angle, real, imag, clip, fabs, sign, absolute,
        i0, sinc, nan_to_num, frexp, modf, divide, frompyfunc)
try:
    from .ufunc import float_power
except ImportError:
    # Absent for NumPy versions prior to 1.12.
    pass
from .reductions import (sum, prod, mean, std, var, any, all, min, max, vnorm,
                         moment,
                         argmin, argmax,
                         nansum, nanmean, nanstd, nanvar, nanmin,
                         nanmax, nanargmin, nanargmax,
                         cumsum, cumprod,
                         topk, argtopk)
from .percentile import percentile
with ignoring(ImportError):
    from .reductions import nanprod, nancumprod, nancumsum
with ignoring(ImportError):
    from . import ma
from . import random, linalg, ghost, learn, fft
from .ghost import map_overlap
from .wrap import ones, zeros, empty, full
from .creation import ones_like, zeros_like, empty_like, full_like
from .rechunk import rechunk
from ..context import set_options
from ..base import compute
from .optimization import optimize
from .creation import (arange, linspace, meshgrid, indices, diag, eye,
                       triu, tril, fromfunction, tile, repeat)
