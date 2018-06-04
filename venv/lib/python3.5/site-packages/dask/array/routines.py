from __future__ import division, print_function, absolute_import

import inspect
import math
import warnings
from collections import Iterable
from distutils.version import LooseVersion
from functools import wraps, partial
from numbers import Number, Real, Integral

import numpy as np
from toolz import concat, merge, sliding_window, interleave

from .. import sharedict
from ..core import flatten
from ..base import tokenize
from . import chunk
from .creation import arange
from .utils import safe_wraps
from .wrap import ones

from .core import (Array, map_blocks, elemwise, from_array, asarray,
                   asanyarray, concatenate, stack, atop, broadcast_shapes,
                   is_scalar_for_elemwise, broadcast_to, tensordot_lookup)

from .einsumfuncs import einsum  # noqa


@wraps(np.array)
def array(x, dtype=None, ndmin=None):
    while ndmin is not None and x.ndim < ndmin:
        x = x[None, :]
    if dtype is not None and x.dtype != dtype:
        x = x.astype(dtype)
    return x


@wraps(np.result_type)
def result_type(*args):
    args = [a if is_scalar_for_elemwise(a) else a.dtype for a in args]
    return np.result_type(*args)


@wraps(np.atleast_3d)
def atleast_3d(*arys):
    new_arys = []
    for x in arys:
        x = asanyarray(x)
        if x.ndim == 0:
            x = x[None, None, None]
        elif x.ndim == 1:
            x = x[None, :, None]
        elif x.ndim == 2:
            x = x[:, :, None]

        new_arys.append(x)

    if len(new_arys) == 1:
        return new_arys[0]
    else:
        return new_arys


@wraps(np.atleast_2d)
def atleast_2d(*arys):
    new_arys = []
    for x in arys:
        x = asanyarray(x)
        if x.ndim == 0:
            x = x[None, None]
        elif x.ndim == 1:
            x = x[None, :]

        new_arys.append(x)

    if len(new_arys) == 1:
        return new_arys[0]
    else:
        return new_arys


@wraps(np.atleast_1d)
def atleast_1d(*arys):
    new_arys = []
    for x in arys:
        x = asanyarray(x)
        if x.ndim == 0:
            x = x[None]

        new_arys.append(x)

    if len(new_arys) == 1:
        return new_arys[0]
    else:
        return new_arys


@wraps(np.vstack)
def vstack(tup):
    tup = tuple(atleast_2d(x) for x in tup)
    return concatenate(tup, axis=0)


@wraps(np.hstack)
def hstack(tup):
    if all(x.ndim == 1 for x in tup):
        return concatenate(tup, axis=0)
    else:
        return concatenate(tup, axis=1)


@wraps(np.dstack)
def dstack(tup):
    tup = tuple(atleast_3d(x) for x in tup)
    return concatenate(tup, axis=2)


@wraps(np.swapaxes)
def swapaxes(a, axis1, axis2):
    if axis1 == axis2:
        return a
    if axis1 < 0:
        axis1 = axis1 + a.ndim
    if axis2 < 0:
        axis2 = axis2 + a.ndim
    ind = list(range(a.ndim))
    out = list(ind)
    out[axis1], out[axis2] = axis2, axis1

    return atop(np.swapaxes, out, a, ind, axis1=axis1, axis2=axis2,
                dtype=a.dtype)


@wraps(np.transpose)
def transpose(a, axes=None):
    if axes:
        if len(axes) != a.ndim:
            raise ValueError("axes don't match array")
    else:
        axes = tuple(range(a.ndim))[::-1]
    axes = tuple(d + a.ndim if d < 0 else d for d in axes)
    return atop(np.transpose, axes, a, tuple(range(a.ndim)),
                dtype=a.dtype, axes=axes)


def flip(m, axis):
    """
    Reverse element order along axis.

    Parameters
    ----------
    axis : int
        Axis to reverse element order of.

    Returns
    -------
    reversed array : ndarray
    """

    m = asanyarray(m)

    sl = m.ndim * [slice(None)]
    try:
        sl[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError(
            "`axis` of %s invalid for %s-D array" % (str(axis), str(m.ndim))
        )
    sl = tuple(sl)

    return m[sl]


@wraps(np.flipud)
def flipud(m):
    return flip(m, 0)


@wraps(np.fliplr)
def fliplr(m):
    return flip(m, 1)


alphabet = 'abcdefghijklmnopqrstuvwxyz'
ALPHABET = alphabet.upper()


def _tensordot(a, b, axes):
    x = max([a, b], key=lambda x: x.__array_priority__)
    tensordot = tensordot_lookup.dispatch(type(x))
    x = tensordot(a, b, axes=axes)
    ind = [slice(None, None)] * x.ndim
    for a in sorted(axes[0]):
        ind.insert(a, None)
    x = x[tuple(ind)]
    return x


@wraps(np.tensordot)
def tensordot(lhs, rhs, axes=2):
    if isinstance(axes, Iterable):
        left_axes, right_axes = axes
    else:
        left_axes = tuple(range(lhs.ndim - 1, lhs.ndim - axes - 1, -1))
        right_axes = tuple(range(0, axes))

    if isinstance(left_axes, int):
        left_axes = (left_axes,)
    if isinstance(right_axes, int):
        right_axes = (right_axes,)
    if isinstance(left_axes, list):
        left_axes = tuple(left_axes)
    if isinstance(right_axes, list):
        right_axes = tuple(right_axes)

    dt = np.promote_types(lhs.dtype, rhs.dtype)

    left_index = list(alphabet[:lhs.ndim])
    right_index = list(ALPHABET[:rhs.ndim])
    out_index = left_index + right_index

    for l, r in zip(left_axes, right_axes):
        out_index.remove(right_index[r])
        right_index[r] = left_index[l]

    intermediate = atop(_tensordot, out_index,
                        lhs, left_index,
                        rhs, right_index, dtype=dt,
                        axes=(left_axes, right_axes))

    result = intermediate.sum(axis=left_axes)
    return result


@wraps(np.dot)
def dot(a, b):
    return tensordot(a, b, axes=((a.ndim - 1,), (b.ndim - 2,)))


@wraps(np.vdot)
def vdot(a, b):
    return dot(a.conj().ravel(), b.ravel())


@wraps(np.matmul)
def matmul(a, b):
    a = asanyarray(a)
    b = asanyarray(b)

    if a.ndim == 0 or b.ndim == 0:
        raise ValueError("`matmul` does not support scalars.")

    a_is_1d = False
    if a.ndim == 1:
        a_is_1d = True
        a = a[np.newaxis, :]

    b_is_1d = False
    if b.ndim == 1:
        b_is_1d = True
        b = b[:, np.newaxis]

    if a.ndim < b.ndim:
        a = a[(b.ndim - a.ndim) * (np.newaxis,)]
    elif a.ndim > b.ndim:
        b = b[(a.ndim - b.ndim) * (np.newaxis,)]

    out = atop(
        np.matmul, tuple(range(1, a.ndim + 1)),
        a, tuple(range(1, a.ndim - 1)) + (a.ndim - 1, 0,),
        b, tuple(range(1, a.ndim - 1)) + (0, a.ndim,),
        dtype=result_type(a, b),
        concatenate=True
    )

    if a_is_1d:
        out = out[..., 0, :]
    if b_is_1d:
        out = out[..., 0]

    return out


def _inner_apply_along_axis(arr,
                            func1d,
                            func1d_axis,
                            func1d_args,
                            func1d_kwargs):
    return np.apply_along_axis(
        func1d, func1d_axis, arr, *func1d_args, **func1d_kwargs
    )


@wraps(np.apply_along_axis)
def apply_along_axis(func1d, axis, arr, *args, **kwargs):
    arr = asarray(arr)

    # Validate and normalize axis.
    arr.shape[axis]
    axis = len(arr.shape[:axis])

    # Test out some data with the function.
    test_data = np.ones((1,), dtype=arr.dtype)
    test_result = np.array(func1d(test_data, *args, **kwargs))

    if (LooseVersion(np.__version__) < LooseVersion("1.13.0") and
            (np.array(test_result.shape) > 1).sum(dtype=int) > 1):
        raise ValueError(
            "No more than one non-trivial dimension allowed in result. "
            "Need NumPy 1.13.0+ for this functionality."
        )

    # Rechunk so that func1d is applied over the full axis.
    arr = arr.rechunk(
        arr.chunks[:axis] + (arr.shape[axis:axis + 1],) + arr.chunks[axis + 1:]
    )

    # Map func1d over the data to get the result
    # Adds other axes as needed.
    result = arr.map_blocks(
        _inner_apply_along_axis,
        token="apply_along_axis",
        dtype=test_result.dtype,
        chunks=(arr.chunks[:axis] + test_result.shape + arr.chunks[axis + 1:]),
        drop_axis=axis,
        new_axis=list(range(axis, axis + test_result.ndim, 1)),
        func1d=func1d,
        func1d_axis=axis,
        func1d_args=args,
        func1d_kwargs=kwargs,
    )

    return result


@wraps(np.apply_over_axes)
def apply_over_axes(func, a, axes):
    # Validate arguments
    a = asarray(a)
    try:
        axes = tuple(axes)
    except TypeError:
        axes = (axes,)

    sl = a.ndim * (slice(None),)

    # Compute using `apply_along_axis`.
    result = a
    for i in axes:
        result = apply_along_axis(func, i, result, 0)

        # Restore original dimensionality or error.
        if result.ndim == (a.ndim - 1):
            result = result[sl[:i] + (None,)]
        elif result.ndim != a.ndim:
            raise ValueError(
                "func must either preserve dimensionality of the input"
                " or reduce it by one."
            )

    return result


@wraps(np.ptp)
def ptp(a, axis=None):
    return a.max(axis=axis) - a.min(axis=axis)


@wraps(np.diff)
def diff(a, n=1, axis=-1):
    a = asarray(a)
    n = int(n)
    axis = int(axis)

    sl_1 = a.ndim * [slice(None)]
    sl_2 = a.ndim * [slice(None)]

    sl_1[axis] = slice(1, None)
    sl_2[axis] = slice(None, -1)

    sl_1 = tuple(sl_1)
    sl_2 = tuple(sl_2)

    r = a
    for i in range(n):
        r = r[sl_1] - r[sl_2]

    return r


@wraps(np.ediff1d)
def ediff1d(ary, to_end=None, to_begin=None):
    ary = asarray(ary)

    aryf = ary.flatten()
    r = aryf[1:] - aryf[:-1]

    r = [r]
    if to_begin is not None:
        r = [asarray(to_begin).flatten()] + r
    if to_end is not None:
        r = r + [asarray(to_end).flatten()]
    r = concatenate(r)

    return r


def _gradient_kernel(f, grad_varargs, grad_kwargs):
    return np.gradient(f, *grad_varargs, **grad_kwargs)


@wraps(np.gradient)
def gradient(f, *varargs, **kwargs):
    f = asarray(f)

    if not all([isinstance(e, Number) for e in varargs]):
        raise NotImplementedError("Only numeric scalar spacings supported.")

    if varargs == ():
        varargs = (1,)
    if len(varargs) == 1:
        varargs = f.ndim * varargs
    if len(varargs) != f.ndim:
        raise TypeError(
            "Spacing must either be a scalar or a scalar per dimension."
        )

    kwargs["edge_order"] = math.ceil(kwargs.get("edge_order", 1))
    if kwargs["edge_order"] > 2:
        raise ValueError("edge_order must be less than or equal to 2.")

    drop_result_list = False
    axis = kwargs.pop("axis", None)
    if axis is None:
        axis = tuple(range(f.ndim))
    elif isinstance(axis, Integral):
        drop_result_list = True
        axis = (axis,)

    for e in axis:
        if not isinstance(e, Integral):
            raise TypeError("%s, invalid value for axis" % repr(e))
        if not (-f.ndim <= e < f.ndim):
            raise ValueError("axis, %s, is out of bounds" % repr(e))

    if len(axis) != len(set(axis)):
        raise ValueError("duplicate axes not allowed")

    axis = tuple(ax % f.ndim for ax in axis)

    if issubclass(f.dtype.type, (np.bool8, Integral)):
        f = f.astype(float)
    elif issubclass(f.dtype.type, Real) and f.dtype.itemsize < 4:
        f = f.astype(float)

    r = [
        f.map_overlap(
            _gradient_kernel,
            dtype=f.dtype,
            depth={j: 1 if j == ax else 0 for j in range(f.ndim)},
            boundary="none",
            grad_varargs=(varargs[i],),
            grad_kwargs=merge(kwargs, {"axis": ax}),
        )
        for i, ax in enumerate(axis)
    ]
    if drop_result_list:
        r = r[0]

    return r


@wraps(np.bincount)
def bincount(x, weights=None, minlength=None):
    if minlength is None:
        raise TypeError("Must specify minlength argument in da.bincount")
    assert x.ndim == 1
    if weights is not None:
        assert weights.chunks == x.chunks

    # Call np.bincount on each block, possibly with weights
    token = tokenize(x, weights, minlength)
    name = 'bincount-' + token
    if weights is not None:
        dsk = {(name, i): (np.bincount, (x.name, i), (weights.name, i), minlength)
               for i, _ in enumerate(x.__dask_keys__())}
        dtype = np.bincount([1], weights=[1]).dtype
    else:
        dsk = {(name, i): (np.bincount, (x.name, i), None, minlength)
               for i, _ in enumerate(x.__dask_keys__())}
        dtype = np.bincount([]).dtype

    # Sum up all of the intermediate bincounts per block
    name = 'bincount-sum-' + token
    dsk[(name, 0)] = (np.sum, list(dsk), 0)

    chunks = ((minlength,),)

    dsk = sharedict.merge((name, dsk), x.dask)
    if weights is not None:
        dsk.update(weights.dask)

    return Array(dsk, name, chunks, dtype)


@wraps(np.digitize)
def digitize(a, bins, right=False):
    bins = np.asarray(bins)
    dtype = np.digitize([0], bins, right=False).dtype
    return a.map_blocks(np.digitize, dtype=dtype, bins=bins, right=right)


def histogram(a, bins=None, range=None, normed=False, weights=None, density=None):
    """
    Blocked variant of :func:`numpy.histogram`.

    Follows the signature of :func:`numpy.histogram` exactly with the following
    exceptions:

    - Either an iterable specifying the ``bins`` or the number of ``bins``
      and a ``range`` argument is required as computing ``min`` and ``max``
      over blocked arrays is an expensive operation that must be performed
      explicitly.

    - ``weights`` must be a dask.array.Array with the same block structure
      as ``a``.

    Examples
    --------
    Using number of bins and range:

    >>> import dask.array as da
    >>> import numpy as np
    >>> x = da.from_array(np.arange(10000), chunks=10)
    >>> h, bins = da.histogram(x, bins=10, range=[0, 10000])
    >>> bins
    array([    0.,  1000.,  2000.,  3000.,  4000.,  5000.,  6000.,  7000.,
            8000.,  9000., 10000.])
    >>> h.compute()
    array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])

    Explicitly specifying the bins:

    >>> h, bins = da.histogram(x, bins=np.array([0, 5000, 10000]))
    >>> bins
    array([    0,  5000, 10000])
    >>> h.compute()
    array([5000, 5000])
    """
    if bins is None or (range is None and bins is None):
        raise ValueError('dask.array.histogram requires either bins '
                         'or bins and range to be defined.')

    if weights is not None and weights.chunks != a.chunks:
        raise ValueError('Input array and weights must have the same '
                         'chunked structure')

    if not np.iterable(bins):
        bin_token = bins
        mn, mx = range
        if mn == mx:
            mn -= 0.5
            mx += 0.5

        bins = np.linspace(mn, mx, bins + 1, endpoint=True)
    else:
        bin_token = bins
    token = tokenize(a, bin_token, range, normed, weights, density)

    nchunks = len(list(flatten(a.__dask_keys__())))
    chunks = ((1,) * nchunks, (len(bins) - 1,))

    name = 'histogram-sum-' + token

    # Map the histogram to all bins
    def block_hist(x, weights=None):
        return np.histogram(x, bins, weights=weights)[0][np.newaxis]

    if weights is None:
        dsk = {(name, i, 0): (block_hist, k)
               for i, k in enumerate(flatten(a.__dask_keys__()))}
        dtype = np.histogram([])[0].dtype
    else:
        a_keys = flatten(a.__dask_keys__())
        w_keys = flatten(weights.__dask_keys__())
        dsk = {(name, i, 0): (block_hist, k, w)
               for i, (k, w) in enumerate(zip(a_keys, w_keys))}
        dtype = weights.dtype

    all_dsk = sharedict.merge(a.dask, (name, dsk))
    if weights is not None:
        all_dsk.update(weights.dask)

    mapped = Array(all_dsk, name, chunks, dtype=dtype)
    n = mapped.sum(axis=0)

    # We need to replicate normed and density options from numpy
    if density is not None:
        if density:
            db = from_array(np.diff(bins).astype(float), chunks=n.chunks)
            return n / db / n.sum(), bins
        else:
            return n, bins
    else:
        # deprecated, will be removed from Numpy 2.0
        if normed:
            db = from_array(np.diff(bins).astype(float), chunks=n.chunks)
            return n / (n * db).sum(), bins
        else:
            return n, bins


@wraps(np.cov)
def cov(m, y=None, rowvar=1, bias=0, ddof=None):
    # This was copied almost verbatim from np.cov
    # See numpy license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory
    if ddof is not None and ddof != int(ddof):
        raise ValueError(
            "ddof must be integer")

    # Handles complex arrays too
    m = asarray(m)
    if y is None:
        dtype = np.result_type(m, np.float64)
    else:
        y = asarray(y)
        dtype = np.result_type(m, y, np.float64)
    X = array(m, ndmin=2, dtype=dtype)

    if X.shape[0] == 1:
        rowvar = 1
    if rowvar:
        N = X.shape[1]
        axis = 0
    else:
        N = X.shape[0]
        axis = 1

    # check ddof
    if ddof is None:
        if bias == 0:
            ddof = 1
        else:
            ddof = 0
    fact = float(N - ddof)
    if fact <= 0:
        warnings.warn("Degrees of freedom <= 0 for slice", RuntimeWarning)
        fact = 0.0

    if y is not None:
        y = array(y, ndmin=2, dtype=dtype)
        X = concatenate((X, y), axis)

    X = X - X.mean(axis=1 - axis, keepdims=True)
    if not rowvar:
        return (dot(X.T, X.conj()) / fact).squeeze()
    else:
        return (dot(X, X.T.conj()) / fact).squeeze()


@wraps(np.corrcoef)
def corrcoef(x, y=None, rowvar=1):

    from .ufunc import sqrt
    from .creation import diag

    c = cov(x, y, rowvar)
    if c.shape == ():
        return c / c
    d = diag(c)
    d = d.reshape((d.shape[0], 1))
    sqr_d = sqrt(d)
    return (c / sqr_d) / sqr_d.T


@wraps(np.round)
def round(a, decimals=0):
    return a.map_blocks(np.round, decimals=decimals, dtype=a.dtype)


def _unique_internal(ar, indices, counts, return_inverse=False):
    """
    Helper/wrapper function for :func:`numpy.unique`.

    Uses :func:`numpy.unique` to find the unique values for the array chunk.
    Given this chunk may not represent the whole array, also take the
    ``indices`` and ``counts`` that are in 1-to-1 correspondence to ``ar``
    and reduce them in the same fashion as ``ar`` is reduced. Namely sum
    any counts that correspond to the same value and take the smallest
    index that corresponds to the same value.

    To handle the inverse mapping from the unique values to the original
    array, simply return a NumPy array created with ``arange`` with enough
    values to correspond 1-to-1 to the unique values. While there is more
    work needed to be done to create the full inverse mapping for the
    original array, this provides enough information to generate the
    inverse mapping in Dask.

    Given Dask likes to have one array returned from functions like
    ``atop``, some formatting is done to stuff all of the resulting arrays
    into one big NumPy structured array. Dask is then able to handle this
    object and can split it apart into the separate results on the Dask
    side, which then can be passed back to this function in concatenated
    chunks for further reduction or can be return to the user to perform
    other forms of analysis.

    By handling the problem in this way, it does not matter where a chunk
    is in a larger array or how big it is. The chunk can still be computed
    on the same way. Also it does not matter if the chunk is the result of
    other chunks being run through this function multiple times. The end
    result will still be just as accurate using this strategy.
    """

    return_index = (indices is not None)
    return_counts = (counts is not None)

    u = np.unique(ar)

    dt = [("values", u.dtype)]
    if return_index:
        dt.append(("indices", np.intp))
    if return_inverse:
        dt.append(("inverse", np.intp))
    if return_counts:
        dt.append(("counts", np.intp))

    r = np.empty(u.shape, dtype=dt)
    r["values"] = u
    if return_inverse:
        r["inverse"] = np.arange(len(r), dtype=np.intp)
    if return_index or return_counts:
        for i, v in enumerate(r["values"]):
            m = (ar == v)
            if return_index:
                indices[m].min(keepdims=True, out=r["indices"][i:i + 1])
            if return_counts:
                counts[m].sum(keepdims=True, out=r["counts"][i:i + 1])

    return r


@wraps(np.unique)
def unique(ar, return_index=False, return_inverse=False, return_counts=False):
    ar = ar.ravel()

    # Run unique on each chunk and collect results in a Dask Array of
    # unknown size.

    args = [ar, "i"]
    out_dtype = [("values", ar.dtype)]
    if return_index:
        args.extend([
            arange(ar.shape[0], dtype=np.intp, chunks=ar.chunks[0]),
            "i"
        ])
        out_dtype.append(("indices", np.intp))
    else:
        args.extend([None, None])
    if return_counts:
        args.extend([
            ones((ar.shape[0],), dtype=np.intp, chunks=ar.chunks[0]),
            "i"
        ])
        out_dtype.append(("counts", np.intp))
    else:
        args.extend([None, None])

    out = atop(
        _unique_internal, "i",
        *args,
        dtype=out_dtype,
        return_inverse=False
    )
    out._chunks = tuple((np.nan,) * len(c) for c in out.chunks)

    # Take the results from the unique chunks and do the following.
    #
    # 1. Collect all results as arguments.
    # 2. Concatenate each result into one big array.
    # 3. Pass all results as arguments to the internal unique again.
    #
    # TODO: This should be replaced with a tree reduction using this strategy.
    # xref: https://github.com/dask/dask/issues/2851

    out_parts = [out["values"]]
    if return_index:
        out_parts.append(out["indices"])
    else:
        out_parts.append(None)
    if return_counts:
        out_parts.append(out["counts"])
    else:
        out_parts.append(None)

    name = 'unique-aggregate-' + out.name
    dsk = {
        (name, 0): (
            (_unique_internal,) +
            tuple(
                (np.concatenate, o. __dask_keys__())
                if hasattr(o, "__dask_keys__") else o
                for o in out_parts
            ) +
            (return_inverse,)
        )
    }
    out_dtype = [("values", ar.dtype)]
    if return_index:
        out_dtype.append(("indices", np.intp))
    if return_inverse:
        out_dtype.append(("inverse", np.intp))
    if return_counts:
        out_dtype.append(("counts", np.intp))

    out = Array(
        sharedict.merge(*(
            [(name, dsk)] +
            [o.dask for o in out_parts if hasattr(o, "__dask_keys__")]
        )),
        name,
        ((np.nan,),),
        out_dtype
    )

    # Split out all results to return to the user.

    result = [out["values"]]
    if return_index:
        result.append(out["indices"])
    if return_inverse:
        # Using the returned unique values and arange of unknown length, find
        # each value matching a unique value and replace it with its
        # corresponding index or `0`. There should be only one entry for this
        # index in axis `1` (the one of unknown length). Reduce axis `1`
        # through summing to get an array with known dimensionality and the
        # mapping of the original values.
        mtches = (ar[:, None] == out["values"][None, :]).astype(np.intp)
        result.append((mtches * out["inverse"]).sum(axis=1))
    if return_counts:
        result.append(out["counts"])

    if len(result) == 1:
        result = result[0]
    else:
        result = tuple(result)

    return result


def _isin_kernel(element, test_elements, assume_unique=False):
    values = np.in1d(element.ravel(), test_elements,
                     assume_unique=assume_unique)
    return values.reshape(element.shape + (1,) * test_elements.ndim)


@safe_wraps(getattr(np, 'isin', None))
def isin(element, test_elements, assume_unique=False, invert=False):
    element = asarray(element)
    test_elements = asarray(test_elements)
    element_axes = tuple(range(element.ndim))
    test_axes = tuple(i + element.ndim for i in range(test_elements.ndim))
    mapped = atop(_isin_kernel, element_axes + test_axes,
                  element, element_axes,
                  test_elements, test_axes,
                  adjust_chunks={axis: lambda _: 1 for axis in test_axes},
                  dtype=bool,
                  assume_unique=assume_unique)
    result = mapped.any(axis=test_axes)
    if invert:
        result = ~result
    return result


@wraps(np.roll)
def roll(array, shift, axis=None):
    result = array

    if axis is None:
        result = ravel(result)

        if not isinstance(shift, Integral):
            raise TypeError(
                "Expect `shift` to be an instance of Integral"
                " when `axis` is None."
            )

        shift = (shift,)
        axis = (0,)
    else:
        try:
            len(shift)
        except TypeError:
            shift = (shift,)
        try:
            len(axis)
        except TypeError:
            axis = (axis,)

    if len(shift) != len(axis):
        raise ValueError("Must have the same number of shifts as axes.")

    for i, s in zip(axis, shift):
        s = -s
        s %= result.shape[i]

        sl1 = result.ndim * [slice(None)]
        sl2 = result.ndim * [slice(None)]

        sl1[i] = slice(s, None)
        sl2[i] = slice(None, s)

        sl1 = tuple(sl1)
        sl2 = tuple(sl2)

        result = concatenate([result[sl1], result[sl2]], axis=i)

    result = result.reshape(array.shape)

    return result


@wraps(np.ravel)
def ravel(array):
    return array.reshape((-1,))


@wraps(np.squeeze)
def squeeze(a, axis=None):
    if axis is None:
        axis = tuple(i for i, d in enumerate(a.shape) if d == 1)
    elif not isinstance(axis, tuple):
        axis = (axis,)

    if any(a.shape[i] != 1 for i in axis):
        raise ValueError("cannot squeeze axis with size other than one")

    for i in axis:
        if not (-a.ndim <= i < a.ndim):
            raise ValueError("%i out of bounds for %i-D array" % (i, a.ndim))

    axis = tuple(i % a.ndim for i in axis)

    sl = tuple(0 if i in axis else slice(None) for i, s in enumerate(a.shape))

    return a[sl]


@wraps(np.compress)
def compress(condition, a, axis=None):
    if axis is None:
        a = a.ravel()
        axis = 0
    if not -a.ndim <= axis < a.ndim:
        raise ValueError('axis=(%s) out of bounds' % axis)
    if axis < 0:
        axis += a.ndim

    # Only coerce non-lazy values to numpy arrays
    if not isinstance(condition, Array):
        condition = np.array(condition, dtype=bool)
    if condition.ndim != 1:
        raise ValueError("Condition must be one dimensional")

    if isinstance(condition, Array):
        if len(condition) < a.shape[axis]:
            a = a[tuple(slice(None, len(condition))
                        if i == axis else slice(None)
                        for i in range(a.ndim))]
        inds = tuple(range(a.ndim))
        out = atop(np.compress, inds, condition, (inds[axis],), a, inds,
                   axis=axis, dtype=a.dtype)
        out._chunks = tuple((np.NaN,) * len(c) if i == axis else c
                            for i, c in enumerate(out.chunks))
        return out
    else:
        # Optimized case when condition is known
        if len(condition) < a.shape[axis]:
            condition = condition.copy()
            condition.resize(a.shape[axis])

        slc = ((slice(None),) * axis + (condition, ) +
               (slice(None),) * (a.ndim - axis - 1))
        return a[slc]


@wraps(np.extract)
def extract(condition, arr):
    if not isinstance(condition, Array):
        condition = np.array(condition, dtype=bool)
    return compress(condition.ravel(), arr.ravel())


@wraps(np.take)
def take(a, indices, axis=0):
    if not -a.ndim <= axis < a.ndim:
        raise ValueError('axis=(%s) out of bounds' % axis)
    if axis < 0:
        axis += a.ndim
    if isinstance(a, np.ndarray) and isinstance(indices, Array):
        return _take_dask_array_from_numpy(a, indices, axis)
    else:
        return a[(slice(None),) * axis + (indices,)]


def _take_dask_array_from_numpy(a, indices, axis):
    assert isinstance(a, np.ndarray)
    assert isinstance(indices, Array)

    return indices.map_blocks(lambda block: np.take(a, block, axis),
                              chunks=indices.chunks,
                              dtype=a.dtype)


@wraps(np.around)
def around(x, decimals=0):
    return map_blocks(partial(np.around, decimals=decimals), x, dtype=x.dtype)


def isnull(values):
    """ pandas.isnull for dask arrays """
    import pandas as pd
    return elemwise(pd.isnull, values, dtype='bool')


def notnull(values):
    """ pandas.notnull for dask arrays """
    return ~isnull(values)


@wraps(np.isclose)
def isclose(arr1, arr2, rtol=1e-5, atol=1e-8, equal_nan=False):
    func = partial(np.isclose, rtol=rtol, atol=atol, equal_nan=equal_nan)
    return elemwise(func, arr1, arr2, dtype='bool')


@wraps(np.allclose)
def allclose(arr1, arr2, rtol=1e-5, atol=1e-8, equal_nan=False):
    return isclose(arr1, arr2, rtol=rtol, atol=atol, equal_nan=equal_nan).all()


def variadic_choose(a, *choices):
    return np.choose(a, choices)


@wraps(np.choose)
def choose(a, choices):
    return elemwise(variadic_choose, a, *choices)


def _isnonzero_vec(v):
    return bool(np.count_nonzero(v))


_isnonzero_vec = np.vectorize(_isnonzero_vec, otypes=[bool])


def isnonzero(a):
    try:
        np.zeros(tuple(), dtype=a.dtype).astype(bool)
    except ValueError:
        ######################################################
        # Handle special cases where conversion to bool does #
        # not work correctly.                                #
        #                                                    #
        # xref: https://github.com/numpy/numpy/issues/9479   #
        ######################################################
        return a.map_blocks(_isnonzero_vec, dtype=bool)
    else:
        return a.astype(bool)


@wraps(np.argwhere)
def argwhere(a):
    from .creation import indices

    a = asarray(a)

    nz = isnonzero(a).flatten()

    ind = indices(a.shape, dtype=np.intp, chunks=a.chunks)
    if ind.ndim > 1:
        ind = stack([ind[i].ravel() for i in range(len(ind))], axis=1)
    ind = compress(nz, ind, axis=0)

    return ind


@wraps(np.where)
def where(condition, x=None, y=None):
    if (x is None) != (y is None):
        raise ValueError("either both or neither of x and y should be given")
    if (x is None) and (y is None):
        return nonzero(condition)

    if np.isscalar(condition):
        dtype = result_type(x, y)
        x = asarray(x)
        y = asarray(y)

        shape = broadcast_shapes(x.shape, y.shape)
        out = x if condition else y

        return broadcast_to(out, shape).astype(dtype)
    else:
        return elemwise(np.where, condition, x, y)


@wraps(np.count_nonzero)
def count_nonzero(a, axis=None):
    return isnonzero(asarray(a)).astype(np.intp).sum(axis=axis)


@wraps(np.flatnonzero)
def flatnonzero(a):
    return argwhere(asarray(a).ravel())[:, 0]


@wraps(np.nonzero)
def nonzero(a):
    ind = argwhere(a)
    if ind.ndim > 1:
        return tuple(ind[:, i] for i in range(ind.shape[1]))
    else:
        return (ind,)


def _int_piecewise(x, *condlist, **kwargs):
    return np.piecewise(
        x, list(condlist), kwargs["funclist"],
        *kwargs["func_args"], **kwargs["func_kw"]
    )


@wraps(np.piecewise)
def piecewise(x, condlist, funclist, *args, **kw):
    return map_blocks(
        _int_piecewise,
        x, *condlist,
        dtype=x.dtype,
        token="piecewise",
        funclist=funclist, func_args=args, func_kw=kw
    )


@wraps(chunk.coarsen)
def coarsen(reduction, x, axes, trim_excess=False):
    if (not trim_excess and
        not all(bd % div == 0 for i, div in axes.items()
                for bd in x.chunks[i])):
        msg = "Coarsening factor does not align with block dimensions"
        raise ValueError(msg)

    if 'dask' in inspect.getfile(reduction):
        reduction = getattr(np, reduction.__name__)

    name = 'coarsen-' + tokenize(reduction, x, axes, trim_excess)
    dsk = {(name,) + key[1:]: (chunk.coarsen, reduction, key, axes, trim_excess)
           for key in flatten(x.__dask_keys__())}
    chunks = tuple(tuple(int(bd // axes.get(i, 1)) for bd in bds)
                   for i, bds in enumerate(x.chunks))

    dt = reduction(np.empty((1,) * x.ndim, dtype=x.dtype)).dtype
    return Array(sharedict.merge(x.dask, (name, dsk)), name, chunks, dtype=dt)


def split_at_breaks(array, breaks, axis=0):
    """ Split an array into a list of arrays (using slices) at the given breaks

    >>> split_at_breaks(np.arange(6), [3, 5])
    [array([0, 1, 2]), array([3, 4]), array([5])]
    """
    padded_breaks = concat([[None], breaks, [None]])
    slices = [slice(i, j) for i, j in sliding_window(2, padded_breaks)]
    preslice = (slice(None),) * axis
    split_array = [array[preslice + (s,)] for s in slices]
    return split_array


@wraps(np.insert)
def insert(arr, obj, values, axis):
    # axis is a required argument here to avoid needing to deal with the numpy
    # default case (which reshapes the array to make it flat)
    if not -arr.ndim <= axis < arr.ndim:
        raise IndexError('axis %r is out of bounds for an array of dimension '
                         '%s' % (axis, arr.ndim))
    if axis < 0:
        axis += arr.ndim

    if isinstance(obj, slice):
        obj = np.arange(*obj.indices(arr.shape[axis]))
    obj = np.asarray(obj)
    scalar_obj = obj.ndim == 0
    if scalar_obj:
        obj = np.atleast_1d(obj)

    obj = np.where(obj < 0, obj + arr.shape[axis], obj)
    if (np.diff(obj) < 0).any():
        raise NotImplementedError(
            'da.insert only implemented for monotonic ``obj`` argument')

    split_arr = split_at_breaks(arr, np.unique(obj), axis)

    if getattr(values, 'ndim', 0) == 0:
        # we need to turn values into a dask array
        name = 'values-' + tokenize(values)
        dtype = getattr(values, 'dtype', type(values))
        values = Array({(name,): values}, name, chunks=(), dtype=dtype)

        values_shape = tuple(len(obj) if axis == n else s
                             for n, s in enumerate(arr.shape))
        values = broadcast_to(values, values_shape)
    elif scalar_obj:
        values = values[(slice(None),) * axis + (None,)]

    values_chunks = tuple(values_bd if axis == n else arr_bd
                          for n, (arr_bd, values_bd)
                          in enumerate(zip(arr.chunks,
                                           values.chunks)))
    values = values.rechunk(values_chunks)

    counts = np.bincount(obj)[:-1]
    values_breaks = np.cumsum(counts[counts > 0])
    split_values = split_at_breaks(values, values_breaks, axis)

    interleaved = list(interleave([split_arr, split_values]))
    interleaved = [i for i in interleaved if i.nbytes]
    return concatenate(interleaved, axis=axis)
