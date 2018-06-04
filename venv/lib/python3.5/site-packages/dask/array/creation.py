from __future__ import absolute_import, division, print_function

from functools import partial, wraps
from itertools import product
from operator import add
from numbers import Integral

import numpy as np
from toolz import accumulate, sliding_window

from .. import sharedict
from ..base import tokenize
from ..utils import ignoring
from . import chunk
from .core import (Array, asarray, normalize_chunks,
                   stack, concatenate,
                   broadcast_arrays)
from .wrap import empty, ones, zeros, full


def empty_like(a, dtype=None, chunks=None):
    """
    Return a new array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of the
        returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    out : ndarray
        Array of uninitialized (arbitrary) data with the same
        shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    zeros_like : Return an array of zeros with shape and type of input.
    empty : Return a new uninitialized array.
    ones : Return a new array setting values to one.
    zeros : Return a new array setting values to zero.

    Notes
    -----
    This function does *not* initialize the returned array; to do that use
    `zeros_like` or `ones_like` instead.  It may be marginally faster than
    the functions that do set the array values.
    """

    a = asarray(a)
    return empty(
        a.shape, dtype=(dtype or a.dtype), chunks=(chunks or a.chunks)
    )


def ones_like(a, dtype=None, chunks=None):
    """
    Return an array of ones with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    out : ndarray
        Array of ones with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a)
    return ones(
        a.shape, dtype=(dtype or a.dtype), chunks=(chunks or a.chunks)
    )


def zeros_like(a, dtype=None, chunks=None):
    """
    Return an array of zeros with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    dtype : data-type, optional
        Overrides the data type of the result.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    out : ndarray
        Array of zeros with the same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    """

    a = asarray(a)
    return zeros(
        a.shape, dtype=(dtype or a.dtype), chunks=(chunks or a.chunks)
    )


def full_like(a, fill_value, dtype=None, chunks=None):
    """
    Return a full array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define these same attributes of
        the returned array.
    fill_value : scalar
        Fill value.
    dtype : data-type, optional
        Overrides the data type of the result.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    out : ndarray
        Array of `fill_value` with the same shape and type as `a`.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.
    full : Fill a new array.
    """

    a = asarray(a)
    return full(
        a.shape,
        fill_value,
        dtype=(dtype or a.dtype),
        chunks=(chunks or a.chunks)
    )


def linspace(start, stop, num=50, chunks=None, dtype=None):
    """
    Return `num` evenly spaced values over the closed interval [`start`,
    `stop`].

    TODO: implement the `endpoint`, `restep`, and `dtype` keyword args

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The last value of the sequence.
    num : int, optional
        Number of samples to include in the returned dask array, including the
        endpoints.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if `num % blocksize != 0`

    Returns
    -------
    samples : dask array

    See Also
    --------
    dask.array.arange
    """
    num = int(num)

    if chunks is None:
        raise ValueError("Must supply a chunks= keyword argument")

    chunks = normalize_chunks(chunks, (num,))

    range_ = stop - start

    space = float(range_) / (num - 1)

    if dtype is None:
        dtype = np.linspace(0, 1, 1).dtype

    name = 'linspace-' + tokenize((start, stop, num, chunks, dtype))

    dsk = {}
    blockstart = start

    for i, bs in enumerate(chunks[0]):
        blockstop = blockstart + ((bs - 1) * space)
        task = (partial(np.linspace, dtype=dtype), blockstart, blockstop, bs)
        blockstart = blockstart + (space * bs)
        dsk[(name, i)] = task

    return Array(dsk, name, chunks, dtype=dtype)


def arange(*args, **kwargs):
    """
    Return evenly spaced values from `start` to `stop` with step size `step`.

    The values are half-open [start, stop), so including start and excluding
    stop. This is basically the same as python's range function but for dask
    arrays.

    When using a non-integer step, such as 0.1, the results will often not be
    consistent. It is better to use linspace for these cases.

    Parameters
    ----------
    start : int, optional
        The starting value of the sequence. The default is 0.
    stop : int
        The end of the interval, this value is excluded from the interval.
    step : int, optional
        The spacing between the values. The default is 1 when not specified.
        The last value of the sequence.
    chunks :  int
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    samples : dask array

    See Also
    --------
    dask.array.linspace
    """
    if len(args) == 1:
        start = 0
        stop = args[0]
        step = 1
    elif len(args) == 2:
        start = args[0]
        stop = args[1]
        step = 1
    elif len(args) == 3:
        start, stop, step = args
    else:
        raise TypeError('''
        arange takes 3 positional arguments: arange([start], stop, [step])
        ''')

    if 'chunks' not in kwargs:
        raise ValueError("Must supply a chunks= keyword argument")
    chunks = kwargs['chunks']

    dtype = kwargs.get('dtype', None)
    if dtype is None:
        dtype = np.arange(0, 1, step).dtype

    num = max(np.ceil((stop - start) / step), 0)
    chunks = normalize_chunks(chunks, (num,))

    name = 'arange-' + tokenize((start, stop, step, chunks, num))
    dsk = {}
    elem_count = 0

    for i, bs in enumerate(chunks[0]):
        blockstart = start + (elem_count * step)
        blockstop = start + ((elem_count + bs) * step)
        task = (chunk.arange, blockstart, blockstop, step, bs, dtype)
        dsk[(name, i)] = task
        elem_count += bs

    return Array(dsk, name, chunks, dtype=dtype)


@wraps(np.meshgrid)
def meshgrid(*xi, **kwargs):
    indexing = kwargs.pop("indexing", "xy")
    sparse = bool(kwargs.pop("sparse", False))

    if "copy" in kwargs:
        raise NotImplementedError("`copy` not supported")

    if kwargs:
        raise TypeError("unsupported keyword argument(s) provided")

    if indexing not in ("ij", "xy"):
        raise ValueError("`indexing` must be `'ij'` or `'xy'`")

    xi = [asarray(e) for e in xi]
    xi = [e.flatten() for e in xi]

    if indexing == "xy" and len(xi) > 1:
        xi[0], xi[1] = xi[1], xi[0]

    grid = []
    for i in range(len(xi)):
        s = len(xi) * [None]
        s[i] = slice(None)
        s = tuple(s)

        r = xi[i][s]

        grid.append(r)

    if not sparse:
        grid = broadcast_arrays(*grid)

    if indexing == "xy" and len(xi) > 1:
        grid[0], grid[1] = grid[1], grid[0]

    return grid


def indices(dimensions, dtype=int, chunks=None):
    """
    Implements NumPy's ``indices`` for Dask Arrays.

    Generates a grid of indices covering the dimensions provided.

    The final array has the shape ``(len(dimensions), *dimensions)``. The
    chunks are used to specify the chunking for axis 1 up to
    ``len(dimensions)``. The 0th axis always has chunks of length 1.

    Parameters
    ----------
    dimensions : sequence of ints
        The shape of the index grid.
    dtype : dtype, optional
        Type to use for the array. Default is ``int``.
    chunks : sequence of ints
        The number of samples on each block. Note that the last block will have
        fewer samples if ``len(array) % chunks != 0``.

    Returns
    -------
    grid : dask array
    """
    if chunks is None:
        raise ValueError("Must supply a chunks= keyword argument")

    dimensions = tuple(dimensions)
    dtype = np.dtype(dtype)
    chunks = tuple(chunks)

    if len(dimensions) != len(chunks):
        raise ValueError("Need same number of chunks as dimensions.")

    xi = []
    for i in range(len(dimensions)):
        xi.append(arange(dimensions[i], dtype=dtype, chunks=(chunks[i],)))

    grid = []
    if np.prod(dimensions):
        grid = meshgrid(*xi, indexing="ij")

    if grid:
        grid = stack(grid)
    else:
        grid = empty(
            (len(dimensions),) + dimensions, dtype=dtype, chunks=(1,) + chunks
        )

    return grid


def eye(N, chunks, M=None, k=0, dtype=float):
    """
    Return a 2-D Array with ones on the diagonal and zeros elsewhere.

    Parameters
    ----------
    N : int
      Number of rows in the output.
    chunks: int
        chunk size of resulting blocks
    M : int, optional
      Number of columns in the output. If None, defaults to `N`.
    k : int, optional
      Index of the diagonal: 0 (the default) refers to the main diagonal,
      a positive value refers to an upper diagonal, and a negative value
      to a lower diagonal.
    dtype : data-type, optional
      Data-type of the returned array.

    Returns
    -------
    I : Array of shape (N,M)
      An array where all elements are equal to zero, except for the `k`-th
      diagonal, whose values are equal to one.
    """
    if not isinstance(chunks, int):
        raise ValueError('chunks must be an int')

    token = tokenize(N, chunk, M, k, dtype)
    name_eye = 'eye-' + token

    eye = {}
    if M is None:
        M = N

    vchunks = [chunks] * (N // chunks)
    if N % chunks != 0:
        vchunks.append(N % chunks)
    hchunks = [chunks] * (M // chunks)
    if M % chunks != 0:
        hchunks.append(M % chunks)

    for i, vchunk in enumerate(vchunks):
        for j, hchunk in enumerate(hchunks):
            if (j - i - 1) * chunks <= k <= (j - i + 1) * chunks:
                eye[name_eye, i, j] = (np.eye, vchunk, hchunk, k - (j - i) * chunks, dtype)
            else:
                eye[name_eye, i, j] = (np.zeros, (vchunk, hchunk), dtype)
    return Array(eye, name_eye, shape=(N, M),
                 chunks=(chunks, chunks), dtype=dtype)


@wraps(np.diag)
def diag(v):
    name = 'diag-' + tokenize(v)
    if isinstance(v, np.ndarray):
        if v.ndim == 1:
            chunks = ((v.shape[0],), (v.shape[0],))
            dsk = {(name, 0, 0): (np.diag, v)}
        elif v.ndim == 2:
            chunks = ((min(v.shape),),)
            dsk = {(name, 0): (np.diag, v)}
        else:
            raise ValueError("Array must be 1d or 2d only")
        return Array(dsk, name, chunks, dtype=v.dtype)
    if not isinstance(v, Array):
        raise TypeError("v must be a dask array or numpy array, "
                        "got {0}".format(type(v)))
    if v.ndim != 1:
        if v.chunks[0] == v.chunks[1]:
            dsk = {(name, i): (np.diag, row[i])
                   for i, row in enumerate(v.__dask_keys__())}
            return Array(sharedict.merge(v.dask, (name, dsk)), name, (v.chunks[0],), dtype=v.dtype)
        else:
            raise NotImplementedError("Extracting diagonals from non-square "
                                      "chunked arrays")
    chunks_1d = v.chunks[0]
    blocks = v.__dask_keys__()
    dsk = {}
    for i, m in enumerate(chunks_1d):
        for j, n in enumerate(chunks_1d):
            key = (name, i, j)
            if i == j:
                dsk[key] = (np.diag, blocks[i])
            else:
                dsk[key] = (np.zeros, (m, n))

    return Array(sharedict.merge(v.dask, (name, dsk)), name, (chunks_1d, chunks_1d),
                 dtype=v.dtype)


def triu(m, k=0):
    """
    Upper triangle of an array with elements above the `k`-th diagonal zeroed.

    Parameters
    ----------
    m : array_like, shape (M, N)
        Input array.
    k : int, optional
        Diagonal above which to zero elements.  `k = 0` (the default) is the
        main diagonal, `k < 0` is below it and `k > 0` is above.

    Returns
    -------
    triu : ndarray, shape (M, N)
        Upper triangle of `m`, of same shape and data-type as `m`.

    See Also
    --------
    tril : lower triangle of an array
    """
    if m.ndim != 2:
        raise ValueError('input must be 2 dimensional')
    if m.shape[0] != m.shape[1]:
        raise NotImplementedError('input must be a square matrix')
    if m.chunks[0][0] != m.chunks[1][0]:
        msg = ('chunks must be a square. '
               'Use .rechunk method to change the size of chunks.')
        raise NotImplementedError(msg)

    rdim = len(m.chunks[0])
    hdim = len(m.chunks[1])
    chunk = m.chunks[0][0]

    token = tokenize(m, k)
    name = 'triu-' + token

    dsk = {}
    for i in range(rdim):
        for j in range(hdim):
            if chunk * (j - i + 1) < k:
                dsk[(name, i, j)] = (np.zeros, (m.chunks[0][i], m.chunks[1][j]))
            elif chunk * (j - i - 1) < k <= chunk * (j - i + 1):
                dsk[(name, i, j)] = (np.triu, (m.name, i, j), k - (chunk * (j - i)))
            else:
                dsk[(name, i, j)] = (m.name, i, j)
    return Array(sharedict.merge((name, dsk), m.dask), name,
                 shape=m.shape, chunks=m.chunks, dtype=m.dtype)


def tril(m, k=0):
    """
    Lower triangle of an array with elements above the `k`-th diagonal zeroed.

    Parameters
    ----------
    m : array_like, shape (M, M)
        Input array.
    k : int, optional
        Diagonal above which to zero elements.  `k = 0` (the default) is the
        main diagonal, `k < 0` is below it and `k > 0` is above.

    Returns
    -------
    tril : ndarray, shape (M, M)
        Lower triangle of `m`, of same shape and data-type as `m`.

    See Also
    --------
    triu : upper triangle of an array
    """
    if m.ndim != 2:
        raise ValueError('input must be 2 dimensional')
    if m.shape[0] != m.shape[1]:
        raise NotImplementedError('input must be a square matrix')
    if not len(set(m.chunks[0] + m.chunks[1])) == 1:
        msg = ('All chunks must be a square matrix to perform lu decomposition. '
               'Use .rechunk method to change the size of chunks.')
        raise ValueError(msg)

    rdim = len(m.chunks[0])
    hdim = len(m.chunks[1])
    chunk = m.chunks[0][0]

    token = tokenize(m, k)
    name = 'tril-' + token

    dsk = {}
    for i in range(rdim):
        for j in range(hdim):
            if chunk * (j - i + 1) < k:
                dsk[(name, i, j)] = (m.name, i, j)
            elif chunk * (j - i - 1) < k <= chunk * (j - i + 1):
                dsk[(name, i, j)] = (np.tril, (m.name, i, j), k - (chunk * (j - i)))
            else:
                dsk[(name, i, j)] = (np.zeros, (m.chunks[0][i], m.chunks[1][j]))
    dsk = sharedict.merge(m.dask, (name, dsk))
    return Array(dsk, name, shape=m.shape, chunks=m.chunks, dtype=m.dtype)


def offset_func(func, offset, *args):
    """  Offsets inputs by offset

    >>> double = lambda x: x * 2
    >>> f = offset_func(double, (10,))
    >>> f(1)
    22
    >>> f(300)
    620
    """
    def _offset(*args):
        args2 = list(map(add, args, offset))
        return func(*args2)

    with ignoring(Exception):
        _offset.__name__ = 'offset_' + func.__name__

    return _offset


@wraps(np.fromfunction)
def fromfunction(func, chunks=None, shape=None, dtype=None):
    if chunks:
        chunks = normalize_chunks(chunks, shape)
    name = 'fromfunction-' + tokenize(func, chunks, shape, dtype)
    keys = list(product([name], *[range(len(bd)) for bd in chunks]))
    aggdims = [list(accumulate(add, (0,) + bd[:-1])) for bd in chunks]
    offsets = list(product(*aggdims))
    shapes = list(product(*chunks))

    values = [(np.fromfunction, offset_func(func, offset), shp)
              for offset, shp in zip(offsets, shapes)]

    dsk = dict(zip(keys, values))

    return Array(dsk, name, chunks, dtype=dtype)


@wraps(np.repeat)
def repeat(a, repeats, axis=None):
    if axis is None:
        if a.ndim == 1:
            axis = 0
        else:
            raise NotImplementedError("Must supply an integer axis value")

    if not isinstance(repeats, Integral):
        raise NotImplementedError("Only integer valued repeats supported")

    if -a.ndim <= axis < 0:
        axis += a.ndim
    elif not 0 <= axis <= a.ndim - 1:
        raise ValueError("axis(=%d) out of bounds" % axis)

    if repeats == 1:
        return a

    cchunks = np.cumsum((0,) + a.chunks[axis])
    slices = []
    for c_start, c_stop in sliding_window(2, cchunks):
        ls = np.linspace(c_start, c_stop, repeats).round(0)
        for ls_start, ls_stop in sliding_window(2, ls):
            if ls_start != ls_stop:
                slices.append(slice(ls_start, ls_stop))

    all_slice = slice(None, None, None)
    slices = [(all_slice,) * axis + (s,) + (all_slice,) * (a.ndim - axis - 1)
              for s in slices]

    slabs = [a[slc] for slc in slices]

    out = []
    for slab in slabs:
        chunks = list(slab.chunks)
        assert len(chunks[axis]) == 1
        chunks[axis] = (chunks[axis][0] * repeats,)
        chunks = tuple(chunks)
        result = slab.map_blocks(np.repeat, repeats, axis=axis, chunks=chunks,
                                 dtype=slab.dtype)
        out.append(result)

    return concatenate(out, axis=axis)


@wraps(np.tile)
def tile(A, reps):
    if not isinstance(reps, Integral):
        raise NotImplementedError("Only integer valued `reps` supported.")

    if reps < 0:
        raise ValueError("Negative `reps` are not allowed.")
    elif reps == 0:
        return A[..., :0]
    elif reps == 1:
        return A

    return concatenate(reps * [A], axis=-1)
