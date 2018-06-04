from __future__ import absolute_import, division, print_function

from operator import getitem
from itertools import product
from numbers import Integral

from toolz import merge, pipe, concat, partial
from toolz.curried import map

from . import chunk, wrap
from .core import Array, map_blocks, concatenate, concatenate3, reshapelist
from .. import sharedict
from ..base import tokenize
from ..core import flatten
from ..utils import concrete


def fractional_slice(task, axes):
    """

    >>> fractional_slice(('x', 5.1), {0: 2})  # doctest: +SKIP
    (getitem, ('x', 6), (slice(0, 2),))

    >>> fractional_slice(('x', 3, 5.1), {0: 2, 1: 3})  # doctest: +SKIP
    (getitem, ('x', 3, 5), (slice(None, None, None), slice(-3, None)))

    >>> fractional_slice(('x', 2.9, 5.1), {0: 2, 1: 3})  # doctest: +SKIP
    (getitem, ('x', 3, 5), (slice(0, 2), slice(-3, None)))
    """
    rounded = (task[0],) + tuple(int(round(i)) for i in task[1:])

    index = []
    for i, (t, r) in enumerate(zip(task[1:], rounded[1:])):
        depth = axes.get(i, 0)
        if t == r:
            index.append(slice(None, None, None))
        elif t < r:
            index.append(slice(0, depth))
        elif t > r and depth == 0:
            index.append(slice(0, 0))
        else:
            index.append(slice(-depth, None))

    index = tuple(index)

    if all(ind == slice(None, None, None) for ind in index):
        return task
    else:
        return (getitem, rounded, index)


def expand_key(k, dims):
    """ Get all neighboring keys around center

    >>> expand_key(('x', 2, 3), dims=[5, 5])  # doctest: +NORMALIZE_WHITESPACE
    [[('x', 1.1, 2.1), ('x', 1.1, 3), ('x', 1.1, 3.9)],
     [('x',   2, 2.1), ('x',   2, 3), ('x',   2, 3.9)],
     [('x', 2.9, 2.1), ('x', 2.9, 3), ('x', 2.9, 3.9)]]

    >>> expand_key(('x', 0, 4), dims=[5, 5])  # doctest: +NORMALIZE_WHITESPACE
    [[('x',   0, 3.1), ('x',   0,   4)],
     [('x', 0.9, 3.1), ('x', 0.9,   4)]]
    """
    def inds(i, ind):
        rv = []
        if ind - 0.9 > 0:
            rv.append(ind - 0.9)
        rv.append(ind)
        if ind + 0.9 < dims[i] - 1:
            rv.append(ind + 0.9)
        return rv

    shape = []
    for i, ind in enumerate(k[1:]):
        num = 1
        if ind > 0:
            num += 1
        if ind < dims[i] - 1:
            num += 1
        shape.append(num)

    seq = list(product([k[0]], *[inds(i, ind)
                                 for i, ind in enumerate(k[1:])]))
    return reshapelist(shape, seq)


def ghost_internal(x, axes):
    """ Share boundaries between neighboring blocks

    Parameters
    ----------

    x: da.Array
        A dask array
    axes: dict
        The size of the shared boundary per axis

    The axes input informs how many cells to overlap between neighboring blocks
    {0: 2, 2: 5} means share two cells in 0 axis, 5 cells in 2 axis
    """
    dims = list(map(len, x.chunks))
    expand_key2 = partial(expand_key, dims=dims)
    interior_keys = pipe(x.__dask_keys__(), flatten, map(expand_key2),
                         map(flatten), concat, list)

    token = tokenize(x, axes)
    name = 'ghost-' + token
    interior_slices = {}
    ghost_blocks = {}
    for k in interior_keys:
        frac_slice = fractional_slice(k, axes)
        if k != frac_slice:
            interior_slices[k] = frac_slice

        else:
            ghost_blocks[(name,) + k[1:]] = (concatenate3,
                                             (concrete, expand_key2(k)))

    chunks = []
    for i, bds in enumerate(x.chunks):
        if len(bds) == 1:
            chunks.append(bds)
        else:
            left = [bds[0] + axes.get(i, 0)]
            right = [bds[-1] + axes.get(i, 0)]
            mid = []
            for bd in bds[1:-1]:
                mid.append(bd + axes.get(i, 0) * 2)
            chunks.append(left + mid + right)

    dsk = merge(interior_slices, ghost_blocks)
    dsk = sharedict.merge(x.dask, (name, dsk))

    return Array(dsk, name, chunks, dtype=x.dtype)


def trim_internal(x, axes):
    """ Trim sides from each block

    This couples well with the ghost operation, which may leave excess data on
    each block

    See also
    --------
    dask.array.chunk.trim
    dask.array.map_blocks
    """
    olist = []
    for i, bd in enumerate(x.chunks):
        ilist = []
        for d in bd:
            ilist.append(d - axes.get(i, 0) * 2)
        olist.append(tuple(ilist))

    chunks = tuple(olist)

    return map_blocks(partial(chunk.trim, axes=axes), x, chunks=chunks,
                      dtype=x.dtype)


def periodic(x, axis, depth):
    """ Copy a slice of an array around to its other side

    Useful to create periodic boundary conditions for ghost
    """

    left = ((slice(None, None, None),) * axis +
            (slice(0, depth),) +
            (slice(None, None, None),) * (x.ndim - axis - 1))
    right = ((slice(None, None, None),) * axis +
             (slice(-depth, None),) +
             (slice(None, None, None),) * (x.ndim - axis - 1))
    l = x[left]
    r = x[right]

    l, r = _remove_ghost_boundaries(l, r, axis, depth)

    return concatenate([r, x, l], axis=axis)


def reflect(x, axis, depth):
    """ Reflect boundaries of array on the same side

    This is the converse of ``periodic``
    """
    if depth == 1:
        left = ((slice(None, None, None),) * axis +
                (slice(0, 1),) +
                (slice(None, None, None),) * (x.ndim - axis - 1))
    else:
        left = ((slice(None, None, None),) * axis +
                (slice(depth - 1, None, -1),) +
                (slice(None, None, None),) * (x.ndim - axis - 1))
    right = ((slice(None, None, None),) * axis +
             (slice(-1, -depth - 1, -1),) +
             (slice(None, None, None),) * (x.ndim - axis - 1))
    l = x[left]
    r = x[right]

    l, r = _remove_ghost_boundaries(l, r, axis, depth)

    return concatenate([l, x, r], axis=axis)


def nearest(x, axis, depth):
    """ Each reflect each boundary value outwards

    This mimics what the skimage.filters.gaussian_filter(... mode="nearest")
    does.
    """
    left = ((slice(None, None, None),) * axis +
            (slice(0, 1),) +
            (slice(None, None, None),) * (x.ndim - axis - 1))
    right = ((slice(None, None, None),) * axis +
             (slice(-1, -2, -1),) +
             (slice(None, None, None),) * (x.ndim - axis - 1))

    l = concatenate([x[left]] * depth, axis=axis)
    r = concatenate([x[right]] * depth, axis=axis)

    l, r = _remove_ghost_boundaries(l, r, axis, depth)

    return concatenate([l, x, r], axis=axis)


def constant(x, axis, depth, value):
    """ Add constant slice to either side of array """
    chunks = list(x.chunks)
    chunks[axis] = (depth,)

    c = wrap.full(tuple(map(sum, chunks)), value,
                  chunks=tuple(chunks), dtype=x.dtype)

    return concatenate([c, x, c], axis=axis)


def _remove_ghost_boundaries(l, r, axis, depth):
    lchunks = list(l.chunks)
    lchunks[axis] = (depth,)
    rchunks = list(r.chunks)
    rchunks[axis] = (depth,)

    l = l.rechunk(tuple(lchunks))
    r = r.rechunk(tuple(rchunks))
    return l, r


def boundaries(x, depth=None, kind=None):
    """ Add boundary conditions to an array before ghosting

    See Also
    --------
    periodic
    constant
    """
    if not isinstance(kind, dict):
        kind = dict((i, kind) for i in range(x.ndim))
    if not isinstance(depth, dict):
        depth = dict((i, depth) for i in range(x.ndim))

    for i in range(x.ndim):
        d = depth.get(i, 0)
        if d == 0:
            continue

        this_kind = kind.get(i, 'none')
        if this_kind == 'none':
            continue
        elif this_kind == 'periodic':
            x = periodic(x, i, d)
        elif this_kind == 'reflect':
            x = reflect(x, i, d)
        elif this_kind == 'nearest':
            x = nearest(x, i, d)
        elif i in kind:
            x = constant(x, i, d, kind[i])

    return x


def ghost(x, depth, boundary):
    """ Share boundaries between neighboring blocks

    Parameters
    ----------

    x: da.Array
        A dask array
    depth: dict
        The size of the shared boundary per axis
    boundary: dict
        The boundary condition on each axis. Options are 'reflect', 'periodic',
        'nearest', 'none', or an array value.  Such a value will fill the
        boundary with that value.

    The depth input informs how many cells to overlap between neighboring
    blocks ``{0: 2, 2: 5}`` means share two cells in 0 axis, 5 cells in 2 axis.
    Axes missing from this input will not be overlapped.

    Examples
    --------
    >>> import numpy as np
    >>> import dask.array as da

    >>> x = np.arange(64).reshape((8, 8))
    >>> d = da.from_array(x, chunks=(4, 4))
    >>> d.chunks
    ((4, 4), (4, 4))

    >>> g = da.ghost.ghost(d, depth={0: 2, 1: 1},
    ...                       boundary={0: 100, 1: 'reflect'})
    >>> g.chunks
    ((8, 8), (6, 6))

    >>> np.array(g)
    array([[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [  0,   0,   1,   2,   3,   4,   3,   4,   5,   6,   7,   7],
           [  8,   8,   9,  10,  11,  12,  11,  12,  13,  14,  15,  15],
           [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
           [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
           [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
           [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
           [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
           [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
           [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
           [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
           [ 48,  48,  49,  50,  51,  52,  51,  52,  53,  54,  55,  55],
           [ 56,  56,  57,  58,  59,  60,  59,  60,  61,  62,  63,  63],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
           [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]])
    """
    depth2 = coerce_depth(x.ndim, depth)
    boundary2 = coerce_boundary(x.ndim, boundary)

    # is depth larger than chunk size?
    depth_values = [depth2.get(i, 0) for i in range(x.ndim)]
    for d, c in zip(depth_values, x.chunks):
        if d > min(c):
            raise ValueError("The overlapping depth %d is larger than your\n"
                             "smallest chunk size %d. Rechunk your array\n"
                             "with a larger chunk size or a chunk size that\n"
                             "more evenly divides the shape of your array." %
                             (d, min(c)))
    x2 = boundaries(x, depth2, boundary2)
    x3 = ghost_internal(x2, depth2)
    trim = dict((k, v * 2 if boundary2.get(k, 'none') != 'none' else 0)
                for k, v in depth2.items())
    x4 = chunk.trim(x3, trim)
    return x4


def add_dummy_padding(x, depth, boundary):
    """
    Pads an array which has 'none' as the boundary type.
    Used to simplify trimming arrays which use 'none'.

    >>> import dask.array as da
    >>> x = da.arange(6, chunks=3)
    >>> add_dummy_padding(x, {0: 1}, {0: 'none'}).compute()  # doctest: +NORMALIZE_WHITESPACE
    array([..., 0, 1, 2, 3, 4, 5, ...])
    """
    for k, v in boundary.items():
        d = depth.get(k, 0)
        if v == 'none' and d > 0:
            empty_shape = list(x.shape)
            empty_shape[k] = d

            empty_chunks = list(x.chunks)
            empty_chunks[k] = (d,)

            empty = wrap.empty(empty_shape, chunks=empty_chunks, dtype=x.dtype)

            out_chunks = list(x.chunks)
            ax_chunks = list(out_chunks[k])
            ax_chunks[0] += d
            ax_chunks[-1] += d
            out_chunks[k] = tuple(ax_chunks)

            x = concatenate([empty, x, empty], axis=k)
            x = x.rechunk(out_chunks)
    return x


def map_overlap(x, func, depth, boundary=None, trim=True, **kwargs):
    """ Map a function over blocks of the array with some overlap

    We share neighboring zones between blocks of the array, then map a
    function, then trim away the neighboring strips.

    Parameters
    ----------
    func: function
        The function to apply to each extended block
    depth: int, tuple, or dict
        The number of elements that each block should share with its neighbors
        If a tuple or dict then this can be different per axis
    boundary: str, tuple, dict
        How to handle the boundaries.
        Values include 'reflect', 'periodic', 'nearest', 'none',
        or any constant value like 0 or np.nan
    trim: bool
        Whether or not to trim ``depth`` elements from each block after
        calling the map function.
        Set this to False if your mapping function already does this for you
    **kwargs:
        Other keyword arguments valid in ``map_blocks``

    Examples
    --------
    >>> import numpy as np
    >>> import dask.array as da

    >>> x = np.array([1, 1, 2, 3, 3, 3, 2, 1, 1])
    >>> x = da.from_array(x, chunks=5)
    >>> def derivative(x):
    ...     return x - np.roll(x, 1)

    >>> y = x.map_overlap(derivative, depth=1, boundary=0)
    >>> y.compute()
    array([ 1,  0,  1,  1,  0,  0, -1, -1,  0])

    >>> x = np.arange(16).reshape((4, 4))
    >>> d = da.from_array(x, chunks=(2, 2))
    >>> d.map_overlap(lambda x: x + x.size, depth=1).compute()
    array([[16, 17, 18, 19],
           [20, 21, 22, 23],
           [24, 25, 26, 27],
           [28, 29, 30, 31]])

    >>> func = lambda x: x + x.size
    >>> depth = {0: 1, 1: 1}
    >>> boundary = {0: 'reflect', 1: 'none'}
    >>> d.map_overlap(func, depth, boundary).compute()  # doctest: +NORMALIZE_WHITESPACE
    array([[12,  13,  14,  15],
           [16,  17,  18,  19],
           [20,  21,  22,  23],
           [24,  25,  26,  27]])
    """
    depth2 = coerce_depth(x.ndim, depth)
    boundary2 = coerce_boundary(x.ndim, boundary)

    g = ghost(x, depth=depth2, boundary=boundary2)
    g2 = g.map_blocks(func, **kwargs)
    if trim:
        g3 = add_dummy_padding(g2, depth2, boundary2)
        return trim_internal(g3, depth2)
    else:
        return g2


def coerce_depth(ndim, depth):
    if isinstance(depth, Integral):
        depth = (depth,) * ndim
    if isinstance(depth, tuple):
        depth = dict(zip(range(ndim), depth))

    return depth


def coerce_boundary(ndim, boundary):
    if boundary is None:
        boundary = 'reflect'
    if not isinstance(boundary, (tuple, dict)):
        boundary = (boundary,) * ndim
    if isinstance(boundary, tuple):
        boundary = dict(zip(range(ndim), boundary))

    return boundary
