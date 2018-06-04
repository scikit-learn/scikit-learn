from __future__ import absolute_import, division, print_function

from itertools import product
import math
from numbers import Integral, Number
from operator import getitem, itemgetter

import numpy as np
from toolz import memoize, merge, pluck, concat

from .. import core
from .. import sharedict
from ..base import tokenize, is_dask_collection

colon = slice(None, None, None)


def _sanitize_index_element(ind):
    """Sanitize a one-element index."""
    if isinstance(ind, Number):
        ind2 = int(ind)
        if ind2 != ind:
            raise IndexError("Bad index.  Must be integer-like: %s" % ind)
        else:
            return ind2
    elif ind is None:
        return None
    else:
        raise TypeError("Invalid index type", type(ind), ind)


def sanitize_index(ind):
    """ Sanitize the elements for indexing along one axis

    >>> sanitize_index([2, 3, 5])
    array([2, 3, 5])
    >>> sanitize_index([True, False, True, False])
    array([0, 2])
    >>> sanitize_index(np.array([1, 2, 3]))
    array([1, 2, 3])
    >>> sanitize_index(np.array([False, True, True]))
    array([1, 2])
    >>> type(sanitize_index(np.int32(0)))
    <type 'int'>
    >>> sanitize_index(1.0)
    1
    >>> sanitize_index(0.5)
    Traceback (most recent call last):
    ...
    IndexError: Bad index.  Must be integer-like: 0.5
    """
    if ind is None:
        return None
    elif isinstance(ind, slice):
        return slice(_sanitize_index_element(ind.start),
                     _sanitize_index_element(ind.stop),
                     _sanitize_index_element(ind.step))
    elif isinstance(ind, Number):
        return _sanitize_index_element(ind)
    elif is_dask_collection(ind):
        return ind
    index_array = np.asanyarray(ind)
    if index_array.dtype == bool:
        nonzero = np.nonzero(index_array)
        if len(nonzero) == 1:
            # If a 1-element tuple, unwrap the element
            nonzero = nonzero[0]
        return np.asanyarray(nonzero)
    elif np.issubdtype(index_array.dtype, np.integer):
        return index_array
    elif np.issubdtype(index_array.dtype, float):
        int_index = index_array.astype(np.intp)
        if np.allclose(index_array, int_index):
            return int_index
        else:
            check_int = np.isclose(index_array, int_index)
            first_err = index_array.ravel(
            )[np.flatnonzero(~check_int)[0]]
            raise IndexError("Bad index.  Must be integer-like: %s" %
                             first_err)
    else:
        raise TypeError("Invalid index type", type(ind), ind)


def slice_array(out_name, in_name, blockdims, index):
    """
    Master function for array slicing

    This function makes a new dask that slices blocks along every
    dimension and aggregates (via cartesian product) each dimension's
    slices so that the resulting block slices give the same results
    as the original slice on the original structure

    Index must be a tuple.  It may contain the following types

        int, slice, list (at most one list), None

    Parameters
    ----------
    in_name - string
      This is the dask variable name that will be used as input
    out_name - string
      This is the dask variable output name
    blockshape - iterable of integers
    index - iterable of integers, slices, lists, or None

    Returns
    -------
    Dict where the keys are tuples of

        (out_name, dim_index[, dim_index[, ...]])

    and the values are

        (function, (in_name, dim_index, dim_index, ...),
                   (slice(...), [slice()[,...]])

    Also new blockdims with shapes of each block

        ((10, 10, 10, 10), (20, 20))

    Examples
    --------
    >>> dsk, blockdims = slice_array('y', 'x', [(20, 20, 20, 20, 20)],
    ...                              (slice(10, 35),))  #  doctest: +SKIP
    >>> dsk  # doctest: +SKIP
    {('y', 0): (getitem, ('x', 0), (slice(10, 20),)),
     ('y', 1): (getitem, ('x', 1), (slice(0, 15),))}
    >>> blockdims  # doctest: +SKIP
    ((10, 15),)

    See Also
    --------
    This function works by successively unwrapping cases and passing down
    through a sequence of functions.

    slice_with_newaxis - handle None/newaxis case
    slice_wrap_lists - handle fancy indexing with lists
    slice_slices_and_integers - handle everything else
    """
    blockdims = tuple(map(tuple, blockdims))

    # x[:, :, :] - Punt and return old value
    if all(isinstance(index, slice) and index == slice(None, None, None)
           for index in index):
        suffixes = product(*[range(len(bd)) for bd in blockdims])
        dsk = dict(((out_name,) + s, (in_name,) + s)
                   for s in suffixes)
        return dsk, blockdims

    # Add in missing colons at the end as needed.  x[5] -> x[5, :, :]
    not_none_count = sum(i is not None for i in index)
    missing = len(blockdims) - not_none_count
    index += (slice(None, None, None),) * missing

    # Pass down to next function
    dsk_out, bd_out = slice_with_newaxes(out_name, in_name, blockdims, index)

    bd_out = tuple(map(tuple, bd_out))
    return dsk_out, bd_out


def slice_with_newaxes(out_name, in_name, blockdims, index):
    """
    Handle indexing with Nones

    Strips out Nones then hands off to slice_wrap_lists
    """
    # Strip Nones from index
    index2 = tuple([ind for ind in index if ind is not None])
    where_none = [i for i, ind in enumerate(index) if ind is None]
    where_none_orig = list(where_none)
    for i, x in enumerate(where_none):
        n = sum(isinstance(ind, int) for ind in index[:x])
        if n:
            where_none[i] -= n

    # Pass down and do work
    dsk, blockdims2 = slice_wrap_lists(out_name, in_name, blockdims, index2)

    if where_none:
        expand = expander(where_none)
        expand_orig = expander(where_none_orig)

        # Insert ",0" into the key:  ('x', 2, 3) -> ('x', 0, 2, 0, 3)
        dsk2 = {(out_name,) + expand(k[1:], 0):
                (v[:2] + (expand_orig(v[2], None),))
                for k, v in dsk.items()
                if k[0] == out_name}

        # Add back intermediate parts of the dask that weren't the output
        dsk3 = merge(dsk2, {k: v for k, v in dsk.items() if k[0] != out_name})

        # Insert (1,) into blockdims:  ((2, 2), (3, 3)) -> ((2, 2), (1,), (3, 3))
        blockdims3 = expand(blockdims2, (1,))

        return dsk3, blockdims3

    else:
        return dsk, blockdims2


def slice_wrap_lists(out_name, in_name, blockdims, index):
    """
    Fancy indexing along blocked array dasks

    Handles index of type list.  Calls slice_slices_and_integers for the rest

    See Also
    --------

    take - handle slicing with lists ("fancy" indexing)
    slice_slices_and_integers - handle slicing with slices and integers
    """
    assert all(isinstance(i, (slice, list, Integral, np.ndarray))
               for i in index)
    if not len(blockdims) == len(index):
        raise IndexError("Too many indices for array")

    # Do we have more than one list in the index?
    where_list = [i for i, ind in enumerate(index)
                  if isinstance(ind, np.ndarray) and ind.ndim > 0]
    if len(where_list) > 1:
        raise NotImplementedError("Don't yet support nd fancy indexing")
    # Is the single list an empty list? In this case just treat it as a zero
    # length slice
    if where_list and not index[where_list[0]].size:
        index = list(index)
        index[where_list.pop()] = slice(0, 0, 1)
        index = tuple(index)

    # No lists, hooray! just use slice_slices_and_integers
    if not where_list:
        return slice_slices_and_integers(out_name, in_name, blockdims, index)

    # Replace all lists with full slices  [3, 1, 0] -> slice(None, None, None)
    index_without_list = tuple(slice(None, None, None)
                               if isinstance(i, np.ndarray) else i
                               for i in index)

    # lists and full slices.  Just use take
    if all(isinstance(i, np.ndarray) or i == slice(None, None, None)
            for i in index):
        axis = where_list[0]
        blockdims2, dsk3 = take(out_name, in_name, blockdims,
                                index[where_list[0]], axis=axis)
    # Mixed case. Both slices/integers and lists. slice/integer then take
    else:
        # Do first pass without lists
        tmp = 'slice-' + tokenize((out_name, in_name, blockdims, index))
        dsk, blockdims2 = slice_slices_and_integers(tmp, in_name, blockdims, index_without_list)

        # After collapsing some axes due to int indices, adjust axis parameter
        axis = where_list[0]
        axis2 = axis - sum(1 for i, ind in enumerate(index)
                           if i < axis and isinstance(ind, Integral))

        # Do work
        blockdims2, dsk2 = take(out_name, tmp, blockdims2, index[axis],
                                axis=axis2)
        dsk3 = merge(dsk, dsk2)

    return dsk3, blockdims2


def slice_slices_and_integers(out_name, in_name, blockdims, index):
    """
    Dask array indexing with slices and integers

    See Also
    --------

    _slice_1d
    """
    shape = tuple(map(sum, blockdims))

    for dim, ind in zip(shape, index):
        if np.isnan(dim) and ind != slice(None, None, None):
            raise ValueError("Arrays chunk sizes are unknown: %s", shape)

    assert all(isinstance(ind, (slice, Integral)) for ind in index)
    assert len(index) == len(blockdims)

    # Get a list (for each dimension) of dicts{blocknum: slice()}
    block_slices = list(map(_slice_1d, shape, blockdims, index))
    sorted_block_slices = [sorted(i.items()) for i in block_slices]

    # (in_name, 1, 1, 2), (in_name, 1, 1, 4), (in_name, 2, 1, 2), ...
    in_names = list(product([in_name], *[pluck(0, s) for s in sorted_block_slices]))

    # (out_name, 0, 0, 0), (out_name, 0, 0, 1), (out_name, 0, 1, 0), ...
    out_names = list(product([out_name],
                             *[range(len(d))[::-1] if i.step and i.step < 0 else range(len(d))
                               for d, i in zip(block_slices, index)
                               if not isinstance(i, Integral)]))

    all_slices = list(product(*[pluck(1, s) for s in sorted_block_slices]))

    dsk_out = {out_name: (getitem, in_name, slices)
               for out_name, in_name, slices
               in zip(out_names, in_names, all_slices)}

    new_blockdims = [new_blockdim(d, db, i)
                     for d, i, db in zip(shape, index, blockdims)
                     if not isinstance(i, Integral)]

    return dsk_out, new_blockdims


def _slice_1d(dim_shape, lengths, index):
    """Returns a dict of {blocknum: slice}

    This function figures out where each slice should start in each
    block for a single dimension. If the slice won't return any elements
    in the block, that block will not be in the output.

    Parameters
    ----------

    dim_shape - the number of elements in this dimension.
      This should be a positive, non-zero integer
    blocksize - the number of elements per block in this dimension
      This should be a positive, non-zero integer
    index - a description of the elements in this dimension that we want
      This might be an integer, a slice(), or an Ellipsis

    Returns
    -------

    dictionary where the keys are the integer index of the blocks that
      should be sliced and the values are the slices

    Examples
    --------

    Trivial slicing

    >>> _slice_1d(100, [60, 40], slice(None, None, None))
    {0: slice(None, None, None), 1: slice(None, None, None)}

    100 length array cut into length 20 pieces, slice 0:35

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(0, 35))
    {0: slice(None, None, None), 1: slice(0, 15, 1)}

    Support irregular blocks and various slices

    >>> _slice_1d(100, [20, 10, 10, 10, 25, 25], slice(10, 35))
    {0: slice(10, 20, 1), 1: slice(None, None, None), 2: slice(0, 5, 1)}

    Support step sizes

    >>> _slice_1d(100, [15, 14, 13], slice(10, 41, 3))
    {0: slice(10, 15, 3), 1: slice(1, 14, 3), 2: slice(2, 12, 3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(0, 100, 40))  # step > blocksize
    {0: slice(0, 20, 40), 2: slice(0, 20, 40), 4: slice(0, 20, 40)}

    Also support indexing single elements

    >>> _slice_1d(100, [20, 20, 20, 20, 20], 25)
    {1: 5}

    And negative slicing

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, 0, -3))
    {0: slice(-2, -20, -3), 1: slice(-1, -21, -3), 2: slice(-3, -21, -3), 3: slice(-2, -21, -3), 4: slice(-1, -21, -3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, 12, -3))
    {0: slice(-2, -8, -3), 1: slice(-1, -21, -3), 2: slice(-3, -21, -3), 3: slice(-2, -21, -3), 4: slice(-1, -21, -3)}

    >>> _slice_1d(100, [20, 20, 20, 20, 20], slice(100, -12, -3))
    {4: slice(-1, -12, -3)}
    """
    chunk_boundaries = np.cumsum(lengths)

    if isinstance(index, Integral):
        # use right-side search to be consistent with previous result
        i = chunk_boundaries.searchsorted(index, side='right')
        if i > 0:
            # the very first chunk has no relative shift
            ind = index - chunk_boundaries[i - 1]
        else:
            ind = index
        return {i: ind}

    assert isinstance(index, slice)

    if index == colon:
        return {k: colon for k in range(len(lengths))}

    step = index.step or 1
    if step > 0:
        start = index.start or 0
        stop = index.stop if index.stop is not None else dim_shape
    else:
        start = index.start if index.start is not None else dim_shape - 1
        start = dim_shape - 1 if start >= dim_shape else start
        stop = -(dim_shape + 1) if index.stop is None else index.stop

    # posify start and stop
    if start < 0:
        start += dim_shape
    if stop < 0:
        stop += dim_shape

    d = dict()
    if step > 0:
        istart = chunk_boundaries.searchsorted(start, side='right')
        istop = chunk_boundaries.searchsorted(stop, side='left')

        # the bound is not exactly tight; make it tighter?
        istop = min(istop + 1, len(lengths))

        # jump directly to istart
        if istart > 0:
            start = start - chunk_boundaries[istart - 1]
            stop = stop - chunk_boundaries[istart - 1]

        for i in range(istart, istop):
            length = lengths[i]
            if start < length and stop > 0:
                d[i] = slice(start, min(stop, length), step)
                start = (start - length) % step
            else:
                start = start - length
            stop -= length
    else:
        rstart = start  # running start

        istart = chunk_boundaries.searchsorted(start, side='left')
        istop = chunk_boundaries.searchsorted(stop, side='right')

        # the bound is not exactly tight; make it tighter?
        istart = min(istart + 1, len(chunk_boundaries) - 1)
        istop = max(istop - 1, -1)

        for i in range(istart, istop, -1):
            chunk_stop = chunk_boundaries[i]
            # create a chunk start and stop
            if i == 0:
                chunk_start = 0
            else:
                chunk_start = chunk_boundaries[i - 1]

            # if our slice is in this chunk
            if (chunk_start <= rstart < chunk_stop) and (rstart > stop):
                d[i] = slice(rstart - chunk_stop,
                             max(chunk_start - chunk_stop - 1,
                                 stop - chunk_stop),
                             step)

                # compute the next running start point,
                offset = (rstart - (chunk_start - 1)) % step
                rstart = chunk_start + offset - 1

    # replace 0:20:1 with : if appropriate
    for k, v in d.items():
        if v == slice(0, lengths[k], 1):
            d[k] = slice(None, None, None)

    if not d:  # special case x[:0]
        d[0] = slice(0, 0, 1)

    return d


def partition_by_size(sizes, seq):
    """

    >>> partition_by_size([10, 20, 10], [1, 5, 9, 12, 29, 35])
    [array([1, 5, 9]), array([ 2, 19]), array([5])]
    """
    seq = np.asanyarray(seq)
    left = np.empty(len(sizes) + 1, dtype=int)
    left[0] = 0

    right = np.cumsum(sizes, out=left[1:])
    locations = np.empty(len(sizes) + 1, dtype=int)
    locations[0] = 0
    locations[1:] = np.searchsorted(seq, right)
    return [(seq[j:k] - l)
            for j, k, l in zip(locations[:-1], locations[1:], left)]


def issorted(seq):
    """ Is sequence sorted?

    >>> issorted([1, 2, 3])
    True
    >>> issorted([3, 1, 2])
    False
    """
    if len(seq) == 0:
        return True
    return np.all(seq[:-1] <= seq[1:])


def take_sorted(outname, inname, blockdims, index, axis=0):
    """ Index array with sorted list index

    Forms a dask for the following case

        x[:, [1, 3, 5, 10], ...]

    where the index, ``[1, 3, 5, 10]`` is sorted in non-decreasing order.

    >>> blockdims, dsk = take('y', 'x', [(20, 20, 20, 20)], [1, 3, 5, 47], axis=0)
    >>> blockdims
    ((3, 1),)
    >>> dsk  # doctest: +SKIP
    {('y', 0): (getitem, ('x', 0), ([1, 3, 5],)),
     ('y', 1): (getitem, ('x', 2), ([7],))}

    See Also
    --------
    take - calls this function
    """
    sizes = blockdims[axis]  # the blocksizes on the axis that we care about

    index_lists = partition_by_size(sizes, index)
    where_index = [i for i, il in enumerate(index_lists) if len(il)]
    index_lists = [il for il in index_lists if len(il)]

    dims = [range(len(bd)) for bd in blockdims]

    indims = list(dims)
    indims[axis] = list(range(len(where_index)))
    keys = list(product([outname], *indims))

    outdims = list(dims)
    outdims[axis] = where_index
    slices = [[colon] * len(bd) for bd in blockdims]
    slices[axis] = index_lists
    slices = list(product(*slices))
    inkeys = list(product([inname], *outdims))
    values = [(getitem, inkey, slc) for inkey, slc in zip(inkeys, slices)]

    blockdims2 = list(blockdims)
    blockdims2[axis] = tuple(map(len, index_lists))

    return tuple(blockdims2), dict(zip(keys, values))


def take(outname, inname, blockdims, index, axis=0):
    """ Index array with an iterable of index

    Handles a single index by a single list

    Mimics ``np.take``

    >>> blockdims, dsk = take('y', 'x', [(20, 20, 20, 20)], [5, 1, 47, 3], axis=0)
    >>> blockdims
    ((4,),)
    >>> dsk  # doctest: +SKIP
    {('y', 0): (getitem, (np.concatenate, [(getitem, ('x', 0), ([1, 3, 5],)),
                                           (getitem, ('x', 2), ([7],))],
                                          0),
                         (2, 0, 4, 1))}

    When list is sorted we retain original block structure

    >>> blockdims, dsk = take('y', 'x', [(20, 20, 20, 20)], [1, 3, 5, 47], axis=0)
    >>> blockdims
    ((3, 1),)
    >>> dsk  # doctest: +SKIP
    {('y', 0): (getitem, ('x', 0), ([1, 3, 5],)),
     ('y', 2): (getitem, ('x', 2), ([7],))}
    """

    index = np.asanyarray(index)

    if issorted(index):
        return take_sorted(outname, inname, blockdims, index, axis)

    if isinstance(index, np.ndarray) and index.ndim > 0:
        sorted_idx = np.sort(index)
    else:
        sorted_idx = index

    n = len(blockdims)
    sizes = blockdims[axis]  # the blocksizes on the axis that we care about

    index_lists = partition_by_size(sizes, sorted_idx)

    dims = [[0] if axis == i else list(range(len(bd)))
            for i, bd in enumerate(blockdims)]
    keys = list(product([outname], *dims))

    rev_index = np.searchsorted(sorted_idx, index)
    vals = [(getitem, (np.concatenate,
                       [(getitem, ((inname, ) + d[:axis] + (i, ) + d[axis + 1:]),
                         ((colon, ) * axis + (IL, ) + (colon, ) * (n - axis - 1)))
                        for i, IL in enumerate(index_lists) if len(IL)], axis),
             ((colon, ) * axis + (rev_index, ) + (colon, ) * (n - axis - 1)))
            for d in product(*dims)]

    blockdims2 = list(blockdims)
    blockdims2[axis] = (len(index), )

    return tuple(blockdims2), dict(zip(keys, vals))


def posify_index(shape, ind):
    """ Flip negative indices around to positive ones

    >>> posify_index(10, 3)
    3
    >>> posify_index(10, -3)
    7
    >>> posify_index(10, [3, -3])
    array([3, 7])

    >>> posify_index((10, 20), (3, -3))
    (3, 17)
    >>> posify_index((10, 20), (3, [3, 4, -3]))  # doctest: +NORMALIZE_WHITESPACE
    (3, array([ 3,  4, 17]))
    """
    if isinstance(ind, tuple):
        return tuple(map(posify_index, shape, ind))
    if isinstance(ind, Integral):
        if ind < 0 and not math.isnan(shape):
            return ind + shape
        else:
            return ind
    if isinstance(ind, (np.ndarray, list)) and not math.isnan(shape):
        ind = np.asanyarray(ind)
        return np.where(ind < 0, ind + shape, ind)
    return ind


@memoize
def _expander(where):
    if not where:
        def expand(seq, val):
            return seq
        return expand
    else:
        decl = """def expand(seq, val):
            return ({left}) + tuple({right})
        """
        left = []
        j = 0
        for i in range(max(where) + 1):
            if i in where:
                left.append("val, ")
            else:
                left.append("seq[%d], " % j)
                j += 1
        right = "seq[%d:]" % j
        left = "".join(left)
        decl = decl.format(**locals())
        ns = {}
        exec(compile(decl, "<dynamic>", "exec"), ns, ns)
        return ns['expand']


def expander(where):
    """Create a function to insert value at many locations in sequence.

    >>> expander([0, 2])(['a', 'b', 'c'], 'z')
    ('z', 'a', 'z', 'b', 'c')
    """
    return _expander(tuple(where))


def new_blockdim(dim_shape, lengths, index):
    """

    >>> new_blockdim(100, [20, 10, 20, 10, 40], slice(0, 90, 2))
    [10, 5, 10, 5, 15]

    >>> new_blockdim(100, [20, 10, 20, 10, 40], [5, 1, 30, 22])
    [4]

    >>> new_blockdim(100, [20, 10, 20, 10, 40], slice(90, 10, -2))
    [16, 5, 10, 5, 4]
    """
    if index == slice(None, None, None):
        return lengths
    if isinstance(index, list):
        return [len(index)]
    assert not isinstance(index, Integral)
    pairs = sorted(_slice_1d(dim_shape, lengths, index).items(),
                   key=itemgetter(0))
    slices = [slice(0, lengths[i], 1) if slc == slice(None, None, None) else slc
              for i, slc in pairs]
    if isinstance(index, slice) and index.step and index.step < 0:
        slices = slices[::-1]
    return [int(math.ceil((1. * slc.stop - slc.start) / slc.step)) for slc in slices]


def replace_ellipsis(n, index):
    """ Replace ... with slices, :, : ,:

    >>> replace_ellipsis(4, (3, Ellipsis, 2))
    (3, slice(None, None, None), slice(None, None, None), 2)

    >>> replace_ellipsis(2, (Ellipsis, None))
    (slice(None, None, None), slice(None, None, None), None)
    """
    # Careful about using in or index because index may contain arrays
    isellipsis = [i for i, ind in enumerate(index) if ind is Ellipsis]
    if not isellipsis:
        return index
    else:
        loc = isellipsis[0]
    extra_dimensions = n - (len(index) - sum(i is None for i in index) - 1)
    return (index[:loc] + (slice(None, None, None),) * extra_dimensions +
            index[loc + 1:])


def normalize_slice(idx, dim):
    """ Normalize slices to canonical form

    Parameters
    ----------
    idx: slice or other index
    dim: dimension length

    Examples
    --------
    >>> normalize_slice(slice(0, 10, 1), 10)
    slice(None, None, None)
    """

    if isinstance(idx, slice):
        start, stop, step = idx.start, idx.stop, idx.step
        if start is not None:
            if start < 0 and not math.isnan(dim):
                start = max(0, start + dim)
            elif start > dim:
                start = dim
        if stop is not None:
            if stop < 0 and not math.isnan(dim):
                stop = max(0, stop + dim)
            elif stop > dim:
                stop = dim
        if start == 0:
            start = None
        if stop == dim:
            stop = None
        if step == 1:
            step = None
        return slice(start, stop, step)
    return idx


def normalize_index(idx, shape):
    """ Normalize slicing indexes

    1.  Replaces ellipses with many full slices
    2.  Adds full slices to end of index
    3.  Checks bounding conditions
    4.  Replaces numpy arrays with lists
    5.  Posify's integers and lists
    6.  Normalizes slices to canonical form

    Examples
    --------
    >>> normalize_index(1, (10,))
    (1,)
    >>> normalize_index(-1, (10,))
    (9,)
    >>> normalize_index([-1], (10,))
    (array([9]),)
    >>> normalize_index(slice(-3, 10, 1), (10,))
    (slice(7, None, None),)
    >>> normalize_index((Ellipsis, None), (10,))
    (slice(None, None, None), None)
    """
    if not isinstance(idx, tuple):
        idx = (idx,)
    idx = replace_ellipsis(len(shape), idx)
    n_sliced_dims = 0
    for i in idx:
        if hasattr(i, 'ndim') and i.ndim >= 1:
            n_sliced_dims += i.ndim
        elif i is None:
            continue
        else:
            n_sliced_dims += 1
    idx = idx + (slice(None),) * (len(shape) - n_sliced_dims)
    if len([i for i in idx if i is not None]) > len(shape):
        raise IndexError("Too many indices for array")

    none_shape = []
    i = 0
    for ind in idx:
        if ind is not None:
            none_shape.append(shape[i])
            i += 1
        else:
            none_shape.append(None)

    for i, d in zip(idx, none_shape):
        if d is not None:
            check_index(i, d)
    idx = tuple(map(sanitize_index, idx))
    idx = tuple(map(normalize_slice, idx, none_shape))
    idx = posify_index(none_shape, idx)
    return idx


def check_index(ind, dimension):
    """ Check validity of index for a given dimension

    Examples
    --------
    >>> check_index(3, 5)
    >>> check_index(5, 5)
    Traceback (most recent call last):
    ...
    IndexError: Index is not smaller than dimension 5 >= 5

    >>> check_index(6, 5)
    Traceback (most recent call last):
    ...
    IndexError: Index is not smaller than dimension 6 >= 5

    >>> check_index(-1, 5)
    >>> check_index(-6, 5)
    Traceback (most recent call last):
    ...
    IndexError: Negative index is not greater than negative dimension -6 <= -5

    >>> check_index([1, 2], 5)
    >>> check_index([6, 3], 5)
    Traceback (most recent call last):
    ...
    IndexError: Index out of bounds 5

    >>> check_index(slice(0, 3), 5)
    """
    # unknown dimension, assumed to be in bounds
    if np.isnan(dimension):
        return
    elif isinstance(ind, (list, np.ndarray)):
        x = np.asanyarray(ind)
        if (x >= dimension).any() or (x < -dimension).any():
            raise IndexError("Index out of bounds %s" % dimension)
    elif isinstance(ind, slice):
        return
    elif is_dask_collection(ind):
        return
    elif ind is None:
        return

    elif ind >= dimension:
        raise IndexError("Index is not smaller than dimension %d >= %d" %
                         (ind, dimension))

    elif ind < -dimension:
        msg = "Negative index is not greater than negative dimension %d <= -%d"
        raise IndexError(msg % (ind, dimension))


def slice_with_dask_array(x, index):
    from .core import Array, atop, elemwise

    out_index = [slice(None)
                 if isinstance(ind, Array) and ind.dtype == bool
                 else ind
                 for ind in index]

    if len(index) == 1 and index[0].ndim == x.ndim:
        y = elemwise(getitem, x, *index, dtype=x.dtype)
        name = 'getitem-' + tokenize(x, index)
        dsk = {(name, i): k for i, k in enumerate(core.flatten(y.__dask_keys__()))}
        chunks = ((np.nan,) * y.npartitions,)
        return (Array(sharedict.merge(y.dask, (name, dsk)), name, chunks, x.dtype),
                out_index)

    if any(isinstance(ind, Array) and ind.dtype == bool and ind.ndim != 1
            for ind in index):
        raise NotImplementedError("Slicing with dask.array only permitted when "
                                  "the indexer has only one dimension or when "
                                  "it has the same dimension as the sliced "
                                  "array")
    indexes = [ind
               if isinstance(ind, Array) and ind.dtype == bool
               else slice(None)
               for ind in index]

    arginds = []
    i = 0
    for ind in indexes:
        if isinstance(ind, Array) and ind.dtype == bool:
            new = (ind, tuple(range(i, i + ind.ndim)))
            i += x.ndim
        else:
            new = (slice(None), None)
            i += 1
        arginds.append(new)

    arginds = list(concat(arginds))

    out = atop(getitem_variadic, tuple(range(x.ndim)), x, tuple(range(x.ndim)), *arginds, dtype=x.dtype)

    chunks = []
    for ind, chunk in zip(index, out.chunks):
        if isinstance(ind, Array) and ind.dtype == bool:
            chunks.append((np.nan,) * len(chunk))
        else:
            chunks.append(chunk)
    out._chunks = tuple(chunks)
    return out, tuple(out_index)


def getitem_variadic(x, *index):
    return x[index]
