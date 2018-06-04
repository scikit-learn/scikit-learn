from __future__ import absolute_import, division, print_function

from itertools import product
from operator import mul

import numpy as np

from .core import Array
from ..base import tokenize
from ..core import flatten
from ..compatibility import reduce
from ..utils import M
from .. import sharedict


def reshape_rechunk(inshape, outshape, inchunks):
    assert all(isinstance(c, tuple) for c in inchunks)
    ii = len(inshape) - 1
    oi = len(outshape) - 1
    result_inchunks = [None for i in range(len(inshape))]
    result_outchunks = [None for i in range(len(outshape))]

    while ii >= 0 or oi >= 0:
        if inshape[ii] == outshape[oi]:
            result_inchunks[ii] = inchunks[ii]
            result_outchunks[oi] = inchunks[ii]
            ii -= 1
            oi -= 1
            continue
        din = inshape[ii]
        dout = outshape[oi]
        if din == 1:
            result_inchunks[ii] = (1,)
            ii -= 1
        elif dout == 1:
            result_outchunks[oi] = (1,)
            oi -= 1
        elif din < dout:  # (4, 4, 4) -> (64,)
            ileft = ii - 1
            while ileft >= 0 and reduce(mul, inshape[ileft:ii + 1]) < dout: # 4 < 64, 4*4 < 64, 4*4*4 == 64
                ileft -= 1
            if reduce(mul, inshape[ileft:ii + 1]) != dout:
                raise ValueError("Shapes not compatible")

            for i in range(ileft + 1, ii + 1):  # need single-shape dimensions
                result_inchunks[i] = (inshape[i],)  # chunks[i] = (4,)

            chunk_reduction = reduce(mul, map(len, inchunks[ileft + 1:ii + 1]))
            result_inchunks[ileft] = expand_tuple(inchunks[ileft], chunk_reduction)

            prod = reduce(mul, inshape[ileft + 1: ii + 1])  # 16
            result_outchunks[oi] = tuple(prod * c for c in result_inchunks[ileft]) # (1, 1, 1, 1) .* 16

            oi -= 1
            ii = ileft - 1
        elif din > dout:  # (64,) -> (4, 4, 4)
            oleft = oi - 1
            while oleft >= 0 and reduce(mul, outshape[oleft:oi + 1]) < din:
                oleft -= 1
            if reduce(mul, outshape[oleft:oi + 1]) != din:
                raise ValueError("Shapes not compatible")

            # TODO: don't coalesce shapes unnecessarily
            cs = reduce(mul, outshape[oleft + 1: oi + 1])

            result_inchunks[ii] = contract_tuple(inchunks[ii], cs) # (16, 16, 16, 16)

            for i in range(oleft + 1, oi + 1):
                result_outchunks[i] = (outshape[i],)

            result_outchunks[oleft] = tuple(c // cs for c in result_inchunks[ii])

            oi = oleft - 1
            ii -= 1

    return tuple(result_inchunks), tuple(result_outchunks)


def expand_tuple(chunks, factor):
    """

    >>> expand_tuple((2, 4), 2)
    (1, 1, 2, 2)

    >>> expand_tuple((2, 4), 3)
    (1, 1, 1, 1, 2)

    >>> expand_tuple((3, 4), 2)
    (1, 2, 2, 2)

    >>> expand_tuple((7, 4), 3)
    (2, 2, 3, 1, 1, 2)
    """
    if factor == 1:
        return chunks

    out = []
    for c in chunks:
        x = c
        part = max(x / factor, 1)
        while x >= 2 * part:
            out.append(int(part))
            x -= int(part)
        if x:
            out.append(x)
    assert sum(chunks) == sum(out)
    return tuple(out)


def contract_tuple(chunks, factor):
    """ Return simple chunks tuple such that factor divides all elements

    Examples
    --------

    >>> contract_tuple((2, 2, 8, 4), 4)
    (4, 8, 4)
    """
    assert sum(chunks) % factor == 0

    out = []
    residual = 0
    for chunk in chunks:
        chunk += residual
        div = chunk // factor
        residual = chunk % factor
        good = factor * div
        if good:
            out.append(good)
    return tuple(out)


def reshape(x, shape):
    """ Reshape array to new shape

    This is a parallelized version of the ``np.reshape`` function with the
    following limitations:

    1.  It assumes that the array is stored in C-order
    2.  It only allows for reshapings that collapse or merge dimensions like
        ``(1, 2, 3, 4) -> (1, 6, 4)`` or ``(64,) -> (4, 4, 4)``

    When communication is necessary this algorithm depends on the logic within
    rechunk.  It endeavors to keep chunk sizes roughly the same when possible.

    See Also
    --------
    dask.array.rechunk
    numpy.reshape
    """
    # Sanitize inputs, look for -1 in shape
    from .slicing import sanitize_index
    shape = tuple(map(sanitize_index, shape))
    known_sizes = [s for s in shape if s != -1]
    if len(known_sizes) < len(shape):
        if len(known_sizes) - len(shape) > 1:
            raise ValueError('can only specify one unknown dimension')
        missing_size = sanitize_index(x.size / reduce(mul, known_sizes, 1))
        shape = tuple(missing_size if s == -1 else s for s in shape)

    if np.isnan(sum(x.shape)):
        raise ValueError("Array chunk size or shape is unknown. shape: %s", x.shape)

    if reduce(mul, shape, 1) != x.size:
        raise ValueError('total size of new array must be unchanged')

    if x.shape == shape:
        return x

    name = 'reshape-' + tokenize(x, shape)

    if x.npartitions == 1:
        key = next(flatten(x.__dask_keys__()))
        dsk = {(name,) + (0,) * len(shape): (M.reshape, key, shape)}
        chunks = tuple((d,) for d in shape)
        return Array(sharedict.merge((name, dsk), x.dask), name, chunks,
                     dtype=x.dtype)

    # Logic for how to rechunk
    inchunks, outchunks = reshape_rechunk(x.shape, shape, x.chunks)
    x2 = x.rechunk(inchunks)

    # Construct graph
    in_keys = list(product([x2.name], *[range(len(c)) for c in inchunks]))
    out_keys = list(product([name], *[range(len(c)) for c in outchunks]))
    shapes = list(product(*outchunks))
    dsk = {a: (M.reshape, b, shape) for a, b, shape in zip(out_keys, in_keys, shapes)}

    return Array(sharedict.merge((name, dsk), x2.dask), name, outchunks,
                 dtype=x.dtype)
