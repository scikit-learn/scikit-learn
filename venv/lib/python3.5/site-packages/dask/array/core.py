from __future__ import absolute_import, division, print_function

from bisect import bisect
from collections import Iterable, Mapping
from collections import Iterator
from functools import partial, wraps
from itertools import product
import math
from numbers import Number
import operator
from operator import add, getitem, mul
import os
import sys
import traceback
import pickle
from threading import Lock
import uuid
import warnings

try:
    from cytoolz import (partition, concat, join, first,
                         groupby, valmap, accumulate, assoc)
    from cytoolz.curried import filter, pluck

except ImportError:
    from toolz import (partition, concat, join, first,
                       groupby, valmap, accumulate, assoc)
    from toolz.curried import filter, pluck
from toolz import pipe, map, reduce
import numpy as np

from . import chunk
from .numpy_compat import _make_sliced_dtype
from .slicing import slice_array, replace_ellipsis
from ..base import (DaskMethodsMixin, tokenize, dont_optimize,
                    compute_as_if_collection, persist, is_dask_collection)
from ..context import _globals, globalmethod
from ..utils import (homogeneous_deepmap, ndeepmap, ignoring, concrete,
                     is_integer, IndexCallable, funcname, derived_from,
                     SerializableLock, ensure_dict, Dispatch)
from ..compatibility import unicode, long, getargspec, zip_longest, apply
from ..core import quote
from ..delayed import Delayed, to_task_dask
from .. import threaded, core
from .. import sharedict
from ..sharedict import ShareDict
from .numpy_compat import _Recurser


concatenate_lookup = Dispatch('concatenate')
tensordot_lookup = Dispatch('tensordot')
concatenate_lookup.register((object, np.ndarray), np.concatenate)
tensordot_lookup.register((object, np.ndarray), np.tensordot)


@tensordot_lookup.register_lazy('sparse')
@concatenate_lookup.register_lazy('sparse')
def register_sparse():
    import sparse
    concatenate_lookup.register(sparse.COO, sparse.concatenate)
    tensordot_lookup.register(sparse.COO, sparse.tensordot)


def getter(a, b, asarray=True, lock=None):
    if isinstance(b, tuple) and any(x is None for x in b):
        b2 = tuple(x for x in b if x is not None)
        b3 = tuple(None if x is None else slice(None, None)
                   for x in b if not isinstance(x, (int, long)))
        return getter(a, b2, asarray=asarray, lock=lock)[b3]

    if lock:
        lock.acquire()
    try:
        c = a[b]
        if asarray:
            c = np.asarray(c)
    finally:
        if lock:
            lock.release()
    return c


def getter_nofancy(a, b, asarray=True, lock=None):
    """ A simple wrapper around ``getter``.

    Used to indicate to the optimization passes that the backend doesn't
    support fancy indexing.
    """
    return getter(a, b, asarray=asarray, lock=lock)


def getter_inline(a, b, asarray=True, lock=None):
    """ A getter function that optimizations feel comfortable inlining

    Slicing operations with this function may be inlined into a graph, such as
    in the following rewrite

    **Before**

    >>> a = x[:10]  # doctest: +SKIP
    >>> b = a + 1  # doctest: +SKIP
    >>> c = a * 2  # doctest: +SKIP

    **After**

    >>> b = x[:10] + 1  # doctest: +SKIP
    >>> c = x[:10] * 2  # doctest: +SKIP

    This inlining can be relevant to operations when running off of disk.
    """
    return getter(a, b, asarray=asarray, lock=lock)


from .optimization import optimize, fuse_slice


def slices_from_chunks(chunks):
    """ Translate chunks tuple to a set of slices in product order

    >>> slices_from_chunks(((2, 2), (3, 3, 3)))  # doctest: +NORMALIZE_WHITESPACE
     [(slice(0, 2, None), slice(0, 3, None)),
      (slice(0, 2, None), slice(3, 6, None)),
      (slice(0, 2, None), slice(6, 9, None)),
      (slice(2, 4, None), slice(0, 3, None)),
      (slice(2, 4, None), slice(3, 6, None)),
      (slice(2, 4, None), slice(6, 9, None))]
    """
    cumdims = [list(accumulate(add, (0,) + bds[:-1])) for bds in chunks]
    shapes = product(*chunks)
    starts = product(*cumdims)
    return [tuple(slice(s, s + dim) for s, dim in zip(start, shape))
            for start, shape in zip(starts, shapes)]


def getem(arr, chunks, getitem=getter, shape=None, out_name=None, lock=False,
          asarray=True):
    """ Dask getting various chunks from an array-like

    >>> getem('X', chunks=(2, 3), shape=(4, 6))  # doctest: +SKIP
    {('X', 0, 0): (getter, 'X', (slice(0, 2), slice(0, 3))),
     ('X', 1, 0): (getter, 'X', (slice(2, 4), slice(0, 3))),
     ('X', 1, 1): (getter, 'X', (slice(2, 4), slice(3, 6))),
     ('X', 0, 1): (getter, 'X', (slice(0, 2), slice(3, 6)))}

    >>> getem('X', chunks=((2, 2), (3, 3)))  # doctest: +SKIP
    {('X', 0, 0): (getter, 'X', (slice(0, 2), slice(0, 3))),
     ('X', 1, 0): (getter, 'X', (slice(2, 4), slice(0, 3))),
     ('X', 1, 1): (getter, 'X', (slice(2, 4), slice(3, 6))),
     ('X', 0, 1): (getter, 'X', (slice(0, 2), slice(3, 6)))}
    """
    out_name = out_name or arr
    chunks = normalize_chunks(chunks, shape)

    keys = list(product([out_name], *[range(len(bds)) for bds in chunks]))
    slices = slices_from_chunks(chunks)

    if not asarray or lock:
        values = [(getitem, arr, x, asarray, lock) for x in slices]
    else:
        # Common case, drop extra parameters
        values = [(getitem, arr, x) for x in slices]

    return dict(zip(keys, values))


def dotmany(A, B, leftfunc=None, rightfunc=None, **kwargs):
    """ Dot product of many aligned chunks

    >>> x = np.array([[1, 2], [1, 2]])
    >>> y = np.array([[10, 20], [10, 20]])
    >>> dotmany([x, x, x], [y, y, y])
    array([[ 90, 180],
           [ 90, 180]])

    Optionally pass in functions to apply to the left and right chunks

    >>> dotmany([x, x, x], [y, y, y], rightfunc=np.transpose)
    array([[150, 150],
           [150, 150]])
    """
    if leftfunc:
        A = map(leftfunc, A)
    if rightfunc:
        B = map(rightfunc, B)
    return sum(map(partial(np.dot, **kwargs), A, B))


def lol_tuples(head, ind, values, dummies):
    """ List of list of tuple keys

    Parameters
    ----------

    head : tuple
        The known tuple so far
    ind : Iterable
        An iterable of indices not yet covered
    values : dict
        Known values for non-dummy indices
    dummies : dict
        Ranges of values for dummy indices

    Examples
    --------

    >>> lol_tuples(('x',), 'ij', {'i': 1, 'j': 0}, {})
    ('x', 1, 0)

    >>> lol_tuples(('x',), 'ij', {'i': 1}, {'j': range(3)})
    [('x', 1, 0), ('x', 1, 1), ('x', 1, 2)]

    >>> lol_tuples(('x',), 'ij', {'i': 1}, {'j': range(3)})
    [('x', 1, 0), ('x', 1, 1), ('x', 1, 2)]

    >>> lol_tuples(('x',), 'ijk', {'i': 1}, {'j': [0, 1, 2], 'k': [0, 1]}) # doctest: +NORMALIZE_WHITESPACE
    [[('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 1, 0), ('x', 1, 1, 1)],
     [('x', 1, 2, 0), ('x', 1, 2, 1)]]
    """
    if not ind:
        return head
    if ind[0] not in dummies:
        return lol_tuples(head + (values[ind[0]],), ind[1:], values, dummies)
    else:
        return [lol_tuples(head + (v,), ind[1:], values, dummies)
                for v in dummies[ind[0]]]


def zero_broadcast_dimensions(lol, nblocks):
    """

    >>> lol = [('x', 1, 0), ('x', 1, 1), ('x', 1, 2)]
    >>> nblocks = (4, 1, 2)  # note singleton dimension in second place
    >>> lol = [[('x', 1, 0, 0), ('x', 1, 0, 1)],
    ...        [('x', 1, 1, 0), ('x', 1, 1, 1)],
    ...        [('x', 1, 2, 0), ('x', 1, 2, 1)]]

    >>> zero_broadcast_dimensions(lol, nblocks)  # doctest: +NORMALIZE_WHITESPACE
    [[('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 0, 0), ('x', 1, 0, 1)]]

    See Also
    --------
    lol_tuples
    """
    f = lambda t: (t[0],) + tuple(0 if d == 1 else i for i, d in zip(t[1:], nblocks))
    return homogeneous_deepmap(f, lol)


def broadcast_dimensions(argpairs, numblocks, sentinels=(1, (1,)),
                         consolidate=None):
    """ Find block dimensions from arguments

    Parameters
    ----------
    argpairs: iterable
        name, ijk index pairs
    numblocks: dict
        maps {name: number of blocks}
    sentinels: iterable (optional)
        values for singleton dimensions
    consolidate: func (optional)
        use this to reduce each set of common blocks into a smaller set

    Examples
    --------
    >>> argpairs = [('x', 'ij'), ('y', 'ji')]
    >>> numblocks = {'x': (2, 3), 'y': (3, 2)}
    >>> broadcast_dimensions(argpairs, numblocks)
    {'i': 2, 'j': 3}

    Supports numpy broadcasting rules

    >>> argpairs = [('x', 'ij'), ('y', 'ij')]
    >>> numblocks = {'x': (2, 1), 'y': (1, 3)}
    >>> broadcast_dimensions(argpairs, numblocks)
    {'i': 2, 'j': 3}

    Works in other contexts too

    >>> argpairs = [('x', 'ij'), ('y', 'ij')]
    >>> d = {'x': ('Hello', 1), 'y': (1, (2, 3))}
    >>> broadcast_dimensions(argpairs, d)
    {'i': 'Hello', 'j': (2, 3)}
    """
    # List like [('i', 2), ('j', 1), ('i', 1), ('j', 2)]
    argpairs2 = [(a, ind) for a, ind in argpairs if ind is not None]
    L = concat([zip(inds, dims) for (x, inds), (x, dims)
                in join(first, argpairs2, first, numblocks.items())])

    g = groupby(0, L)
    g = dict((k, set([d for i, d in v])) for k, v in g.items())

    g2 = dict((k, v - set(sentinels) if len(v) > 1 else v) for k, v in g.items())

    if consolidate:
        return valmap(consolidate, g2)

    if g2 and not set(map(len, g2.values())) == set([1]):
        raise ValueError("Shapes do not align %s" % g)

    return valmap(first, g2)


def top(func, output, out_indices, *arrind_pairs, **kwargs):
    """ Tensor operation

    Applies a function, ``func``, across blocks from many different input
    dasks.  We arrange the pattern with which those blocks interact with sets
    of matching indices.  E.g.::

        top(func, 'z', 'i', 'x', 'i', 'y', 'i')

    yield an embarrassingly parallel communication pattern and is read as

        $$ z_i = func(x_i, y_i) $$

    More complex patterns may emerge, including multiple indices::

        top(func, 'z', 'ij', 'x', 'ij', 'y', 'ji')

        $$ z_{ij} = func(x_{ij}, y_{ji}) $$

    Indices missing in the output but present in the inputs results in many
    inputs being sent to one function (see examples).

    Examples
    --------

    Simple embarrassing map operation

    >>> inc = lambda x: x + 1
    >>> top(inc, 'z', 'ij', 'x', 'ij', numblocks={'x': (2, 2)})  # doctest: +SKIP
    {('z', 0, 0): (inc, ('x', 0, 0)),
     ('z', 0, 1): (inc, ('x', 0, 1)),
     ('z', 1, 0): (inc, ('x', 1, 0)),
     ('z', 1, 1): (inc, ('x', 1, 1))}

    Simple operation on two datasets

    >>> add = lambda x, y: x + y
    >>> top(add, 'z', 'ij', 'x', 'ij', 'y', 'ij', numblocks={'x': (2, 2),
    ...                                                      'y': (2, 2)})  # doctest: +SKIP
    {('z', 0, 0): (add, ('x', 0, 0), ('y', 0, 0)),
     ('z', 0, 1): (add, ('x', 0, 1), ('y', 0, 1)),
     ('z', 1, 0): (add, ('x', 1, 0), ('y', 1, 0)),
     ('z', 1, 1): (add, ('x', 1, 1), ('y', 1, 1))}

    Operation that flips one of the datasets

    >>> addT = lambda x, y: x + y.T  # Transpose each chunk
    >>> #                                        z_ij ~ x_ij y_ji
    >>> #               ..         ..         .. notice swap
    >>> top(addT, 'z', 'ij', 'x', 'ij', 'y', 'ji', numblocks={'x': (2, 2),
    ...                                                       'y': (2, 2)})  # doctest: +SKIP
    {('z', 0, 0): (add, ('x', 0, 0), ('y', 0, 0)),
     ('z', 0, 1): (add, ('x', 0, 1), ('y', 1, 0)),
     ('z', 1, 0): (add, ('x', 1, 0), ('y', 0, 1)),
     ('z', 1, 1): (add, ('x', 1, 1), ('y', 1, 1))}

    Dot product with contraction over ``j`` index.  Yields list arguments

    >>> top(dotmany, 'z', 'ik', 'x', 'ij', 'y', 'jk', numblocks={'x': (2, 2),
    ...                                                          'y': (2, 2)})  # doctest: +SKIP
    {('z', 0, 0): (dotmany, [('x', 0, 0), ('x', 0, 1)],
                            [('y', 0, 0), ('y', 1, 0)]),
     ('z', 0, 1): (dotmany, [('x', 0, 0), ('x', 0, 1)],
                            [('y', 0, 1), ('y', 1, 1)]),
     ('z', 1, 0): (dotmany, [('x', 1, 0), ('x', 1, 1)],
                            [('y', 0, 0), ('y', 1, 0)]),
     ('z', 1, 1): (dotmany, [('x', 1, 0), ('x', 1, 1)],
                            [('y', 0, 1), ('y', 1, 1)])}

    Pass ``concatenate=True`` to concatenate arrays ahead of time

    >>> top(f, 'z', 'i', 'x', 'ij', 'y', 'ij', concatenate=True,
    ...     numblocks={'x': (2, 2), 'y': (2, 2,)})  # doctest: +SKIP
    {('z', 0): (f, (concatenate_axes, [('x', 0, 0), ('x', 0, 1)], (1,)),
                   (concatenate_axes, [('y', 0, 0), ('y', 0, 1)], (1,)))
     ('z', 1): (f, (concatenate_axes, [('x', 1, 0), ('x', 1, 1)], (1,)),
                   (concatenate_axes, [('y', 1, 0), ('y', 1, 1)], (1,)))}

    Supports Broadcasting rules

    >>> top(add, 'z', 'ij', 'x', 'ij', 'y', 'ij', numblocks={'x': (1, 2),
    ...                                                      'y': (2, 2)})  # doctest: +SKIP
    {('z', 0, 0): (add, ('x', 0, 0), ('y', 0, 0)),
     ('z', 0, 1): (add, ('x', 0, 1), ('y', 0, 1)),
     ('z', 1, 0): (add, ('x', 0, 0), ('y', 1, 0)),
     ('z', 1, 1): (add, ('x', 0, 1), ('y', 1, 1))}

    Support keyword arguments with apply

    >>> def f(a, b=0): return a + b
    >>> top(f, 'z', 'i', 'x', 'i', numblocks={'x': (2,)}, b=10)  # doctest: +SKIP
    {('z', 0): (apply, f, [('x', 0)], {'b': 10}),
     ('z', 1): (apply, f, [('x', 1)], {'b': 10})}

    Include literals by indexing with ``None``

    >>> top(add, 'z', 'i', 'x', 'i', 100, None,  numblocks={'x': (2,)})  # doctest: +SKIP
    {('z', 0): (add, ('x', 0), 100),
     ('z', 1): (add, ('x', 1), 100)}


    See Also
    --------
    atop
    """
    numblocks = kwargs.pop('numblocks')
    concatenate = kwargs.pop('concatenate', None)
    new_axes = kwargs.pop('new_axes', {})
    argpairs = list(partition(2, arrind_pairs))

    assert set(numblocks) == {name for name, ind in argpairs if ind is not None}

    all_indices = pipe(argpairs, pluck(1), filter(None), concat, set)
    dummy_indices = all_indices - set(out_indices)

    # Dictionary mapping {i: 3, j: 4, ...} for i, j, ... the dimensions
    dims = broadcast_dimensions(argpairs, numblocks)
    for k in new_axes:
        dims[k] = 1

    # (0, 0), (0, 1), (0, 2), (1, 0), ...
    keytups = list(product(*[range(dims[i]) for i in out_indices]))
    # {i: 0, j: 0}, {i: 0, j: 1}, ...
    keydicts = [dict(zip(out_indices, tup)) for tup in keytups]

    # {j: [1, 2, 3], ...}  For j a dummy index of dimension 3
    dummies = dict((i, list(range(dims[i]))) for i in dummy_indices)

    dsk = {}

    # Unpack dask values in non-array arguments
    for i, (arg, ind) in enumerate(argpairs):
        if ind is None:
            arg2, dsk2 = to_task_dask(arg)
            if dsk2:
                dsk.update(ensure_dict(dsk2))
                argpairs[i] = (arg2, ind)

    # Create argument lists
    valtups = []
    for kd in keydicts:
        args = []
        for arg, ind in argpairs:
            if ind is None:
                args.append(arg)
            else:
                tups = lol_tuples((arg,), ind, kd, dummies)
                if any(nb == 1 for nb in numblocks[arg]):
                    tups2 = zero_broadcast_dimensions(tups, numblocks[arg])
                else:
                    tups2 = tups
                if concatenate and isinstance(tups2, list):
                    axes = [n for n, i in enumerate(ind) if i in dummies]
                    tups2 = (concatenate_axes, tups2, axes)
                args.append(tups2)
        valtups.append(args)

    if not kwargs:  # will not be used in an apply, should be a tuple
        valtups = [tuple(vt) for vt in valtups]

    # Add heads to tuples
    keys = [(output,) + kt for kt in keytups]

    # Unpack delayed objects in kwargs
    if kwargs:
        task, dsk2 = to_task_dask(kwargs)
        if dsk2:
            dsk.update(ensure_dict(dsk2))
            kwargs2 = task
        else:
            kwargs2 = kwargs
        vals = [(apply, func, vt, kwargs2) for vt in valtups]
    else:
        vals = [(func,) + vt for vt in valtups]

    dsk.update(dict(zip(keys, vals)))

    return dsk


def _concatenate2(arrays, axes=[]):
    """ Recursively Concatenate nested lists of arrays along axes

    Each entry in axes corresponds to each level of the nested list.  The
    length of axes should correspond to the level of nesting of arrays.

    >>> x = np.array([[1, 2], [3, 4]])
    >>> _concatenate2([x, x], axes=[0])
    array([[1, 2],
           [3, 4],
           [1, 2],
           [3, 4]])

    >>> _concatenate2([x, x], axes=[1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])

    >>> _concatenate2([[x, x], [x, x]], axes=[0, 1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4],
           [1, 2, 1, 2],
           [3, 4, 3, 4]])

    Supports Iterators
    >>> _concatenate2(iter([x, x]), axes=[1])
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])
    """
    if isinstance(arrays, Iterator):
        arrays = list(arrays)
    if not isinstance(arrays, (list, tuple)):
        return arrays
    if len(axes) > 1:
        arrays = [_concatenate2(a, axes=axes[1:]) for a in arrays]
    concatenate = concatenate_lookup.dispatch(type(max(arrays, key=lambda x: x.__array_priority__)))
    return concatenate(arrays, axis=axes[0])


def apply_infer_dtype(func, args, kwargs, funcname, suggest_dtype=True):
    args = [np.ones((1,) * x.ndim, dtype=x.dtype)
            if isinstance(x, Array) else x for x in args]
    try:
        with np.errstate(all='ignore'):
            o = func(*args, **kwargs)
    except Exception as e:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = ''.join(traceback.format_tb(exc_traceback))
        suggest = ("Please specify the dtype explicitly using the "
                   "`dtype` kwarg.\n\n") if suggest_dtype else ""
        msg = ("`dtype` inference failed in `{0}`.\n\n"
               "{1}"
               "Original error is below:\n"
               "------------------------\n"
               "{2}\n\n"
               "Traceback:\n"
               "---------\n"
               "{3}").format(funcname, suggest, repr(e), tb)
    else:
        msg = None
    if msg is not None:
        raise ValueError(msg)
    return o.dtype


def map_blocks(func, *args, **kwargs):
    """ Map a function across all blocks of a dask array.

    Parameters
    ----------
    func : callable
        Function to apply to every block in the array.
    args : dask arrays or constants
    dtype : np.dtype, optional
        The ``dtype`` of the output array. It is recommended to provide this.
        If not provided, will be inferred by applying the function to a small
        set of fake data.
    chunks : tuple, optional
        Chunk shape of resulting blocks if the function does not preserve
        shape. If not provided, the resulting array is assumed to have the same
        block structure as the first input array.
    drop_axis : number or iterable, optional
        Dimensions lost by the function.
    new_axis : number or iterable, optional
        New dimensions created by the function. Note that these are applied
        after ``drop_axis`` (if present).
    token : string, optional
        The key prefix to use for the output array. If not provided, will be
        determined from the function name.
    name : string, optional
        The key name to use for the output array. Note that this fully
        specifies the output key name, and must be unique. If not provided,
        will be determined by a hash of the arguments.
    **kwargs :
        Other keyword arguments to pass to function. Values must be constants
        (not dask.arrays)

    Examples
    --------
    >>> import dask.array as da
    >>> x = da.arange(6, chunks=3)

    >>> x.map_blocks(lambda x: x * 2).compute()
    array([ 0,  2,  4,  6,  8, 10])

    The ``da.map_blocks`` function can also accept multiple arrays.

    >>> d = da.arange(5, chunks=2)
    >>> e = da.arange(5, chunks=2)

    >>> f = map_blocks(lambda a, b: a + b**2, d, e)
    >>> f.compute()
    array([ 0,  2,  6, 12, 20])

    If the function changes shape of the blocks then you must provide chunks
    explicitly.

    >>> y = x.map_blocks(lambda x: x[::2], chunks=((2, 2),))

    You have a bit of freedom in specifying chunks.  If all of the output chunk
    sizes are the same, you can provide just that chunk size as a single tuple.

    >>> a = da.arange(18, chunks=(6,))
    >>> b = a.map_blocks(lambda x: x[:3], chunks=(3,))

    If the function changes the dimension of the blocks you must specify the
    created or destroyed dimensions.

    >>> b = a.map_blocks(lambda x: x[None, :, None], chunks=(1, 6, 1),
    ...                  new_axis=[0, 2])

    Map_blocks aligns blocks by block positions without regard to shape. In the
    following example we have two arrays with the same number of blocks but
    with different shape and chunk sizes.

    >>> x = da.arange(1000, chunks=(100,))
    >>> y = da.arange(100, chunks=(10,))

    The relevant attribute to match is numblocks.

    >>> x.numblocks
    (10,)
    >>> y.numblocks
    (10,)

    If these match (up to broadcasting rules) then we can map arbitrary
    functions across blocks

    >>> def func(a, b):
    ...     return np.array([a.max(), b.max()])

    >>> da.map_blocks(func, x, y, chunks=(2,), dtype='i8')
    dask.array<func, shape=(20,), dtype=int64, chunksize=(2,)>

    >>> _.compute()
    array([ 99,   9, 199,  19, 299,  29, 399,  39, 499,  49, 599,  59, 699,
            69, 799,  79, 899,  89, 999,  99])

    Your block function can learn where in the array it is if it supports a
    ``block_id`` keyword argument.  This will receive entries like (2, 0, 1),
    the position of the block in the dask array.

    >>> def func(block, block_id=None):
    ...     pass

    You may specify the key name prefix of the resulting task in the graph with
    the optional ``token`` keyword argument.

    >>> x.map_blocks(lambda x: x + 1, token='increment')  # doctest: +SKIP
    dask.array<increment, shape=(100,), dtype=int64, chunksize=(10,)>
    """
    if not callable(func):
        msg = ("First argument must be callable function, not %s\n"
               "Usage:   da.map_blocks(function, x)\n"
               "   or:   da.map_blocks(function, x, y, z)")
        raise TypeError(msg % type(func).__name__)
    name = kwargs.pop('name', None)
    token = kwargs.pop('token', None)
    if not name:
        name = '%s-%s' % (token or funcname(func),
                          tokenize(token or func, args, **kwargs))
    dtype = kwargs.pop('dtype', None)
    chunks = kwargs.pop('chunks', None)
    drop_axis = kwargs.pop('drop_axis', [])
    new_axis = kwargs.pop('new_axis', [])
    if isinstance(drop_axis, Number):
        drop_axis = [drop_axis]
    if isinstance(new_axis, Number):
        new_axis = [new_axis]

    arrs = [a for a in args if isinstance(a, Array)]

    argpairs = [(a.name, tuple(range(a.ndim))[::-1])
                if isinstance(a, Array)
                else (a, None)
                for a in args]
    numblocks = {a.name: a.numblocks for a in arrs}
    arginds = list(concat(argpairs))
    out_ind = tuple(range(max(a.ndim for a in arrs)))[::-1]

    try:
        spec = getargspec(func)
        block_id = ('block_id' in spec.args or
                    'block_id' in getattr(spec, 'kwonly_args', ()))
    except Exception:
        block_id = False

    if block_id:
        kwargs['block_id'] = '__dummy__'

    dsk = top(func, name, out_ind, *arginds, numblocks=numblocks,
              **kwargs)

    # If func has block_id as an argument, add it to the kwargs for each call
    if block_id:
        for k in dsk.keys():
            dsk[k] = dsk[k][:-1] + (assoc(dsk[k][-1], 'block_id', k[1:]),)

    if dtype is None:
        if block_id:
            kwargs2 = assoc(kwargs, 'block_id', first(dsk.keys())[1:])
        else:
            kwargs2 = kwargs
        dtype = apply_infer_dtype(func, args, kwargs2, 'map_blocks')

    if len(arrs) == 1:
        numblocks = list(arrs[0].numblocks)
    else:
        dims = broadcast_dimensions(argpairs, numblocks)
        numblocks = [b for (_, b) in sorted(dims.items(), reverse=True)]

    if drop_axis:
        if any(numblocks[i] > 1 for i in drop_axis):
            raise ValueError("Can't drop an axis with more than 1 block. "
                             "Please use `atop` instead.")
        dsk = dict((tuple(k for i, k in enumerate(k)
                          if i - 1 not in drop_axis), v)
                   for k, v in dsk.items())
        numblocks = [n for i, n in enumerate(numblocks) if i not in drop_axis]
    if new_axis:
        new_axis = sorted(new_axis)
        for i in new_axis:
            if not 0 <= i <= len(numblocks):
                ndim = len(numblocks)
                raise ValueError("Can't add axis %d when current "
                                 "axis are %r. Missing axis: "
                                 "%r" % (i, list(range(ndim)),
                                         list(range(ndim, i))))
            numblocks.insert(i, 1)
        dsk, old_dsk = dict(), dsk
        for key in old_dsk:
            new_key = list(key)
            for i in new_axis:
                new_key.insert(i + 1, 0)
            dsk[tuple(new_key)] = old_dsk[key]

    if chunks:
        if len(chunks) != len(numblocks):
            raise ValueError("Provided chunks have {0} dims, expected {1} "
                             "dims.".format(len(chunks), len(numblocks)))
        chunks2 = []
        for i, (c, nb) in enumerate(zip(chunks, numblocks)):
            if isinstance(c, tuple):
                if not len(c) == nb:
                    raise ValueError("Dimension {0} has {1} blocks, "
                                     "chunks specified with "
                                     "{2} blocks".format(i, nb, len(c)))
                chunks2.append(c)
            else:
                chunks2.append(nb * (c,))
    else:
        if len(arrs) == 1:
            chunks2 = list(arrs[0].chunks)
        else:
            try:
                chunks2 = list(broadcast_chunks(*[a.chunks for a in arrs]))
            except Exception:
                raise ValueError("Arrays in `map_blocks` don't align, can't "
                                 "infer output chunks. Please provide "
                                 "`chunks` kwarg.")
        if drop_axis:
            chunks2 = [c for (i, c) in enumerate(chunks2) if i not in drop_axis]
        if new_axis:
            for i in sorted(new_axis):
                chunks2.insert(i, (1,))

    chunks = tuple(chunks2)

    return Array(sharedict.merge((name, dsk), *[a.dask for a in arrs]),
                 name, chunks, dtype)


def broadcast_chunks(*chunkss):
    """ Construct a chunks tuple that broadcasts many chunks tuples

    >>> a = ((5, 5),)
    >>> b = ((5, 5),)
    >>> broadcast_chunks(a, b)
    ((5, 5),)

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((5, 5),)
    >>> broadcast_chunks(a, b)
    ((10, 10, 10), (5, 5))

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((1,), (5, 5),)
    >>> broadcast_chunks(a, b)
    ((10, 10, 10), (5, 5))

    >>> a = ((10, 10, 10), (5, 5),)
    >>> b = ((3, 3,), (5, 5),)
    >>> broadcast_chunks(a, b)
    Traceback (most recent call last):
        ...
    ValueError: Chunks do not align: [(10, 10, 10), (3, 3)]
    """
    if not chunkss:
        return ()
    elif len(chunkss) == 1:
        return chunkss[0]
    n = max(map(len, chunkss))
    chunkss2 = [((1,),) * (n - len(c)) + c for c in chunkss]
    result = []
    for i in range(n):
        step1 = [c[i] for c in chunkss2]
        if all(c == (1,) for c in step1):
            step2 = step1
        else:
            step2 = [c for c in step1 if c != (1,)]
        if len(set(step2)) != 1:
            raise ValueError("Chunks do not align: %s" % str(step2))
        result.append(step2[0])
    return tuple(result)


def store(sources, targets, lock=True, regions=None, compute=True,
          return_stored=False, **kwargs):
    """ Store dask arrays in array-like objects, overwrite data in target

    This stores dask arrays into object that supports numpy-style setitem
    indexing.  It stores values chunk by chunk so that it does not have to
    fill up memory.  For best performance you can align the block size of
    the storage target with the block size of your array.

    If your data fits in memory then you may prefer calling
    ``np.array(myarray)`` instead.

    Parameters
    ----------

    sources: Array or iterable of Arrays
    targets: array-like or Delayed or iterable of array-likes and/or Delayeds
        These should support setitem syntax ``target[10:20] = ...``
    lock: boolean or threading.Lock, optional
        Whether or not to lock the data stores while storing.
        Pass True (lock each file individually), False (don't lock) or a
        particular ``threading.Lock`` object to be shared among all writes.
    regions: tuple of slices or iterable of tuple of slices
        Each ``region`` tuple in ``regions`` should be such that
        ``target[region].shape = source.shape``
        for the corresponding source and target in sources and targets, respectively.
    compute: boolean, optional
        If true compute immediately, return ``dask.delayed.Delayed`` otherwise
    return_stored: boolean, optional
        Optionally return the stored result (default False).

    Examples
    --------
    >>> x = ...  # doctest: +SKIP

    >>> import h5py  # doctest: +SKIP
    >>> f = h5py.File('myfile.hdf5')  # doctest: +SKIP
    >>> dset = f.create_dataset('/data', shape=x.shape,
    ...                                  chunks=x.chunks,
    ...                                  dtype='f8')  # doctest: +SKIP

    >>> store(x, dset)  # doctest: +SKIP

    Alternatively store many arrays at the same time

    >>> store([x, y, z], [dset1, dset2, dset3])  # doctest: +SKIP
    """

    if isinstance(sources, Array):
        sources = [sources]
        targets = [targets]

    if any(not isinstance(s, Array) for s in sources):
        raise ValueError("All sources must be dask array objects")

    if len(sources) != len(targets):
        raise ValueError("Different number of sources [%d] and targets [%d]"
                         % (len(sources), len(targets)))

    if isinstance(regions, tuple) or regions is None:
        regions = [regions]

    if len(sources) > 1 and len(regions) == 1:
        regions *= len(sources)

    if len(sources) != len(regions):
        raise ValueError("Different number of sources [%d] and targets [%d] than regions [%d]"
                         % (len(sources), len(targets), len(regions)))

    # Optimize all sources together
    sources_dsk = sharedict.merge(*[e.__dask_graph__() for e in sources])
    sources_dsk = Array.__dask_optimize__(
        sources_dsk,
        list(core.flatten([e.__dask_keys__() for e in sources]))
    )
    sources2 = [Array(sources_dsk, e.name, e.chunks, e.dtype) for e in sources]

    # Optimize all targets together
    targets2 = []
    targets_keys = []
    targets_dsk = []
    for e in targets:
        if isinstance(e, Delayed):
            targets2.append(e.key)
            targets_keys.extend(e.__dask_keys__())
            targets_dsk.append(e.__dask_graph__())
        elif is_dask_collection(e):
            raise TypeError(
                "Targets must be either Delayed objects or array-likes"
            )
        else:
            targets2.append(e)

    targets_dsk = sharedict.merge(*targets_dsk)
    targets_dsk = Delayed.__dask_optimize__(targets_dsk, targets_keys)

    load_stored = (return_stored and not compute)
    store_dsk = sharedict.merge(*[
        insert_to_ooc(s, t, lock, r, return_stored, load_stored)
        for s, t, r in zip(sources2, targets2, regions)
    ])
    store_keys = list(store_dsk.keys())

    store_dsk = sharedict.merge(store_dsk, targets_dsk, sources_dsk)

    if return_stored:
        load_store_dsk = store_dsk
        if compute:
            store_dlyds = [Delayed(k, store_dsk) for k in store_keys]
            store_dlyds = persist(*store_dlyds, **kwargs)
            store_dsk_2 = sharedict.merge(*[e.dask for e in store_dlyds])

            load_store_dsk = retrieve_from_ooc(
                store_keys, store_dsk, store_dsk_2
            )

        result = tuple(
            Array(load_store_dsk, 'load-store-%s' % s.name, s.chunks, s.dtype)
            for s in sources
        )

        return result
    else:
        name = 'store-' + tokenize(*store_keys)
        dsk = sharedict.merge({name: store_keys}, store_dsk)
        result = Delayed(name, dsk)

        if compute:
            result.compute(**kwargs)
            return None
        else:
            return result


def blockdims_from_blockshape(shape, chunks):
    """

    >>> blockdims_from_blockshape((10, 10), (4, 3))
    ((4, 4, 2), (3, 3, 3, 1))
    >>> blockdims_from_blockshape((10, 0), (4, 0))
    ((4, 4, 2), (0,))
    """
    if chunks is None:
        raise TypeError("Must supply chunks= keyword argument")
    if shape is None:
        raise TypeError("Must supply shape= keyword argument")
    if np.isnan(sum(shape)) or np.isnan(sum(chunks)):
        raise ValueError("Array chunk sizes are unknown. shape: %s, chunks: %s"
                         % (shape, chunks))
    if not all(map(is_integer, chunks)):
        raise ValueError("chunks can only contain integers.")
    if not all(map(is_integer, shape)):
        raise ValueError("shape can only contain integers.")
    shape = tuple(map(int, shape))
    chunks = tuple(map(int, chunks))
    return tuple(((bd,) * (d // bd) + ((d % bd,) if d % bd else ())
                 if d else (0,))
                 for d, bd in zip(shape, chunks))


def finalize(results):
    if not results:
        return concatenate3(results)
    results2 = results
    while isinstance(results2, (tuple, list)):
        if len(results2) > 1:
            return concatenate3(results)
        else:
            results2 = results2[0]
    return unpack_singleton(results)


CHUNKS_NONE_ERROR_MESSAGE = """
You must specify a chunks= keyword argument.
This specifies the chunksize of your array blocks.

See the following documentation page for details:
  http://dask.pydata.org/en/latest/array-creation.html#chunks
""".strip()


class Array(DaskMethodsMixin):
    """ Parallel Dask Array

    A parallel nd-array comprised of many numpy arrays arranged in a grid.

    This constructor is for advanced uses only.  For normal use see the
    ``da.from_array`` function.

    Parameters
    ----------

    dask : dict
        Task dependency graph
    name : string
        Name of array in dask
    shape : tuple of ints
        Shape of the entire array
    chunks: iterable of tuples
        block sizes along each dimension

    See Also
    --------
    dask.array.from_array
    """
    __slots__ = 'dask', '_name', '_cached_keys', '_chunks', 'dtype'

    def __new__(cls, dask, name, chunks, dtype, shape=None):
        self = super(Array, cls).__new__(cls)
        assert isinstance(dask, Mapping)
        if not isinstance(dask, ShareDict):
            s = ShareDict()
            s.update_with_key(dask, key=name)
            dask = s
        self.dask = dask
        self.name = name
        self._chunks = normalize_chunks(chunks, shape)
        if self._chunks is None:
            raise ValueError(CHUNKS_NONE_ERROR_MESSAGE)
        if dtype is None:
            raise ValueError("You must specify the dtype of the array")
        self.dtype = np.dtype(dtype)

        for plugin in _globals.get('array_plugins', ()):
            result = plugin(self)
            if result is not None:
                self = result

        return self

    def __reduce__(self):
        return (Array, (self.dask, self.name, self.chunks, self.dtype))

    def __dask_graph__(self):
        return self.dask

    def __dask_keys__(self):
        if self._cached_keys is not None:
            return self._cached_keys

        name, chunks, numblocks = self.name, self.chunks, self.numblocks

        def keys(*args):
            if not chunks:
                return [(name,)]
            ind = len(args)
            if ind + 1 == len(numblocks):
                result = [(name,) + args + (i,) for i in range(numblocks[ind])]
            else:
                result = [keys(*(args + (i,))) for i in range(numblocks[ind])]
            return result

        self._cached_keys = result = keys()
        return result

    def __dask_tokenize__(self):
        return self.name

    __dask_optimize__ = globalmethod(optimize, key='array_optimize',
                                     falsey=dont_optimize)
    __dask_scheduler__ = staticmethod(threaded.get)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return Array, (self.name, self.chunks, self.dtype)

    @property
    def numblocks(self):
        return tuple(map(len, self.chunks))

    @property
    def npartitions(self):
        return reduce(mul, self.numblocks, 1)

    @property
    def shape(self):
        return tuple(map(sum, self.chunks))

    @property
    def _meta(self):
        return np.empty(shape=(), dtype=self.dtype)

    def _get_chunks(self):
        return self._chunks

    def _set_chunks(self, chunks):
        raise TypeError("Can not set chunks directly\n\n"
                        "Please use the rechunk method instead:\n"
                        "  x.rechunk(%s)" % str(chunks))

    chunks = property(_get_chunks, _set_chunks, "chunks property")

    def __len__(self):
        if not self.chunks:
            raise TypeError("len() of unsized object")
        return sum(self.chunks[0])

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get('out', ())
        for x in inputs + out:
            if not isinstance(x, (np.ndarray, Number, Array)):
                return NotImplemented

        if method == '__call__':
            if numpy_ufunc.signature is not None:
                return NotImplemented
            if numpy_ufunc.nout > 1:
                from . import ufunc
                try:
                    da_ufunc = getattr(ufunc, numpy_ufunc.__name__)
                except AttributeError:
                    return NotImplemented
                return da_ufunc(*inputs, **kwargs)
            else:
                return elemwise(numpy_ufunc, *inputs, **kwargs)
        elif method == 'outer':
            from . import ufunc
            try:
                da_ufunc = getattr(ufunc, numpy_ufunc.__name__)
            except AttributeError:
                return NotImplemented
            return da_ufunc.outer(*inputs, **kwargs)
        else:
            return NotImplemented

    def __repr__(self):
        """

        >>> import dask.array as da
        >>> da.ones((10, 10), chunks=(5, 5), dtype='i4')
        dask.array<..., shape=(10, 10), dtype=int32, chunksize=(5, 5)>
        """
        chunksize = str(tuple(c[0] for c in self.chunks))
        name = self.name.rsplit('-', 1)[0]
        return ("dask.array<%s, shape=%s, dtype=%s, chunksize=%s>" %
                (name, self.shape, self.dtype, chunksize))

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        """ Number of elements in array """
        return reduce(mul, self.shape, 1)

    @property
    def nbytes(self):
        """ Number of bytes in array """
        return self.size * self.dtype.itemsize

    @property
    def itemsize(self):
        """ Length of one array element in bytes """
        return self.dtype.itemsize

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        self._name = val
        # Clear the key cache when the name is reset
        self._cached_keys = None

    __array_priority__ = 11  # higher than numpy.ndarray and numpy.matrix

    def __array__(self, dtype=None, **kwargs):
        x = self.compute()
        if dtype and x.dtype != dtype:
            x = x.astype(dtype)
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        return x

    @property
    def _elemwise(self):
        return elemwise

    @wraps(store)
    def store(self, target, **kwargs):
        r = store([self], [target], **kwargs)

        if kwargs.get("return_stored", False):
            r = r[0]

        return r

    def to_hdf5(self, filename, datapath, **kwargs):
        """ Store array in HDF5 file

        >>> x.to_hdf5('myfile.hdf5', '/x')  # doctest: +SKIP

        Optionally provide arguments as though to ``h5py.File.create_dataset``

        >>> x.to_hdf5('myfile.hdf5', '/x', compression='lzf', shuffle=True)  # doctest: +SKIP

        See Also
        --------
        da.store
        h5py.File.create_dataset
        """
        return to_hdf5(filename, datapath, self, **kwargs)

    def to_dask_dataframe(self, columns=None):
        """ Convert dask Array to dask Dataframe

        Parameters
        ----------
        columns: list or string
            list of column names if DataFrame, single string if Series

        See Also
        --------
        dask.dataframe.from_dask_array
        """
        from ..dataframe import from_dask_array
        return from_dask_array(self, columns=columns)

    def __bool__(self):
        if self.size > 1:
            raise ValueError("The truth value of a {0} is ambiguous. "
                             "Use a.any() or a.all()."
                             .format(self.__class__.__name__))
        else:
            return bool(self.compute())

    __nonzero__ = __bool__  # python 2

    def _scalarfunc(self, cast_type):
        if self.size > 1:
            raise TypeError("Only length-1 arrays can be converted "
                            "to Python scalars")
        else:
            return cast_type(self.compute())

    def __int__(self):
        return self._scalarfunc(int)

    __long__ = __int__  # python 2

    def __float__(self):
        return self._scalarfunc(float)

    def __complex__(self):
        return self._scalarfunc(complex)

    def __setitem__(self, key, value):
        from .routines import where
        if isinstance(key, Array):
            if isinstance(value, Array) and value.ndim > 1:
                raise ValueError('boolean index array should have 1 dimension')
            y = where(key, value, self)
            self.dtype = y.dtype
            self.dask = y.dask
            self.name = y.name
            return self
        else:
            raise NotImplementedError("Item assignment with %s not supported"
                                      % type(key))

    def __getitem__(self, index):
        # Field access, e.g. x['a'] or x[['a', 'b']]
        if (isinstance(index, (str, unicode)) or
                (isinstance(index, list) and index and
                 all(isinstance(i, (str, unicode)) for i in index))):
            out = 'getitem-' + tokenize(self, index)
            if isinstance(index, (str, unicode)):
                dt = self.dtype[index]
            else:
                dt = _make_sliced_dtype(self.dtype, index)

            if dt.shape:
                new_axis = list(range(self.ndim, self.ndim + len(dt.shape)))
                chunks = self.chunks + tuple((i,) for i in dt.shape)
                return self.map_blocks(getitem, index, dtype=dt.base, name=out,
                                       chunks=chunks, new_axis=new_axis)
            else:
                return self.map_blocks(getitem, index, dtype=dt, name=out)

        if not isinstance(index, tuple):
            index = (index,)

        from .slicing import normalize_index, slice_with_dask_array
        index2 = normalize_index(index, self.shape)

        if any(isinstance(i, Array) for i in index2):
            self, index2 = slice_with_dask_array(self, index2)

        if all(isinstance(i, slice) and i == slice(None) for i in index2):
            return self

        out = 'getitem-' + tokenize(self, index2)
        dsk, chunks = slice_array(out, self.name, self.chunks, index2)

        dsk2 = sharedict.merge(self.dask, (out, dsk))

        return Array(dsk2, out, chunks, dtype=self.dtype)

    def _vindex(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if any(k is None for k in key):
            raise IndexError(
                "vindex does not support indexing with None (np.newaxis), "
                "got {}".format(key))
        if all(isinstance(k, slice) for k in key):
            raise IndexError(
                "vindex requires at least one non-slice to vectorize over. "
                "Use normal slicing instead when only using slices. Got: {}"
                .format(key))
        return _vindex(self, *key)

    @property
    def vindex(self):
        """Vectorized indexing with broadcasting.

        This is equivalent to numpy's advanced indexing, using arrays that are
        broadcast against each other. This allows for pointwise indexing:

        >>> x = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> x = from_array(x, chunks=2)
        >>> x.vindex[[0, 1, 2], [0, 1, 2]].compute()
        array([1, 5, 9])

        Mixed basic/advanced indexing with slices/arrays is also supported. The
        order of dimensions in the result follows those proposed for
        ndarray.vindex [1]_: the subspace spanned by arrays is followed by all
        slices.

        Note: ``vindex`` provides more general functionality than standard
        indexing, but it also has fewer optimizations and can be significantly
        slower.

        _[1]: https://github.com/numpy/numpy/pull/6256
        """
        return IndexCallable(self._vindex)

    @derived_from(np.ndarray)
    def dot(self, other):
        from .routines import tensordot
        return tensordot(self, other,
                         axes=((self.ndim - 1,), (other.ndim - 2,)))

    @property
    def A(self):
        return self

    @property
    def T(self):
        return self.transpose()

    @derived_from(np.ndarray)
    def transpose(self, *axes):
        from .routines import transpose
        if not axes:
            axes = None
        elif len(axes) == 1 and isinstance(axes[0], Iterable):
            axes = axes[0]
        return transpose(self, axes=axes)

    @derived_from(np.ndarray)
    def ravel(self):
        from .routines import ravel
        return ravel(self)

    flatten = ravel

    @derived_from(np.ndarray)
    def choose(self, choices):
        from .routines import choose
        return choose(self, choices)

    @derived_from(np.ndarray)
    def reshape(self, *shape):
        from .reshape import reshape
        if len(shape) == 1 and not isinstance(shape[0], Number):
            shape = shape[0]
        return reshape(self, shape)

    def topk(self, k, axis=-1, split_every=None):
        """The top k elements of an array.

        See ``da.topk`` for docstring"""
        from .reductions import topk
        return topk(self, k, axis=axis, split_every=split_every)

    def argtopk(self, k, axis=-1, split_every=None):
        """The indices of the top k elements of an array.

        See ``da.argtopk`` for docstring"""
        from .reductions import argtopk
        return argtopk(self, k, axis=axis, split_every=split_every)

    def astype(self, dtype, **kwargs):
        """Copy of the array, cast to a specified type.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        casting : {'no', 'equiv', 'safe', 'same_kind', 'unsafe'}, optional
            Controls what kind of data casting may occur. Defaults to 'unsafe'
            for backwards compatibility.

            * 'no' means the data types should not be cast at all.
            * 'equiv' means only byte-order changes are allowed.
            * 'safe' means only casts which can preserve values are allowed.
            * 'same_kind' means only safe casts or casts within a kind,
                like float64 to float32, are allowed.
            * 'unsafe' means any data conversions may be done.
        copy : bool, optional
            By default, astype always returns a newly allocated array. If this
            is set to False and the `dtype` requirement is satisfied, the input
            array is returned instead of a copy.
        """
        # Scalars don't take `casting` or `copy` kwargs - as such we only pass
        # them to `map_blocks` if specified by user (different than defaults).
        extra = set(kwargs) - {'casting', 'copy'}
        if extra:
            raise TypeError("astype does not take the following keyword "
                            "arguments: {0!s}".format(list(extra)))
        casting = kwargs.get('casting', 'unsafe')
        copy = kwargs.get('copy', True)
        dtype = np.dtype(dtype)
        if self.dtype == dtype:
            return self
        elif not np.can_cast(self.dtype, dtype, casting=casting):
            raise TypeError("Cannot cast array from {0!r} to {1!r}"
                            " according to the rule "
                            "{2!r}".format(self.dtype, dtype, casting))
        name = 'astype-' + tokenize(self, dtype, casting, copy)
        return self.map_blocks(chunk.astype, dtype=dtype, name=name,
                               astype_dtype=dtype, **kwargs)

    def __abs__(self):
        return elemwise(operator.abs, self)

    def __add__(self, other):
        return elemwise(operator.add, self, other)

    def __radd__(self, other):
        return elemwise(operator.add, other, self)

    def __and__(self, other):
        return elemwise(operator.and_, self, other)

    def __rand__(self, other):
        return elemwise(operator.and_, other, self)

    def __div__(self, other):
        return elemwise(operator.div, self, other)

    def __rdiv__(self, other):
        return elemwise(operator.div, other, self)

    def __eq__(self, other):
        return elemwise(operator.eq, self, other)

    def __gt__(self, other):
        return elemwise(operator.gt, self, other)

    def __ge__(self, other):
        return elemwise(operator.ge, self, other)

    def __invert__(self):
        return elemwise(operator.invert, self)

    def __lshift__(self, other):
        return elemwise(operator.lshift, self, other)

    def __rlshift__(self, other):
        return elemwise(operator.lshift, other, self)

    def __lt__(self, other):
        return elemwise(operator.lt, self, other)

    def __le__(self, other):
        return elemwise(operator.le, self, other)

    def __mod__(self, other):
        return elemwise(operator.mod, self, other)

    def __rmod__(self, other):
        return elemwise(operator.mod, other, self)

    def __mul__(self, other):
        return elemwise(operator.mul, self, other)

    def __rmul__(self, other):
        return elemwise(operator.mul, other, self)

    def __ne__(self, other):
        return elemwise(operator.ne, self, other)

    def __neg__(self):
        return elemwise(operator.neg, self)

    def __or__(self, other):
        return elemwise(operator.or_, self, other)

    def __pos__(self):
        return self

    def __ror__(self, other):
        return elemwise(operator.or_, other, self)

    def __pow__(self, other):
        return elemwise(operator.pow, self, other)

    def __rpow__(self, other):
        return elemwise(operator.pow, other, self)

    def __rshift__(self, other):
        return elemwise(operator.rshift, self, other)

    def __rrshift__(self, other):
        return elemwise(operator.rshift, other, self)

    def __sub__(self, other):
        return elemwise(operator.sub, self, other)

    def __rsub__(self, other):
        return elemwise(operator.sub, other, self)

    def __truediv__(self, other):
        return elemwise(operator.truediv, self, other)

    def __rtruediv__(self, other):
        return elemwise(operator.truediv, other, self)

    def __floordiv__(self, other):
        return elemwise(operator.floordiv, self, other)

    def __rfloordiv__(self, other):
        return elemwise(operator.floordiv, other, self)

    def __xor__(self, other):
        return elemwise(operator.xor, self, other)

    def __rxor__(self, other):
        return elemwise(operator.xor, other, self)

    def __matmul__(self, other):
        from .routines import matmul
        return matmul(self, other)

    def __rmatmul__(self, other):
        from .routines import matmul
        return matmul(other, self)

    @derived_from(np.ndarray)
    def any(self, axis=None, keepdims=False, split_every=None, out=None):
        from .reductions import any
        return any(self, axis=axis, keepdims=keepdims, split_every=split_every,
                   out=out)

    @derived_from(np.ndarray)
    def all(self, axis=None, keepdims=False, split_every=None, out=None):
        from .reductions import all
        return all(self, axis=axis, keepdims=keepdims, split_every=split_every,
                   out=out)

    @derived_from(np.ndarray)
    def min(self, axis=None, keepdims=False, split_every=None, out=None):
        from .reductions import min
        return min(self, axis=axis, keepdims=keepdims, split_every=split_every,
                   out=out)

    @derived_from(np.ndarray)
    def max(self, axis=None, keepdims=False, split_every=None, out=None):
        from .reductions import max
        return max(self, axis=axis, keepdims=keepdims, split_every=split_every,
                   out=out)

    @derived_from(np.ndarray)
    def argmin(self, axis=None, split_every=None, out=None):
        from .reductions import argmin
        return argmin(self, axis=axis, split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def argmax(self, axis=None, split_every=None, out=None):
        from .reductions import argmax
        return argmax(self, axis=axis, split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def sum(self, axis=None, dtype=None, keepdims=False, split_every=None,
            out=None):
        from .reductions import sum
        return sum(self, axis=axis, dtype=dtype, keepdims=keepdims,
                   split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def prod(self, axis=None, dtype=None, keepdims=False, split_every=None,
             out=None):
        from .reductions import prod
        return prod(self, axis=axis, dtype=dtype, keepdims=keepdims,
                    split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def mean(self, axis=None, dtype=None, keepdims=False, split_every=None,
             out=None):
        from .reductions import mean
        return mean(self, axis=axis, dtype=dtype, keepdims=keepdims,
                    split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def std(self, axis=None, dtype=None, keepdims=False, ddof=0,
            split_every=None, out=None):
        from .reductions import std
        return std(self, axis=axis, dtype=dtype, keepdims=keepdims, ddof=ddof,
                   split_every=split_every, out=out)

    @derived_from(np.ndarray)
    def var(self, axis=None, dtype=None, keepdims=False, ddof=0,
            split_every=None, out=None):
        from .reductions import var
        return var(self, axis=axis, dtype=dtype, keepdims=keepdims, ddof=ddof,
                   split_every=split_every, out=out)

    def moment(self, order, axis=None, dtype=None, keepdims=False, ddof=0,
               split_every=None, out=None):
        """Calculate the nth centralized moment.

        Parameters
        ----------
        order : int
            Order of the moment that is returned, must be >= 2.
        axis : int, optional
            Axis along which the central moment is computed. The default is to
            compute the moment of the flattened array.
        dtype : data-type, optional
            Type to use in computing the moment. For arrays of integer type the
            default is float64; for arrays of float types it is the same as the
            array type.
        keepdims : bool, optional
            If this is set to True, the axes which are reduced are left in the
            result as dimensions with size one. With this option, the result
            will broadcast correctly against the original array.
        ddof : int, optional
            "Delta Degrees of Freedom": the divisor used in the calculation is
            N - ddof, where N represents the number of elements. By default
            ddof is zero.

        Returns
        -------
        moment : ndarray

        References
        ----------
        .. [1] Pebay, Philippe (2008), "Formulas for Robust, One-Pass Parallel
        Computation of Covariances and Arbitrary-Order Statistical Moments"
        (PDF), Technical Report SAND2008-6212, Sandia National Laboratories

        """

        from .reductions import moment
        return moment(self, order, axis=axis, dtype=dtype, keepdims=keepdims,
                      ddof=ddof, split_every=split_every, out=out)

    def vnorm(self, ord=None, axis=None, keepdims=False, split_every=None,
              out=None):
        """ Vector norm """
        from .reductions import vnorm
        return vnorm(self, ord=ord, axis=axis, keepdims=keepdims,
                     split_every=split_every, out=out)

    @wraps(map_blocks)
    def map_blocks(self, func, *args, **kwargs):
        return map_blocks(func, self, *args, **kwargs)

    def map_overlap(self, func, depth, boundary=None, trim=True, **kwargs):
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
        >>> x = np.array([1, 1, 2, 3, 3, 3, 2, 1, 1])
        >>> x = from_array(x, chunks=5)
        >>> def derivative(x):
        ...     return x - np.roll(x, 1)

        >>> y = x.map_overlap(derivative, depth=1, boundary=0)
        >>> y.compute()
        array([ 1,  0,  1,  1,  0,  0, -1, -1,  0])

        >>> import dask.array as da
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
        from .ghost import map_overlap
        return map_overlap(self, func, depth, boundary, trim, **kwargs)

    def cumsum(self, axis, dtype=None, out=None):
        """ See da.cumsum for docstring """
        from .reductions import cumsum
        return cumsum(self, axis, dtype, out=out)

    def cumprod(self, axis, dtype=None, out=None):
        """ See da.cumprod for docstring """
        from .reductions import cumprod
        return cumprod(self, axis, dtype, out=out)

    @derived_from(np.ndarray)
    def squeeze(self, axis=None):
        from .routines import squeeze
        return squeeze(self, axis)

    def rechunk(self, chunks, threshold=None, block_size_limit=None):
        """ See da.rechunk for docstring """
        from . import rechunk   # avoid circular import
        return rechunk(self, chunks, threshold, block_size_limit)

    @property
    def real(self):
        from .ufunc import real
        return real(self)

    @property
    def imag(self):
        from .ufunc import imag
        return imag(self)

    def conj(self):
        from .ufunc import conj
        return conj(self)

    @derived_from(np.ndarray)
    def clip(self, min=None, max=None):
        from .ufunc import clip
        return clip(self, min, max)

    def view(self, dtype, order='C'):
        """ Get a view of the array as a new data type

        Parameters
        ----------
        dtype:
            The dtype by which to view the array
        order: string
            'C' or 'F' (Fortran) ordering

        This reinterprets the bytes of the array under a new dtype.  If that
        dtype does not have the same size as the original array then the shape
        will change.

        Beware that both numpy and dask.array can behave oddly when taking
        shape-changing views of arrays under Fortran ordering.  Under some
        versions of NumPy this function will fail when taking shape-changing
        views of Fortran ordered arrays if the first dimension has chunks of
        size one.
        """
        dtype = np.dtype(dtype)
        mult = self.dtype.itemsize / dtype.itemsize

        if order == 'C':
            chunks = self.chunks[:-1] + (tuple(ensure_int(c * mult)
                                         for c in self.chunks[-1]),)
        elif order == 'F':
            chunks = ((tuple(ensure_int(c * mult) for c in self.chunks[0]), ) +
                      self.chunks[1:])
        else:
            raise ValueError("Order must be one of 'C' or 'F'")

        return self.map_blocks(chunk.view, dtype, order=order,
                               dtype=dtype, chunks=chunks)

    @derived_from(np.ndarray)
    def swapaxes(self, axis1, axis2):
        from .routines import swapaxes
        return swapaxes(self, axis1, axis2)

    @derived_from(np.ndarray)
    def round(self, decimals=0):
        from .routines import round
        return round(self, decimals=decimals)

    def copy(self):
        """
        Copy array.  This is a no-op for dask.arrays, which are immutable
        """
        return Array(self.dask, self.name, self.chunks, self.dtype)

    def __deepcopy__(self, memo):
        c = self.copy()
        memo[id(self)] = c
        return c

    def to_delayed(self, optimize_graph=True):
        """Convert into an array of ``dask.delayed`` objects, one per chunk.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.

        See Also
        --------
        dask.array.from_delayed
        """
        from ..delayed import Delayed
        keys = self.__dask_keys__()
        dsk = self.__dask_graph__()
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, keys)
        L = ndeepmap(self.ndim, lambda k: Delayed(k, dsk), keys)
        return np.array(L, dtype=object)

    @derived_from(np.ndarray)
    def repeat(self, repeats, axis=None):
        from .creation import repeat
        return repeat(self, repeats, axis=axis)

    @derived_from(np.ndarray)
    def nonzero(self):
        from .routines import nonzero
        return nonzero(self)


def ensure_int(f):
    i = int(f)
    if i != f:
        raise ValueError("Could not coerce %f to integer" % f)
    return i


def normalize_chunks(chunks, shape=None):
    """ Normalize chunks to tuple of tuples

    >>> normalize_chunks((2, 2), shape=(5, 6))
    ((2, 2, 1), (2, 2, 2))

    >>> normalize_chunks(((2, 2, 1), (2, 2, 2)), shape=(5, 6))  # Idempotent
    ((2, 2, 1), (2, 2, 2))

    >>> normalize_chunks([[2, 2], [3, 3]])  # Cleans up lists to tuples
    ((2, 2), (3, 3))

    >>> normalize_chunks(10, shape=(30, 5))  # Supports integer inputs
    ((10, 10, 10), (5,))

    >>> normalize_chunks((-1,), shape=(10,))  # -1 gets mapped to full size
    ((10,),)

    >>> normalize_chunks((), shape=(0, 0))  #  respects null dimensions
    ((0,), (0,))
    """
    if chunks is None:
        raise ValueError(CHUNKS_NONE_ERROR_MESSAGE)
    if type(chunks) is not tuple:
        if type(chunks) is list:
            chunks = tuple(chunks)
        if isinstance(chunks, Number):
            chunks = (chunks,) * len(shape)
    if not chunks and shape and all(s == 0 for s in shape):
        chunks = ((0,),) * len(shape)

    if shape and len(chunks) != len(shape):
        if not (len(shape) == 1 and sum(chunks) == shape[0]):
            raise ValueError(
                "Chunks and shape must be of the same length/dimension. "
                "Got chunks=%s, shape=%s" % (chunks, shape))

    if shape is not None:
        chunks = tuple(c if c not in {None, -1} else s
                       for c, s in zip(chunks, shape))

    if chunks and shape is not None:
        chunks = sum((blockdims_from_blockshape((s,), (c,))
                      if not isinstance(c, (tuple, list)) else (c,)
                      for s, c in zip(shape, chunks)), ())
    for c in chunks:
        if not c:
            raise ValueError("Empty tuples are not allowed in chunks. Express "
                             "zero length dimensions with 0(s) in chunks")

    if shape is not None:
        if len(chunks) != len(shape):
            raise ValueError("Input array has %d dimensions but the supplied "
                             "chunks has only %d dimensions" %
                             (len(shape), len(chunks)))
        if not all(c == s or (math.isnan(c) and math.isnan(s))
                   for c, s in zip(map(sum, chunks), shape)):
            raise ValueError("Chunks do not add up to shape. "
                             "Got chunks=%s, shape=%s" % (chunks, shape))

    return tuple(tuple(int(x) if not math.isnan(x) else x for x in c) for c in chunks)


def from_array(x, chunks, name=None, lock=False, asarray=True, fancy=True,
               getitem=None):
    """ Create dask array from something that looks like an array

    Input must have a ``.shape`` and support numpy-style slicing.

    Parameters
    ----------
    x : array_like
    chunks : int, tuple
        How to chunk the array. Must be one of the following forms:
        -   A blocksize like 1000.
        -   A blockshape like (1000, 1000).
        -   Explicit sizes of all blocks along all dimensions like
            ((1000, 1000, 500), (400, 400)).

        -1 as a blocksize indicates the size of the corresponding dimension.
    name : str, optional
        The key name to use for the array. Defaults to a hash of ``x``.
        Use ``name=False`` to generate a random name instead of hashing (fast)
    lock : bool or Lock, optional
        If ``x`` doesn't support concurrent reads then provide a lock here, or
        pass in True to have dask.array create one for you.
    asarray : bool, optional
        If True (default), then chunks will be converted to instances of
        ``ndarray``. Set to False to pass passed chunks through unchanged.
    fancy : bool, optional
        If ``x`` doesn't support fancy indexing (e.g. indexing with lists or
        arrays) then set to False. Default is True.

    Examples
    --------

    >>> x = h5py.File('...')['/data/path']  # doctest: +SKIP
    >>> a = da.from_array(x, chunks=(1000, 1000))  # doctest: +SKIP

    If your underlying datastore does not support concurrent reads then include
    the ``lock=True`` keyword argument or ``lock=mylock`` if you want multiple
    arrays to coordinate around the same lock.

    >>> a = da.from_array(x, chunks=(1000, 1000), lock=True)  # doctest: +SKIP
    """
    chunks = normalize_chunks(chunks, x.shape)
    if name in (None, True):
        token = tokenize(x, chunks)
        original_name = 'array-original-' + token
        name = name or 'array-' + token
    elif name is False:
        original_name = name = 'array-' + str(uuid.uuid1())
    else:
        original_name = name
    if lock is True:
        lock = SerializableLock()

    if getitem is None:
        getitem = getter if fancy else getter_nofancy

    dsk = getem(original_name, chunks, getitem=getitem, shape=x.shape,
                out_name=name, lock=lock, asarray=asarray)
    dsk[original_name] = x

    return Array(dsk, name, chunks, dtype=x.dtype)


def from_delayed(value, shape, dtype, name=None):
    """ Create a dask array from a dask delayed value

    This routine is useful for constructing dask arrays in an ad-hoc fashion
    using dask delayed, particularly when combined with stack and concatenate.

    The dask array will consist of a single chunk.

    Examples
    --------
    >>> from dask import delayed
    >>> value = delayed(np.ones)(5)
    >>> array = from_delayed(value, (5,), float)
    >>> array
    dask.array<from-value, shape=(5,), dtype=float64, chunksize=(5,)>
    >>> array.compute()
    array([1., 1., 1., 1., 1.])
    """
    from dask.delayed import delayed, Delayed
    if not isinstance(value, Delayed) and hasattr(value, 'key'):
        value = delayed(value)
    name = name or 'from-value-' + tokenize(value, shape, dtype)
    dsk = {(name,) + (0,) * len(shape): value.key}
    chunks = tuple((d,) for d in shape)
    return Array(sharedict.merge(value.dask, (name, dsk)), name, chunks, dtype)


def from_func(func, shape, dtype=None, name=None, args=(), kwargs={}):
    """ Create dask array in a single block by calling a function

    Calling the provided function with func(*args, **kwargs) should return a
    NumPy array of the indicated shape and dtype.

    Examples
    --------

    >>> a = from_func(np.arange, (3,), dtype='i8', args=(3,))
    >>> a.compute()
    array([0, 1, 2])

    This works particularly well when coupled with dask.array functions like
    concatenate and stack:

    >>> arrays = [from_func(np.array, (), dtype='i8', args=(n,)) for n in range(5)]
    >>> stack(arrays).compute()
    array([0, 1, 2, 3, 4])
    """
    name = name or 'from_func-' + tokenize(func, shape, dtype, args, kwargs)
    if args or kwargs:
        func = partial(func, *args, **kwargs)
    dsk = {(name,) + (0,) * len(shape): (func,)}
    chunks = tuple((i,) for i in shape)
    return Array(dsk, name, chunks, dtype)


def common_blockdim(blockdims):
    """ Find the common block dimensions from the list of block dimensions

    Currently only implements the simplest possible heuristic: the common
    block-dimension is the only one that does not span fully span a dimension.
    This is a conservative choice that allows us to avoid potentially very
    expensive rechunking.

    Assumes that each element of the input block dimensions has all the same
    sum (i.e., that they correspond to dimensions of the same size).

    Examples
    --------
    >>> common_blockdim([(3,), (2, 1)])
    (2, 1)
    >>> common_blockdim([(1, 2), (2, 1)])
    (1, 1, 1)
    >>> common_blockdim([(2, 2), (3, 1)])  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    ValueError: Chunks do not align
    """
    if not any(blockdims):
        return ()
    non_trivial_dims = set([d for d in blockdims if len(d) > 1])
    if len(non_trivial_dims) == 1:
        return first(non_trivial_dims)
    if len(non_trivial_dims) == 0:
        return max(blockdims, key=first)

    if np.isnan(sum(map(sum, blockdims))):
        raise ValueError("Arrays chunk sizes are unknown: %s", blockdims)

    if len(set(map(sum, non_trivial_dims))) > 1:
        raise ValueError("Chunks do not add up to same value", blockdims)

    # We have multiple non-trivial chunks on this axis
    # e.g. (5, 2) and (4, 3)

    # We create a single chunk tuple with the same total length
    # that evenly divides both, e.g. (4, 1, 2)

    # To accomplish this we walk down all chunk tuples together, finding the
    # smallest element, adding it to the output, and subtracting it from all
    # other elements and remove the element itself.  We stop once we have
    # burned through all of the chunk tuples.
    # For efficiency's sake we reverse the lists so that we can pop off the end
    rchunks = [list(ntd)[::-1] for ntd in non_trivial_dims]
    total = sum(first(non_trivial_dims))
    i = 0

    out = []
    while i < total:
        m = min(c[-1] for c in rchunks)
        out.append(m)
        for c in rchunks:
            c[-1] -= m
            if c[-1] == 0:
                c.pop()
        i += m

    return tuple(out)


def unify_chunks(*args, **kwargs):
    """
    Unify chunks across a sequence of arrays

    Parameters
    ----------
    *args: sequence of Array, index pairs
        Sequence like (x, 'ij', y, 'jk', z, 'i')

    Examples
    --------
    >>> import dask.array as da
    >>> x = da.ones(10, chunks=((5, 2, 3),))
    >>> y = da.ones(10, chunks=((2, 3, 5),))
    >>> chunkss, arrays = unify_chunks(x, 'i', y, 'i')
    >>> chunkss
    {'i': (2, 3, 2, 3)}

    >>> x = da.ones((100, 10), chunks=(20, 5))
    >>> y = da.ones((10, 100), chunks=(4, 50))
    >>> chunkss, arrays = unify_chunks(x, 'ij', y, 'jk')
    >>> chunkss  # doctest: +SKIP
    {'k': (50, 50), 'i': (20, 20, 20, 20, 20), 'j': (4, 1, 3, 2)}

    Returns
    -------
    chunkss : dict
        Map like {index: chunks}.
    arrays : list
        List of rechunked arrays.

    See Also
    --------
    common_blockdim
    """
    arginds = [(asarray(a) if ind is not None else a, ind)
               for a, ind in partition(2, args)]  # [x, ij, y, jk]
    args = list(concat(arginds))  # [(x, ij), (y, jk)]
    warn = kwargs.get('warn', True)

    arrays, inds = zip(*arginds)
    if all(ind == inds[0] for ind in inds) and all(a.chunks == arrays[0].chunks for a in arrays):
        return dict(zip(inds[0], arrays[0].chunks)), arrays

    nameinds = [(a.name if i is not None else a, i) for a, i in arginds]
    blockdim_dict = {a.name: a.chunks
                     for a, ind in arginds
                     if ind is not None}

    chunkss = broadcast_dimensions(nameinds, blockdim_dict,
                                   consolidate=common_blockdim)
    max_parts = max(arg.npartitions for arg, ind in arginds if ind is not None)
    nparts = np.prod(list(map(len, chunkss.values())))

    if warn and nparts and nparts >= max_parts * 10:
        warnings.warn("Increasing number of chunks by factor of %d" %
                      (nparts / max_parts))

    arrays = []
    for a, i in arginds:
        if i is None:
            arrays.append(a)
        else:
            chunks = tuple(chunkss[j] if a.shape[n] > 1 else a.shape[n]
                           if not np.isnan(sum(chunkss[j])) else None
                           for n, j in enumerate(i))
            if chunks != a.chunks and all(a.chunks):
                arrays.append(a.rechunk(chunks))
            else:
                arrays.append(a)
    return chunkss, arrays


def atop(func, out_ind, *args, **kwargs):
    """ Tensor operation: Generalized inner and outer products

    A broad class of blocked algorithms and patterns can be specified with a
    concise multi-index notation.  The ``atop`` function applies an in-memory
    function across multiple blocks of multiple inputs in a variety of ways.
    Many dask.array operations are special cases of atop including elementwise,
    broadcasting, reductions, tensordot, and transpose.

    Parameters
    ----------
    func : callable
        Function to apply to individual tuples of blocks
    out_ind : iterable
        Block pattern of the output, something like 'ijk' or (1, 2, 3)
    *args : sequence of Array, index pairs
        Sequence like (x, 'ij', y, 'jk', z, 'i')
    **kwargs : dict
        Extra keyword arguments to pass to function
    dtype : np.dtype
        Datatype of resulting array.
    concatenate : bool, keyword only
        If true concatenate arrays along dummy indices, else provide lists
    adjust_chunks : dict
        Dictionary mapping index to function to be applied to chunk sizes
    new_axes : dict, keyword only
        New indexes and their dimension lengths

    Examples
    --------
    2D embarrassingly parallel operation from two arrays, x, and y.

    >>> z = atop(operator.add, 'ij', x, 'ij', y, 'ij', dtype='f8')  # z = x + y  # doctest: +SKIP

    Outer product multiplying x by y, two 1-d vectors

    >>> z = atop(operator.mul, 'ij', x, 'i', y, 'j', dtype='f8')  # doctest: +SKIP

    z = x.T

    >>> z = atop(np.transpose, 'ji', x, 'ij', dtype=x.dtype)  # doctest: +SKIP

    The transpose case above is illustrative because it does same transposition
    both on each in-memory block by calling ``np.transpose`` and on the order
    of the blocks themselves, by switching the order of the index ``ij -> ji``.

    We can compose these same patterns with more variables and more complex
    in-memory functions

    z = X + Y.T

    >>> z = atop(lambda x, y: x + y.T, 'ij', x, 'ij', y, 'ji', dtype='f8')  # doctest: +SKIP

    Any index, like ``i`` missing from the output index is interpreted as a
    contraction (note that this differs from Einstein convention; repeated
    indices do not imply contraction.)  In the case of a contraction the passed
    function should expect an iterable of blocks on any array that holds that
    index.  To receive arrays concatenated along contracted dimensions instead
    pass ``concatenate=True``.

    Inner product multiplying x by y, two 1-d vectors

    >>> def sequence_dot(x_blocks, y_blocks):
    ...     result = 0
    ...     for x, y in zip(x_blocks, y_blocks):
    ...         result += x.dot(y)
    ...     return result

    >>> z = atop(sequence_dot, '', x, 'i', y, 'i', dtype='f8')  # doctest: +SKIP

    Add new single-chunk dimensions with the ``new_axes=`` keyword, including
    the length of the new dimension.  New dimensions will always be in a single
    chunk.

    >>> def f(x):
    ...     return x[:, None] * np.ones((1, 5))

    >>> z = atop(f, 'az', x, 'a', new_axes={'z': 5}, dtype=x.dtype)  # doctest: +SKIP

    If the applied function changes the size of each chunk you can specify this
    with a ``adjust_chunks={...}`` dictionary holding a function for each index
    that modifies the dimension size in that index.

    >>> def double(x):
    ...     return np.concatenate([x, x])

    >>> y = atop(double, 'ij', x, 'ij',
    ...          adjust_chunks={'i': lambda n: 2 * n}, dtype=x.dtype)  # doctest: +SKIP

    Include literals by indexing with None

    >>> y = atop(add, 'ij', x, 'ij', 1234, None, dtype=x.dtype)  # doctest: +SKIP

    See Also
    --------
    top - dict formulation of this function, contains most logic
    """
    out = kwargs.pop('name', None)      # May be None at this point
    token = kwargs.pop('token', None)
    dtype = kwargs.pop('dtype', None)
    adjust_chunks = kwargs.pop('adjust_chunks', None)
    new_axes = kwargs.get('new_axes', {})

    if dtype is None:
        raise ValueError("Must specify dtype of output array")

    chunkss, arrays = unify_chunks(*args)
    for k, v in new_axes.items():
        chunkss[k] = (v,)
    arginds = list(zip(arrays, args[1::2]))

    for arg, ind in arginds:
        if hasattr(arg, 'ndim') and hasattr(ind, '__len__') and arg.ndim != len(ind):
            raise ValueError("Index string %s does not match array dimension %d"
                             % (ind, arg.ndim))

    numblocks = {a.name: a.numblocks for a, ind in arginds if ind is not None}
    argindsstr = list(concat([(a if ind is None else a.name, ind) for a, ind in arginds]))
    # Finish up the name
    if not out:
        out = '%s-%s' % (token or funcname(func).strip('_'),
                         tokenize(func, out_ind, argindsstr, dtype, **kwargs))

    dsk = top(func, out, out_ind, *argindsstr, numblocks=numblocks, **kwargs)
    dsks = [a.dask for a, ind in arginds if ind is not None]

    chunks = [chunkss[i] for i in out_ind]
    if adjust_chunks:
        for i, ind in enumerate(out_ind):
            if ind in adjust_chunks:
                if callable(adjust_chunks[ind]):
                    chunks[i] = tuple(map(adjust_chunks[ind], chunks[i]))
                elif isinstance(adjust_chunks[ind], int):
                    chunks[i] = tuple(adjust_chunks[ind] for _ in chunks[i])
                elif isinstance(adjust_chunks[ind], (tuple, list)):
                    chunks[i] = tuple(adjust_chunks[ind])
                else:
                    raise NotImplementedError(
                        "adjust_chunks values must be callable, int, or tuple")
    chunks = tuple(chunks)

    return Array(sharedict.merge((out, dsk), *dsks), out, chunks, dtype=dtype)


def unpack_singleton(x):
    """

    >>> unpack_singleton([[[[1]]]])
    1
    >>> unpack_singleton(np.array(np.datetime64('2000-01-01')))
    array('2000-01-01', dtype='datetime64[D]')
    """
    while isinstance(x, (list, tuple)):
        try:
            x = x[0]
        except (IndexError, TypeError, KeyError):
            break
    return x


def block(arrays, allow_unknown_chunksizes=False):
    """
    Assemble an nd-array from nested lists of blocks.

    Blocks in the innermost lists are concatenated along the last
    dimension (-1), then these are concatenated along the second-last
    dimension (-2), and so on until the outermost list is reached

    Blocks can be of any dimension, but will not be broadcasted using the normal
    rules. Instead, leading axes of size 1 are inserted, to make ``block.ndim``
    the same for all blocks. This is primarily useful for working with scalars,
    and means that code like ``block([v, 1])`` is valid, where
    ``v.ndim == 1``.

    When the nested list is two levels deep, this allows block matrices to be
    constructed from their components.

    Parameters
    ----------
    arrays : nested list of array_like or scalars (but not tuples)
        If passed a single ndarray or scalar (a nested list of depth 0), this
        is returned unmodified (and not copied).

        Elements shapes must match along the appropriate axes (without
        broadcasting), but leading 1s will be prepended to the shape as
        necessary to make the dimensions match.

    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.

    Returns
    -------
    block_array : ndarray
        The array assembled from the given blocks.

        The dimensionality of the output is equal to the greatest of:
        * the dimensionality of all the inputs
        * the depth to which the input list is nested

    Raises
    ------
    ValueError
        * If list depths are mismatched - for instance, ``[[a, b], c]`` is
          illegal, and should be spelt ``[[a, b], [c]]``
        * If lists are empty - for instance, ``[[a, b], []]``

    See Also
    --------
    concatenate : Join a sequence of arrays together.
    stack : Stack arrays in sequence along a new dimension.
    hstack : Stack arrays in sequence horizontally (column wise).
    vstack : Stack arrays in sequence vertically (row wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).
    vsplit : Split array into a list of multiple sub-arrays vertically.

    Notes
    -----

    When called with only scalars, ``block`` is equivalent to an ndarray
    call. So ``block([[1, 2], [3, 4]])`` is equivalent to
    ``array([[1, 2], [3, 4]])``.

    This function does not enforce that the blocks lie on a fixed grid.
    ``block([[a, b], [c, d]])`` is not restricted to arrays of the form::

        AAAbb
        AAAbb
        cccDD

    But is also allowed to produce, for some ``a, b, c, d``::

        AAAbb
        AAAbb
        cDDDD

    Since concatenation happens along the last axis first, `block` is _not_
    capable of producing the following directly::

        AAAbb
        cccbb
        cccDD

    Matlab's "square bracket stacking", ``[A, B, ...; p, q, ...]``, is
    equivalent to ``block([[A, B, ...], [p, q, ...]])``.
    """

    # This was copied almost verbatim from numpy.core.shape_base.block
    # See numpy license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory

    def atleast_nd(x, ndim):
        x = asanyarray(x)
        diff = max(ndim - x.ndim, 0)
        return x[(None,) * diff + (Ellipsis,)]

    def format_index(index):
        return 'arrays' + ''.join('[{}]'.format(i) for i in index)

    rec = _Recurser(recurse_if=lambda x: type(x) is list)

    # ensure that the lists are all matched in depth
    list_ndim = None
    any_empty = False
    for index, value, entering in rec.walk(arrays):
        if type(value) is tuple:
            # not strictly necessary, but saves us from:
            #  - more than one way to do things - no point treating tuples like
            #    lists
            #  - horribly confusing behaviour that results when tuples are
            #    treated like ndarray
            raise TypeError(
                '{} is a tuple. '
                'Only lists can be used to arrange blocks, and np.block does '
                'not allow implicit conversion from tuple to ndarray.'.format(
                    format_index(index)
                )
            )
        if not entering:
            curr_depth = len(index)
        elif len(value) == 0:
            curr_depth = len(index) + 1
            any_empty = True
        else:
            continue

        if list_ndim is not None and list_ndim != curr_depth:
            raise ValueError(
                "List depths are mismatched. First element was at depth {}, "
                "but there is an element at depth {} ({})".format(
                    list_ndim,
                    curr_depth,
                    format_index(index)
                )
            )
        list_ndim = curr_depth

    # do this here so we catch depth mismatches first
    if any_empty:
        raise ValueError('Lists cannot be empty')

    # convert all the arrays to ndarrays
    arrays = rec.map_reduce(
        arrays,
        f_map=asanyarray,
        f_reduce=list
    )

    # determine the maximum dimension of the elements
    elem_ndim = rec.map_reduce(
        arrays,
        f_map=lambda xi: xi.ndim,
        f_reduce=max
    )
    ndim = max(list_ndim, elem_ndim)

    # first axis to concatenate along
    first_axis = ndim - list_ndim

    # Make all the elements the same dimension
    arrays = rec.map_reduce(
        arrays,
        f_map=lambda xi: atleast_nd(xi, ndim),
        f_reduce=list
    )

    # concatenate innermost lists on the right, outermost on the left
    return rec.map_reduce(
        arrays,
        f_reduce=lambda xs, axis: concatenate(
            list(xs),
            axis=axis,
            allow_unknown_chunksizes=allow_unknown_chunksizes
        ),
        f_kwargs=lambda axis: dict(axis=(axis + 1)),
        axis=first_axis
    )


def concatenate(seq, axis=0, allow_unknown_chunksizes=False):
    """
    Concatenate arrays along an existing axis

    Given a sequence of dask Arrays form a new dask Array by stacking them
    along an existing dimension (axis=0 by default)

    Parameters
    ----------
    seq: list of dask.arrays
    axis: int
        Dimension along which to align all of the arrays
    allow_unknown_chunksizes: bool
        Allow unknown chunksizes, such as come from converting from dask
        dataframes.  Dask.array is unable to verify that chunks line up.  If
        data comes from differently aligned sources then this can cause
        unexpected results.

    Examples
    --------

    Create slices

    >>> import dask.array as da
    >>> import numpy as np

    >>> data = [from_array(np.ones((4, 4)), chunks=(2, 2))
    ...          for i in range(3)]

    >>> x = da.concatenate(data, axis=0)
    >>> x.shape
    (12, 4)

    >>> da.concatenate(data, axis=1).shape
    (4, 12)

    Result is a new dask Array

    See Also
    --------
    stack
    """
    n = len(seq)
    ndim = len(seq[0].shape)

    if axis < 0:
        axis = ndim + axis
    if axis >= ndim:
        msg = ("Axis must be less than than number of dimensions"
               "\nData has %d dimensions, but got axis=%d")
        raise ValueError(msg % (ndim, axis))

    if n == 1:
        return seq[0]

    if (not allow_unknown_chunksizes and
        not all(i == axis or all(x.shape[i] == seq[0].shape[i] for x in seq)
                for i in range(ndim))):
        if any(map(np.isnan, seq[0].shape)):
            raise ValueError("Tried to concatenate arrays with unknown"
                             " shape %s.  To force concatenation pass"
                             " allow_unknown_chunksizes=True."
                             % str(seq[0].shape))
        raise ValueError("Shapes do not align: %s", [x.shape for x in seq])

    inds = [list(range(ndim)) for i in range(n)]
    for i, ind in enumerate(inds):
        ind[axis] = -(i + 1)

    uc_args = list(concat(zip(seq, inds)))
    _, seq = unify_chunks(*uc_args, warn=False)

    bds = [a.chunks for a in seq]

    chunks = (seq[0].chunks[:axis] + (sum([bd[axis] for bd in bds], ()), ) +
              seq[0].chunks[axis + 1:])

    cum_dims = [0] + list(accumulate(add, [len(a.chunks[axis]) for a in seq]))

    seq_dtypes = [a.dtype for a in seq]
    if len(set(seq_dtypes)) > 1:
        dt = reduce(np.promote_types, seq_dtypes)
        seq = [x.astype(dt) for x in seq]
    else:
        dt = seq_dtypes[0]

    names = [a.name for a in seq]

    name = 'concatenate-' + tokenize(names, axis)
    keys = list(product([name], *[range(len(bd)) for bd in chunks]))

    values = [(names[bisect(cum_dims, key[axis + 1]) - 1],) + key[1:axis + 1] +
              (key[axis + 1] - cum_dims[bisect(cum_dims, key[axis + 1]) - 1], ) +
              key[axis + 2:] for key in keys]

    dsk = dict(zip(keys, values))
    dsk2 = sharedict.merge((name, dsk), * [a.dask for a in seq])

    return Array(dsk2, name, chunks, dtype=dt)


def load_store_chunk(x, out, index, lock, return_stored, load_stored):
    """
    A function inserted in a Dask graph for storing a chunk.

    Parameters
    ----------
    x: array-like
        An array (potentially a NumPy one)
    out: array-like
        Where to store results too.
    index: slice-like
        Where to store result from ``x`` in ``out``.
    lock: Lock-like or False
        Lock to use before writing to ``out``.
    return_stored: bool
        Whether to return ``out``.
    load_stored: bool
        Whether to return the array stored in ``out``.
        Ignored if ``return_stored`` is not ``True``.

    Examples
    --------

    >>> a = np.ones((5, 6))
    >>> b = np.empty(a.shape)
    >>> load_store_chunk(a, b, (slice(None), slice(None)), False, False, False)
    """

    result = None
    if return_stored and not load_stored:
        result = out

    if lock:
        lock.acquire()
    try:
        if x is not None:
            out[index] = np.asanyarray(x)
        if return_stored and load_stored:
            result = out[index]
    finally:
        if lock:
            lock.release()

    return result


def store_chunk(x, out, index, lock, return_stored):
    return load_store_chunk(x, out, index, lock, return_stored, False)


def load_chunk(out, index, lock):
    return load_store_chunk(None, out, index, lock, True, True)


def insert_to_ooc(arr, out, lock=True, region=None,
                  return_stored=False, load_stored=False):
    """
    Creates a Dask graph for storing chunks from ``arr`` in ``out``.

    Parameters
    ----------
    arr: da.Array
        A dask array
    out: array-like
        Where to store results too.
    lock: Lock-like or bool, optional
        Whether to lock or with what (default is ``True``,
        which means a ``threading.Lock`` instance).
    region: slice-like, optional
        Where in ``out`` to store ``arr``'s results
        (default is ``None``, meaning all of ``out``).
    return_stored: bool, optional
        Whether to return ``out``
        (default is ``False``, meaning ``None`` is returned).
    load_stored: bool, optional
        Whether to handling loading from ``out`` at the same time.
        Ignored if ``return_stored`` is not ``True``.
        (default is ``False``, meaning defer to ``return_stored``).

    Examples
    --------
    >>> import dask.array as da
    >>> d = da.ones((5, 6), chunks=(2, 3))
    >>> a = np.empty(d.shape)
    >>> insert_to_ooc(d, a)  # doctest: +SKIP
    """

    if lock is True:
        lock = Lock()

    slices = slices_from_chunks(arr.chunks)
    if region:
        slices = [fuse_slice(region, slc) for slc in slices]

    name = 'store-%s' % arr.name
    func = store_chunk
    args = ()
    if return_stored and load_stored:
        name = 'load-%s' % name
        func = load_store_chunk
        args = args + (load_stored,)

    dsk = {
        (name,) + t[1:]: (func, t, out, slc, lock, return_stored) + args
        for t, slc in zip(core.flatten(arr.__dask_keys__()), slices)
    }

    return dsk


def retrieve_from_ooc(keys, dsk_pre, dsk_post=None):
    """
    Creates a Dask graph for loading stored ``keys`` from ``dsk``.

    Parameters
    ----------
    keys: Sequence
        A sequence containing Dask graph keys to load
    dsk_pre: Mapping
        A Dask graph corresponding to a Dask Array before computation
    dsk_post: Mapping, optional
        A Dask graph corresponding to a Dask Array after computation

    Examples
    --------
    >>> import dask.array as da
    >>> d = da.ones((5, 6), chunks=(2, 3))
    >>> a = np.empty(d.shape)
    >>> g = insert_to_ooc(d, a)
    >>> retrieve_from_ooc(g.keys(), g)  # doctest: +SKIP
    """

    if not dsk_post:
        dsk_post = {k: k for k in keys}

    load_dsk = {
        ('load-' + k[0],) + k[1:]: (load_chunk, dsk_post[k]) + dsk_pre[k][3:-1]
        for k in keys
    }

    return load_dsk


def asarray(a):
    """Convert the input to a dask array.

    Parameters
    ----------
    a : array-like
        Input data, in any form that can be converted to a dask array.

    Returns
    -------
    out : dask array
        Dask array interpretation of a.

    Examples
    --------
    >>> import dask.array as da
    >>> import numpy as np
    >>> x = np.arange(3)
    >>> da.asarray(x)
    dask.array<array, shape=(3,), dtype=int64, chunksize=(3,)>

    >>> y = [[1, 2, 3], [4, 5, 6]]
    >>> da.asarray(y)
    dask.array<array, shape=(2, 3), dtype=int64, chunksize=(2, 3)>
    """
    if isinstance(a, Array):
        return a
    if isinstance(a, (list, tuple)) and any(isinstance(i, Array) for i in a):
        a = stack(a)
    elif not isinstance(getattr(a, 'shape', None), Iterable):
        a = np.asarray(a)
    return from_array(a, chunks=a.shape, getitem=getter_inline)


def asanyarray(a):
    """Convert the input to a dask array.

    Subclasses of ``np.ndarray`` will be passed through as chunks unchanged.

    Parameters
    ----------
    a : array-like
        Input data, in any form that can be converted to a dask array.

    Returns
    -------
    out : dask array
        Dask array interpretation of a.

    Examples
    --------
    >>> import dask.array as da
    >>> import numpy as np
    >>> x = np.arange(3)
    >>> da.asanyarray(x)
    dask.array<array, shape=(3,), dtype=int64, chunksize=(3,)>

    >>> y = [[1, 2, 3], [4, 5, 6]]
    >>> da.asanyarray(y)
    dask.array<array, shape=(2, 3), dtype=int64, chunksize=(2, 3)>
    """
    if isinstance(a, Array):
        return a
    if isinstance(a, (list, tuple)) and any(isinstance(i, Array) for i in a):
        a = stack(a)
    elif not isinstance(getattr(a, 'shape', None), Iterable):
        a = np.asanyarray(a)
    return from_array(a, chunks=a.shape, getitem=getter_inline,
                      asarray=False)


def is_scalar_for_elemwise(arg):
    """

    >>> is_scalar_for_elemwise(42)
    True
    >>> is_scalar_for_elemwise('foo')
    True
    >>> is_scalar_for_elemwise(True)
    True
    >>> is_scalar_for_elemwise(np.array(42))
    True
    >>> is_scalar_for_elemwise([1, 2, 3])
    True
    >>> is_scalar_for_elemwise(np.array([1, 2, 3]))
    False
    >>> is_scalar_for_elemwise(from_array(np.array(0), chunks=()))
    False
    >>> is_scalar_for_elemwise(np.dtype('i4'))
    True
    """
    return (np.isscalar(arg) or
            not isinstance(getattr(arg, 'shape', None), Iterable) or
            isinstance(arg, np.dtype) or
            (isinstance(arg, np.ndarray) and arg.ndim == 0))


def broadcast_shapes(*shapes):
    """
    Determines output shape from broadcasting arrays.

    Parameters
    ----------
    shapes : tuples
        The shapes of the arguments.

    Returns
    -------
    output_shape : tuple

    Raises
    ------
    ValueError
        If the input shapes cannot be successfully broadcast together.
    """
    if len(shapes) == 1:
        return shapes[0]
    out = []
    for sizes in zip_longest(*map(reversed, shapes), fillvalue=-1):
        dim = 0 if 0 in sizes else np.max(sizes)
        if any(i not in [-1, 0, 1, dim] and not np.isnan(i) for i in sizes):
            raise ValueError("operands could not be broadcast together with "
                             "shapes {0}".format(' '.join(map(str, shapes))))
        out.append(dim)
    return tuple(reversed(out))


def elemwise(op, *args, **kwargs):
    """ Apply elementwise function across arguments

    Respects broadcasting rules

    Examples
    --------
    >>> elemwise(add, x, y)  # doctest: +SKIP
    >>> elemwise(sin, x)  # doctest: +SKIP

    See Also
    --------
    atop
    """
    out = kwargs.pop('out', None)
    if not set(['name', 'dtype']).issuperset(kwargs):
        msg = "%s does not take the following keyword arguments %s"
        raise TypeError(msg % (op.__name__, str(sorted(set(kwargs) - set(['name', 'dtype'])))))

    args = [np.asarray(a) if isinstance(a, (list, tuple)) else a for a in args]

    shapes = [getattr(arg, 'shape', ()) for arg in args]
    shapes = [s if isinstance(s, Iterable) else () for s in shapes]
    out_ndim = len(broadcast_shapes(*shapes))   # Raises ValueError if dimensions mismatch
    expr_inds = tuple(range(out_ndim))[::-1]

    need_enforce_dtype = False
    if 'dtype' in kwargs:
        dt = kwargs['dtype']
    else:
        # We follow NumPy's rules for dtype promotion, which special cases
        # scalars and 0d ndarrays (which it considers equivalent) by using
        # their values to compute the result dtype:
        # https://github.com/numpy/numpy/issues/6240
        # We don't inspect the values of 0d dask arrays, because these could
        # hold potentially very expensive calculations. Instead, we treat
        # them just like other arrays, and if necessary cast the result of op
        # to match.
        vals = [np.empty((1,) * max(1, a.ndim), dtype=a.dtype)
                if not is_scalar_for_elemwise(a) else a
                for a in args]
        try:
            dt = apply_infer_dtype(op, vals, {}, 'elemwise', suggest_dtype=False)
        except Exception:
            return NotImplemented
        need_enforce_dtype = any(not is_scalar_for_elemwise(a) and a.ndim == 0 for a in args)

    name = kwargs.get('name', None) or '%s-%s' % (funcname(op),
                                                  tokenize(op, dt, *args))

    atop_kwargs = dict(dtype=dt, name=name, token=funcname(op).strip('_'))
    if need_enforce_dtype:
        atop_kwargs['enforce_dtype'] = dt
        atop_kwargs['enforce_dtype_function'] = op
        op = _enforce_dtype
    result = atop(op, expr_inds,
                  *concat((a, tuple(range(a.ndim)[::-1])
                           if not is_scalar_for_elemwise(a)
                           else None) for a in args),
                  **atop_kwargs)

    return handle_out(out, result)


def handle_out(out, result):
    """ Handle out parameters

    If out is a dask.array then this overwrites the contents of that array with
    the result
    """
    if isinstance(out, tuple):
        if len(out) == 1:
            out = out[0]
        elif len(out) > 1:
            raise NotImplementedError("The out parameter is not fully supported")
        else:
            out = None
    if isinstance(out, Array):
        if out.shape != result.shape:
            raise ValueError(
                "Mismatched shapes between result and out parameter. "
                "out=%s, result=%s" % (str(out.shape), str(result.shape)))
        out._chunks = result.chunks
        out.dask = result.dask
        out.dtype = result.dtype
        out.name = result.name
    elif out is not None:
        msg = ("The out parameter is not fully supported."
               " Received type %s, expected Dask Array" % type(out).__name__)
        raise NotImplementedError(msg)
    else:
        return result


def _enforce_dtype(*args, **kwargs):
    """Calls a function and converts its result to the given dtype.

    The parameters have deliberately been given unwieldy names to avoid
    clashes with keyword arguments consumed by atop

    A dtype of `object` is treated as a special case and not enforced,
    because it is used as a dummy value in some places when the result will
    not be a block in an Array.

    Parameters
    ----------
    enforce_dtype : dtype
        Result dtype
    enforce_dtype_function : callable
        The wrapped function, which will be passed the remaining arguments
    """
    dtype = kwargs.pop('enforce_dtype')
    function = kwargs.pop('enforce_dtype_function')

    result = function(*args, **kwargs)
    if dtype != result.dtype and dtype != object:
        if not np.can_cast(result, dtype, casting='same_kind'):
            raise ValueError("Inferred dtype from function %r was %r "
                             "but got %r, which can't be cast using "
                             "casting='same_kind'" %
                             (funcname(function), str(dtype), str(result.dtype)))
        if np.isscalar(result):
            # scalar astype method doesn't take the keyword arguments, so
            # have to convert via 0-dimensional array and back.
            result = result.astype(dtype)
        else:
            try:
                result = result.astype(dtype, copy=False)
            except TypeError:
                # Missing copy kwarg
                result = result.astype(dtype)
    return result


def broadcast_to(x, shape, chunks=None):
    """Broadcast an array to a new shape.

    Parameters
    ----------
    x : array_like
        The array to broadcast.
    shape : tuple
        The shape of the desired array.
    chunks : tuple, optional
        If provided, then the result will use these chunks instead of the same
        chunks as the source array. Setting chunks explicitly as part of
        broadcast_to is more efficient than rechunking afterwards. Chunks are
        only allowed to differ from the original shape along dimensions that
        are new on the result or have size 1 the input array.

    Returns
    -------
    broadcast : dask array

    See Also
    --------
    :func:`numpy.broadcast_to`
    """
    x = asarray(x)
    shape = tuple(shape)

    if x.shape == shape and (chunks is None or chunks == x.chunks):
        return x

    ndim_new = len(shape) - x.ndim
    if ndim_new < 0 or any(new != old
                           for new, old in zip(shape[ndim_new:], x.shape)
                           if old != 1):
        raise ValueError('cannot broadcast shape %s to shape %s'
                         % (x.shape, shape))

    if chunks is None:
        chunks = (tuple((s,) for s in shape[:ndim_new]) +
                  tuple(bd if old > 1 else (new,)
                  for bd, old, new in zip(x.chunks, x.shape, shape[ndim_new:])))
    else:
        chunks = normalize_chunks(chunks, shape)
        for old_bd, new_bd in zip(x.chunks, chunks[ndim_new:]):
            if old_bd != new_bd and old_bd != (1,):
                raise ValueError('cannot broadcast chunks %s to chunks %s: '
                                 'new chunks must either be along a new '
                                 'dimension or a dimension of size 1'
                                 % (x.chunks, chunks))

    name = 'broadcast_to-' + tokenize(x, shape, chunks)
    dsk = {}

    enumerated_chunks = product(*(enumerate(bds) for bds in chunks))
    for new_index, chunk_shape in (zip(*ec) for ec in enumerated_chunks):
        old_index = tuple(0 if bd == (1,) else i
                          for bd, i in zip(x.chunks, new_index[ndim_new:]))
        old_key = (x.name,) + old_index
        new_key = (name,) + new_index
        dsk[new_key] = (chunk.broadcast_to, old_key, quote(chunk_shape))

    return Array(sharedict.merge((name, dsk), x.dask), name, chunks,
                 dtype=x.dtype)


@wraps(np.broadcast_arrays)
def broadcast_arrays(*args, **kwargs):
    subok = bool(kwargs.pop("subok", False))

    to_array = asanyarray if subok else asarray
    args = tuple(to_array(e) for e in args)

    if kwargs:
        raise TypeError("unsupported keyword argument(s) provided")

    shape = broadcast_shapes(*(e.shape for e in args))
    chunks = broadcast_chunks(*(e.chunks for e in args))

    result = [broadcast_to(e, shape=shape, chunks=chunks) for e in args]

    return result


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


def chunks_from_arrays(arrays):
    """ Chunks tuple from nested list of arrays

    >>> x = np.array([1, 2])
    >>> chunks_from_arrays([x, x])
    ((2, 2),)

    >>> x = np.array([[1, 2]])
    >>> chunks_from_arrays([[x], [x]])
    ((1, 1), (2,))

    >>> x = np.array([[1, 2]])
    >>> chunks_from_arrays([[x, x]])
    ((1,), (2, 2))

    >>> chunks_from_arrays([1, 1])
    ((1, 1),)
    """
    if not arrays:
        return ()
    result = []
    dim = 0

    def shape(x):
        try:
            return x.shape
        except AttributeError:
            return (1,)

    while isinstance(arrays, (list, tuple)):
        result.append(tuple([shape(deepfirst(a))[dim] for a in arrays]))
        arrays = arrays[0]
        dim += 1
    return tuple(result)


def deepfirst(seq):
    """ First element in a nested list

    >>> deepfirst([[[1, 2], [3, 4]], [5, 6], [7, 8]])
    1
    """
    if not isinstance(seq, (list, tuple)):
        return seq
    else:
        return deepfirst(seq[0])


def ndimlist(seq):
    if not isinstance(seq, (list, tuple)):
        return 0
    elif not seq:
        return 1
    else:
        return 1 + ndimlist(seq[0])


def shapelist(a):
    """ Get the shape of nested list """
    if type(a) is list:
        return tuple([len(a)] + list(shapelist(a[0])))
    else:
        return ()


def reshapelist(shape, seq):
    """ Reshape iterator to nested shape

    >>> reshapelist((2, 3), range(6))
    [[0, 1, 2], [3, 4, 5]]
    """
    if len(shape) == 1:
        return list(seq)
    else:
        n = int(len(seq) / shape[0])
        return [reshapelist(shape[1:], part) for part in partition(n, seq)]


def transposelist(arrays, axes, extradims=0):
    """ Permute axes of nested list

    >>> transposelist([[1,1,1],[1,1,1]], [2,1])
    [[[1, 1], [1, 1], [1, 1]]]

    >>> transposelist([[1,1,1],[1,1,1]], [2,1], extradims=1)
    [[[[1], [1]], [[1], [1]], [[1], [1]]]]
    """
    if len(axes) != ndimlist(arrays):
        raise ValueError("Length of axes should equal depth of nested arrays")
    if extradims < 0:
        raise ValueError("`newdims` should be positive")
    if len(axes) > len(set(axes)):
        raise ValueError("`axes` should be unique")

    ndim = max(axes) + 1
    shape = shapelist(arrays)
    newshape = [shape[axes.index(i)] if i in axes else 1 for i in range(ndim + extradims)]

    result = list(core.flatten(arrays))
    return reshapelist(newshape, result)


def stack(seq, axis=0):
    """
    Stack arrays along a new axis

    Given a sequence of dask Arrays form a new dask Array by stacking them
    along a new dimension (axis=0 by default)

    Examples
    --------

    Create slices

    >>> import dask.array as da
    >>> import numpy as np

    >>> data = [from_array(np.ones((4, 4)), chunks=(2, 2))
    ...          for i in range(3)]

    >>> x = da.stack(data, axis=0)
    >>> x.shape
    (3, 4, 4)

    >>> da.stack(data, axis=1).shape
    (4, 3, 4)

    >>> da.stack(data, axis=-1).shape
    (4, 4, 3)

    Result is a new dask Array

    See Also
    --------
    concatenate
    """
    n = len(seq)
    ndim = len(seq[0].shape)
    if axis < 0:
        axis = ndim + axis + 1
    if axis > ndim:
        raise ValueError("Axis must not be greater than number of dimensions"
                         "\nData has %d dimensions, but got axis=%d" %
                         (ndim, axis))
    if not all(x.shape == seq[0].shape for x in seq):
        raise ValueError("Stacked arrays must have the same shape. Got %s",
                         [x.shape for x in seq])

    ind = list(range(ndim))
    uc_args = list(concat((x, ind) for x in seq))
    _, seq = unify_chunks(*uc_args)

    dt = reduce(np.promote_types, [a.dtype for a in seq])
    seq = [x.astype(dt) for x in seq]

    assert len(set(a.chunks for a in seq)) == 1  # same chunks
    chunks = (seq[0].chunks[:axis] + ((1,) * n,) + seq[0].chunks[axis:])

    names = [a.name for a in seq]
    name = 'stack-' + tokenize(names, axis)
    keys = list(product([name], *[range(len(bd)) for bd in chunks]))

    inputs = [(names[key[axis + 1]], ) + key[1:axis + 1] + key[axis + 2:]
              for key in keys]
    values = [(getitem, inp, (slice(None, None, None),) * axis +
              (None, ) + (slice(None, None, None), ) * (ndim - axis))
              for inp in inputs]

    dsk = dict(zip(keys, values))
    dsk2 = sharedict.merge((name, dsk), *[a.dask for a in seq])

    return Array(dsk2, name, chunks, dtype=dt)


def concatenate3(arrays):
    """ Recursive np.concatenate

    Input should be a nested list of numpy arrays arranged in the order they
    should appear in the array itself.  Each array should have the same number
    of dimensions as the desired output and the nesting of the lists.

    >>> x = np.array([[1, 2]])
    >>> concatenate3([[x, x, x], [x, x, x]])
    array([[1, 2, 1, 2, 1, 2],
           [1, 2, 1, 2, 1, 2]])

    >>> concatenate3([[x, x], [x, x], [x, x]])
    array([[1, 2, 1, 2],
           [1, 2, 1, 2],
           [1, 2, 1, 2]])
    """
    arrays = concrete(arrays)
    if not arrays:
        return np.empty(0)

    advanced = max(core.flatten(arrays, container=(list, tuple)),
                   key=lambda x: getattr(x, '__array_priority__', 0))
    if concatenate_lookup.dispatch(type(advanced)) is not np.concatenate:
        x = unpack_singleton(arrays)
        return _concatenate2(arrays, axes=list(range(x.ndim)))

    ndim = ndimlist(arrays)
    if not ndim:
        return arrays
    chunks = chunks_from_arrays(arrays)
    shape = tuple(map(sum, chunks))

    def dtype(x):
        try:
            return x.dtype
        except AttributeError:
            return type(x)

    result = np.empty(shape=shape, dtype=dtype(deepfirst(arrays)))

    for (idx, arr) in zip(slices_from_chunks(chunks), core.flatten(arrays)):
        if hasattr(arr, 'ndim'):
            while arr.ndim < ndim:
                arr = arr[None, ...]
        result[idx] = arr

    return result


def concatenate_axes(arrays, axes):
    """ Recursively call np.concatenate along axes """
    if len(axes) != ndimlist(arrays):
        raise ValueError("Length of axes should equal depth of nested arrays")

    extradims = max(0, deepfirst(arrays).ndim - (max(axes) + 1))
    return concatenate3(transposelist(arrays, axes, extradims=extradims))


def to_hdf5(filename, *args, **kwargs):
    """ Store arrays in HDF5 file

    This saves several dask arrays into several datapaths in an HDF5 file.
    It creates the necessary datasets and handles clean file opening/closing.

    >>> da.to_hdf5('myfile.hdf5', '/x', x)  # doctest: +SKIP

    or

    >>> da.to_hdf5('myfile.hdf5', {'/x': x, '/y': y})  # doctest: +SKIP

    Optionally provide arguments as though to ``h5py.File.create_dataset``

    >>> da.to_hdf5('myfile.hdf5', '/x', x, compression='lzf', shuffle=True)  # doctest: +SKIP

    This can also be used as a method on a single Array

    >>> x.to_hdf5('myfile.hdf5', '/x')  # doctest: +SKIP

    See Also
    --------
    da.store
    h5py.File.create_dataset
    """
    if len(args) == 1 and isinstance(args[0], dict):
        data = args[0]
    elif (len(args) == 2 and
          isinstance(args[0], str) and
          isinstance(args[1], Array)):
        data = {args[0]: args[1]}
    else:
        raise ValueError("Please provide {'/data/path': array} dictionary")

    chunks = kwargs.pop('chunks', True)

    import h5py
    with h5py.File(filename) as f:
        dsets = [f.require_dataset(dp, shape=x.shape, dtype=x.dtype,
                                   chunks=tuple([c[0] for c in x.chunks])
                                   if chunks is True else chunks, **kwargs)
                 for dp, x in data.items()]
        store(list(data.values()), dsets)


def interleave_none(a, b):
    """

    >>> interleave_none([0, None, 2, None], [1, 3])
    (0, 1, 2, 3)
    """
    result = []
    i = j = 0
    n = len(a) + len(b)
    while i + j < n:
        if a[i] is not None:
            result.append(a[i])
            i += 1
        else:
            result.append(b[j])
            i += 1
            j += 1
    return tuple(result)


def keyname(name, i, okey):
    """

    >>> keyname('x', 3, [None, None, 0, 2])
    ('x', 3, 0, 2)
    """
    return (name, i) + tuple(k for k in okey if k is not None)


def _vindex(x, *indexes):
    """Point wise indexing with broadcasting.

    >>> x = np.arange(56).reshape((7, 8))
    >>> x
    array([[ 0,  1,  2,  3,  4,  5,  6,  7],
           [ 8,  9, 10, 11, 12, 13, 14, 15],
           [16, 17, 18, 19, 20, 21, 22, 23],
           [24, 25, 26, 27, 28, 29, 30, 31],
           [32, 33, 34, 35, 36, 37, 38, 39],
           [40, 41, 42, 43, 44, 45, 46, 47],
           [48, 49, 50, 51, 52, 53, 54, 55]])

    >>> d = from_array(x, chunks=(3, 4))
    >>> result = _vindex(d, [0, 1, 6, 0], [0, 1, 0, 7])
    >>> result.compute()
    array([ 0,  9, 48,  7])
    """
    indexes = replace_ellipsis(x.ndim, indexes)

    nonfancy_indexes = []
    reduced_indexes = []
    for i, ind in enumerate(indexes):
        if isinstance(ind, Number):
            nonfancy_indexes.append(ind)
        elif isinstance(ind, slice):
            nonfancy_indexes.append(ind)
            reduced_indexes.append(slice(None))
        else:
            nonfancy_indexes.append(slice(None))
            reduced_indexes.append(ind)

    nonfancy_indexes = tuple(nonfancy_indexes)
    reduced_indexes = tuple(reduced_indexes)

    x = x[nonfancy_indexes]

    array_indexes = {}
    for i, (ind, size) in enumerate(zip(reduced_indexes, x.shape)):
        if not isinstance(ind, slice):
            ind = np.array(ind, copy=True)
            if ind.dtype.kind == 'b':
                raise IndexError('vindex does not support indexing with '
                                 'boolean arrays')
            if ((ind >= size) | (ind < -size)).any():
                raise IndexError('vindex key has entries out of bounds for '
                                 'indexing along axis %s of size %s: %r'
                                 % (i, size, ind))
            ind %= size
            array_indexes[i] = ind

    if array_indexes:
        x = _vindex_array(x, array_indexes)

    return x


def _vindex_array(x, dict_indexes):
    """Point wise indexing with only NumPy Arrays."""

    try:
        broadcast_indexes = np.broadcast_arrays(*dict_indexes.values())
    except ValueError:
        # note: error message exactly matches numpy
        shapes_str = ' '.join(str(a.shape) for a in dict_indexes.values())
        raise IndexError('shape mismatch: indexing arrays could not be '
                         'broadcast together with shapes ' + shapes_str)
    broadcast_shape = broadcast_indexes[0].shape

    lookup = dict(zip(dict_indexes, broadcast_indexes))
    flat_indexes = [lookup[i].ravel().tolist() if i in lookup else None
                    for i in range(x.ndim)]
    flat_indexes.extend([None] * (x.ndim - len(flat_indexes)))

    flat_indexes = [
        list(index) if index is not None else index for index in flat_indexes
    ]
    bounds = [list(accumulate(add, (0,) + c)) for c in x.chunks]
    bounds2 = [
        b for i, b in zip(flat_indexes, bounds) if i is not None
    ]
    axis = _get_axis(flat_indexes)
    token = tokenize(x, flat_indexes)
    out_name = 'vindex-merge-' + token

    points = list()
    for i, idx in enumerate(zip(*[i for i in flat_indexes if i is not None])):
        block_idx = [np.searchsorted(b, ind, 'right') - 1
                     for b, ind in zip(bounds2, idx)]
        inblock_idx = [ind - bounds2[k][j]
                       for k, (ind, j) in enumerate(zip(idx, block_idx))]
        points.append((i, tuple(block_idx), tuple(inblock_idx)))

    chunks = [c for i, c in zip(flat_indexes, x.chunks) if i is None]
    chunks.insert(0, (len(points),) if points else (0,))
    chunks = tuple(chunks)

    if points:
        per_block = groupby(1, points)
        per_block = dict((k, v) for k, v in per_block.items() if v)

        other_blocks = list(product(*[list(range(len(c))) if i is None else [None]
                                    for i, c in zip(flat_indexes, x.chunks)]))

        full_slices = [
            slice(None, None) if i is None else None for i in flat_indexes
        ]

        name = 'vindex-slice-' + token
        dsk = dict((keyname(name, i, okey),
                    (_vindex_transpose,
                    (_vindex_slice, (x.name,) + interleave_none(okey, key),
                     interleave_none(full_slices, list(zip(*pluck(2, per_block[key]))))),
                     axis))
                   for i, key in enumerate(per_block)
                   for okey in other_blocks)

        dsk.update((keyname('vindex-merge-' + token, 0, okey),
                   (_vindex_merge,
                    [list(pluck(0, per_block[key])) for key in per_block],
                    [keyname(name, i, okey) for i in range(len(per_block))]))
                   for okey in other_blocks)

        result_1d = Array(
            sharedict.merge(x.dask, (out_name, dsk)), out_name, chunks, x.dtype
        )
        return result_1d.reshape(broadcast_shape + result_1d.shape[1:])

    # output has a zero dimension, just create a new zero-shape array with the
    # same dtype
    from .wrap import empty
    result_1d = empty(
        tuple(map(sum, chunks)), chunks=chunks, dtype=x.dtype, name=out_name
    )
    return result_1d.reshape(broadcast_shape + result_1d.shape[1:])


def _get_axis(indexes):
    """ Get axis along which point-wise slicing results lie

    This is mostly a hack because I can't figure out NumPy's rule on this and
    can't be bothered to go reading.

    >>> _get_axis([[1, 2], None, [1, 2], None])
    0
    >>> _get_axis([None, [1, 2], [1, 2], None])
    1
    >>> _get_axis([None, None, [1, 2], [1, 2]])
    2
    """
    ndim = len(indexes)
    indexes = [slice(None, None) if i is None else [0] for i in indexes]
    x = np.empty((2,) * ndim)
    x2 = x[tuple(indexes)]
    return x2.shape.index(1)


def _vindex_slice(block, points):
    """ Pull out point-wise slices from block """
    points = [p if isinstance(p, slice) else list(p) for p in points]
    return block[tuple(points)]


def _vindex_transpose(block, axis):
    """ Rotate block so that points are on the first dimension """
    axes = [axis] + list(range(axis)) + list(range(axis + 1, block.ndim))
    return block.transpose(axes)


def _vindex_merge(locations, values):
    """

    >>> locations = [0], [2, 1]
    >>> values = [np.array([[1, 2, 3]]),
    ...           np.array([[10, 20, 30], [40, 50, 60]])]

    >>> _vindex_merge(locations, values)
    array([[ 1,  2,  3],
           [40, 50, 60],
           [10, 20, 30]])
    """
    locations = list(map(list, locations))
    values = list(values)

    n = sum(map(len, locations))

    shape = list(values[0].shape)
    shape[0] = n
    shape = tuple(shape)

    dtype = values[0].dtype

    x = np.empty(shape, dtype=dtype)

    ind = [slice(None, None) for i in range(x.ndim)]
    for loc, val in zip(locations, values):
        ind[0] = loc
        x[tuple(ind)] = val

    return x


def to_npy_stack(dirname, x, axis=0):
    """ Write dask array to a stack of .npy files

    This partitions the dask.array along one axis and stores each block along
    that axis as a single .npy file in the specified directory

    Examples
    --------
    >>> x = da.ones((5, 10, 10), chunks=(2, 4, 4))  # doctest: +SKIP
    >>> da.to_npy_stack('data/', x, axis=0)  # doctest: +SKIP

        $ tree data/
        data/
        |-- 0.npy
        |-- 1.npy
        |-- 2.npy
        |-- info

    The ``.npy`` files store numpy arrays for ``x[0:2], x[2:4], and x[4:5]``
    respectively, as is specified by the chunk size along the zeroth axis.  The
    info file stores the dtype, chunks, and axis information of the array.

    You can load these stacks with the ``da.from_npy_stack`` function.

    >>> y = da.from_npy_stack('data/')  # doctest: +SKIP

    See Also
    --------
    from_npy_stack
    """

    chunks = tuple((c if i == axis else (sum(c),))
                   for i, c in enumerate(x.chunks))
    xx = x.rechunk(chunks)

    if not os.path.exists(dirname):
        os.mkdir(dirname)

    meta = {'chunks': chunks, 'dtype': x.dtype, 'axis': axis}

    with open(os.path.join(dirname, 'info'), 'wb') as f:
        pickle.dump(meta, f)

    name = 'to-npy-stack-' + str(uuid.uuid1())
    dsk = {(name, i): (np.save, os.path.join(dirname, '%d.npy' % i), key)
           for i, key in enumerate(core.flatten(xx.__dask_keys__()))}

    compute_as_if_collection(Array, sharedict.merge(dsk, xx.dask), list(dsk))


def from_npy_stack(dirname, mmap_mode='r'):
    """ Load dask array from stack of npy files

    See ``da.to_npy_stack`` for docstring

    Parameters
    ----------
    dirname: string
        Directory of .npy files
    mmap_mode: (None or 'r')
        Read data in memory map mode
    """
    with open(os.path.join(dirname, 'info'), 'rb') as f:
        info = pickle.load(f)

    dtype = info['dtype']
    chunks = info['chunks']
    axis = info['axis']

    name = 'from-npy-stack-%s' % dirname
    keys = list(product([name], *[range(len(c)) for c in chunks]))
    values = [(np.load, os.path.join(dirname, '%d.npy' % i), mmap_mode)
              for i in range(len(chunks[axis]))]
    dsk = dict(zip(keys, values))

    return Array(dsk, name, chunks, dtype)
