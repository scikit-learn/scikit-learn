from __future__ import absolute_import, division, print_function
import copy

import pytest
np = pytest.importorskip('numpy')

import os
import sys
import time
from distutils.version import LooseVersion
import operator
from operator import add, sub, getitem
from threading import Lock
import warnings

from toolz import merge, countby, concat
from toolz.curried import identity

import dask
import dask.array as da
from dask.base import tokenize, compute_as_if_collection
from dask.delayed import Delayed, delayed
from dask.local import get_sync
from dask.utils import ignoring, tmpfile, tmpdir
from dask.utils_test import inc

from dask.array import chunk

from dask.array.core import (getem, getter, top, dotmany, concatenate3,
                             broadcast_dimensions, Array, stack, concatenate,
                             from_array, elemwise, broadcast_shapes,
                             broadcast_to, blockdims_from_blockshape, store,
                             optimize, from_func, normalize_chunks,
                             broadcast_chunks, atop, from_delayed,
                             concatenate_axes, common_blockdim)
from dask.array.utils import assert_eq, same_keys

# temporary until numpy functions migrated
try:
    from numpy import nancumsum, nancumprod
except ImportError:  # pragma: no cover
    import dask.array.numpy_compat as npcompat
    nancumsum = npcompat.nancumsum
    nancumprod = npcompat.nancumprod


def test_getem():
    sol = {('X', 0, 0): (getter, 'X', (slice(0, 2), slice(0, 3))),
           ('X', 1, 0): (getter, 'X', (slice(2, 4), slice(0, 3))),
           ('X', 1, 1): (getter, 'X', (slice(2, 4), slice(3, 6))),
           ('X', 0, 1): (getter, 'X', (slice(0, 2), slice(3, 6)))}
    assert getem('X', (2, 3), shape=(4, 6)) == sol


def test_top():
    assert top(inc, 'z', 'ij', 'x', 'ij', numblocks={'x': (2, 2)}) == \
        {('z', 0, 0): (inc, ('x', 0, 0)),
         ('z', 0, 1): (inc, ('x', 0, 1)),
         ('z', 1, 0): (inc, ('x', 1, 0)),
         ('z', 1, 1): (inc, ('x', 1, 1))}

    assert top(add, 'z', 'ij', 'x', 'ij', 'y', 'ij',
               numblocks={'x': (2, 2), 'y': (2, 2)}) == \
        {('z', 0, 0): (add, ('x', 0, 0), ('y', 0, 0)),
         ('z', 0, 1): (add, ('x', 0, 1), ('y', 0, 1)),
         ('z', 1, 0): (add, ('x', 1, 0), ('y', 1, 0)),
         ('z', 1, 1): (add, ('x', 1, 1), ('y', 1, 1))}

    assert top(dotmany, 'z', 'ik', 'x', 'ij', 'y', 'jk',
               numblocks={'x': (2, 2), 'y': (2, 2)}) == \
        {('z', 0, 0): (dotmany, [('x', 0, 0), ('x', 0, 1)],
                                [('y', 0, 0), ('y', 1, 0)]),
         ('z', 0, 1): (dotmany, [('x', 0, 0), ('x', 0, 1)],
                                [('y', 0, 1), ('y', 1, 1)]),
         ('z', 1, 0): (dotmany, [('x', 1, 0), ('x', 1, 1)],
                                [('y', 0, 0), ('y', 1, 0)]),
         ('z', 1, 1): (dotmany, [('x', 1, 0), ('x', 1, 1)],
                                [('y', 0, 1), ('y', 1, 1)])}

    assert top(identity, 'z', '', 'x', 'ij', numblocks={'x': (2, 2)}) ==\
        {('z',): (identity, [[('x', 0, 0), ('x', 0, 1)],
                             [('x', 1, 0), ('x', 1, 1)]])}


def test_top_supports_broadcasting_rules():
    assert top(add, 'z', 'ij', 'x', 'ij', 'y', 'ij',
               numblocks={'x': (1, 2), 'y': (2, 1)}) == \
        {('z', 0, 0): (add, ('x', 0, 0), ('y', 0, 0)),
         ('z', 0, 1): (add, ('x', 0, 1), ('y', 0, 0)),
         ('z', 1, 0): (add, ('x', 0, 0), ('y', 1, 0)),
         ('z', 1, 1): (add, ('x', 0, 1), ('y', 1, 0))}


def test_top_literals():
    assert top(add, 'z', 'ij', 'x', 'ij', 123, None, numblocks={'x': (2, 2)}) == \
        {('z', 0, 0): (add, ('x', 0, 0), 123),
         ('z', 0, 1): (add, ('x', 0, 1), 123),
         ('z', 1, 0): (add, ('x', 1, 0), 123),
         ('z', 1, 1): (add, ('x', 1, 1), 123)}


def test_atop_literals():
    x = da.ones((10, 10), chunks=(5, 5))
    z = atop(add, 'ij', x, 'ij', 100, None, dtype=x.dtype)
    assert_eq(z, x + 100)

    z = atop(lambda x, y, z: x * y + z, 'ij', 2, None,  x, 'ij', 100, None, dtype=x.dtype)
    assert_eq(z, 2 * x + 100)

    z = atop(getitem, 'ij', x, 'ij', slice(None), None, dtype=x.dtype)
    assert_eq(z, x)


def test_concatenate3_on_scalars():
    assert_eq(concatenate3([1, 2]), np.array([1, 2]))


def test_chunked_dot_product():
    x = np.arange(400).reshape((20, 20))
    o = np.ones((20, 20))

    d = {'x': x, 'o': o}

    getx = getem('x', (5, 5), shape=(20, 20))
    geto = getem('o', (5, 5), shape=(20, 20))

    result = top(dotmany, 'out', 'ik', 'x', 'ij', 'o', 'jk',
                 numblocks={'x': (4, 4), 'o': (4, 4)})

    dsk = merge(d, getx, geto, result)
    out = dask.get(dsk, [[('out', i, j) for j in range(4)] for i in range(4)])

    assert_eq(np.dot(x, o), concatenate3(out))


def test_chunked_transpose_plus_one():
    x = np.arange(400).reshape((20, 20))

    d = {'x': x}

    getx = getem('x', (5, 5), shape=(20, 20))

    f = lambda x: x.T + 1
    comp = top(f, 'out', 'ij', 'x', 'ji', numblocks={'x': (4, 4)})

    dsk = merge(d, getx, comp)
    out = dask.get(dsk, [[('out', i, j) for j in range(4)] for i in range(4)])

    assert_eq(concatenate3(out), x.T + 1)


def test_broadcast_dimensions_works_with_singleton_dimensions():
    argpairs = [('x', 'i')]
    numblocks = {'x': ((1,),)}
    assert broadcast_dimensions(argpairs, numblocks) == {'i': (1,)}


def test_broadcast_dimensions():
    argpairs = [('x', 'ij'), ('y', 'ij')]
    d = {'x': ('Hello', 1), 'y': (1, (2, 3))}
    assert broadcast_dimensions(argpairs, d) == {'i': 'Hello', 'j': (2, 3)}


def test_Array():
    shape = (1000, 1000)
    chunks = (100, 100)
    name = 'x'
    dsk = merge({name: 'some-array'}, getem(name, chunks, shape=shape))
    a = Array(dsk, name, chunks, shape=shape, dtype='f8')

    assert a.numblocks == (10, 10)

    assert a.__dask_keys__() == [[('x', i, j) for j in range(10)]
                                 for i in range(10)]

    assert a.chunks == ((100,) * 10, (100,) * 10)

    assert a.shape == shape

    assert len(a) == shape[0]


def test_uneven_chunks():
    a = Array({}, 'x', chunks=(3, 3), shape=(10, 10), dtype='f8')
    assert a.chunks == ((3, 3, 3, 1), (3, 3, 3, 1))


def test_numblocks_suppoorts_singleton_block_dims():
    shape = (100, 10)
    chunks = (10, 10)
    name = 'x'
    dsk = merge({name: 'some-array'}, getem(name, shape=shape, chunks=chunks))
    a = Array(dsk, name, chunks, shape=shape, dtype='f8')

    assert set(concat(a.__dask_keys__())) == {('x', i, 0) for i in range(10)}


def test_keys():
    dsk = dict((('x', i, j), ()) for i in range(5) for j in range(6))
    dx = Array(dsk, 'x', chunks=(10, 10), shape=(50, 60), dtype='f8')
    assert dx.__dask_keys__() == [[(dx.name, i, j) for j in range(6)]
                                  for i in range(5)]
    # Cache works
    assert dx.__dask_keys__() is dx.__dask_keys__()
    # Test mutating names clears key cache
    dx.dask = {('y', i, j): () for i in range(5) for j in range(6)}
    dx.name = 'y'
    assert dx.__dask_keys__() == [[(dx.name, i, j) for j in range(6)]
                                  for i in range(5)]
    d = Array({}, 'x', (), shape=(), dtype='f8')
    assert d.__dask_keys__() == [('x',)]


def test_Array_computation():
    a = Array({('x', 0, 0): np.eye(3)}, 'x', shape=(3, 3), chunks=(3, 3), dtype='f8')
    assert_eq(np.array(a), np.eye(3))
    assert isinstance(a.compute(), np.ndarray)
    assert float(a[0, 0]) == 1


def test_stack():
    a, b, c = [Array(getem(name, chunks=(2, 3), shape=(4, 6)),
                     name, chunks=(2, 3), dtype='f8', shape=(4, 6))
               for name in 'ABC']

    s = stack([a, b, c], axis=0)

    colon = slice(None, None, None)

    assert s.shape == (3, 4, 6)
    assert s.chunks == ((1, 1, 1), (2, 2), (3, 3))
    assert s.dask[(s.name, 0, 1, 0)] == (getitem, ('A', 1, 0),
                                         (None, colon, colon))
    assert s.dask[(s.name, 2, 1, 0)] == (getitem, ('C', 1, 0),
                                         (None, colon, colon))
    assert same_keys(s, stack([a, b, c], axis=0))

    s2 = stack([a, b, c], axis=1)
    assert s2.shape == (4, 3, 6)
    assert s2.chunks == ((2, 2), (1, 1, 1), (3, 3))
    assert s2.dask[(s2.name, 0, 1, 0)] == (getitem, ('B', 0, 0),
                                           (colon, None, colon))
    assert s2.dask[(s2.name, 1, 1, 0)] == (getitem, ('B', 1, 0),
                                           (colon, None, colon))
    assert same_keys(s2, stack([a, b, c], axis=1))

    s2 = stack([a, b, c], axis=2)
    assert s2.shape == (4, 6, 3)
    assert s2.chunks == ((2, 2), (3, 3), (1, 1, 1))
    assert s2.dask[(s2.name, 0, 1, 0)] == (getitem, ('A', 0, 1),
                                           (colon, colon, None))
    assert s2.dask[(s2.name, 1, 1, 2)] == (getitem, ('C', 1, 1),
                                           (colon, colon, None))
    assert same_keys(s2, stack([a, b, c], axis=2))

    pytest.raises(ValueError, lambda: stack([a, b, c], axis=3))

    assert set(b.dask.keys()).issubset(s2.dask.keys())

    assert stack([a, b, c], axis=-1).chunks == stack([a, b, c], axis=2).chunks


def test_short_stack():
    x = np.array([1])
    d = da.from_array(x, chunks=(1,))
    s = da.stack([d])
    assert s.shape == (1, 1)
    chunks = compute_as_if_collection(Array, s.dask, s.__dask_keys__())
    assert chunks[0][0].shape == (1, 1)


def test_stack_scalars():
    d = da.arange(4, chunks=2)

    s = da.stack([d.mean(), d.sum()])

    assert s.compute().tolist() == [np.arange(4).mean(), np.arange(4).sum()]


def test_stack_promote_type():
    i = np.arange(10, dtype='i4')
    f = np.arange(10, dtype='f4')
    di = da.from_array(i, chunks=5)
    df = da.from_array(f, chunks=5)
    res = da.stack([di, df])
    assert_eq(res, np.stack([i, f]))


def test_stack_rechunk():
    x = da.random.random(10, chunks=5)
    y = da.random.random(10, chunks=4)

    z = da.stack([x, y], axis=0)
    assert z.shape == (2, 10)
    assert z.chunks == ((1, 1), (4, 1, 3, 2))

    assert_eq(z, np.stack([x.compute(), y.compute()], axis=0))


def test_concatenate():
    a, b, c = [Array(getem(name, chunks=(2, 3), shape=(4, 6)),
                     name, chunks=(2, 3), dtype='f8', shape=(4, 6))
               for name in 'ABC']

    x = concatenate([a, b, c], axis=0)

    assert x.shape == (12, 6)
    assert x.chunks == ((2, 2, 2, 2, 2, 2), (3, 3))
    assert x.dask[(x.name, 0, 1)] == ('A', 0, 1)
    assert x.dask[(x.name, 5, 0)] == ('C', 1, 0)
    assert same_keys(x, concatenate([a, b, c], axis=0))

    y = concatenate([a, b, c], axis=1)

    assert y.shape == (4, 18)
    assert y.chunks == ((2, 2), (3, 3, 3, 3, 3, 3))
    assert y.dask[(y.name, 1, 0)] == ('A', 1, 0)
    assert y.dask[(y.name, 1, 5)] == ('C', 1, 1)
    assert same_keys(y, concatenate([a, b, c], axis=1))

    assert set(b.dask.keys()).issubset(y.dask.keys())

    z = concatenate([a], axis=0)

    assert z.shape == a.shape
    assert z.chunks == a.chunks
    assert z.dask == a.dask
    assert z is a

    assert (concatenate([a, b, c], axis=-1).chunks ==
            concatenate([a, b, c], axis=1).chunks)

    pytest.raises(ValueError, lambda: concatenate([a, b, c], axis=2))


@pytest.mark.parametrize('dtypes', [(('>f8', '>f8'), '>f8'),
                                    (('<f4', '<f8'), '<f8')])
def test_concatenate_types(dtypes):
    dts_in, dt_out = dtypes
    arrs = [np.zeros(4, dtype=dt) for dt in dts_in]
    darrs = [from_array(arr, chunks=(2,)) for arr in arrs]

    x = concatenate(darrs, axis=0)
    assert x.dtype == dt_out


def test_concatenate_unknown_axes():
    dd = pytest.importorskip('dask.dataframe')
    pd = pytest.importorskip('pandas')

    a_df = pd.DataFrame({'x': np.arange(12)})
    b_df = pd.DataFrame({'y': np.arange(12) * 10})

    a_ddf = dd.from_pandas(a_df, sort=False, npartitions=3)
    b_ddf = dd.from_pandas(b_df, sort=False, npartitions=3)

    a_x = a_ddf.values
    b_x = b_ddf.values

    assert np.isnan(a_x.shape[0])
    assert np.isnan(b_x.shape[0])

    da.concatenate([a_x, b_x], axis=0)  # works fine

    with pytest.raises(ValueError) as exc_info:
        da.concatenate([a_x, b_x], axis=1)  # unknown chunks

    assert 'nan' in str(exc_info.value)
    assert 'allow_unknown_chunksize' in str(exc_info.value)

    c_x = da.concatenate([a_x, b_x], axis=1, allow_unknown_chunksizes=True)  # unknown chunks

    assert_eq(c_x, np.concatenate([a_df.values, b_df.values], axis=1))


def test_concatenate_rechunk():
    x = da.random.random((6, 6), chunks=(3, 3))
    y = da.random.random((6, 6), chunks=(2, 2))

    z = da.concatenate([x, y], axis=0)
    assert z.shape == (12, 6)
    assert z.chunks == ((3, 3, 2, 2, 2), (2, 1, 1, 2))
    assert_eq(z, np.concatenate([x.compute(), y.compute()], axis=0))

    z = da.concatenate([x, y], axis=1)
    assert z.shape == (6, 12)
    assert z.chunks == ((2, 1, 1, 2), (3, 3, 2, 2, 2))
    assert_eq(z, np.concatenate([x.compute(), y.compute()], axis=1))


def test_concatenate_fixlen_strings():
    x = np.array(['a', 'b', 'c'])
    y = np.array(['aa', 'bb', 'cc'])

    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))

    assert_eq(np.concatenate([x, y]),
              da.concatenate([a, b]))


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_simple_row_wise():
    a1 = np.ones((2, 2))
    a2 = 2 * a1

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([a1, a2])
    result = da.block([d1, d2])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_simple_column_wise():
    a1 = np.ones((2, 2))
    a2 = 2 * a1

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1], [a2]])
    result = da.block([[d1], [d2]])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_with_1d_arrays_row_wise():
    # # # 1-D vectors are treated as row arrays
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([a1, a2])
    result = da.block([d1, d2])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_with_1d_arrays_multiple_rows():
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1, a2], [a1, a2]])
    result = da.block([[d1, d2], [d1, d2]])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_with_1d_arrays_column_wise():
    # # # 1-D vectors are treated as row arrays
    a1 = np.array([1, 2, 3])
    a2 = np.array([2, 3, 4])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[a1], [a2]])
    result = da.block([[d1], [d2]])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_mixed_1d_and_2d():
    a1 = np.ones((2, 2))
    a2 = np.array([2, 2])

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)

    expected = np.block([[d1], [d2]])
    result = da.block([[a1], [a2]])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_complicated():
    # a bit more complicated
    a1 = np.array([[1, 1, 1]])
    a2 = np.array([[2, 2, 2]])
    a3 = np.array([[3, 3, 3, 3, 3, 3]])
    a4 = np.array([4, 4, 4, 4, 4, 4])
    a5 = np.array(5)
    a6 = np.array([6, 6, 6, 6, 6])
    a7 = np.zeros((2, 6))

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)
    d3 = da.asarray(a3)
    d4 = da.asarray(a4)
    d5 = da.asarray(a5)
    d6 = da.asarray(a6)
    d7 = da.asarray(a7)

    expected = np.block([[a1, a2],
                         [a3],
                         [a4],
                         [a5, a6],
                         [a7]])
    result = da.block([[d1, d2],
                       [d3],
                       [d4],
                       [d5, d6],
                       [d7]])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_nested():
    a1 = np.array([1, 1, 1])
    a2 = np.array([[2, 2, 2], [2, 2, 2], [2, 2, 2]])
    a3 = np.array([3, 3, 3])
    a4 = np.array([4, 4, 4])
    a5 = np.array(5)
    a6 = np.array([6, 6, 6, 6, 6])
    a7 = np.zeros((2, 6))

    d1 = da.asarray(a1)
    d2 = da.asarray(a2)
    d3 = da.asarray(a3)
    d4 = da.asarray(a4)
    d5 = da.asarray(a5)
    d6 = da.asarray(a6)
    d7 = da.asarray(a7)

    expected = np.block([
        [
            np.block([
                [a1],
                [a3],
                [a4]
            ]),
            a2
        ],
        [a5, a6],
        [a7]
    ])
    result = da.block([
        [
            da.block([
                [d1],
                [d3],
                [d4]
            ]),
            d2
        ],
        [d5, d6],
        [d7]
    ])

    assert_eq(expected, result)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_3d():
    a000 = np.ones((2, 2, 2), int) * 1

    a100 = np.ones((3, 2, 2), int) * 2
    a010 = np.ones((2, 3, 2), int) * 3
    a001 = np.ones((2, 2, 3), int) * 4

    a011 = np.ones((2, 3, 3), int) * 5
    a101 = np.ones((3, 2, 3), int) * 6
    a110 = np.ones((3, 3, 2), int) * 7

    a111 = np.ones((3, 3, 3), int) * 8

    d000 = da.asarray(a000)

    d100 = da.asarray(a100)
    d010 = da.asarray(a010)
    d001 = da.asarray(a001)

    d011 = da.asarray(a011)
    d101 = da.asarray(a101)
    d110 = da.asarray(a110)

    d111 = da.asarray(a111)

    expected = np.block([
        [
            [a000, a001],
            [a010, a011],
        ],
        [
            [a100, a101],
            [a110, a111],
        ]
    ])
    result = da.block([
        [
            [d000, d001],
            [d010, d011],
        ],
        [
            [d100, d101],
            [d110, d111],
        ]
    ])

    assert_eq(expected, result)


def test_block_with_mismatched_shape():
    a = np.array([0, 0])
    b = np.eye(2)

    for arrays in [[a, b],
                   [b, a]]:
        with pytest.raises(ValueError):
            da.block(arrays)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="NumPy doesn't support `block` yet")
def test_block_no_lists():
    assert_eq(da.block(1),         np.block(1))
    assert_eq(da.block(np.eye(3)), np.block(np.eye(3)))


def test_block_invalid_nesting():
    for arrays in [
            [1, [2]],
            [1, []],
            [[1], 2],
            [[], 2],
            [
                [[1], [2]],
                [[3, 4]],
                [5]  # missing brackets
            ],
    ]:
        with pytest.raises(ValueError) as e:
            da.block(arrays)
        e.match(r'depths are mismatched')


def test_block_empty_lists():
    for arrays in [
            [],
            [[]],
            [[1], []],
    ]:
        with pytest.raises(ValueError) as e:
            da.block(arrays)
        e.match(r'empty')


def test_block_tuple():
    for arrays in [
            ([1, 2], [3, 4]),
            [(1, 2), (3, 4)],
    ]:
        with pytest.raises(TypeError) as e:
            da.block(arrays)
        e.match(r'tuple')


def test_binops():
    a = Array(dict((('a', i), np.array([0])) for i in range(3)),
              'a', chunks=((1, 1, 1),), dtype='i8')
    b = Array(dict((('b', i), np.array([0])) for i in range(3)),
              'b', chunks=((1, 1, 1),), dtype='i8')

    result = elemwise(add, a, b, name='c')
    assert result.dask == merge(a.dask, b.dask,
                                dict((('c', i), (add, ('a', i), ('b', i)))
                                     for i in range(3)))

    result = elemwise(pow, a, 2, name='c')
    assert "'a', 0" in str(result.dask[('c', 0)])
    assert "2" in str(result.dask[('c', 0)])


def test_broadcast_shapes():
    assert () == broadcast_shapes()
    assert (2, 5) == broadcast_shapes((2, 5))
    assert (0, 5) == broadcast_shapes((0, 1), (1, 5))
    assert np.allclose(
        (2, np.nan), broadcast_shapes((1, np.nan), (2, 1)), equal_nan=True
    )
    assert np.allclose(
        (2, np.nan), broadcast_shapes((2, 1), (1, np.nan)), equal_nan=True
    )
    assert (3, 4, 5) == broadcast_shapes((3, 4, 5), (4, 1), ())
    assert (3, 4) == broadcast_shapes((3, 1), (1, 4), (4,))
    assert (5, 6, 7, 3, 4) == broadcast_shapes((3, 1), (), (5, 6, 7, 1, 4))
    pytest.raises(ValueError, lambda: broadcast_shapes((3,), (3, 4)))
    pytest.raises(ValueError, lambda: broadcast_shapes((2, 3), (2, 3, 1)))
    pytest.raises(ValueError, lambda: broadcast_shapes((2, 3), (1, np.nan)))


def test_elemwise_on_scalars():
    x = np.arange(10, dtype=np.int64)
    a = from_array(x, chunks=(5,))
    assert len(a.__dask_keys__()) == 2
    assert_eq(a.sum()**2, x.sum()**2)

    y = np.arange(10, dtype=np.int32)
    b = from_array(y, chunks=(5,))
    result = a.sum() * b
    # Dask 0-d arrays do not behave like numpy scalars for type promotion
    assert result.dtype == np.int64
    assert result.compute().dtype == np.int64
    assert (x.sum() * y).dtype == np.int32
    assert_eq((x.sum() * y).astype(np.int64), result)


def test_elemwise_with_ndarrays():
    x = np.arange(3)
    y = np.arange(12).reshape(4, 3)
    a = from_array(x, chunks=(3,))
    b = from_array(y, chunks=(2, 3))

    assert_eq(x + a, 2 * x)
    assert_eq(a + x, 2 * x)

    assert_eq(x + b, x + y)
    assert_eq(b + x, x + y)
    assert_eq(a + y, x + y)
    assert_eq(y + a, x + y)
    # Error on shape mismatch
    pytest.raises(ValueError, lambda: a + y.T)
    pytest.raises(ValueError, lambda: a + np.arange(2))


def test_elemwise_differently_chunked():
    x = np.arange(3)
    y = np.arange(12).reshape(4, 3)
    a = from_array(x, chunks=(3,))
    b = from_array(y, chunks=(2, 2))

    assert_eq(a + b, x + y)
    assert_eq(b + a, x + y)


def test_elemwise_dtype():
    values = [
        da.from_array(np.ones(5, np.float32), chunks=3),
        da.from_array(np.ones(5, np.int16), chunks=3),
        da.from_array(np.ones(5, np.int64), chunks=3),
        da.from_array(np.ones((), np.float64), chunks=()) * 1e200,
        np.ones(5, np.float32),
        1, 1.0, 1e200, np.int64(1), np.ones((), np.int64),
    ]
    for x in values:
        for y in values:
            assert da.maximum(x, y).dtype == da.result_type(x, y)


def test_operators():
    x = np.arange(10)
    y = np.arange(10).reshape((10, 1))
    a = from_array(x, chunks=(5,))
    b = from_array(y, chunks=(5, 1))

    c = a + 1
    assert_eq(c, x + 1)

    c = a + b
    assert_eq(c, x + x.reshape((10, 1)))

    expr = (3 / a * b)**2 > 5
    with pytest.warns(None):  # ZeroDivisionWarning
        assert_eq(expr, (3 / x * y)**2 > 5)

    with pytest.warns(None):  # OverflowWarning
        c = da.exp(a)
    assert_eq(c, np.exp(x))

    assert_eq(abs(-a), a)
    assert_eq(a, +x)


def test_operator_dtype_promotion():
    x = np.arange(10, dtype=np.float32)
    y = np.array([1])
    a = from_array(x, chunks=(5,))

    assert_eq(x + 1, a + 1)  # still float32
    assert_eq(x + 1e50, a + 1e50)  # now float64
    assert_eq(x + y, a + y)  # also float64


def test_field_access():
    x = np.array([(1, 1.0), (2, 2.0)], dtype=[('a', 'i4'), ('b', 'f4')])
    y = from_array(x, chunks=(1,))
    assert_eq(y['a'], x['a'])
    assert_eq(y[['b', 'a']], x[['b', 'a']])
    assert same_keys(y[['b', 'a']], y[['b', 'a']])


def test_field_access_with_shape():
    dtype = [('col1', ('f4', (3, 2))), ('col2', ('f4', 3))]
    data = np.ones((100, 50), dtype=dtype)
    x = da.from_array(data, 10)
    assert_eq(x['col1'], data['col1'])
    assert_eq(x[['col1']], data[['col1']])
    assert_eq(x['col2'], data['col2'])
    assert_eq(x[['col1', 'col2']], data[['col1', 'col2']])


@pytest.mark.skipif(sys.version_info < (3, 5),
                    reason="Matrix multiplication operator only after Py3.5")
def test_matmul():
    x = np.random.random((5, 5))
    y = np.random.random((5, 2))
    a = from_array(x, chunks=(1, 5))
    b = from_array(y, chunks=(5, 1))
    assert_eq(operator.matmul(a, b), a.dot(b))
    assert_eq(operator.matmul(a, b), operator.matmul(x, y))
    assert_eq(operator.matmul(a, y), operator.matmul(x, b))
    list_vec = list(range(1, 6))
    assert_eq(operator.matmul(list_vec, b), operator.matmul(list_vec, y))
    assert_eq(operator.matmul(x, list_vec), operator.matmul(a, list_vec))
    z = np.random.random((5, 5, 5))
    c = from_array(z, chunks=(1, 5, 1))
    assert_eq(operator.matmul(a, z), operator.matmul(x, c))
    assert_eq(operator.matmul(z, a), operator.matmul(c, x))


def test_T():
    x = np.arange(400).reshape((20, 20))
    a = from_array(x, chunks=(5, 5))

    assert_eq(x.T, a.T)


def test_norm():
    a = np.arange(200, dtype='f8').reshape((20, 10))
    a = a + (a.max() - a) * 1j
    b = from_array(a, chunks=(5, 5))

    # TODO: Deprecated method, remove test when method removed
    with pytest.warns(UserWarning):
        assert_eq(b.vnorm(), np.linalg.norm(a))
        assert_eq(b.vnorm(ord=1), np.linalg.norm(a.flatten(), ord=1))
        assert_eq(b.vnorm(ord=4, axis=0), np.linalg.norm(a, ord=4, axis=0))
        assert b.vnorm(ord=4, axis=0, keepdims=True).ndim == b.ndim
        split_every = {0: 3, 1: 3}
        assert_eq(b.vnorm(ord=1, axis=0, split_every=split_every),
                  np.linalg.norm(a, ord=1, axis=0))
        assert_eq(b.vnorm(ord=np.inf, axis=0, split_every=split_every),
                  np.linalg.norm(a, ord=np.inf, axis=0))
        assert_eq(b.vnorm(ord=np.inf, split_every=split_every),
                  np.linalg.norm(a.flatten(), ord=np.inf))


def test_broadcast_to():
    x = np.random.randint(10, size=(5, 1, 6))
    a = from_array(x, chunks=(3, 1, 3))

    for shape in [a.shape, (5, 0, 6), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        xb = chunk.broadcast_to(x, shape)
        ab = broadcast_to(a, shape)

        assert_eq(xb, ab)

        if a.shape == ab.shape:
            assert a is ab

    pytest.raises(ValueError, lambda: broadcast_to(a, (2, 1, 6)))
    pytest.raises(ValueError, lambda: broadcast_to(a, (3,)))


def test_broadcast_to_array():
    x = np.random.randint(10, size=(5, 1, 6))

    for shape in [(5, 0, 6), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        a = np.broadcast_to(x, shape)
        d = broadcast_to(x, shape)

        assert_eq(a, d)


def test_broadcast_to_scalar():
    x = 5

    for shape in [tuple(), (0,), (2, 3), (5, 4, 6), (2, 5, 1, 6), (3, 4, 5, 4, 6)]:
        a = np.broadcast_to(x, shape)
        d = broadcast_to(x, shape)

        assert_eq(a, d)


def test_broadcast_to_chunks():
    x = np.random.randint(10, size=(5, 1, 6))
    a = from_array(x, chunks=(3, 1, 3))

    for shape, chunks, expected_chunks in [
            ((5, 3, 6), (3, -1, 3), ((3, 2), (3,), (3, 3))),
            ((5, 3, 6), (3, 1, 3), ((3, 2), (1, 1, 1,), (3, 3))),
            ((2, 5, 3, 6), (1, 3, 1, 3), ((1, 1), (3, 2), (1, 1, 1,), (3, 3)))]:
        xb = chunk.broadcast_to(x, shape)
        ab = broadcast_to(a, shape, chunks=chunks)
        assert_eq(xb, ab)
        assert ab.chunks == expected_chunks

    with pytest.raises(ValueError):
        broadcast_to(a, a.shape, chunks=((2, 3), (1,), (3, 3)))
    with pytest.raises(ValueError):
        broadcast_to(a, a.shape, chunks=((3, 2), (3,), (3, 3)))
    with pytest.raises(ValueError):
        broadcast_to(a, (5, 2, 6), chunks=((3, 2), (3,), (3, 3)))


def test_broadcast_arrays():
    # Calling `broadcast_arrays` with no arguments only works in NumPy 1.13.0+.
    if LooseVersion(np.__version__) >= LooseVersion("1.13.0"):
        assert np.broadcast_arrays() == da.broadcast_arrays()

    a = np.arange(4)
    d_a = da.from_array(a, chunks=tuple(s // 2 for s in a.shape))

    a_0 = np.arange(4)[None, :]
    a_1 = np.arange(4)[:, None]

    d_a_0 = d_a[None, :]
    d_a_1 = d_a[:, None]

    a_r = np.broadcast_arrays(a_0, a_1)
    d_r = da.broadcast_arrays(d_a_0, d_a_1)

    assert isinstance(d_r, list)
    assert len(a_r) == len(d_r)

    for e_a_r, e_d_r in zip(a_r, d_r):
        assert_eq(e_a_r, e_d_r)


@pytest.mark.parametrize('u_shape, v_shape', [
    [tuple(), (2, 3)],
    [(1,), (2, 3)],
    [(1, 1), (2, 3)],
    [(0, 3), (1, 3)],
    [(2, 0), (2, 1)],
    [(1, 0), (2, 1)],
    [(0, 1), (1, 3)],
])
def test_broadcast_operator(u_shape, v_shape):
    u = np.random.random(u_shape)
    v = np.random.random(v_shape)

    d_u = from_array(u, chunks=1)
    d_v = from_array(v, chunks=1)

    w = u * v
    d_w = d_u * d_v

    assert_eq(w, d_w)


@pytest.mark.parametrize('original_shape,new_shape,chunks', [
    ((10,), (10,), (3, 3, 4)),
    ((10,), (10, 1, 1), 5),
    ((10,), (1, 10,), 5),
    ((24,), (2, 3, 4), 12),
    ((1, 24,), (2, 3, 4), 12),
    ((2, 3, 4), (24,), (1, 3, 4)),
    ((2, 3, 4), (24,), 4),
    ((2, 3, 4), (24, 1), 4),
    ((2, 3, 4), (1, 24), 4),
    ((4, 4, 1), (4, 4), 2),
    ((4, 4), (4, 4, 1), 2),
    ((1, 4, 4), (4, 4), 2),
    ((1, 4, 4), (4, 4, 1), 2),
    ((1, 4, 4), (1, 1, 4, 4), 2),
    ((4, 4), (1, 4, 4, 1), 2),
    ((4, 4), (1, 4, 4), 2),
    ((2, 3), (2, 3), (1, 2)),
    ((2, 3), (3, 2), 3),
    ((4, 2, 3), (4, 6), 4),
    ((3, 4, 5, 6), (3, 4, 5, 6), (2, 3, 4, 5)),
    ((), (1,), 1),
    ((1,), (), 1),
    ((24,), (3, 8), 24),
    ((24,), (4, 6), 6),
    ((24,), (4, 3, 2), 6),
    ((24,), (4, 6, 1), 6),
    ((24,), (4, 6), (6, 12, 6)),
    ((64, 4), (8, 8, 4), (16, 2)),
    ((4, 64), (4, 8, 4, 2), (2, 16)),
    ((4, 8, 4, 2), (2, 1, 2, 32, 2), (2, 4, 2, 2)),
    ((4, 1, 4), (4, 4), (2, 1, 2)),
    ((0, 10), (0, 5, 2), (5, 5)),
    ((5, 0, 2), (0, 10), (5, 2, 2)),
    ((0,), (2, 0, 2), (4,)),
    ((2, 0, 2), (0,), (4, 4, 4)),
])
def test_reshape(original_shape, new_shape, chunks):
    x = np.random.randint(10, size=original_shape)
    a = from_array(x, chunks=chunks)

    xr = x.reshape(new_shape)
    ar = a.reshape(new_shape)

    if a.shape == new_shape:
        assert a is ar

    assert_eq(xr, ar)


def test_reshape_exceptions():
    x = np.random.randint(10, size=(5,))
    a = from_array(x, chunks=(2,))
    with pytest.raises(ValueError):
        da.reshape(a, (100,))


def test_reshape_splat():
    x = da.ones((5, 5), chunks=(2, 2))
    assert_eq(x.reshape((25,)), x.reshape(25))


def test_reshape_fails_for_dask_only():
    cases = [
        ((3, 4), (4, 3), 2),
    ]
    for original_shape, new_shape, chunks in cases:
        x = np.random.randint(10, size=original_shape)
        a = from_array(x, chunks=chunks)
        assert x.reshape(new_shape).shape == new_shape
        with pytest.raises(ValueError):
            da.reshape(a, new_shape)


def test_reshape_unknown_dimensions():
    for original_shape in [(24,), (2, 12), (2, 3, 4)]:
        for new_shape in [(-1,), (2, -1), (-1, 3, 4)]:
            x = np.random.randint(10, size=original_shape)
            a = from_array(x, 24)
            assert_eq(x.reshape(new_shape), a.reshape(new_shape))

    pytest.raises(ValueError, lambda: da.reshape(a, (-1, -1)))


def test_full():
    d = da.full((3, 4), 2, chunks=((2, 1), (2, 2)))
    assert d.chunks == ((2, 1), (2, 2))
    assert_eq(d, np.full((3, 4), 2))


def test_map_blocks():
    x = np.arange(400).reshape((20, 20))
    d = from_array(x, chunks=(7, 7))

    e = d.map_blocks(inc, dtype=d.dtype)

    assert d.chunks == e.chunks
    assert_eq(e, x + 1)

    e = d.map_blocks(inc, name='increment')
    assert e.name == 'increment'

    e = d.map_blocks(inc, token='increment')
    assert e.name != 'increment'
    assert e.name.startswith('increment')

    d = from_array(x, chunks=(10, 10))
    e = d.map_blocks(lambda x: x[::2, ::2], chunks=(5, 5), dtype=d.dtype)

    assert e.chunks == ((5, 5), (5, 5))
    assert_eq(e, x[::2, ::2])

    d = from_array(x, chunks=(8, 8))
    e = d.map_blocks(lambda x: x[::2, ::2], chunks=((4, 4, 2), (4, 4, 2)),
                     dtype=d.dtype)

    assert_eq(e, x[::2, ::2])


def test_map_blocks2():
    x = np.arange(10, dtype='i8')
    d = from_array(x, chunks=(2,))

    def func(block, block_id=None, c=0):
        return np.ones_like(block) * sum(block_id) + c

    out = d.map_blocks(func, dtype='i8')
    expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4], dtype='i8')

    assert_eq(out, expected)
    assert same_keys(d.map_blocks(func, dtype='i8'), out)

    out = d.map_blocks(func, dtype='i8', c=1)
    expected = expected + 1

    assert_eq(out, expected)
    assert same_keys(d.map_blocks(func, dtype='i8', c=1), out)


def test_map_blocks_with_constants():
    d = da.arange(10, chunks=3)
    e = d.map_blocks(add, 100, dtype=d.dtype)

    assert_eq(e, np.arange(10) + 100)

    assert_eq(da.map_blocks(sub, d, 10, dtype=d.dtype),
              np.arange(10) - 10)
    assert_eq(da.map_blocks(sub, 10, d, dtype=d.dtype),
              10 - np.arange(10))


def test_map_blocks_with_kwargs():
    d = da.arange(10, chunks=5)

    result = d.map_blocks(np.max, axis=0, keepdims=True, dtype=d.dtype,
                          chunks=(1,))

    assert_eq(result, np.array([4, 9]))


def test_map_blocks_with_chunks():
    dx = da.ones((5, 3), chunks=(2, 2))
    dy = da.ones((5, 3), chunks=(2, 2))
    dz = da.map_blocks(np.add, dx, dy, chunks=dx.chunks)
    assert_eq(dz, np.ones((5, 3)) * 2)


def test_map_blocks_dtype_inference():
    x = np.arange(50).reshape((5, 10))
    y = np.arange(10)
    dx = da.from_array(x, chunks=5)
    dy = da.from_array(y, chunks=5)

    def foo(x, *args, **kwargs):
        cast = kwargs.pop('cast', 'i8')
        return (x + sum(args)).astype(cast)

    assert_eq(dx.map_blocks(foo, dy, 1), foo(dx, dy, 1))
    assert_eq(dx.map_blocks(foo, dy, 1, cast='f8'), foo(dx, dy, 1, cast='f8'))
    assert_eq(dx.map_blocks(foo, dy, 1, cast='f8', dtype='f8'),
              foo(dx, dy, 1, cast='f8', dtype='f8'))

    def foo(x):
        raise RuntimeError("Woops")

    try:
        dx.map_blocks(foo)
    except Exception as e:
        assert e.args[0].startswith("`dtype` inference failed")
        assert "Please specify the dtype explicitly" in e.args[0]
        assert 'RuntimeError' in e.args[0]
    else:
        assert False, "Should have errored"


def test_from_function_requires_block_args():
    x = np.arange(10)
    pytest.raises(Exception, lambda: from_array(x))


def test_repr():
    d = da.ones((4, 4), chunks=(2, 2))
    assert d.name[:5] in repr(d)
    assert str(d.shape) in repr(d)
    assert str(d.dtype) in repr(d)
    d = da.ones((4000, 4), chunks=(4, 2))
    assert len(str(d)) < 1000


def test_slicing_with_ellipsis():
    x = np.arange(256).reshape((4, 4, 4, 4))
    d = da.from_array(x, chunks=((2, 2, 2, 2)))

    assert_eq(d[..., 1], x[..., 1])
    assert_eq(d[0, ..., 1], x[0, ..., 1])


def test_slicing_with_ndarray():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=((4, 4)))

    assert_eq(d[np.arange(8)], x)
    assert_eq(d[np.ones(8, dtype=bool)], x)
    assert_eq(d[np.array([1])], x[[1]])
    assert_eq(d[np.array([True])], x[[0]])


def test_dtype():
    d = da.ones((4, 4), chunks=(2, 2))

    assert d.dtype == d.compute().dtype
    assert (d * 1.0).dtype == (d + 1.0).compute().dtype
    assert d.sum().dtype == d.sum().compute().dtype  # no shape


def test_blockdims_from_blockshape():
    assert blockdims_from_blockshape((10, 10), (4, 3)) == ((4, 4, 2), (3, 3, 3, 1))
    pytest.raises(TypeError, lambda: blockdims_from_blockshape((10,), None))
    assert blockdims_from_blockshape((1e2, 3), [1e1, 3]) == ((10, ) * 10, (3, ))
    assert blockdims_from_blockshape((np.int8(10), ), (5, )) == ((5, 5), )


def test_coerce():
    d0 = da.from_array(np.array(1), chunks=(1,))
    d1 = da.from_array(np.array([1]), chunks=(1,))
    with dask.set_options(get=dask.get):
        for d in d0, d1:
            assert bool(d) is True
            assert int(d) == 1
            assert float(d) == 1.0
            assert complex(d) == complex(1)

    a2 = np.arange(2)
    d2 = da.from_array(a2, chunks=(2,))
    for func in (int, float, complex):
        pytest.raises(TypeError, lambda :func(d2))


def test_bool():
    arr = np.arange(100).reshape((10,10))
    darr = da.from_array(arr, chunks=(10,10))
    with pytest.raises(ValueError):
        bool(darr)
        bool(darr == darr)


def test_store_kwargs():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    called = [False]

    def get_func(*args, **kwargs):
        assert kwargs.pop("foo") == "test kwarg"
        r = dask.get(*args, **kwargs)
        called[0] = True
        return r

    called[0] = False
    at = np.zeros(shape=(10, 10))
    store([a], [at], get=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = np.zeros(shape=(10, 10))
    a.store(at, get=get_func, foo="test kwarg")
    assert called[0]

    called[0] = False
    at = np.zeros(shape=(10, 10))
    store([a], [at], get=get_func, return_store=True, foo="test kwarg")
    assert called[0]


def test_store_delayed_target():
    from dask.delayed import delayed
    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    # empty buffers to be used as targets
    targs = {}

    def make_target(key):
        a = np.empty((4, 4))
        targs[key] = a
        return a

    # delayed calls to these targets
    atd = delayed(make_target)('at')
    btd = delayed(make_target)('bt')

    # test not keeping result
    st = store([a, b], [atd, btd])

    at = targs['at']
    bt = targs['bt']

    assert st is None
    assert_eq(at, a)
    assert_eq(bt, b)

    # test keeping result
    for st_compute in [False, True]:
        targs.clear()

        st = store([a, b], [atd, btd], return_stored=True, compute=st_compute)
        if st_compute:
            assert all(
                not any(dask.core.get_deps(e.dask)[0].values()) for e in st
            )

        st = dask.compute(*st)

        at = targs['at']
        bt = targs['bt']

        assert st is not None
        assert isinstance(st, tuple)
        assert all([isinstance(v, np.ndarray) for v in st])
        assert_eq(at, a)
        assert_eq(bt, b)
        assert_eq(st[0], a)
        assert_eq(st[1], b)

        pytest.raises(ValueError, lambda: store([a], [at, bt]))
        pytest.raises(ValueError, lambda: store(at, at))
        pytest.raises(ValueError, lambda: store([at, bt], [at, bt]))


def test_store():
    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.empty(shape=(4, 4))
    bt = np.empty(shape=(4, 4))

    st = store([a, b], [at, bt])
    assert st is None
    assert (at == 2).all()
    assert (bt == 3).all()

    pytest.raises(ValueError, lambda: store([a], [at, bt]))
    pytest.raises(ValueError, lambda: store(at, at))
    pytest.raises(ValueError, lambda: store([at, bt], [at, bt]))


def test_store_regions():
    d = da.ones((4, 4, 4), dtype=int, chunks=(2, 2, 2))
    a, b = d + 1, d + 2
    a = a[:, 1:, :].astype(float)

    region = (slice(None,None,2), slice(None), [1, 2, 4, 5])

    # Single region:
    at = np.zeros(shape=(8, 3, 6))
    bt = np.zeros(shape=(8, 4, 6))
    v = store([a, b], [at, bt], regions=region, compute=False)
    assert isinstance(v, Delayed)
    assert (at == 0).all() and (bt[region] == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at[region] == 2).all() and (bt[region] == 3).all()
    assert not (bt == 3).all() and not ( bt == 0 ).all()
    assert not (at == 2).all() and not ( at == 0 ).all()

    # Multiple regions:
    at = np.zeros(shape=(8, 3, 6))
    bt = np.zeros(shape=(8, 4, 6))
    v = store([a, b], [at, bt], regions=[region, region], compute=False)
    assert isinstance(v, Delayed)
    assert (at == 0).all() and (bt[region] == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at[region] == 2).all() and (bt[region] == 3).all()
    assert not (bt == 3).all() and not ( bt == 0 ).all()
    assert not (at == 2).all() and not ( at == 0 ).all()

    # Single region (keep result):
    for st_compute in [False, True]:
        at = np.zeros(shape=(8, 3, 6))
        bt = np.zeros(shape=(8, 4, 6))
        v = store(
            [a, b], [at, bt], regions=region,
            compute=st_compute, return_stored=True
        )
        assert isinstance(v, tuple)
        assert all([isinstance(e, da.Array) for e in v])
        if st_compute:
            assert all(
                not any(dask.core.get_deps(e.dask)[0].values()) for e in v
            )
        else:
            assert (at == 0).all() and (bt[region] == 0).all()

        ar, br = v
        assert ar.dtype == a.dtype
        assert br.dtype == b.dtype
        assert ar.shape == a.shape
        assert br.shape == b.shape
        assert ar.chunks == a.chunks
        assert br.chunks == b.chunks

        ar, br = da.compute(ar, br)
        assert (at[region] == 2).all() and (bt[region] == 3).all()
        assert not (bt == 3).all() and not ( bt == 0 ).all()
        assert not (at == 2).all() and not ( at == 0 ).all()
        assert (br == 3).all()
        assert (ar == 2).all()

    # Multiple regions (keep result):
    for st_compute in [False, True]:
        at = np.zeros(shape=(8, 3, 6))
        bt = np.zeros(shape=(8, 4, 6))
        v = store(
            [a, b], [at, bt], regions=[region, region],
            compute=st_compute, return_stored=True
        )
        assert isinstance(v, tuple)
        assert all([isinstance(e, da.Array) for e in v])
        if st_compute:
            assert all(
                not any(dask.core.get_deps(e.dask)[0].values()) for e in v
            )
        else:
            assert (at == 0).all() and (bt[region] == 0).all()

        ar, br = v
        assert ar.dtype == a.dtype
        assert br.dtype == b.dtype
        assert ar.shape == a.shape
        assert br.shape == b.shape
        assert ar.chunks == a.chunks
        assert br.chunks == b.chunks

        ar, br = da.compute(ar, br)
        assert (at[region] == 2).all() and (bt[region] == 3).all()
        assert not (bt == 3).all() and not ( bt == 0 ).all()
        assert not (at == 2).all() and not ( at == 0 ).all()
        assert (br == 3).all()
        assert (ar == 2).all()


def test_store_compute_false():
    d = da.ones((4, 4), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.zeros(shape=(4, 4))
    bt = np.zeros(shape=(4, 4))

    v = store([a, b], [at, bt], compute=False)
    assert isinstance(v, Delayed)
    assert (at == 0).all() and (bt == 0).all()
    assert all([ev is None for ev in v.compute()])
    assert (at == 2).all() and (bt == 3).all()

    at = np.zeros(shape=(4, 4))
    bt = np.zeros(shape=(4, 4))

    dat, dbt = store([a, b], [at, bt], compute=False, return_stored=True)
    assert isinstance(dat, Array) and isinstance(dbt, Array)
    assert (at == 0).all() and (bt == 0).all()
    assert (dat.compute() == at).all() and (dbt.compute() == bt).all()
    assert (at == 2).all() and (bt == 3).all()


class ThreadSafetyError(Exception):
    pass


class NonthreadSafeStore(object):
    def __init__(self):
        self.in_use = False

    def __setitem__(self, key, value):
        if self.in_use:
            raise ThreadSafetyError()
        self.in_use = True
        time.sleep(0.001)
        self.in_use = False


class ThreadSafeStore(object):
    def __init__(self):
        self.concurrent_uses = 0
        self.max_concurrent_uses = 0

    def __setitem__(self, key, value):
        self.concurrent_uses += 1
        self.max_concurrent_uses = max(self.concurrent_uses, self.max_concurrent_uses)
        time.sleep(0.01)
        self.concurrent_uses -= 1


class CounterLock(object):
    def __init__(self, *args, **kwargs):
        self.lock = Lock(*args, **kwargs)

        self.acquire_count = 0
        self.release_count = 0

    def acquire(self, *args, **kwargs):
        self.acquire_count += 1
        return self.lock.acquire(*args, **kwargs)

    def release(self, *args, **kwargs):
        self.release_count += 1
        return self.lock.release(*args, **kwargs)


def test_store_locks():
    _Lock = type(Lock())
    d = da.ones((10, 10), chunks=(2, 2))
    a, b = d + 1, d + 2

    at = np.zeros(shape=(10, 10))
    bt = np.zeros(shape=(10, 10))

    lock = Lock()
    v = store([a, b], [at, bt], compute=False, lock=lock)
    assert isinstance(v, Delayed)
    dsk = v.dask
    locks = set(vv for v in dsk.values() for vv in v if isinstance(vv, _Lock))
    assert locks == set([lock])

    # Ensure same lock applies over multiple stores
    at = NonthreadSafeStore()
    v = store([a, b], [at, at], lock=lock,
              get=dask.threaded.get, num_workers=10)
    assert v is None

    # Don't assume thread safety by default
    at = NonthreadSafeStore()
    assert store(a, at, get=dask.threaded.get, num_workers=10) is None
    assert a.store(at, get=dask.threaded.get, num_workers=10) is None

    # Ensure locks can be removed
    at = ThreadSafeStore()
    for i in range(10):
        st = a.store(at, lock=False, get=dask.threaded.get, num_workers=10)
        assert st is None
        if at.max_concurrent_uses > 1:
            break
        if i == 9:
            assert False

    # Verify number of lock calls
    nchunks = np.sum([np.prod([len(c) for c in e.chunks]) for e in [a, b]])
    for c in (False, True):
        at = np.zeros(shape=(10, 10))
        bt = np.zeros(shape=(10, 10))
        lock = CounterLock()

        v = store([a, b], [at, bt], lock=lock, compute=c, return_stored=True)
        assert all(isinstance(e, Array) for e in v)

        da.compute(v)

        # When `return_stored=True` and `compute=False`,
        # the lock should be acquired only once for store and load steps
        # as they are fused together into one step.
        assert lock.acquire_count == lock.release_count
        if c:
            assert lock.acquire_count == 2 * nchunks
        else:
            assert lock.acquire_count == nchunks


def test_store_method_return():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    for compute in [False, True]:
        for return_stored in [False, True]:
            at = np.zeros(shape=(10, 10))
            r = a.store(
                at, get=dask.threaded.get,
                compute=compute, return_stored=return_stored
            )

            if return_stored:
                assert isinstance(r, Array)
            elif compute:
                assert r is None
            else:
                assert isinstance(r, Delayed)


@pytest.mark.xfail(reason="can't lock with multiprocessing")
def test_store_multiprocessing_lock():
    d = da.ones((10, 10), chunks=(2, 2))
    a = d + 1

    at = np.zeros(shape=(10, 10))
    st = a.store(at, get=dask.multiprocessing.get, num_workers=10)
    assert st is None


def test_to_hdf5():
    h5py = pytest.importorskip('h5py')
    x = da.ones((4, 4), chunks=(2, 2))
    y = da.ones(4, chunks=2, dtype='i4')

    with tmpfile('.hdf5') as fn:
        x.to_hdf5(fn, '/x')
        with h5py.File(fn) as f:
            d = f['/x']

            assert_eq(d[:], x)
            assert d.chunks == (2, 2)

    with tmpfile('.hdf5') as fn:
        x.to_hdf5(fn, '/x', chunks=None)
        with h5py.File(fn) as f:
            d = f['/x']

            assert_eq(d[:], x)
            assert d.chunks is None

    with tmpfile('.hdf5') as fn:
        x.to_hdf5(fn, '/x', chunks=(1, 1))
        with h5py.File(fn) as f:
            d = f['/x']

            assert_eq(d[:], x)
            assert d.chunks == (1, 1)

    with tmpfile('.hdf5') as fn:
        da.to_hdf5(fn, {'/x': x, '/y': y})

        with h5py.File(fn) as f:
            assert_eq(f['/x'][:], x)
            assert f['/x'].chunks == (2, 2)
            assert_eq(f['/y'][:], y)
            assert f['/y'].chunks == (2,)


def test_to_dask_dataframe():
    dd = pytest.importorskip('dask.dataframe')
    a = da.ones((4,), chunks=(2,))
    d = a.to_dask_dataframe()
    assert isinstance(d, dd.Series)

    a = da.ones((4, 4), chunks=(2, 2))
    d = a.to_dask_dataframe()
    assert isinstance(d, dd.DataFrame)


def test_np_array_with_zero_dimensions():
    d = da.ones((4, 4), chunks=(2, 2))
    assert_eq(np.array(d.sum()), np.array(d.compute().sum()))


def test_dtype_complex():
    x = np.arange(24).reshape((4, 6)).astype('f4')
    y = np.arange(24).reshape((4, 6)).astype('i8')
    z = np.arange(24).reshape((4, 6)).astype('i2')

    a = da.from_array(x, chunks=(2, 3))
    b = da.from_array(y, chunks=(2, 3))
    c = da.from_array(z, chunks=(2, 3))

    def assert_eq(a, b):
        return (isinstance(a, np.dtype) and
                isinstance(b, np.dtype) and
                str(a) == str(b))

    assert_eq(a.dtype, x.dtype)
    assert_eq(b.dtype, y.dtype)

    assert_eq((a + 1).dtype, (x + 1).dtype)
    assert_eq((a + b).dtype, (x + y).dtype)
    assert_eq(a.T.dtype, x.T.dtype)
    assert_eq(a[:3].dtype, x[:3].dtype)
    assert_eq((a.dot(b.T)).dtype, (x.dot(y.T)).dtype)

    assert_eq(stack([a, b]).dtype, np.vstack([x, y]).dtype)
    assert_eq(concatenate([a, b]).dtype, np.concatenate([x, y]).dtype)

    assert_eq(b.std().dtype, y.std().dtype)
    assert_eq(c.sum().dtype, z.sum().dtype)
    assert_eq(a.min().dtype, a.min().dtype)
    assert_eq(b.std().dtype, b.std().dtype)
    assert_eq(a.argmin(axis=0).dtype, a.argmin(axis=0).dtype)

    assert_eq(da.sin(c).dtype, np.sin(z).dtype)
    assert_eq(da.exp(b).dtype, np.exp(y).dtype)
    assert_eq(da.floor(a).dtype, np.floor(x).dtype)
    assert_eq(da.isnan(b).dtype, np.isnan(y).dtype)
    with ignoring(ImportError):
        assert da.isnull(b).dtype == 'bool'
        assert da.notnull(b).dtype == 'bool'

    x = np.array([('a', 1)], dtype=[('text', 'S1'), ('numbers', 'i4')])
    d = da.from_array(x, chunks=(1,))

    assert_eq(d['text'].dtype, x['text'].dtype)
    assert_eq(d[['numbers', 'text']].dtype, x[['numbers', 'text']].dtype)


def test_astype():
    x = np.ones((5, 5), dtype='f8')
    d = da.from_array(x, chunks=(2,2))

    assert d.astype('i8').dtype == 'i8'
    assert_eq(d.astype('i8'), x.astype('i8'))
    assert same_keys(d.astype('i8'), d.astype('i8'))

    with pytest.raises(TypeError):
        d.astype('i8', casting='safe')

    with pytest.raises(TypeError):
        d.astype('i8', not_a_real_kwarg='foo')

    # smoketest with kwargs
    assert_eq(d.astype('i8', copy=False), x.astype('i8', copy=False))

    # Check it's a noop
    assert d.astype('f8') is d


def test_arithmetic():
    x = np.arange(5).astype('f4') + 2
    y = np.arange(5).astype('i8') + 2
    z = np.arange(5).astype('i4') + 2
    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))
    c = da.from_array(z, chunks=(2,))
    assert_eq(a + b, x + y)
    assert_eq(a * b, x * y)
    assert_eq(a - b, x - y)
    assert_eq(a / b, x / y)
    assert_eq(b & b, y & y)
    assert_eq(b | b, y | y)
    assert_eq(b ^ b, y ^ y)
    assert_eq(a // b, x // y)
    assert_eq(a ** b, x ** y)
    assert_eq(a % b, x % y)
    assert_eq(a > b, x > y)
    assert_eq(a < b, x < y)
    assert_eq(a >= b, x >= y)
    assert_eq(a <= b, x <= y)
    assert_eq(a == b, x == y)
    assert_eq(a != b, x != y)

    assert_eq(a + 2, x + 2)
    assert_eq(a * 2, x * 2)
    assert_eq(a - 2, x - 2)
    assert_eq(a / 2, x / 2)
    assert_eq(b & True, y & True)
    assert_eq(b | True, y | True)
    assert_eq(b ^ True, y ^ True)
    assert_eq(a // 2, x // 2)
    assert_eq(a ** 2, x ** 2)
    assert_eq(a % 2, x % 2)
    assert_eq(a > 2, x > 2)
    assert_eq(a < 2, x < 2)
    assert_eq(a >= 2, x >= 2)
    assert_eq(a <= 2, x <= 2)
    assert_eq(a == 2, x == 2)
    assert_eq(a != 2, x != 2)

    assert_eq(2 + b, 2 + y)
    assert_eq(2 * b, 2 * y)
    assert_eq(2 - b, 2 - y)
    assert_eq(2 / b, 2 / y)
    assert_eq(True & b, True & y)
    assert_eq(True | b, True | y)
    assert_eq(True ^ b, True ^ y)
    assert_eq(2 // b, 2 // y)
    assert_eq(2 ** b, 2 ** y)
    assert_eq(2 % b, 2 % y)
    assert_eq(2 > b, 2 > y)
    assert_eq(2 < b, 2 < y)
    assert_eq(2 >= b, 2 >= y)
    assert_eq(2 <= b, 2 <= y)
    assert_eq(2 == b, 2 == y)
    assert_eq(2 != b, 2 != y)

    assert_eq(-a, -x)
    assert_eq(abs(a), abs(x))
    assert_eq(~(a == b), ~(x == y))
    assert_eq(~(a == b), ~(x == y))

    assert_eq(da.logaddexp(a, b), np.logaddexp(x, y))
    assert_eq(da.logaddexp2(a, b), np.logaddexp2(x, y))
    with pytest.warns(None):  # Overflow warning
        assert_eq(da.exp(b), np.exp(y))
    assert_eq(da.log(a), np.log(x))
    assert_eq(da.log10(a), np.log10(x))
    assert_eq(da.log1p(a), np.log1p(x))
    with pytest.warns(None):  # Overflow warning
        assert_eq(da.expm1(b), np.expm1(y))
    assert_eq(da.sqrt(a), np.sqrt(x))
    assert_eq(da.square(a), np.square(x))

    assert_eq(da.sin(a), np.sin(x))
    assert_eq(da.cos(b), np.cos(y))
    assert_eq(da.tan(a), np.tan(x))
    assert_eq(da.arcsin(b / 10), np.arcsin(y / 10))
    assert_eq(da.arccos(b / 10), np.arccos(y / 10))
    assert_eq(da.arctan(b / 10), np.arctan(y / 10))
    assert_eq(da.arctan2(b * 10, a), np.arctan2(y * 10, x))
    assert_eq(da.hypot(b, a), np.hypot(y, x))
    assert_eq(da.sinh(a), np.sinh(x))
    with pytest.warns(None):  # Overflow warning
        assert_eq(da.cosh(b), np.cosh(y))
    assert_eq(da.tanh(a), np.tanh(x))
    assert_eq(da.arcsinh(b * 10), np.arcsinh(y * 10))
    assert_eq(da.arccosh(b * 10), np.arccosh(y * 10))
    assert_eq(da.arctanh(b / 10), np.arctanh(y / 10))
    assert_eq(da.deg2rad(a), np.deg2rad(x))
    assert_eq(da.rad2deg(a), np.rad2deg(x))

    assert_eq(da.logical_and(a < 1, b < 4), np.logical_and(x < 1, y < 4))
    assert_eq(da.logical_or(a < 1, b < 4), np.logical_or(x < 1, y < 4))
    assert_eq(da.logical_xor(a < 1, b < 4), np.logical_xor(x < 1, y < 4))
    assert_eq(da.logical_not(a < 1), np.logical_not(x < 1))
    assert_eq(da.maximum(a, 5 - a), np.maximum(a, 5 - a))
    assert_eq(da.minimum(a, 5 - a), np.minimum(a, 5 - a))
    assert_eq(da.fmax(a, 5 - a), np.fmax(a, 5 - a))
    assert_eq(da.fmin(a, 5 - a), np.fmin(a, 5 - a))

    assert_eq(da.isreal(a + 1j * b), np.isreal(x + 1j * y))
    assert_eq(da.iscomplex(a + 1j * b), np.iscomplex(x + 1j * y))
    assert_eq(da.isfinite(a), np.isfinite(x))
    assert_eq(da.isinf(a), np.isinf(x))
    assert_eq(da.isnan(a), np.isnan(x))
    assert_eq(da.signbit(a - 3), np.signbit(x - 3))
    assert_eq(da.copysign(a - 3, b), np.copysign(x - 3, y))
    assert_eq(da.nextafter(a - 3, b), np.nextafter(x - 3, y))
    with pytest.warns(None):  # overflow warning
        assert_eq(da.ldexp(c, c), np.ldexp(z, z))
    assert_eq(da.fmod(a * 12, b), np.fmod(x * 12, y))
    assert_eq(da.floor(a * 0.5), np.floor(x * 0.5))
    assert_eq(da.ceil(a), np.ceil(x))
    assert_eq(da.trunc(a / 2), np.trunc(x / 2))

    assert_eq(da.degrees(b), np.degrees(y))
    assert_eq(da.radians(a), np.radians(x))

    assert_eq(da.rint(a + 0.3), np.rint(x + 0.3))
    assert_eq(da.fix(a - 2.5), np.fix(x - 2.5))

    assert_eq(da.angle(a + 1j), np.angle(x + 1j))
    assert_eq(da.real(a + 1j), np.real(x + 1j))
    assert_eq((a + 1j).real, np.real(x + 1j))
    assert_eq(da.imag(a + 1j), np.imag(x + 1j))
    assert_eq((a + 1j).imag, np.imag(x + 1j))
    assert_eq(da.conj(a + 1j * b), np.conj(x + 1j * y))
    assert_eq((a + 1j * b).conj(), (x + 1j * y).conj())

    assert_eq(da.clip(b, 1, 4), np.clip(y, 1, 4))
    assert_eq(b.clip(1, 4), y.clip(1, 4))
    assert_eq(da.fabs(b), np.fabs(y))
    assert_eq(da.sign(b - 2), np.sign(y - 2))
    assert_eq(da.absolute(b - 2), np.absolute(y - 2))
    assert_eq(da.absolute(b - 2 + 1j), np.absolute(y - 2 + 1j))

    l1, l2 = da.frexp(a)
    r1, r2 = np.frexp(x)
    assert_eq(l1, r1)
    assert_eq(l2, r2)

    l1, l2 = da.modf(a)
    r1, r2 = np.modf(x)
    assert_eq(l1, r1)
    assert_eq(l2, r2)

    assert_eq(da.around(a, -1), np.around(x, -1))


def test_elemwise_consistent_names():
    a = da.from_array(np.arange(5, dtype='f4'), chunks=(2,))
    b = da.from_array(np.arange(5, dtype='f4'), chunks=(2,))
    assert same_keys(a + b, a + b)
    assert same_keys(a + 2, a + 2)
    assert same_keys(da.exp(a), da.exp(a))
    assert same_keys(da.exp(a, dtype='f8'), da.exp(a, dtype='f8'))
    assert same_keys(da.maximum(a, b), da.maximum(a, b))


def test_optimize():
    x = np.arange(5).astype('f4')
    a = da.from_array(x, chunks=(2,))
    expr = a[1:4] + 1
    result = optimize(expr.dask, expr.__dask_keys__())
    assert isinstance(result, dict)
    assert all(key in result for key in expr.__dask_keys__())


def test_slicing_with_non_ndarrays():
    class ARangeSlice(object):
        def __init__(self, start, stop):
            self.start = start
            self.stop = stop

        def __array__(self):
            return np.arange(self.start, self.stop)

    class ARangeSlicable(object):
        dtype = 'i8'

        def __init__(self, n):
            self.n = n

        @property
        def shape(self):
            return (self.n,)

        def __getitem__(self, key):
            return ARangeSlice(key[0].start, key[0].stop)

    x = da.from_array(ARangeSlicable(10), chunks=(4,))

    assert_eq((x + 1).sum(), (np.arange(10, dtype=x.dtype) + 1).sum())


def test_getter():
    assert type(getter(np.matrix([[1]]), 0)) is np.ndarray
    assert type(getter(np.matrix([[1]]), 0, asarray=False)) is np.matrix
    assert_eq(getter([1, 2, 3, 4, 5], slice(1, 4)), np.array([2, 3, 4]))

    assert_eq(getter(np.arange(5), (None, slice(None, None))),
              np.arange(5)[None, :])


def test_size():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.size == np.array(x).size
    assert isinstance(x.size, int)


def test_nbytes():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.nbytes == np.array(x).nbytes


def test_itemsize():
    x = da.ones((10, 2), chunks=(3, 1))
    assert x.itemsize == 8


def test_Array_normalizes_dtype():
    x = da.ones((3,), chunks=(1,), dtype=int)
    assert isinstance(x.dtype, np.dtype)


def test_from_array_with_lock():
    x = np.arange(10)
    d = da.from_array(x, chunks=5, lock=True)

    tasks = [v for k, v in d.dask.items() if k[0] == d.name]

    assert hasattr(tasks[0][4], 'acquire')
    assert len(set(task[4] for task in tasks)) == 1

    assert_eq(d, x)

    lock = Lock()
    e = da.from_array(x, chunks=5, lock=lock)
    f = da.from_array(x, chunks=5, lock=lock)

    assert_eq(e + f, x + x)


class MyArray(object):
    def __init__(self, x):
        self.x = x
        self.dtype = x.dtype
        self.shape = x.shape

    def __getitem__(self, i):
        return self.x[i]


def test_from_array_tasks_always_call_getter():
    x1 = np.arange(25).reshape((5, 5))
    x2 = np.array([[1]])
    x3 = np.array(1)[()]

    dx1a = da.from_array(MyArray(x1), chunks=(5, 5), asarray=False)
    dx1b = da.from_array(MyArray(x1), chunks=x1.shape, asarray=False)
    dx2 = da.from_array(MyArray(x2), chunks=1, asarray=False)
    dx3 = da.from_array(MyArray(x3), chunks=1, asarray=False)

    for res, sol in [(dx1a, x1), (dx1b, x1), (dx2, x2), (dx3, x3)]:
        assert_eq(res, sol)


def test_from_array_no_asarray():

    def assert_chunks_are_of_type(x, cls):
        chunks = compute_as_if_collection(Array, x.dask, x.__dask_keys__())
        for c in concat(chunks):
            assert type(c) is cls

    x = np.matrix(np.arange(100).reshape((10, 10)))

    for asarray, cls in [(True, np.ndarray), (False, np.matrix)]:
        dx = da.from_array(x, chunks=(5, 5), asarray=asarray)
        assert_chunks_are_of_type(dx, cls)
        assert_chunks_are_of_type(dx[0:5], cls)
        assert_chunks_are_of_type(dx[0:5][:, 0], cls)


def test_from_array_getitem():
    x = np.arange(10)

    def my_getitem(x, ind):
        return x[ind]

    y = da.from_array(x, chunks=(5,), getitem=my_getitem)

    for k, v in y.dask.items():
        if isinstance(v, tuple):
            assert v[0] is my_getitem

    assert_eq(x, y)


def test_from_array_minus_one():
    x = np.arange(10)
    y = da.from_array(x, -1)
    assert y.chunks == ((10,),)
    assert_eq(x, y)


def test_asarray():
    assert_eq(da.asarray([1, 2, 3]), np.asarray([1, 2, 3]))

    x = da.asarray([1, 2, 3])
    assert da.asarray(x) is x


def test_asarray_h5py():
    h5py = pytest.importorskip('h5py')

    with tmpfile('.hdf5') as fn:
        with h5py.File(fn) as f:
            d = f.create_dataset('/x', shape=(2, 2), dtype=float)
            x = da.asarray(d)
            assert d in x.dask.values()
            assert not any(isinstance(v, np.ndarray) for v in x.dask.values())


def test_asanyarray():
    x = np.matrix([1, 2, 3])
    dx = da.asanyarray(x)
    assert dx.numblocks == (1, 1)
    chunks = compute_as_if_collection(Array, dx.dask, dx.__dask_keys__())
    assert isinstance(chunks[0][0], np.matrix)
    assert da.asanyarray(dx) is dx


def test_from_func():
    x = np.arange(10)
    f = lambda n: n * x
    d = from_func(f, (10,), x.dtype, kwargs={'n': 2})

    assert d.shape == x.shape
    assert d.dtype == x.dtype
    assert_eq(d.compute(), 2 * x)
    assert same_keys(d, from_func(f, (10,), x.dtype, kwargs={'n': 2}))


def test_concatenate3_2():
    x = np.array([1, 2])
    assert_eq(concatenate3([x, x, x]),
              np.array([1, 2, 1, 2, 1, 2]))

    x = np.array([[1, 2]])
    assert (concatenate3([[x, x, x], [x, x, x]]) ==
            np.array([[1, 2, 1, 2, 1, 2],
                      [1, 2, 1, 2, 1, 2]])).all()

    assert (concatenate3([[x, x], [x, x], [x, x]]) ==
            np.array([[1, 2, 1, 2],
                      [1, 2, 1, 2],
                      [1, 2, 1, 2]])).all()

    x = np.arange(12).reshape((2, 2, 3))
    assert_eq(concatenate3([[[x, x, x], [x, x, x]],
                           [[x, x, x], [x, x, x]]]),
              np.array([[[ 0,  1,  2,  0,  1,  2,  0,  1,  2],
                         [ 3,  4,  5,  3,  4,  5,  3,  4,  5],
                         [ 0,  1,  2,  0,  1,  2,  0,  1,  2],
                         [ 3,  4,  5,  3,  4,  5,  3,  4,  5]],

                        [[ 6,  7,  8,  6,  7,  8,  6,  7,  8],
                         [ 9, 10, 11,  9, 10, 11,  9, 10, 11],
                         [ 6,  7,  8,  6,  7,  8,  6,  7,  8],
                         [ 9, 10, 11,  9, 10, 11,  9, 10, 11]],

                        [[ 0,  1,  2,  0,  1,  2,  0,  1,  2],
                         [ 3,  4,  5,  3,  4,  5,  3,  4,  5],
                         [ 0,  1,  2,  0,  1,  2,  0,  1,  2],
                         [ 3,  4,  5,  3,  4,  5,  3,  4,  5]],

                        [[ 6,  7,  8,  6,  7,  8,  6,  7,  8],
                         [ 9, 10, 11,  9, 10, 11,  9, 10, 11],
                         [ 6,  7,  8,  6,  7,  8,  6,  7,  8],
                         [ 9, 10, 11,  9, 10, 11,  9, 10, 11]]]))


def test_map_blocks3():
    x = np.arange(10)
    y = np.arange(10) * 2

    d = da.from_array(x, chunks=5)
    e = da.from_array(y, chunks=5)

    assert_eq(da.core.map_blocks(lambda a, b: a + 2 * b, d, e, dtype=d.dtype),
              x + 2 * y)

    z = np.arange(100).reshape((10, 10))
    f = da.from_array(z, chunks=5)

    func = lambda a, b: a + 2 * b
    res = da.core.map_blocks(func, d, f, dtype=d.dtype)
    assert_eq(res, x + 2 * z)
    assert same_keys(da.core.map_blocks(func, d, f, dtype=d.dtype), res)

    assert_eq(da.map_blocks(func, f, d, dtype=d.dtype), z + 2 * x)


def test_from_array_with_missing_chunks():
    x = np.random.randn(2, 4, 3)
    d = da.from_array(x, chunks=(None, 2, None))
    assert d.chunks == da.from_array(x, chunks=(2, 2, 3)).chunks


def test_normalize_chunks():
    assert normalize_chunks(3, (4, 6)) == ((3, 1), (3, 3))


def test_raise_on_no_chunks():
    x = da.ones(6, chunks=3)
    try:
        Array(x.dask, x.name, chunks=None, dtype=x.dtype, shape=None)
        assert False
    except ValueError as e:
        assert "dask.pydata.org" in str(e)

    pytest.raises(ValueError, lambda: da.ones(6))


def test_chunks_is_immutable():
    x = da.ones(6, chunks=3)
    try:
        x.chunks = 2
        assert False
    except TypeError as e:
        assert 'rechunk(2)' in str(e)


def test_raise_on_bad_kwargs():
    x = da.ones(5, chunks=3)
    try:
        da.minimum(x, foo=None)
    except TypeError as e:
        assert 'minimum' in str(e)
        assert 'foo' in str(e)


def test_long_slice():
    x = np.arange(10000)
    d = da.from_array(x, chunks=1)

    assert_eq(d[8000:8200], x[8000:8200])


def test_h5py_newaxis():
    h5py = pytest.importorskip('h5py')

    with tmpfile('h5') as fn:
        with h5py.File(fn) as f:
            x = f.create_dataset('/x', shape=(10, 10), dtype='f8')
            d = da.from_array(x, chunks=(5, 5))
            assert d[None, :, :].compute(get=get_sync).shape == (1, 10, 10)
            assert d[:, None, :].compute(get=get_sync).shape == (10, 1, 10)
            assert d[:, :, None].compute(get=get_sync).shape == (10, 10, 1)
            assert same_keys(d[:, :, None], d[:, :, None])


def test_ellipsis_slicing():
    assert_eq(da.ones(4, chunks=2)[...], np.ones(4))


def test_point_slicing():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    result = d.vindex[[1, 2, 5, 5], [3, 1, 6, 1]]
    assert_eq(result, x[[1, 2, 5, 5], [3, 1, 6, 1]])

    result = d.vindex[[0, 1, 6, 0], [0, 1, 0, 7]]
    assert_eq(result, x[[0, 1, 6, 0], [0, 1, 0, 7]])
    assert same_keys(result, d.vindex[[0, 1, 6, 0], [0, 1, 0, 7]])


def test_point_slicing_with_full_slice():
    from dask.array.core import _vindex_transpose, _get_axis
    x = np.arange(4 * 5 * 6 * 7).reshape((4, 5, 6, 7))
    d = da.from_array(x, chunks=(2, 3, 3, 4))

    inds = [[[1, 2, 3], None, [3, 2, 1], [5, 3, 4]],
            [[1, 2, 3], None, [4, 3, 2], None],
            [[1, 2, 3], [3, 2, 1]],
            [[1, 2, 3], [3, 2, 1], [3, 2, 1], [5, 3, 4]],
            [[], [], [], None],
            [np.array([1, 2, 3]), None, np.array([4, 3, 2]), None],
            [None, None, [1, 2, 3], [4, 3, 2]],
            [None, [0, 2, 3], None, [0, 3, 2]]]

    for ind in inds:
        slc = [i if isinstance(i, (np.ndarray, list)) else slice(None, None)
               for i in ind]
        result = d.vindex[tuple(slc)]

        # Rotate the expected result accordingly
        axis = _get_axis(ind)
        expected = _vindex_transpose(x[tuple(slc)], axis)

        assert_eq(result, expected)

        # Always have the first axis be the length of the points
        k = len(next(i for i in ind if isinstance(i, (np.ndarray, list))))
        assert result.shape[0] == k


def test_slice_with_floats():
    d = da.ones((5,), chunks=(3,))
    with pytest.raises(IndexError):
        d[1.5]
    with pytest.raises(IndexError):
        d[0:1.5]
    with pytest.raises(IndexError):
        d[[1, 1.5]]


def test_slice_with_integer_types():
    x = np.arange(10)
    dx = da.from_array(x, chunks=5)
    inds = np.array([0, 3, 6], dtype='u8')
    assert_eq(dx[inds], x[inds])
    assert_eq(dx[inds.astype('u4')], x[inds.astype('u4')])

    inds = np.array([0, 3, 6], dtype=np.int64)
    assert_eq(dx[inds], x[inds])
    assert_eq(dx[inds.astype('u4')], x[inds.astype('u4')])


def test_index_with_integer_types():
    x = np.arange(10)
    dx = da.from_array(x, chunks=5)
    inds = int(3)
    assert_eq(dx[inds], x[inds])

    inds = np.int64(3)
    assert_eq(dx[inds], x[inds])


def test_vindex_basic():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    # cases where basic and advanced indexing coincide
    result = d.vindex[0]
    assert_eq(result, x[0])

    result = d.vindex[0, 1]
    assert_eq(result, x[0, 1])

    result = d.vindex[[0, 1], ::-1]  # slices last
    assert_eq(result, x[:2, ::-1])


def test_vindex_nd():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(3, 4))

    result = d.vindex[[[0, 1], [6, 0]], [[0, 1], [0, 7]]]
    assert_eq(result, x[[[0, 1], [6, 0]], [[0, 1], [0, 7]]])

    result = d.vindex[np.arange(7)[:, None], np.arange(8)[None, :]]
    assert_eq(result, x)

    result = d.vindex[np.arange(7)[None, :], np.arange(8)[:, None]]
    assert_eq(result, x.T)


def test_vindex_negative():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    result = d.vindex[np.array([0, -1])]
    assert_eq(result, x[np.array([0, -1])])


def test_vindex_errors():
    d = da.ones((5, 5, 5), chunks=(3, 3, 3))
    pytest.raises(IndexError, lambda: d.vindex[np.newaxis])
    pytest.raises(IndexError, lambda: d.vindex[:5])
    pytest.raises(IndexError, lambda: d.vindex[[1, 2], [1, 2, 3]])
    pytest.raises(IndexError, lambda: d.vindex[[True] * 5])
    pytest.raises(IndexError, lambda: d.vindex[[0], [5]])
    pytest.raises(IndexError, lambda: d.vindex[[0], [-6]])


def test_vindex_merge():
    from dask.array.core import _vindex_merge
    locations = [1], [2, 0]
    values = [np.array([[1, 2, 3]]),
              np.array([[10, 20, 30], [40, 50, 60]])]

    assert (_vindex_merge(locations, values) == np.array([[40, 50, 60],
                                                          [1, 2, 3],
                                                          [10, 20, 30]])).all()


def test_empty_array():
    assert_eq(np.arange(0), da.arange(0, chunks=5))


def test_memmap():
    with tmpfile('npy') as fn_1:
        with tmpfile('npy') as fn_2:
            try:
                x = da.arange(100, chunks=15)
                target = np.memmap(fn_1, shape=x.shape, mode='w+', dtype=x.dtype)

                x.store(target)

                assert_eq(target, x)

                np.save(fn_2, target)

                assert_eq(np.load(fn_2, mmap_mode='r'), x)
            finally:
                target._mmap.close()


def test_to_npy_stack():
    x = np.arange(5 * 10 * 10).reshape((5, 10, 10))
    d = da.from_array(x, chunks=(2, 4, 4))

    with tmpdir() as dirname:
        stackdir = os.path.join(dirname, 'test')
        da.to_npy_stack(stackdir, d, axis=0)
        assert os.path.exists(os.path.join(stackdir, '0.npy'))
        assert (np.load(os.path.join(stackdir, '1.npy')) == x[2:4]).all()

        e = da.from_npy_stack(stackdir)
        assert_eq(d, e)


def test_view():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(2, 3))

    assert_eq(x.view('i4'), d.view('i4'))
    assert_eq(x.view('i2'), d.view('i2'))
    assert all(isinstance(s, int) for s in d.shape)

    x = np.arange(8, dtype='i1')
    d = da.from_array(x, chunks=(4,))
    assert_eq(x.view('i4'), d.view('i4'))

    with pytest.raises(ValueError):
        x = np.arange(8, dtype='i1')
        d = da.from_array(x, chunks=(3,))
        d.view('i4')

    with pytest.raises(ValueError):
        d.view('i4', order='asdf')


def test_view_fortran():
    x = np.asfortranarray(np.arange(64).reshape((8, 8)))
    d = da.from_array(x, chunks=(2, 3))
    assert_eq(x.T.view('i4').T, d.view('i4', order='F'))
    assert_eq(x.T.view('i2').T, d.view('i2', order='F'))


def test_h5py_tokenize():
    h5py = pytest.importorskip('h5py')
    with tmpfile('hdf5') as fn1:
        with tmpfile('hdf5') as fn2:
            f = h5py.File(fn1)
            g = h5py.File(fn2)

            f['x'] = np.arange(10).astype(float)
            g['x'] = np.ones(10).astype(float)

            x1 = f['x']
            x2 = g['x']

            assert tokenize(x1) != tokenize(x2)


def test_map_blocks_with_changed_dimension():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(7, 4))

    e = d.map_blocks(lambda b: b.sum(axis=0), chunks=(4,), drop_axis=0,
                     dtype=d.dtype)
    assert e.chunks == ((4, 4),)
    assert_eq(e, x.sum(axis=0))

    # Provided chunks have wrong shape
    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=0), chunks=(7, 4), drop_axis=0)

    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=0), chunks=((4, 4, 4),), drop_axis=0)

    # Can't drop axis with more than 1 block
    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b.sum(axis=1), drop_axis=1, dtype=d.dtype)

    # Adding axis with a gap
    with pytest.raises(ValueError):
        d.map_blocks(lambda b: b, new_axis=(3, 4))

    d = da.from_array(x, chunks=(4, 8))
    e = d.map_blocks(lambda b: b.sum(axis=1), drop_axis=1, dtype=d.dtype)
    assert e.chunks == ((4, 3),)
    assert_eq(e, x.sum(axis=1))

    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))
    e = d.map_blocks(lambda b: b[None, :, :, None],
                     chunks=(1, 4, 4, 1), new_axis=[0, 3], dtype=d.dtype)
    assert e.chunks == ((1,), (4, 4), (4, 4), (1,))
    assert_eq(e, x[None, :, :, None])

    e = d.map_blocks(lambda b: b[None, :, :, None],
                     new_axis=[0, 3], dtype=d.dtype)
    assert e.chunks == ((1,), (4, 4), (4, 4), (1,))
    assert_eq(e, x[None, :, :, None])

    # Both new_axis and drop_axis
    d = da.from_array(x, chunks=(8, 4))
    e = d.map_blocks(lambda b: b.sum(axis=0)[:, None, None],
                     drop_axis=0, new_axis=(1, 2), dtype=d.dtype)
    assert e.chunks == ((4, 4), (1,), (1,))
    assert_eq(e, x.sum(axis=0)[:, None, None])

    d = da.from_array(x, chunks=(4, 8))
    e = d.map_blocks(lambda b: b.sum(axis=1)[:, None, None],
                     drop_axis=1, new_axis=(1, 2), dtype=d.dtype)
    assert e.chunks == ((4, 4), (1,), (1,))
    assert_eq(e, x.sum(axis=1)[:, None, None])


def test_broadcast_chunks():
    assert broadcast_chunks() == ()

    assert broadcast_chunks(((2, 3),)) == ((2, 3),)

    assert broadcast_chunks(((5, 5),), ((5, 5),)) == ((5, 5),)

    a = ((10, 10, 10), (5, 5),)
    b = ((5, 5),)
    assert broadcast_chunks(a, b) == ((10, 10, 10), (5, 5),)
    assert broadcast_chunks(b, a) == ((10, 10, 10), (5, 5),)

    a = ((10, 10, 10), (5, 5),)
    b = ((1,), (5, 5),)
    assert broadcast_chunks(a, b) == ((10, 10, 10), (5, 5),)

    a = ((10, 10, 10), (5, 5),)
    b = ((3, 3,), (5, 5),)
    with pytest.raises(ValueError):
        broadcast_chunks(a, b)

    a = ((1,), (5, 5),)
    b = ((1,), (5, 5),)
    assert broadcast_chunks(a, b) == a

    a = ((1,), (np.nan, np.nan, np.nan),)
    b = ((3, 3), (1,),)
    r = broadcast_chunks(a, b)
    assert r[0] == b[0] and np.allclose(r[1], a[1], equal_nan=True)

    a = ((3, 3), (1,),)
    b = ((1,), (np.nan, np.nan, np.nan),)
    r = broadcast_chunks(a, b)
    assert r[0] == a[0] and np.allclose(r[1], b[1], equal_nan=True)

    a = ((3, 3,), (5, 5),)
    b = ((1,), (np.nan, np.nan, np.nan),)
    with pytest.raises(ValueError):
        broadcast_chunks(a, b)


def test_chunks_error():
    x = np.ones((10, 10))
    with pytest.raises(ValueError):
        da.from_array(x, chunks=(5,))


def test_array_compute_forward_kwargs():
    x = da.arange(10, chunks=2).sum()
    x.compute(bogus_keyword=10)


def test_dont_fuse_outputs():
    dsk = {('x', 0): np.array([1, 2]),
           ('x', 1): (inc, ('x', 0))}

    a = da.Array(dsk, 'x', chunks=(2,), shape=(4,), dtype=np.array([1]).dtype)
    assert_eq(a, np.array([1, 2, 2, 3], dtype=a.dtype))


def test_dont_dealias_outputs():
    dsk = {('x', 0, 0): np.ones((2, 2)),
           ('x', 0, 1): np.ones((2, 2)),
           ('x', 1, 0): np.ones((2, 2)),
           ('x', 1, 1): ('x', 0, 0)}

    a = da.Array(dsk, 'x', chunks=(2, 2), shape=(4, 4), dtype=np.ones(1).dtype)
    assert_eq(a, np.ones((4, 4)))


def test_timedelta_op():
    x = np.array([np.timedelta64(10, 'h')])
    y = np.timedelta64(1, 'h')
    a = da.from_array(x, chunks=(1,)) / y
    assert a.compute() == x / y


def test_to_delayed():
    x = da.random.random((4, 4), chunks=(2, 2))
    y = x + 10

    [[a, b], [c, d]] = y.to_delayed()
    assert_eq(a.compute(), y[:2, :2])

    s = 2
    x = da.from_array(np.array(s), chunks=0)
    a = x.to_delayed()[tuple()]
    assert a.compute() == s


def test_to_delayed_optimize_graph():
    x = da.ones((4, 4), chunks=(2, 2))
    y = x[1:][1:][1:][:, 1:][:, 1:][:, 1:]

    # optimizations
    d = y.to_delayed().flatten().tolist()[0]
    assert len([k for k in d.dask if k[0].startswith('getitem')]) == 1

    # no optimizations
    d2 = y.to_delayed(optimize_graph=False).flatten().tolist()[0]
    assert dict(d2.dask) == dict(y.dask)

    assert (d.compute() == d2.compute()).all()


def test_cumulative():
    x = da.arange(20, chunks=5)
    assert_eq(x.cumsum(axis=0), np.arange(20).cumsum())
    assert_eq(x.cumprod(axis=0), np.arange(20).cumprod())

    assert_eq(da.nancumsum(x, axis=0), nancumsum(np.arange(20)))
    assert_eq(da.nancumprod(x, axis=0), nancumprod(np.arange(20)))

    a = np.random.random((20))
    rs = np.random.RandomState(0)
    a[rs.rand(*a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=5)
    assert_eq(da.nancumsum(x, axis=0), nancumsum(a))
    assert_eq(da.nancumprod(x, axis=0), nancumprod(a))

    a = np.random.random((20, 24))
    x = da.from_array(a, chunks=(6, 5))
    assert_eq(x.cumsum(axis=0), a.cumsum(axis=0))
    assert_eq(x.cumsum(axis=1), a.cumsum(axis=1))
    assert_eq(x.cumprod(axis=0), a.cumprod(axis=0))
    assert_eq(x.cumprod(axis=1), a.cumprod(axis=1))

    assert_eq(da.nancumsum(x, axis=0), nancumsum(a, axis=0))
    assert_eq(da.nancumsum(x, axis=1), nancumsum(a, axis=1))
    assert_eq(da.nancumprod(x, axis=0), nancumprod(a, axis=0))
    assert_eq(da.nancumprod(x, axis=1), nancumprod(a, axis=1))

    a = np.random.random((20, 24))
    rs = np.random.RandomState(0)
    a[rs.rand(*a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=(6, 5))
    assert_eq(da.nancumsum(x, axis=0), nancumsum(a, axis=0))
    assert_eq(da.nancumsum(x, axis=1), nancumsum(a, axis=1))
    assert_eq(da.nancumprod(x, axis=0), nancumprod(a, axis=0))
    assert_eq(da.nancumprod(x, axis=1), nancumprod(a, axis=1))

    a = np.random.random((20, 24, 13))
    x = da.from_array(a, chunks=(6, 5, 4))
    for axis in [0, 1, 2, -1, -2, -3]:
        assert_eq(x.cumsum(axis=axis), a.cumsum(axis=axis))
        assert_eq(x.cumprod(axis=axis), a.cumprod(axis=axis))

        assert_eq(da.nancumsum(x, axis=axis), nancumsum(a, axis=axis))
        assert_eq(da.nancumprod(x, axis=axis), nancumprod(a, axis=axis))

    a = np.random.random((20, 24, 13))
    rs = np.random.RandomState(0)
    a[rs.rand(*a.shape) < 0.5] = np.nan
    x = da.from_array(a, chunks=(6, 5, 4))
    for axis in [0, 1, 2, -1, -2, -3]:
        assert_eq(da.nancumsum(x, axis=axis), nancumsum(a, axis=axis))
        assert_eq(da.nancumprod(x, axis=axis), nancumprod(a, axis=axis))

    with pytest.raises(ValueError):
        x.cumsum(axis=3)

    with pytest.raises(ValueError):
        x.cumsum(axis=-4)


def test_atop_names():
    x = da.ones(5, chunks=(2,))
    y = atop(add, 'i', x, 'i', dtype=x.dtype)
    assert y.name.startswith('add')


def test_atop_new_axes():
    def f(x):
        return x[:, None] * np.ones((1, 7))
    x = da.ones(5, chunks=2)
    y = atop(f, 'aq', x, 'a', new_axes={'q': 7}, concatenate=True,
             dtype=x.dtype)
    assert y.chunks == ((2, 2, 1), (7,))
    assert_eq(y, np.ones((5, 7)))

    def f(x):
        return x[None, :] * np.ones((7, 1))
    x = da.ones(5, chunks=2)
    y = atop(f, 'qa', x, 'a', new_axes={'q': 7}, concatenate=True,
             dtype=x.dtype)
    assert y.chunks == ((7,), (2, 2, 1))
    assert_eq(y, np.ones((7, 5)))

    def f(x):
        y = x.sum(axis=1)
        return y[:, None] * np.ones((1, 5))

    x = da.ones((4, 6), chunks=(2, 2))
    y = atop(f, 'aq', x, 'ab', new_axes={'q': 5}, concatenate=True,
             dtype=x.dtype)
    assert y.chunks == ((2, 2), (5,))
    assert_eq(y, np.ones((4, 5)) * 6)


def test_atop_kwargs():
    def f(a, b=0):
        return a + b

    x = da.ones(5, chunks=(2,))
    y = atop(f, 'i', x, 'i', b=10, dtype=x.dtype)
    assert_eq(y, np.ones(5) + 10)


def test_atop_chunks():
    x = da.ones((5, 5), chunks=((2, 1, 2), (3, 2)))

    def double(a, axis=0):
        return np.concatenate([a, a], axis=axis)

    y = atop(double, 'ij', x, 'ij',
             adjust_chunks={'i': lambda n: 2 * n}, axis=0, dtype=x.dtype)
    assert y.chunks == ((4, 2, 4), (3, 2))
    assert_eq(y, np.ones((10, 5)))

    y = atop(double, 'ij', x, 'ij',
             adjust_chunks={'j': lambda n: 2 * n}, axis=1, dtype=x.dtype)
    assert y.chunks == ((2, 1, 2), (6, 4))
    assert_eq(y, np.ones((5, 10)))

    x = da.ones((10, 10), chunks=(5, 5))
    y = atop(double, 'ij', x, 'ij', axis=0,
             adjust_chunks={'i': 10}, dtype=x.dtype)
    assert y.chunks == ((10, 10), (5, 5))
    assert_eq(y, np.ones((20, 10)))

    y = atop(double, 'ij', x, 'ij', axis=0,
             adjust_chunks={'i': (10, 10)}, dtype=x.dtype)
    assert y.chunks == ((10, 10), (5, 5))
    assert_eq(y, np.ones((20, 10)))


def test_atop_raises_on_incorrect_indices():
    x = da.arange(5, chunks=3)
    with pytest.raises(ValueError) as info:
        da.atop(lambda x: x, 'ii', x, 'ii', dtype=int)

    assert 'ii' in str(info.value)
    assert '1' in str(info.value)


def test_from_delayed():
    v = delayed(np.ones)((5, 3))
    x = from_delayed(v, shape=(5, 3), dtype=np.ones(0).dtype)
    assert isinstance(x, Array)
    assert_eq(x, np.ones((5, 3)))


def test_A_property():
    x = da.ones(5, chunks=(2,))
    assert x.A is x


def test_copy_mutate():
    x = da.arange(5, chunks=(2,))
    y = x.copy()
    memo = {}
    y2 = copy.deepcopy(x, memo=memo)
    x[x % 2 == 0] = -1

    xx = np.arange(5)
    xx[xx % 2 == 0] = -1
    assert_eq(x, xx)

    assert_eq(y, np.arange(5))
    assert_eq(y2, np.arange(5))
    assert memo[id(x)] is y2


def test_npartitions():
    assert da.ones(5, chunks=(2,)).npartitions == 3
    assert da.ones((5, 5), chunks=(2, 3)).npartitions == 6


def test_astype_gh1151():
    a = np.arange(5).astype(np.int32)
    b = da.from_array(a, (1,))
    assert_eq(a.astype(np.int16), b.astype(np.int16))


def test_elemwise_name():
    assert (da.ones(5, chunks=2) + 1).name.startswith('add-')


def test_map_blocks_name():
    assert da.ones(5, chunks=2).map_blocks(inc).name.startswith('inc-')


def test_from_array_names():
    pytest.importorskip('distributed')
    from distributed.utils import key_split

    x = np.ones(10)
    d = da.from_array(x, chunks=2)

    names = countby(key_split, d.dask)
    assert set(names.values()) == set([1, 5])


def test_array_picklable():
    from pickle import loads, dumps

    a = da.arange(100, chunks=25)
    a2 = loads(dumps(a))
    assert_eq(a, a2)


def test_from_array_raises_on_bad_chunks():
    x = np.ones(10)

    with pytest.raises(ValueError):
        da.from_array(x, chunks=(5, 5, 5))

    # with pytest.raises(ValueError):
    #      da.from_array(x, chunks=100)

    with pytest.raises(ValueError):
        da.from_array(x, chunks=((5, 5, 5),))


def test_concatenate_axes():
    x = np.ones((2, 2, 2))

    assert_eq(concatenate_axes([x, x], axes=[0]),
              np.ones((4, 2, 2)))
    assert_eq(concatenate_axes([x, x, x], axes=[0]),
              np.ones((6, 2, 2)))
    assert_eq(concatenate_axes([x, x], axes=[1]),
              np.ones((2, 4, 2)))
    assert_eq(concatenate_axes([[x, x], [x, x]], axes=[0, 1]),
              np.ones((4, 4, 2)))
    assert_eq(concatenate_axes([[x, x], [x, x]], axes=[0, 2]),
              np.ones((4, 2, 4)))
    assert_eq(concatenate_axes([[x, x, x], [x, x, x]], axes=[1, 2]),
              np.ones((2, 4, 6)))

    with pytest.raises(ValueError):
        concatenate_axes([[x, x], [x, x]], axes=[0])  # not all nested lists accounted for
    with pytest.raises(ValueError):
        concatenate_axes([x, x], axes=[0, 1, 2, 3])  # too many axes


def test_atop_concatenate():
    x = da.ones((4, 4, 4), chunks=(2, 2, 2))
    y = da.ones((4, 4), chunks=(2, 2))

    def f(a, b):
        assert isinstance(a, np.ndarray)
        assert isinstance(b, np.ndarray)

        assert a.shape == (2, 4, 4)
        assert b.shape == (4, 4)

        return (a + b).sum(axis=(1, 2))

    z = atop(f, 'i', x, 'ijk', y, 'jk', concatenate=True, dtype=x.dtype)
    assert_eq(z, np.ones(4) * 32)

    z = atop(add, 'ij', y, 'ij', y, 'ij', concatenate=True, dtype=x.dtype)
    assert_eq(z, np.ones((4, 4)) * 2)

    def f(a, b, c):
        assert isinstance(a, np.ndarray)
        assert isinstance(b, np.ndarray)
        assert isinstance(c, np.ndarray)

        assert a.shape == (4, 2, 4)
        assert b.shape == (4, 4)
        assert c.shape == (4, 2)

        return np.ones(5)

    z = atop(f, 'j', x, 'ijk', y, 'ki', y, 'ij', concatenate=True,
             dtype=x.dtype)
    assert_eq(z, np.ones(10), check_shape=False)


def test_common_blockdim():
    assert common_blockdim([(5,), (5,)]) == (5,)
    assert common_blockdim([(5,), (2, 3,)]) == (2, 3)
    assert common_blockdim([(5, 5), (2, 3, 5)]) == (2, 3, 5)
    assert common_blockdim([(5, 5), (2, 3, 5)]) == (2, 3, 5)
    assert common_blockdim([(5, 2, 3), (2, 3, 5)]) == (2, 3, 2, 3)

    assert common_blockdim([(1, 2), (2, 1)]) == (1, 1, 1)
    assert common_blockdim([(1, 2, 2), (2, 1, 2), (2, 2, 1)]) == (1, 1, 1, 1, 1)


def test_uneven_chunks_that_fit_neatly():
    x = da.arange(10, chunks=((5, 5),))
    y = da.ones(10, chunks=((5, 2, 3),))

    assert_eq(x + y, np.arange(10) + np.ones(10))

    z = x + y
    assert z.chunks == ((5, 2, 3),)


def test_elemwise_uneven_chunks():
    x = da.arange(10, chunks=((4, 6),))
    y = da.ones(10, chunks=((6, 4),))

    assert_eq(x + y, np.arange(10) + np.ones(10))

    z = x + y
    assert z.chunks == ((4, 2, 4),)

    x = da.random.random((10, 10), chunks=((4, 6), (5, 2, 3)))
    y = da.random.random((4, 10, 10), chunks=((2, 2), (6, 4), (2, 3, 5)))

    z = x + y
    assert_eq(x + y, x.compute() + y.compute())
    assert z.chunks == ((2, 2), (4, 2, 4), (2, 3, 2, 3))


def test_uneven_chunks_atop():
    x = da.random.random((10, 10), chunks=((2, 3, 2, 3), (5, 5)))
    y = da.random.random((10, 10), chunks=((4, 4, 2), (4, 2, 4)))
    z = atop(np.dot, 'ik', x, 'ij', y, 'jk', dtype=x.dtype, concatenate=True)
    assert z.chunks == (x.chunks[0], y.chunks[1])

    assert_eq(z, x.compute().dot(y))


def test_warn_bad_rechunking():
    x = da.ones((20, 20), chunks=(20, 1))
    y = da.ones((20, 20), chunks=(1, 20))

    with warnings.catch_warnings(record=True) as record:
        x + y

    assert record
    assert '20' in record[0].message.args[0]


def test_optimize_fuse_keys():
    x = da.ones(10, chunks=(5,))
    y = x + 1
    z = y + 1

    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__())
    assert not set(y.dask) & set(dsk)

    dsk = z.__dask_optimize__(z.dask, z.__dask_keys__(),
                              fuse_keys=y.__dask_keys__())
    assert all(k in dsk for k in y.__dask_keys__())


def test_concatenate_stack_dont_warn():
    with warnings.catch_warnings(record=True) as record:
        da.concatenate([da.ones(2, chunks=1)] * 62)
    assert not record

    with warnings.catch_warnings(record=True) as record:
        da.stack([da.ones(2, chunks=1)] * 62)
    assert not record


def test_map_blocks_delayed():
    x = da.ones((10, 10), chunks=(5, 5))
    y = np.ones((5, 5))

    z = x.map_blocks(add, y, dtype=x.dtype)

    yy = delayed(y)
    zz = x.map_blocks(add, yy, dtype=x.dtype)

    assert_eq(z, zz)

    assert yy.key in zz.dask


def test_no_chunks():
    X = np.arange(11)
    dsk = {('x', 0): np.arange(5), ('x', 1): np.arange(5, 11)}
    x = Array(dsk, 'x', ((np.nan, np.nan,),), np.arange(1).dtype)
    assert_eq(x + 1, X + 1)
    assert_eq(x.sum(), X.sum())
    assert_eq((x + 1).std(), (X + 1).std())
    assert_eq((x + x).std(), (X + X).std())
    assert_eq((x + x).std(keepdims=True), (X + X).std(keepdims=True))


def test_no_chunks_2d():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((np.nan, np.nan), (np.nan, np.nan, np.nan))

    with pytest.warns(None):  # zero division warning
        assert_eq(da.log(x), np.log(X))
    assert_eq(x.T, X.T)
    assert_eq(x.sum(axis=0, keepdims=True), X.sum(axis=0, keepdims=True))
    assert_eq(x.sum(axis=1, keepdims=True), X.sum(axis=1, keepdims=True))
    assert_eq(x.dot(x.T + 1), X.dot(X.T + 1))


def test_no_chunks_yes_chunks():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    assert (x + 1).chunks == ((2, 2), (np.nan, np.nan, np.nan))
    assert (x.T).chunks == ((np.nan, np.nan, np.nan), (2, 2))
    assert (x.dot(x.T)).chunks == ((2, 2), (2, 2))


def test_raise_informative_errors_no_chunks():
    X = np.arange(10)
    a = da.from_array(X, chunks=(5, 5))
    a._chunks = ((np.nan, np.nan),)

    b = da.from_array(X, chunks=(4, 4, 2))
    b._chunks = ((np.nan, np.nan, np.nan),)

    for op in [lambda: a + b,
               lambda: a[1],
               lambda: a[::2],
               lambda: a[-5],
               lambda: a.rechunk(3),
               lambda: a.reshape(2, 5)]:
        with pytest.raises(ValueError) as e:
            op()
        if 'chunk' not in str(e) or 'unknown' not in str(e):
            op()


def test_no_chunks_slicing_2d():
    X = np.arange(24).reshape((4, 6))
    x = da.from_array(X, chunks=(2, 2))
    x._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    assert_eq(x[0], X[0])

    for op in [lambda: x[:, 4],
               lambda: x[:, ::2],
               lambda: x[0, 2:4]]:
        with pytest.raises(ValueError) as e:
            op()
        assert 'chunk' in str(e) and 'unknown' in str(e)


def test_index_array_with_array_1d():
    x = np.arange(10)
    dx = da.from_array(x, chunks=(5,))
    dx._chunks = ((np.nan, np.nan),)

    assert_eq(x[x > 6], dx[dx > 6])
    assert_eq(x[x % 2 == 0], dx[dx % 2 == 0])

    dy = da.ones(11, chunks=(3,))

    with pytest.raises(ValueError):
        dx[dy > 5]


def test_index_array_with_array_2d():
    x = np.arange(24).reshape((4, 6))
    dx = da.from_array(x, chunks=(2, 2))
    dx._chunks = ((2, 2), (np.nan, np.nan, np.nan))

    assert (sorted(x[x % 2 == 0].tolist()) ==
            sorted(dx[dx % 2 == 0].compute().tolist()))
    assert (sorted(x[x > 6].tolist()) ==
            sorted(dx[dx > 6].compute().tolist()))


@pytest.mark.xfail(reason='Chunking does not align well')
def test_index_array_with_array_3d_2d():
    x = np.arange(4**3).reshape((4, 4, 4))
    dx = da.from_array(x, chunks=(2, 2, 2))

    ind = np.random.random((4, 4)) > 0.5
    ind = np.arange(4 ** 2).reshape((4, 4)) % 2 == 0
    dind = da.from_array(ind, (2, 2))

    assert_eq(x[ind], dx[dind])
    assert_eq(x[:, ind], dx[:, dind])


def test_setitem_1d():
    x = np.arange(10)
    dx = da.from_array(x.copy(), chunks=(5,))

    x[x > 6] = -1
    x[x % 2 == 0] = -2

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2

    assert_eq(x, dx)


def test_setitem_2d():
    x = np.arange(24).reshape((4, 6))
    dx = da.from_array(x.copy(), chunks=(2, 2))

    x[x > 6] = -1
    x[x % 2 == 0] = -2

    dx[dx > 6] = -1
    dx[dx % 2 == 0] = -2

    assert_eq(x, dx)


@pytest.mark.skipif(np.__version__ >= '1.13.0',
                    reason='boolean slicing rules changed')
def test_setitem_mixed_d():
    x = np.arange(24).reshape((4, 6))
    dx = da.from_array(x, chunks=(2, 2))

    x[x[0, None] > 2] = -1
    dx[dx[0, None] > 2] = -1
    assert_eq(x, dx)

    x[x[None, 0] > 2] = -1
    dx[dx[None, 0] > 2] = -1
    assert_eq(x, dx)


def test_setitem_errs():
    x = da.ones((4, 4), chunks=(2, 2))

    with pytest.raises(ValueError):
        x[x > 1] = x


def test_zero_slice_dtypes():
    x = da.arange(5, chunks=1)
    y = x[[]]
    assert y.dtype == x.dtype
    assert y.shape == (0,)
    assert_eq(x[[]], np.arange(5)[[]])


def test_zero_sized_array_rechunk():
    x = da.arange(5, chunks=1)[:0]
    y = da.atop(identity, 'i', x, 'i', dtype=x.dtype)
    assert_eq(x, y)


def test_atop_zero_shape():
    da.atop(lambda x: x, 'i',
            da.arange(10, chunks=10), 'i',
            da.from_array(np.ones((0, 2)), ((0,), 2)), 'ab',
            da.from_array(np.ones((0,)), ((0,),)), 'a',
            dtype='float64')


def test_atop_zero_shape_new_axes():
    da.atop(lambda x: np.ones(42), 'i',
            da.from_array(np.ones((0, 2)), ((0,), 2)), 'ab',
            da.from_array(np.ones((0,)), ((0,),)), 'a',
            dtype='float64', new_axes={'i': 42})


def test_broadcast_against_zero_shape():
    assert_eq(da.arange(1, chunks=1)[:0] + 0,
              np.arange(1)[:0] + 0)
    assert_eq(da.arange(1, chunks=1)[:0] + 0.1,
              np.arange(1)[:0] + 0.1)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:0] + 0,
              np.ones((5, 5))[:0] + 0)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:0] + 0.1,
              np.ones((5, 5))[:0] + 0.1)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:, :0] + 0,
              np.ones((5, 5))[:, :0] + 0)
    assert_eq(da.ones((5, 5), chunks=(2, 3))[:, :0] + 0.1,
              np.ones((5, 5))[:, :0] + 0.1)


def test_from_array_name():
    x = np.array([1, 2, 3, 4, 5])
    chunks = x.shape
    # Default is tokenize the array
    dx = da.from_array(x, chunks=chunks)
    hashed_name = dx.name
    assert da.from_array(x, chunks=chunks).name == hashed_name
    # Specify name directly
    assert da.from_array(x, chunks=chunks, name='x').name == 'x'
    # False gives a random name
    dx2 = da.from_array(x, chunks=chunks, name=False)
    dx3 = da.from_array(x, chunks=chunks, name=False)
    assert dx2.name != hashed_name
    assert dx3.name != hashed_name
    assert dx2.name != dx3.name


def test_concatenate_errs():
    with pytest.raises(ValueError) as e:
        da.concatenate([da.zeros((2, 1), chunks=(2, 1)),
                        da.zeros((2, 3), chunks=(2, 3))])

    assert 'shape' in str(e).lower()
    assert '(2, 1)' in str(e)

    with pytest.raises(ValueError):
        da.concatenate([da.zeros((1, 2), chunks=(1, 2)),
                        da.zeros((3, 2), chunks=(3, 2))], axis=1)


def test_stack_errs():
    with pytest.raises(ValueError) as e:
        da.stack([da.zeros((2), chunks=(2)),
                  da.zeros((3), chunks=(3))])
    assert 'shape' in str(e).lower()
    assert '(2,)' in str(e)


def test_atop_with_numpy_arrays():
    x = np.ones(10)
    y = da.ones(10, chunks=(5,))

    assert_eq(x + y, x + x)

    s = da.sum(x)
    assert any(x is v for v in s.dask.values())


@pytest.mark.parametrize('chunks', (100, 6))
@pytest.mark.parametrize('other', [[0, 0, 1], [2, 1, 3], (0, 0, 1)])
def test_elemwise_with_lists(chunks, other):
    x = np.arange(12).reshape((4, 3))
    d = da.arange(12, chunks=chunks).reshape((4, 3))

    x2 = np.vstack([x[:, 0], x[:, 1], x[:, 2]]).T
    d2 = da.vstack([d[:, 0], d[:, 1], d[:, 2]]).T

    assert_eq(x2, d2)

    x3 = x2 * other
    d3 = d2 * other

    assert_eq(x3, d3)


def test_constructor_plugin():
    L = []
    L2 = []
    with dask.set_options(array_plugins=[L.append, L2.append]):
        x = da.ones(10, chunks=5)
        y = x + 1

    assert L == L2 == [x, y]

    with dask.set_options(array_plugins=[lambda x: x.compute()]):
        x = da.ones(10, chunks=5)
        y = x + 1

    assert isinstance(y, np.ndarray)
    assert len(L) == 2


def test_no_warnings_on_metadata():
    x = da.ones(5, chunks=3)
    with warnings.catch_warnings(record=True) as record:
        da.arccos(x)

    assert not record


def test_delayed_array_key_hygeine():
    a = da.zeros((1,), chunks=(1,))
    d = delayed(identity)(a)
    b = da.from_delayed(d, shape=a.shape, dtype=a.dtype)
    assert_eq(a, b)


def test_empty_chunks_in_array_len():
    x = da.ones((), chunks=())
    with pytest.raises(TypeError) as exc_info:
        len(x)

    err_msg = 'len() of unsized object'
    assert err_msg in str(exc_info.value)


@pytest.mark.parametrize('dtype', [None, [('a', 'f4'), ('b', object)]])
def test_meta(dtype):
    a = da.zeros((1,), chunks=(1,))
    assert a._meta.dtype == a.dtype
    assert isinstance(a._meta, np.ndarray)
    assert a.nbytes < 1000
