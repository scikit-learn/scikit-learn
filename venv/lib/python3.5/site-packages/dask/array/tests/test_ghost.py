import pytest
pytest.importorskip('numpy')

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

import dask.array as da
from dask.array.ghost import (fractional_slice, getitem, trim_internal,
                              ghost_internal, nearest, constant, boundaries,
                              reflect, periodic, ghost)
from dask.core import get
from dask.array.utils import assert_eq, same_keys


def test_fractional_slice():
    assert (fractional_slice(('x', 4.9), {0: 2}) ==
            (getitem, ('x', 5), (slice(0, 2), )))

    assert (fractional_slice(('x', 3, 5.1), {0: 2, 1: 3}) ==
            (getitem, ('x', 3, 5), (slice(None, None, None), slice(-3, None))))

    assert (fractional_slice(('x', 2.9, 5.1), {0: 2, 1: 3}) ==
            (getitem, ('x', 3, 5), (slice(0, 2), slice(-3, None))))

    fs = fractional_slice(('x', 4.9), {0: 2})
    assert isinstance(fs[1][1], int)


def test_ghost_internal():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    g = ghost_internal(d, {0: 2, 1: 1})
    result = g.compute(get=get)
    assert g.chunks == ((6, 6), (5, 5))

    expected = np.array([
        [ 0,  1,  2,  3,  4,    3,  4,  5,  6,  7],
        [ 8,  9, 10, 11, 12,   11, 12, 13, 14, 15],
        [16, 17, 18, 19, 20,   19, 20, 21, 22, 23],
        [24, 25, 26, 27, 28,   27, 28, 29, 30, 31],
        [32, 33, 34, 35, 36,   35, 36, 37, 38, 39],
        [40, 41, 42, 43, 44,   43, 44, 45, 46, 47],

        [16, 17, 18, 19, 20,   19, 20, 21, 22, 23],
        [24, 25, 26, 27, 28,   27, 28, 29, 30, 31],
        [32, 33, 34, 35, 36,   35, 36, 37, 38, 39],
        [40, 41, 42, 43, 44,   43, 44, 45, 46, 47],
        [48, 49, 50, 51, 52,   51, 52, 53, 54, 55],
        [56, 57, 58, 59, 60,   59, 60, 61, 62, 63]])

    assert_eq(result, expected)
    assert same_keys(ghost_internal(d, {0: 2, 1: 1}), g)


def test_trim_internal():
    d = da.ones((40, 60), chunks=(10, 10))
    e = trim_internal(d, axes={0: 1, 1: 2})

    assert e.chunks == ((8, 8, 8, 8), (6, 6, 6, 6, 6, 6))


def test_periodic():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = periodic(d, axis=0, depth=2)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], d[-1, :])
    assert_eq(e[0, :], d[-2, :])


def test_reflect():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    e = reflect(d, axis=0, depth=2)
    expected = np.array([1, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 8])
    assert_eq(e, expected)

    e = reflect(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected)


def test_nearest():
    x = np.arange(10)
    d = da.from_array(x, chunks=(5, 5))

    e = nearest(d, axis=0, depth=2)
    expected = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9])
    assert_eq(e, expected)

    e = nearest(d, axis=0, depth=1)
    expected = np.array([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9])
    assert_eq(e, expected)


def test_constant():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = constant(d, axis=0, depth=2, value=10)
    assert e.shape[0] == d.shape[0] + 4
    assert e.shape[1] == d.shape[1]

    assert_eq(e[1, :], np.ones(8, dtype=x.dtype) * 10)
    assert_eq(e[-1, :], np.ones(8, dtype=x.dtype) * 10)


def test_boundaries():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))

    e = boundaries(d, {0: 2, 1: 1}, {0: 0, 1: 'periodic'})

    expected = np.array(
        [[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [ 7, 0, 1, 2, 3, 4, 5, 6, 7, 0],
         [15, 8, 9,10,11,12,13,14,15, 8],
         [23,16,17,18,19,20,21,22,23,16],
         [31,24,25,26,27,28,29,30,31,24],
         [39,32,33,34,35,36,37,38,39,32],
         [47,40,41,42,43,44,45,46,47,40],
         [55,48,49,50,51,52,53,54,55,48],
         [63,56,57,58,59,60,61,62,63,56],
         [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    assert_eq(e, expected)


def test_ghost():
    x = np.arange(64).reshape((8, 8))
    d = da.from_array(x, chunks=(4, 4))
    g = ghost(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: 'reflect'})
    assert g.chunks == ((8, 8), (6, 6))
    expected = np.array(
        [[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
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
    assert_eq(g, expected)
    assert same_keys(g, ghost(d, depth={0: 2, 1: 1},
                              boundary={0: 100, 1: 'reflect'}))

    g = ghost(d, depth={0: 2, 1: 1}, boundary={0: 100, 1: 'none'})
    expected = np.array(
        [[100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
         [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
         [  0,   1,   2,   3,   4,   3,   4,   5,   6,   7],
         [  8,   9,  10,  11,  12,  11,  12,  13,  14,  15],
         [ 16,  17,  18,  19,  20,  19,  20,  21,  22,  23],
         [ 24,  25,  26,  27,  28,  27,  28,  29,  30,  31],
         [ 32,  33,  34,  35,  36,  35,  36,  37,  38,  39],
         [ 40,  41,  42,  43,  44,  43,  44,  45,  46,  47],
         [ 16,  17,  18,  19,  20,  19,  20,  21,  22,  23],
         [ 24,  25,  26,  27,  28,  27,  28,  29,  30,  31],
         [ 32,  33,  34,  35,  36,  35,  36,  37,  38,  39],
         [ 40,  41,  42,  43,  44,  43,  44,  45,  46,  47],
         [ 48,  49,  50,  51,  52,  51,  52,  53,  54,  55],
         [ 56,  57,  58,  59,  60,  59,  60,  61,  62,  63],
         [100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
         [100, 100, 100, 100, 100, 100, 100, 100, 100, 100]])
    assert_eq(g, expected)
    assert g.chunks == ((8, 8), (5, 5))


def test_map_overlap():
    x = da.arange(10, chunks=5)
    y = x.map_overlap(lambda x: x + len(x), depth=2, dtype=x.dtype)
    assert_eq(y, np.arange(10) + 5 + 2 + 2)

    x = da.arange(10, chunks=5)
    y = x.map_overlap(lambda x: x + len(x), depth=np.int64(2), dtype=x.dtype)
    assert all([(type(s) is int) for s in y.shape])
    assert_eq(y, np.arange(10) + 5 + 2 + 2)

    x = np.arange(16).reshape((4, 4))
    d = da.from_array(x, chunks=(2, 2))
    exp1 = d.map_overlap(lambda x: x + x.size, depth=1, dtype=d.dtype)
    exp2 = d.map_overlap(lambda x: x + x.size, depth={0: 1, 1: 1},
                         boundary={0: 'reflect', 1: 'none'}, dtype=d.dtype)
    exp3 = d.map_overlap(lambda x: x + x.size, depth={1: 1},
                         boundary={1: 'reflect'}, dtype=d.dtype)
    assert_eq(exp1, x + 16)
    assert_eq(exp2, x + 12)
    assert_eq(exp3, x + 8)


@pytest.mark.parametrize("boundary", [
    None,
    "reflect",
    "periodic",
    "nearest",
    "none",
    0
])
def test_map_overlap_no_depth(boundary):
    x = da.arange(10, chunks=5)

    y = x.map_overlap(lambda i: i, depth=0, boundary=boundary, dtype=x.dtype)
    assert_eq(y, x)


def test_nearest_ghost():
    a = np.arange(144).reshape(12, 12).astype(float)

    darr = da.from_array(a, chunks=(6, 6))
    garr = ghost(darr, depth={0: 5, 1: 5},
                 boundary={0: 'nearest', 1: 'nearest'})
    tarr = trim_internal(garr, {0: 5, 1: 5})
    assert_array_almost_equal(tarr, a)


def test_0_depth():
    expected = np.arange(100).reshape(10, 10)
    darr = da.from_array(expected, chunks=(5, 2))

    depth = {0: 0, 1: 0}

    reflected = ghost(darr, depth=depth, boundary='reflect')
    nearest = ghost(darr, depth=depth, boundary='nearest')
    periodic = ghost(darr, depth=depth, boundary='periodic')
    constant = ghost(darr, depth=depth, boundary=42)

    result = trim_internal(reflected, depth)
    assert_array_equal(result, expected)

    result = trim_internal(nearest, depth)
    assert_array_equal(result, expected)

    result = trim_internal(periodic, depth)
    assert_array_equal(result, expected)

    result = trim_internal(constant, depth)
    assert_array_equal(result, expected)


def test_some_0_depth():
    expected = np.arange(100).reshape(10, 10)
    darr = da.from_array(expected, chunks=(5, 5))

    depth = {0: 4, 1: 0}

    reflected = ghost(darr, depth=depth, boundary='reflect')
    nearest = ghost(darr, depth=depth, boundary='nearest')
    periodic = ghost(darr, depth=depth, boundary='periodic')
    constant = ghost(darr, depth=depth, boundary=42)

    result = trim_internal(reflected, depth)
    assert_array_equal(result, expected)

    result = trim_internal(nearest, depth)
    assert_array_equal(result, expected)

    result = trim_internal(periodic, depth)
    assert_array_equal(result, expected)

    result = trim_internal(constant, depth)
    assert_array_equal(result, expected)


def test_one_chunk_along_axis():
    a = np.arange(2 * 9).reshape(2, 9)
    darr = da.from_array(a, chunks=((2,), (2, 2, 2, 3)))
    g = ghost(darr, depth=0, boundary=0)
    assert a.shape == g.shape


def test_constant_boundaries():
    a = np.arange(1 * 9).reshape(1, 9)
    darr = da.from_array(a, chunks=((1,), (2, 2, 2, 3)))
    b = boundaries(darr, {0: 0, 1: 0}, {0: 0, 1: 0})
    assert b.chunks == darr.chunks


def test_depth_equals_boundary_length():
    expected = np.arange(100).reshape(10, 10)
    darr = da.from_array(expected, chunks=(5, 5))

    depth = {0: 5, 1: 5}

    reflected = ghost(darr, depth=depth, boundary='reflect')
    nearest = ghost(darr, depth=depth, boundary='nearest')
    periodic = ghost(darr, depth=depth, boundary='periodic')
    constant = ghost(darr, depth=depth, boundary=42)

    result = trim_internal(reflected, depth)
    assert_array_equal(result, expected)

    result = trim_internal(nearest, depth)
    assert_array_equal(result, expected)

    result = trim_internal(periodic, depth)
    assert_array_equal(result, expected)

    result = trim_internal(constant, depth)
    assert_array_equal(result, expected)


@pytest.mark.xfail
def test_depth_greater_than_boundary_length():
    expected = np.arange(100).reshape(10, 10)
    darr = da.from_array(expected, chunks=(5, 5))

    depth = {0: 8, 1: 7}

    reflected = ghost(darr, depth=depth, boundary='reflect')
    nearest = ghost(darr, depth=depth, boundary='nearest')
    periodic = ghost(darr, depth=depth, boundary='periodic')
    constant = ghost(darr, depth=depth, boundary=42)

    result = trim_internal(reflected, depth)
    assert_array_equal(result, expected)

    result = trim_internal(nearest, depth)
    assert_array_equal(result, expected)

    result = trim_internal(periodic, depth)
    assert_array_equal(result, expected)

    result = trim_internal(constant, depth)
    assert_array_equal(result, expected)


def test_bad_depth_raises():
    expected = np.arange(144).reshape(12, 12)
    darr = da.from_array(expected, chunks=(5, 5))

    depth = {0: 4, 1: 2}

    pytest.raises(ValueError, ghost, darr, depth=depth, boundary=1)


def test_none_boundaries():
    x = da.from_array(np.arange(16).reshape(4, 4), chunks=(2, 2))
    exp = boundaries(x, 2, {0: 'none', 1: 33})
    res = np.array(
        [[33, 33,  0,  1,  2,  3, 33, 33],
         [33, 33,  4,  5,  6,  7, 33, 33],
         [33, 33,  8,  9, 10, 11, 33, 33],
         [33, 33, 12, 13, 14, 15, 33, 33]])
    assert_eq(exp, res)


def test_ghost_small():
    x = da.ones((10, 10), chunks=(5, 5))

    y = x.map_overlap(lambda x: x, depth=1)
    assert len(y.dask) < 200

    y = x.map_overlap(lambda x: x, depth=1, boundary='none')
    assert len(y.dask) < 100
