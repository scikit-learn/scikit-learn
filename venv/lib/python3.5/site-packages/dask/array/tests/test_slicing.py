import itertools
from operator import getitem

import pytest
from toolz import merge

np = pytest.importorskip('numpy')

import dask
import dask.array as da
from dask.array.slicing import (_sanitize_index_element, _slice_1d,
                                new_blockdim, sanitize_index, slice_array,
                                take, normalize_index)
from dask.array.utils import assert_eq, same_keys


def test_slice_1d():
    expected = {0: slice(10, 25, 1), 1: slice(None, None, None), 2: slice(0, 1, 1)}
    result = _slice_1d(100, [25] * 4, slice(10, 51, None))
    assert expected == result

    # x[100:12:-3]
    expected = {0: slice(-2, -8, -3),
                1: slice(-1, -21, -3),
                2: slice(-3, -21, -3),
                3: slice(-2, -21, -3),
                4: slice(-1, -21, -3)}
    result = _slice_1d(100, [20] * 5, slice(100, 12, -3))
    assert expected == result

    # x[102::-3]
    expected = {0: slice(-2, -21, -3),
                1: slice(-1, -21, -3),
                2: slice(-3, -21, -3),
                3: slice(-2, -21, -3),
                4: slice(-1, -21, -3)}
    result = _slice_1d(100, [20] * 5, slice(102, None, -3))
    assert expected == result

    # x[::-4]
    expected = {0: slice(-1, -21, -4),
                1: slice(-1, -21, -4),
                2: slice(-1, -21, -4),
                3: slice(-1, -21, -4),
                4: slice(-1, -21, -4)}
    result = _slice_1d(100, [20] * 5, slice(None, None, -4))
    assert expected == result

    # x[::-7]
    expected = {0: slice(-5, -21, -7),
                1: slice(-4, -21, -7),
                2: slice(-3, -21, -7),
                3: slice(-2, -21, -7),
                4: slice(-1, -21, -7)}
    result = _slice_1d(100, [20] * 5, slice(None, None, -7))
    assert expected == result

    # x=range(115)
    # x[::-7]
    expected = {0: slice(-7, -24, -7),
                1: slice(-2, -24, -7),
                2: slice(-4, -24, -7),
                3: slice(-6, -24, -7),
                4: slice(-1, -24, -7)}
    result = _slice_1d(115, [23] * 5, slice(None, None, -7))
    assert expected == result

    # x[79::-3]
    expected = {0: slice(-1, -21, -3),
                1: slice(-3, -21, -3),
                2: slice(-2, -21, -3),
                3: slice(-1, -21, -3)}
    result = _slice_1d(100, [20] * 5, slice(79, None, -3))
    assert expected == result

    # x[-1:-8:-1]
    expected = {4: slice(-1, -8, -1)}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(-1, 92, -1))
    assert expected == result

    # x[20:0:-1]
    expected = {0: slice(-1, -20, -1),
                1: slice(-20, -21, -1)}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(20, 0, -1))
    assert expected == result

    # x[:0]
    expected = {}
    result = _slice_1d(100, [20, 20, 20, 20, 20], slice(0))
    assert result

    # x=range(99)
    expected = {0: slice(-3, -21, -3),
                1: slice(-2, -21, -3),
                2: slice(-1, -21, -3),
                3: slice(-2, -20, -3),
                4: slice(-1, -21, -3)}
    # This array has non-uniformly sized blocks
    result = _slice_1d(99, [20, 20, 20, 19, 20], slice(100, None, -3))
    assert expected == result

    # x=range(104)
    # x[::-3]
    expected = {0: slice(-1, -21, -3),
                1: slice(-3, -24, -3),
                2: slice(-3, -28, -3),
                3: slice(-1, -14, -3),
                4: slice(-1, -22, -3)}
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(None, None, -3))
    assert expected == result

    # x=range(104)
    # x[:27:-3]
    expected = {1: slice(-3, -16, -3),
                2: slice(-3, -28, -3),
                3: slice(-1, -14, -3),
                4: slice(-1, -22, -3)}
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(None, 27, -3))
    assert expected == result

    # x=range(104)
    # x[100:27:-3]
    expected = {1: slice(-3, -16, -3),
                2: slice(-3, -28, -3),
                3: slice(-1, -14, -3),
                4: slice(-4, -22, -3)}
    # This array has non-uniformly sized blocks
    result = _slice_1d(104, [20, 23, 27, 13, 21], slice(100, 27, -3))
    assert expected == result


def test_slice_singleton_value_on_boundary():
    assert _slice_1d(15, [5, 5, 5], 10) == {2: 0}
    assert _slice_1d(30, (5, 5, 5, 5, 5, 5), 10) == {2: 0}


def test_slice_array_1d():
    #x[24::2]
    expected = {('y', 0): (getitem, ('x', 0), (slice(24, 25, 2),)),
                ('y', 1): (getitem, ('x', 1), (slice(1, 25, 2),)),
                ('y', 2): (getitem, ('x', 2), (slice(0, 25, 2),)),
                ('y', 3): (getitem, ('x', 3), (slice(1, 25, 2),))}
    result, chunks = slice_array('y', 'x', [[25] * 4], [slice(24, None, 2)])

    assert expected == result

    #x[26::2]
    expected = {('y', 0): (getitem, ('x', 1), (slice(1, 25, 2),)),
                ('y', 1): (getitem, ('x', 2), (slice(0, 25, 2),)),
                ('y', 2): (getitem, ('x', 3), (slice(1, 25, 2),))}

    result, chunks = slice_array('y', 'x', [[25] * 4], [slice(26, None, 2)])
    assert expected == result

    #x[24::2]
    expected = {('y', 0): (getitem, ('x', 0), (slice(24, 25, 2),)),
                ('y', 1): (getitem, ('x', 1), (slice(1, 25, 2),)),
                ('y', 2): (getitem, ('x', 2), (slice(0, 25, 2),)),
                ('y', 3): (getitem, ('x', 3), (slice(1, 25, 2),))}
    result, chunks = slice_array('y', 'x', [(25, ) * 4], (slice(24, None, 2), ))

    assert expected == result

    #x[26::2]
    expected = {('y', 0): (getitem, ('x', 1), (slice(1, 25, 2),)),
                ('y', 1): (getitem, ('x', 2), (slice(0, 25, 2),)),
                ('y', 2): (getitem, ('x', 3), (slice(1, 25, 2),))}

    result, chunks = slice_array('y', 'x', [(25, ) * 4], (slice(26, None, 2), ))
    assert expected == result


def test_slice_array_2d():
    #2d slices: x[13::2,10::1]
    expected = {('y', 0, 0): (getitem, ('x', 0, 0),
                              (slice(13, 20, 2), slice(10, 20, 1))),
                ('y', 0, 1): (getitem, ('x', 0, 1),
                              (slice(13, 20, 2), slice(None, None, None))),
                ('y', 0, 2): (getitem, ('x', 0, 2),
                              (slice(13, 20, 2), slice(None, None, None)))}

    result, chunks = slice_array('y', 'x', [[20], [20, 20, 5]],
                                 [slice(13, None, 2), slice(10, None, 1)])

    assert expected == result

    #2d slices with one dimension: x[5,10::1]
    expected = {('y', 0): (getitem, ('x', 0, 0),
                           (5, slice(10, 20, 1))),
                ('y', 1): (getitem, ('x', 0, 1),
                           (5, slice(None, None, None))),
                ('y', 2): (getitem, ('x', 0, 2),
                           (5, slice(None, None, None)))}

    result, chunks = slice_array('y', 'x', ([20], [20, 20, 5]),
                                 [5, slice(10, None, 1)])

    assert expected == result


def test_slice_optimizations():
    #bar[:]
    expected = {('foo', 0): ('bar', 0)}
    result, chunks = slice_array('foo', 'bar', [[100]], (slice(None, None, None),))
    assert expected == result

    #bar[:,:,:]
    expected = {('foo', 0): ('bar', 0),
                ('foo', 1): ('bar', 1),
                ('foo', 2): ('bar', 2)}
    result, chunks = slice_array('foo', 'bar', [(100, 1000, 10000)],
                                 (slice(None, None, None),
                                  slice(None, None, None),
                                  slice(None, None, None)))
    assert expected == result


def test_slicing_with_singleton_indices():
    result, chunks = slice_array('y', 'x', ([5, 5], [5, 5]), (slice(0, 5), 8))

    expected = {('y', 0): (getitem, ('x', 0, 1), (slice(None, None, None), 3))}

    assert expected == result


def test_slicing_with_newaxis():
    result, chunks = slice_array('y', 'x', ([5, 5], [5, 5]),
                                 (slice(0, 3), None, slice(None, None, None)))

    expected = {
        ('y', 0, 0, 0): (getitem, ('x', 0, 0),
                         (slice(0, 3, 1), None, slice(None, None, None))),
        ('y', 0, 0, 1): (getitem, ('x', 0, 1),
                         (slice(0, 3, 1), None, slice(None, None, None)))}

    assert expected == result
    assert chunks == ((3,), (1,), (5, 5))


def test_take():
    chunks, dsk = take('y', 'x', [(20, 20, 20, 20)], [5, 1, 47, 3], axis=0)
    expected = {('y', 0): (getitem, (np.concatenate,
                                     [(getitem, ('x', 0), (np.array([1, 3, 5]),)),
                                      (getitem, ('x', 2), (np.array([7]),))], 0),
                           (np.array([2, 0, 3, 1]), ))}
    np.testing.assert_equal(sorted(dsk.items()), sorted(expected.items()))
    assert chunks == ((4,),)

    chunks, dsk = take('y', 'x', [(20, 20, 20, 20), (20, 20)], [
                       5, 1, 47, 3], axis=0)
    expected = {('y', 0, j): (getitem, (np.concatenate,
                                        [(getitem, ('x', 0, j),
                                          ([1, 3, 5], slice(None, None, None))),
                                         (getitem, ('x', 2, j),
                                          ([7], slice(None, None, None)))], 0),
                              ([2, 0, 3, 1], slice(None, None, None)))
                for j in range(2)}
    np.testing.assert_equal(sorted(dsk.items()), sorted(expected.items()))
    assert chunks == ((4,), (20, 20))

    chunks, dsk = take('y', 'x', [(20, 20, 20, 20), (20, 20)], [
                       5, 1, 37, 3], axis=1)
    expected = {('y', i, 0): (getitem, (np.concatenate,
                                        [(getitem, ('x', i, 0),
                                          (slice(None, None, None), [1, 3, 5])),
                                         (getitem, ('x', i, 1),
                                          (slice(None, None, None), [17]))], 1),
                              (slice(None, None, None), [2, 0, 3, 1]))
                for i in range(4)}
    np.testing.assert_equal(sorted(dsk.items()), sorted(expected.items()))
    assert chunks == ((20, 20, 20, 20), (4,))


def test_take_sorted():
    chunks, dsk = take('y', 'x', [(20, 20, 20, 20)], [1, 3, 5, 47], axis=0)
    expected = {('y', 0): (getitem, ('x', 0), ([1, 3, 5],)),
                ('y', 1): (getitem, ('x', 2), ([7],))}
    np.testing.assert_equal(dsk, expected)
    assert chunks == ((3, 1),)

    chunks, dsk = take('y', 'x', [(20, 20, 20, 20), (20, 20)], [1, 3, 5, 37], axis=1)
    expected = merge(dict((('y', i, 0), (getitem, ('x', i, 0),
                                         (slice(None, None, None), [1, 3, 5])))
                          for i in range(4)),
                     dict((('y', i, 1), (getitem, ('x', i, 1),
                                         (slice(None, None, None), [17])))
                          for i in range(4)))
    np.testing.assert_equal(dsk, expected)
    assert chunks == ((20, 20, 20, 20), (3, 1))


def test_slice_lists():
    y, chunks = slice_array('y', 'x', ((3, 3, 3, 1), (3, 3, 3, 1)),
                            (np.array([2, 1, 9]), slice(None, None, None)))
    exp = {('y', 0, i): (getitem, (np.concatenate,
                                   [(getitem, ('x', 0, i),
                                     ([1, 2], slice(None, None, None))),
                                    (getitem, ('x', 3, i),
                                     ([0], slice(None, None, None)))], 0),
                         ([1, 0, 2], slice(None, None, None)))
           for i in range(4)}
    np.testing.assert_equal(y, exp)
    assert chunks == ((3,), (3, 3, 3, 1))


def test_slicing_chunks():
    result, chunks = slice_array('y', 'x', ([5, 5], [5, 5]),
                                 (1, np.array([2, 0, 3])))
    assert chunks == ((3,), )

    result, chunks = slice_array('y', 'x', ([5, 5], [5, 5]),
                                 (slice(0, 7), np.array([2, 0, 3])))
    assert chunks == ((5, 2), (3, ))

    result, chunks = slice_array('y', 'x', ([5, 5], [5, 5]),
                                 (slice(0, 7), 1))
    assert chunks == ((5, 2), )


def test_slicing_with_numpy_arrays():
    a, bd1 = slice_array('y', 'x', ((3, 3, 3, 1), (3, 3, 3, 1)),
                         (np.array([1, 2, 9]), slice(None, None, None)))
    b, bd2 = slice_array('y', 'x', ((3, 3, 3, 1), (3, 3, 3, 1)),
                         (np.array([1, 2, 9]), slice(None, None, None)))

    assert bd1 == bd2
    np.testing.assert_equal(a, b)

    i = [False, True, True, False, False,
         False, False, False, False, True, False]
    index = (i, slice(None, None, None))
    index = normalize_index(index, (10, 10))
    c, bd3 = slice_array('y', 'x', ((3, 3, 3, 1), (3, 3, 3, 1)), index)
    assert bd1 == bd3
    np.testing.assert_equal(a, c)


def test_slicing_and_chunks():
    o = da.ones((24, 16), chunks=((4, 8, 8, 4), (2, 6, 6, 2)))
    t = o[4:-4, 2:-2]
    assert t.chunks == ((8, 8), (6, 6))


def test_slicing_identities():
    a = da.ones((24, 16), chunks=((4, 8, 8, 4), (2, 6, 6, 2)))

    assert a is a[slice(None)]
    assert a is a[:]
    assert a is a[::]
    assert a is a[...]
    assert a is a[0:]
    assert a is a[0::]
    assert a is a[::1]
    assert a is a[0:len(a)]
    assert a is a[0::1]
    assert a is a[0:len(a):1]


def test_slice_stop_0():
    # from gh-125
    a = da.ones(10, chunks=(10,))[:0].compute()
    b = np.ones(10)[:0]
    assert_eq(a, b)


def test_slice_list_then_None():
    x = da.zeros(shape=(5, 5), chunks=(3, 3))
    y = x[[2, 1]][None]

    assert_eq(y, np.zeros((1, 2, 5)))


class ReturnItem(object):

    def __getitem__(self, key):
        return key


@pytest.mark.skip(reason='really long test')
def test_slicing_exhaustively():
    x = np.random.rand(6, 7, 8)
    a = da.from_array(x, chunks=(3, 3, 3))
    I = ReturnItem()

    # independent indexing along different axes
    indexers = [0, -2, I[:], I[:5], [0, 1], [0, 1, 2], [4, 2], I[::-1], None, I[:0], []]
    for i in indexers:
        assert_eq(x[i], a[i]), i
        for j in indexers:
            assert_eq(x[i][:, j], a[i][:, j]), (i, j)
            assert_eq(x[:, i][j], a[:, i][j]), (i, j)
            for k in indexers:
                assert_eq(x[..., i][:, j][k], a[..., i][:, j][k]), (i, j, k)

    # repeated indexing along the first axis
    first_indexers = [I[:], I[:5], np.arange(5), [3, 1, 4, 5, 0], np.arange(6) < 6]
    second_indexers = [0, -1, 3, I[:], I[:3], I[2:-1], [2, 4], [], I[:0]]
    for i in first_indexers:
        for j in second_indexers:
            assert_eq(x[i][j], a[i][j]), (i, j)


def test_slicing_with_negative_step_flops_keys():
    x = da.arange(10, chunks=5)
    y = x[:1:-1]
    assert (x.name, 1) in y.dask[(y.name, 0)]
    assert (x.name, 0) in y.dask[(y.name, 1)]

    assert_eq(y, np.arange(10)[:1:-1])

    assert y.chunks == ((5, 3),)

    assert y.dask[(y.name, 0)] == (getitem, (x.name, 1),
                                            (slice(-1, -6, -1),))
    assert y.dask[(y.name, 1)] == (getitem, (x.name, 0),
                                            (slice(-1, -4, -1),))


def test_empty_slice():
    x = da.ones((5, 5), chunks=(2, 2), dtype='i4')
    y = x[:0]

    assert_eq(y, np.ones((5, 5), dtype='i4')[:0])


def test_multiple_list_slicing():
    x = np.random.rand(6, 7, 8)
    a = da.from_array(x, chunks=(3, 3, 3))
    assert_eq(x[:, [0, 1, 2]][[0, 1]], a[:, [0, 1, 2]][[0, 1]])


def test_empty_list():
    x = np.ones((5, 5, 5), dtype='i4')
    dx = da.from_array(x, chunks=2)

    assert_eq(dx[[], :3, :2], x[[], :3, :2])
    assert_eq(dx[:3, [], :2], x[:3, [], :2])
    assert_eq(dx[:3, :2, []], x[:3, :2, []])


def test_uneven_chunks():
    assert da.ones(20, chunks=5)[::2].chunks == ((3, 2, 3, 2),)


def test_new_blockdim():
    assert new_blockdim(20, [5, 5, 5, 5], slice(0, None, 2)) == [3, 2, 3, 2]


def test_slicing_consistent_names():
    x = np.arange(100).reshape((10, 10))
    a = da.from_array(x, chunks=(5, 5))
    assert same_keys(a[0], a[0])
    assert same_keys(a[:, [1, 2, 3]], a[:, [1, 2, 3]])
    assert same_keys(a[:, 5:2:-1], a[:, 5:2:-1])
    assert same_keys(a[0, ...], a[0, ...])
    assert same_keys(a[...], a[...])
    assert same_keys(a[[1, 3, 5]], a[[1, 3, 5]])
    assert same_keys(a[-11:11], a[:])
    assert same_keys(a[-11:-9], a[:1])
    assert same_keys(a[-1], a[9])


def test_slicing_consistent_names_after_normalization():
    x = da.zeros(10, chunks=(5,))
    assert same_keys(x[0:], x[:10])
    assert same_keys(x[0:], x[0:10])
    assert same_keys(x[0:], x[0:10:1])
    assert same_keys(x[:], x[0:10:1])


def test_sanitize_index_element():
    with pytest.raises(TypeError):
        _sanitize_index_element('Hello!')


def test_sanitize_index():
    pd = pytest.importorskip('pandas')
    with pytest.raises(TypeError):
        sanitize_index('Hello!')

    np.testing.assert_equal(sanitize_index(pd.Series([1, 2, 3])), [1, 2, 3])
    np.testing.assert_equal(sanitize_index((1, 2, 3)), [1, 2, 3])


def test_uneven_blockdims():
    blockdims = ((31, 28, 31,  30,  31,  30,  31,  31,  30,  31, 30), (100,))
    index = (slice(240, 270), slice(None))
    dsk_out, bd_out = slice_array('in', 'out', blockdims, index)
    sol = {('in', 0, 0): (getitem, ('out', 7, 0), (slice(28, 31, 1), slice(None))),
           ('in', 1, 0): (getitem, ('out', 8, 0), (slice(0, 27, 1), slice(None)))}
    assert dsk_out == sol
    assert bd_out == ((3, 27), (100,))

    blockdims = ((31, 28, 31,  30,  31,  30,  31,  31,  30,  31, 30),) * 2
    index = (slice(240, 270), slice(180, 230))
    dsk_out, bd_out = slice_array('in', 'out', blockdims, index)
    sol = {('in', 0, 0): (getitem, ('out', 7, 5), (slice(28, 31, 1), slice(29, 30, 1))),
           ('in', 0, 1): (getitem, ('out', 7, 6), (slice(28, 31, 1), slice(None))),
           ('in', 0, 2): (getitem, ('out', 7, 7), (slice(28, 31, 1), slice(0, 18, 1))),
           ('in', 1, 0): (getitem, ('out', 8, 5), (slice(0, 27, 1), slice(29, 30, 1))),
           ('in', 1, 1): (getitem, ('out', 8, 6), (slice(0, 27, 1), slice(None))),
           ('in', 1, 2): (getitem, ('out', 8, 7), (slice(0, 27, 1), slice(0, 18, 1)))}
    assert dsk_out == sol
    assert bd_out == ((3, 27), (1, 31, 18))


def test_oob_check():
    x = da.ones(5, chunks=(2,))
    with pytest.raises(IndexError):
        x[6]
    with pytest.raises(IndexError):
        x[[6]]
    with pytest.raises(IndexError):
        x[-10]
    with pytest.raises(IndexError):
        x[[-10]]
    with pytest.raises(IndexError):
        x[0, 0]


def test_index_with_dask_array():
    x = np.arange(36).reshape((6, 6))
    d = da.from_array(x, chunks=(3, 3))
    ind = np.asarray([True, True, False, True, False, False], dtype=bool)
    ind = da.from_array(ind, chunks=2)
    for index in [ind, (slice(1, 9, 2), ind), (ind, slice(2, 8, 1))]:
        x_index = dask.compute(index)[0]
        assert_eq(x[x_index], d[index])


def test_index_with_dask_array_2():
    x = np.random.random((10, 10, 10))
    ind = np.random.random(10) > 0.5

    d = da.from_array(x, chunks=(3, 4, 5))
    dind = da.from_array(ind, chunks=4)

    index = [slice(1, 9, 1), slice(None)]

    for i in range(x.ndim):
        index2 = index[:]
        index2.insert(i, dind)

        index3 = index[:]
        index3.insert(i, ind)

        assert_eq(x[tuple(index3)], d[tuple(index2)])


@pytest.mark.xfail
def test_cull():
    x = da.ones(1000, chunks=(10,))

    for slc in [1, slice(0, 30), slice(0, None, 100)]:
        y = x[slc]
        assert len(y.dask) < len(x.dask)


@pytest.mark.parametrize('shape', [(2,), (2, 3), (2, 3, 5)])
@pytest.mark.parametrize('index', [(Ellipsis,),
                                   (None, Ellipsis),
                                   (Ellipsis, None),
                                   (None, Ellipsis, None)])
def test_slicing_with_Nones(shape, index):
    x = np.random.random(shape)
    d = da.from_array(x, chunks=shape)

    assert_eq(x[index], d[index])


indexers = [Ellipsis, slice(2), 0, 1, -2, -1, slice(-2, None), None]


"""
# We comment this out because it is 4096 tests
@pytest.mark.parametrize('a', indexers)
@pytest.mark.parametrize('b', indexers)
@pytest.mark.parametrize('c', indexers)
@pytest.mark.parametrize('d', indexers)
def test_slicing_none_int_ellipses(a, b, c, d):
    if (a, b, c, d).count(Ellipsis) > 1:
        return
    shape = (2,3,5,7,11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)

    xx = x[a, b, c, d]
    yy = y[a, b, c, d]
    assert_eq(xx, yy)
"""


def test_slicing_integer_no_warnings():
    # https://github.com/dask/dask/pull/2457/
    X = da.random.random((100, 2), (2, 2))
    idx = np.array([0, 0, 1, 1])
    with pytest.warns(None) as rec:
        X[idx].compute()
    assert len(rec) == 0


@pytest.mark.slow
def test_slicing_none_int_ellipes():
    shape = (2, 3, 5, 7, 11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)
    for ind in itertools.product(indexers, indexers, indexers, indexers):
        if ind.count(Ellipsis) > 1:
            continue

        assert_eq(x[ind], y[ind])


def test_None_overlap_int():
    a, b, c, d = (0, slice(None, 2, None), None, Ellipsis)
    shape = (2, 3, 5, 7, 11)
    x = np.arange(np.prod(shape)).reshape(shape)
    y = da.core.asarray(x)

    xx = x[a, b, c, d]
    yy = y[a, b, c, d]
    assert_eq(xx, yy)


def test_negative_n_slicing():
    assert_eq(da.ones(2, chunks=2)[-2], np.ones(2)[-2])


def test_negative_list_slicing():
    x = np.arange(5)
    dx = da.from_array(x, chunks=2)
    assert_eq(dx[[0, -5]], x[[0, -5]])
    assert_eq(dx[[4, -1]], x[[4, -1]])


def test_permit_oob_slices():
    x = np.arange(5)
    dx = da.from_array(x, chunks=2)

    assert_eq(x[-102:], dx[-102:])
    assert_eq(x[102:], dx[102:])
    assert_eq(x[:102], dx[:102])
    assert_eq(x[:-102], dx[:-102])


def test_normalize_index():
    assert normalize_index((Ellipsis, None), (10,)) == (slice(None), None)
    assert normalize_index(5, (np.nan,)) == (5,)
    assert normalize_index(-5, (np.nan,)) == (-5,)
    (result,) = normalize_index([-5, -2, 1], (np.nan,))
    assert result.tolist() == [-5, -2, 1]
    assert normalize_index(slice(-5, -2), (np.nan,)) == (slice(-5, -2),)
