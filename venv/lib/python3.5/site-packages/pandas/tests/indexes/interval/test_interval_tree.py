from __future__ import division

import pytest
import numpy as np
from pandas import compat
from pandas._libs.interval import IntervalTree
import pandas.util.testing as tm


@pytest.fixture(scope='class', params=['left', 'right', 'both', 'neither'])
def closed(request):
    return request.param


@pytest.fixture(
    scope='class', params=['int32', 'int64', 'float32', 'float64', 'uint64'])
def dtype(request):
    return request.param


@pytest.fixture(scope='class')
def tree(dtype):
    left = np.arange(5, dtype=dtype)
    return IntervalTree(left, left + 2)


class TestIntervalTree(object):

    def test_get_loc(self, tree):
        tm.assert_numpy_array_equal(tree.get_loc(1),
                                    np.array([0], dtype='int64'))
        tm.assert_numpy_array_equal(np.sort(tree.get_loc(2)),
                                    np.array([0, 1], dtype='int64'))
        with pytest.raises(KeyError):
            tree.get_loc(-1)

    def test_get_indexer(self, tree):
        tm.assert_numpy_array_equal(
            tree.get_indexer(np.array([1.0, 5.5, 6.5])),
            np.array([0, 4, -1], dtype='int64'))
        with pytest.raises(KeyError):
            tree.get_indexer(np.array([3.0]))

    def test_get_indexer_non_unique(self, tree):
        indexer, missing = tree.get_indexer_non_unique(
            np.array([1.0, 2.0, 6.5]))
        tm.assert_numpy_array_equal(indexer[:1],
                                    np.array([0], dtype='int64'))
        tm.assert_numpy_array_equal(np.sort(indexer[1:3]),
                                    np.array([0, 1], dtype='int64'))
        tm.assert_numpy_array_equal(np.sort(indexer[3:]),
                                    np.array([-1], dtype='int64'))
        tm.assert_numpy_array_equal(missing, np.array([2], dtype='int64'))

    def test_duplicates(self, dtype):
        left = np.array([0, 0, 0], dtype=dtype)
        tree = IntervalTree(left, left + 1)
        tm.assert_numpy_array_equal(np.sort(tree.get_loc(0.5)),
                                    np.array([0, 1, 2], dtype='int64'))

        with pytest.raises(KeyError):
            tree.get_indexer(np.array([0.5]))

        indexer, missing = tree.get_indexer_non_unique(np.array([0.5]))
        tm.assert_numpy_array_equal(np.sort(indexer),
                                    np.array([0, 1, 2], dtype='int64'))
        tm.assert_numpy_array_equal(missing, np.array([], dtype='int64'))

    def test_get_loc_closed(self, closed):
        tree = IntervalTree([0], [1], closed=closed)
        for p, errors in [(0, tree.open_left),
                          (1, tree.open_right)]:
            if errors:
                with pytest.raises(KeyError):
                    tree.get_loc(p)
            else:
                tm.assert_numpy_array_equal(tree.get_loc(p),
                                            np.array([0], dtype='int64'))

    @pytest.mark.skipif(compat.is_platform_32bit(),
                        reason="int type mismatch on 32bit")
    @pytest.mark.parametrize('leaf_size', [1, 10, 100, 10000])
    def test_get_indexer_closed(self, closed, leaf_size):
        x = np.arange(1000, dtype='float64')
        found = x.astype('intp')
        not_found = (-1 * np.ones(1000)).astype('intp')

        tree = IntervalTree(x, x + 0.5, closed=closed, leaf_size=leaf_size)
        tm.assert_numpy_array_equal(found, tree.get_indexer(x + 0.25))

        expected = found if tree.closed_left else not_found
        tm.assert_numpy_array_equal(expected, tree.get_indexer(x + 0.0))

        expected = found if tree.closed_right else not_found
        tm.assert_numpy_array_equal(expected, tree.get_indexer(x + 0.5))
