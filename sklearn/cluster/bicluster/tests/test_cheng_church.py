"""Testing for Spectral Biclustering methods"""

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises

from sklearn.metrics import consensus_score
from sklearn.datasets import make_msr

from sklearn.cluster.bicluster.cheng_church import ChengChurch
from sklearn.cluster.bicluster.cheng_church import EmptyBiclusterException


def test_cheng_church():
    """Test Cheng and Church algorithm on a simple problem."""
    for shape in ((150, 150), (50, 50)):
        for noise in (0, 5):
            for deletion_threshold in (1.5, 2):
                data, rows, cols = make_msr(shape, 3, random_state=0)
                model = ChengChurch(n_clusters=3, max_msr=max(noise * 2, 1),
                                    deletion_threshold=deletion_threshold,
                                    random_state=0)
                model.fit(data)
                assert(consensus_score((rows, cols), model.biclusters_) > 0.9)


def test_empty_biclusters():
    """Cheng and Church should always find at least one bicluster.

    The MSR of a bicluster with one row or one column is zero.

    """
    for i in range(10):
        generator = np.random.RandomState(i)
        data = generator.uniform(0, 1000, (50, 50))
        model = ChengChurch(n_clusters=1, max_msr=0)
        model.fit(data)
        assert_equal(len(model.rows_), 1)
        assert_equal(len(model.columns_), 1)


def test_single_node_deletion():
    generator = np.random.RandomState(0)
    model = ChengChurch()
    arr = generator.uniform(20, 100, (20, 20))
    arr[5:15, 5:15] = 10
    rows = np.arange(20)[np.newaxis].T
    cols = np.arange(20)
    rows, cols = model._single_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(5, 15)[np.newaxis].T)
    assert_array_equal(cols, np.arange(5, 15))


def test_multiple_node_deletion():
    generator = np.random.RandomState(0)
    model = ChengChurch(max_msr=250,
                        deletion_threshold=1.25,
                        row_deletion_cutoff=1,
                        column_deletion_cutoff=1)
    arr = generator.uniform(20, 100, (20, 20))
    arr[3:18, 3:18] = 10
    rows = np.arange(20)[np.newaxis].T
    cols = np.arange(20)
    rows, cols = model._multiple_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(3, 18)[np.newaxis].T)
    assert_array_equal(cols, np.arange(3, 18))


def test_node_deletion_when_unnecessary():
    model = ChengChurch()
    arr = np.zeros((20, 20))
    rows = np.arange(20)[np.newaxis].T
    cols = np.arange(20)
    rows, cols = model._single_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(20)[np.newaxis].T)
    assert_array_equal(cols, np.arange(20))

    rows, cols = model._multiple_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(20)[np.newaxis].T)
    assert_array_equal(cols, np.arange(20))


def test_node_addition():
    generator = np.random.RandomState(0)
    model = ChengChurch()
    arr = np.zeros((20, 20)) + 15
    arr[5:10, 5:10] = generator.uniform(10, 20, (5, 5))
    rows = np.arange(5, 10)[np.newaxis].T
    cols = np.arange(5, 10)
    rows, cols = model._node_addition(rows, cols, arr)
    assert_array_equal(rows, np.arange(20)[np.newaxis].T)
    assert_array_equal(cols, np.arange(20))


def test_row_msr():
    X = np.arange(9).reshape(3, 3)
    rows = np.arange(2)[np.newaxis].T
    cols = np.arange(2)
    model = ChengChurch()
    vec = model._row_msr(rows, cols, X)
    expected = np.zeros(3)
    assert_array_equal(vec, expected)


def test_col_msr():
    X = np.arange(9).reshape(3, 3)
    rows = np.arange(2)[np.newaxis].T
    cols = np.arange(2)
    model = ChengChurch()
    vec = model._col_msr(rows, cols, X)
    expected = np.zeros(3)
    assert_array_equal(vec, expected)


def test_errors():
    data = np.arange(25).reshape((5, 5))

    model = ChengChurch(n_clusters=0)
    assert_raises(ValueError, model.fit, data)

    model = ChengChurch(max_msr=-1)
    assert_raises(ValueError, model.fit, data)

    model = ChengChurch(deletion_threshold=0.9)
    assert_raises(ValueError, model.fit, data)

    model = ChengChurch(row_deletion_cutoff=0)
    assert_raises(ValueError, model.fit, data)

    model = ChengChurch(column_deletion_cutoff=0)
    assert_raises(ValueError, model.fit, data)

    assert_raises(EmptyBiclusterException, model._all_msr,
                  np.array([]), np.array([]), data)

    assert_raises(EmptyBiclusterException, model._row_msr,
                  np.array([]), np.array([]), data)

    assert_raises(EmptyBiclusterException, model._col_msr,
                  np.array([]), np.array([]), data)
