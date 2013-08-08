"""Testing for Spectral Biclustering methods"""

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.metrics import consensus_score
from sklearn.datasets import make_msr

from sklearn.cluster.bicluster.cheng_church import IncrementalMSR
from sklearn.cluster.bicluster.cheng_church import ChengChurch
from sklearn.cluster.bicluster.cheng_church import EmptyBiclusterException


def test_incremental_msr():
    generator = np.random.RandomState(0)
    data = generator.uniform(1, 100, (20, 20))
    rows = np.ones(20, dtype=np.bool)
    cols = np.ones(20, dtype=np.bool)
    inc = IncrementalMSR(rows, cols, data)

    new_rows = rows.copy()
    new_rows[0] = False

    new_cols = cols.copy()
    new_cols[0] = False

    inc.remove_row(0)
    inc.remove_col(0)

    arr = data[new_rows][:, new_cols]
    sr = arr - arr.mean(axis=1, keepdims=True) - arr.mean(axis=0) + arr.mean()
    sr = np.power(sr, 2)

    assert_almost_equal(inc.msr, sr.mean())
    assert_array_almost_equal(inc.row_msr, sr.mean(axis=1))
    assert_array_almost_equal(inc.col_msr, sr.mean(axis=0))


def test_cheng_church():
    """Test Cheng and Church algorithm on a simple problem."""
    for shape in ((150, 150), (50, 50)):
        for noise in (0, 1):
            for deletion_threshold in (1.5, 2):
                data, rows, cols = make_msr(shape, 3, noise=noise,
                                            random_state=0)
                model = ChengChurch(n_clusters=3, max_msr=10,
                                    deletion_threshold=deletion_threshold,
                                    random_state=0)
                model.fit(data)
                assert(consensus_score((rows, cols), model.biclusters_) > 0.7)


def test_inverse_rows():
    data = np.zeros((15, 10))
    data[:10] = np.arange(100).reshape(10, 10)
    data[10:15] = -np.arange(100, 150).reshape(5, 10)
    rows = np.ones(15, dtype=np.bool)
    cols = np.ones(10, dtype=np.bool)

    model = ChengChurch(n_clusters=1, max_msr=1, random_state=0,
                        inverse_rows=False)
    model.fit(data)
    old_score = consensus_score((rows, cols), model.biclusters_)

    model = ChengChurch(n_clusters=1, max_msr=1, random_state=0,
                        inverse_rows=True)
    model.fit(data)
    new_score = consensus_score((rows, cols), model.biclusters_)
    assert(new_score > old_score)


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
    expected_rows = np.zeros(20, dtype=np.bool)
    expected_cols = np.zeros(20, dtype=np.bool)
    expected_rows[5:15] = True
    expected_cols[5:15] = True

    rows = np.ones(20, dtype=np.bool)
    cols = np.ones(20, dtype=np.bool)
    rows, cols = model._single_node_deletion(rows, cols, arr)

    assert_array_equal(rows, expected_rows)
    assert_array_equal(cols, expected_cols)


def test_multiple_node_deletion():
    generator = np.random.RandomState(0)
    model = ChengChurch(max_msr=250,
                        deletion_threshold=1.25,
                        row_deletion_cutoff=1,
                        column_deletion_cutoff=1)
    arr = generator.uniform(20, 100, (20, 20))
    arr[3:18, 3:18] = 10
    expected_rows = np.zeros(20, dtype=np.bool)
    expected_cols = np.zeros(20, dtype=np.bool)
    expected_rows[3:18] = True
    expected_cols[3:18] = True

    rows = np.ones(20, dtype=np.bool)
    cols = np.ones(20, dtype=np.bool)
    rows, cols = model._multiple_node_deletion(rows, cols, arr)
    assert_array_equal(rows, expected_rows)
    assert_array_equal(cols, expected_cols)


def test_node_deletion_when_unnecessary():
    model = ChengChurch()
    arr = np.zeros((20, 20))
    rows = np.ones(20, dtype=np.bool)
    cols = np.ones(20, dtype=np.bool)

    expected_rows = np.ones(20, dtype=np.bool)
    expected_cols = np.ones(20, dtype=np.bool)

    rows, cols = model._single_node_deletion(rows, cols, arr)
    assert_array_equal(rows, expected_rows)
    assert_array_equal(cols, expected_cols)

    rows, cols = model._multiple_node_deletion(rows, cols, arr)
    assert_array_equal(rows, expected_rows)
    assert_array_equal(cols, expected_cols)


def test_node_addition():
    generator = np.random.RandomState(0)
    model = ChengChurch()
    arr = np.zeros((20, 20)) + 15
    arr[5:10, 5:10] = generator.uniform(10, 20, (5, 5))
    rows = np.zeros(20, dtype=np.bool)
    cols = np.zeros(20, dtype=np.bool)
    rows[5:10] = True
    cols[5:10] = True
    expected_rows = np.ones(20, dtype=np.bool)
    expected_cols = np.ones(20, dtype=np.bool)

    rows, cols = model._node_addition(rows, cols, arr)
    assert_array_equal(rows, expected_rows)
    assert_array_equal(cols, expected_cols)


def test_row_and_col_msr():
    X = np.arange(9).reshape(3, 3)
    rows = np.zeros(3, dtype=np.bool)
    cols = np.zeros(3, dtype=np.bool)
    rows[0:2] = True
    cols[0:2] = True
    model = ChengChurch()

    vec = model._row_msr(rows, cols, X)
    expected = np.zeros(3)
    assert_array_equal(vec, expected)

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

    assert_raises(EmptyBiclusterException, model._msr,
                  np.array([]), np.array([]), data)

    assert_raises(EmptyBiclusterException, model._row_msr,
                  np.array([]), np.array([]), data)

    assert_raises(EmptyBiclusterException, model._col_msr,
                  np.array([]), np.array([]), data)


def test_incremental_msr_errors():
    data = np.arange(25).reshape((5, 5))
    rows = np.zeros(5, dtype=np.bool)
    cols = np.zeros(5, dtype=np.bool)
    rows[0] = True
    cols[0] = True

    inc = IncrementalMSR(rows, cols, data)
    assert_raises(EmptyBiclusterException, inc.remove_row, 0)

    inc = IncrementalMSR(rows, cols, data)
    assert_raises(EmptyBiclusterException, inc.remove_col, 0)

    inc = IncrementalMSR(rows, cols, data)
    assert_raises(ValueError, inc.remove_row, 1)

    inc = IncrementalMSR(rows, cols, data)
    assert_raises(ValueError, inc.remove_col, 1)
