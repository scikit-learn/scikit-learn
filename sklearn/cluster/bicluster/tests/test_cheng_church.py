"""Testing for Spectral Biclustering methods"""

import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises

from sklearn.cluster.bicluster import ChengChurch
from sklearn.metrics import consensus_score
from sklearn.datasets import make_msr


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
    generator = np.random.RandomState(0)
    data = generator.uniform(0, 1000, (50, 50))
    model = ChengChurch(n_clusters=1, max_msr=0)
    model.fit(data)
    assert_equal(len(model.rows_), 0)
    assert_equal(len(model.columns_), 0)


def test_single_node_deletion():
    generator = np.random.RandomState(0)
    model = ChengChurch()
    arr = generator.uniform(20, 100, (20, 20))
    arr[5:15, 5:15] = 10
    rows = np.arange(20)
    cols = np.arange(20)
    rows, cols = model._single_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(5, 15)[None].T)
    assert_array_equal(cols, np.arange(5, 15))


def test_multiple_node_deletion():
    generator = np.random.RandomState(0)
    model = ChengChurch(max_msr=250,
                        deletion_threshold=1.25,
                        row_deletion_cutoff=1,
                        column_deletion_cutoff=1)
    arr = generator.uniform(20, 100, (20, 20))
    arr[3:18, 3:18] = 10
    rows = np.arange(20)
    cols = np.arange(20)
    rows, cols = model._multiple_node_deletion(rows, cols, arr)
    assert_array_equal(rows, np.arange(3, 18)[None].T)
    assert_array_equal(cols, np.arange(3, 18))


def test_node_addition():
    generator = np.random.RandomState(0)
    model = ChengChurch()
    arr = np.zeros((20, 20)) + 15
    arr[5:10, 5:10] = generator.uniform(10, 20, (5, 5))
    rows = np.arange(5, 10)
    cols = np.arange(5, 10)
    rows, cols = model._node_addition(rows, cols, arr)
    assert_array_equal(rows, np.arange(20)[None].T)
    assert_array_equal(cols, np.arange(20))


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
