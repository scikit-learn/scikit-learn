"""Testing for BiMax."""

import numpy as np
from sklearn.cluster.bicluster import BiMax
from sklearn.cluster.bicluster.bimax import _precompute_neighbors


def test_neighbors():
    data = np.zeros((20, 20), dtype=np.bool)
    data[0:10, 0:10] = 1
    neighbors = _precompute_neighbors(data)
    for i in range(40):
        if (0 <= i < 10) or (20 <= i < 30):
            assert len(neighbors[i]) == 19
        else:
            assert len(neighbors[i]) == 0


def test_bimax():
    data = np.zeros((20, 20), dtype=np.bool)
    data[0:10, 0:10] = 1
    model = BiMax(random_state=0)
    model.fit(data)
    assert len(model.rows_) == 1
    assert len(model.columns_) == 1
    rows, cols = model.get_indices(0)
    assert set(rows) == set(range(10))
    assert set(cols) == set(range(10))
