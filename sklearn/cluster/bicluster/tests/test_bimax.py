"""Testing for BiMax."""

import numpy as np
from sklearn.cluster.bicluster import BiMax


def test_bimax():
    data = np.zeros((20, 20), dtype=np.int8)
    data[0:10, 0:10] = 1
    model = BiMax()
    model.fit(data)
    assert len(model.rows_) == 1
    assert len(model.columns_) == 1
    rows, cols = model.get_indices(0)
    assert set(rows) == set(range(10))
    assert set(cols) == set(range(10))
