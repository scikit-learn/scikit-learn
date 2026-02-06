# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.utils._sorting import _py_simultaneous_sort


def test_simultaneous_sort_correctness():
    for dist in [
        np.random.rand(10**3),
        np.random.rand(3),
        np.random.rand(10),
        np.random.geometric(0.2, 10**2).astype("float32"),
    ]:
        ind = np.arange(dist.shape[0], dtype=np.intp)

        _py_simultaneous_sort(dist, ind, dist.shape[0])
        assert (dist[:-1] <= dist[1:]).all()


def test_simultaneous_sort_recursion_depth():
    dist = np.zeros(1000000, dtype=np.float64)
    ind = np.arange(dist.shape[0], dtype=np.intp)

    _py_simultaneous_sort(dist, ind, dist.shape[0])
