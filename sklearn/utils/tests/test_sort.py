# SPDX-License-Identifier: BSD-3-Clause

from time import perf_counter

import numpy as np

from sklearn.utils._sorting import _py_simultaneous_sort


def test_simultaneous_sort_correctness():
    for x in [
        np.random.rand(10**3),
        np.random.rand(3),
        np.random.rand(10),
        np.random.geometric(0.2, 10**2).astype("float32"),
        np.random.choice(2, size=10**3).astype("float32"),
    ]:
        n = x.size
        ind = np.arange(n, dtype=np.intp)
        x_sorted = x.copy()
        _py_simultaneous_sort(x_sorted, ind, n)
        assert (x_sorted[:-1] <= x_sorted[1:]).all()
        assert (x[ind] == x_sorted).all()


def test_simultaneous_sort_not_quadratic():
    dist = np.zeros(100000, dtype=np.float64)
    ind = np.arange(dist.shape[0], dtype=np.intp)

    t0 = perf_counter()
    _py_simultaneous_sort(dist, ind, dist.shape[0])
    dt = perf_counter() - t0
    assert dt < 1e-2
