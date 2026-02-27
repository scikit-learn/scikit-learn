# SPDX-License-Identifier: BSD-3-Clause

from time import perf_counter

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._sorting import _py_simultaneous_sort


def test_simultaneous_sort_correctness():
    rng = np.random.default_rng(0)
    for x in [
        rng.uniform(size=3),
        rng.uniform(size=10),
        rng.uniform(size=1000),
        # with duplicates:
        rng.geometric(0.2, size=100).astype("float32"),
        rng.integers(0, 2, size=1000).astype("float32"),
    ]:
        n = x.size
        ind = np.arange(n, dtype=np.intp)
        x_sorted = x.copy()
        _py_simultaneous_sort(x_sorted, ind, n)
        assert (x_sorted[:-1] <= x_sorted[1:]).all()
        assert_array_equal(x[ind], x_sorted)
        assert_array_equal(np.sort(ind), np.arange(n, dtype=np.intp))


def test_simultaneous_sort_not_quadratic():
    n = int(3e5)
    dist = np.zeros(n)
    ind = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_simultaneous_sort(dist, ind, dist.shape[0])
    dt = perf_counter() - t0
    # sorting 300k elements should take less than 1s unless the sort goes quadratic
    # (it should take something like ~10ms)
    assert dt < 1
