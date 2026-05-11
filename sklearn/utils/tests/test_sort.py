# SPDX-License-Identifier: BSD-3-Clause

from time import perf_counter

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._sorting import _py_simultaneous_sort


@pytest.mark.parametrize("kind", ["2-way", "3-way"])
def test_simultaneous_sort_correctness(kind):
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
        _py_simultaneous_sort(x_sorted, ind, n, use_three_way_partition=kind == "3-way")
        assert (x_sorted[:-1] <= x_sorted[1:]).all()
        assert_array_equal(x[ind], x_sorted)
        assert_array_equal(np.sort(ind), np.arange(n, dtype=np.intp))


@pytest.mark.parametrize("kind", ["2-way", "3-way"])
def test_simultaneous_sort_not_quadratic(kind):
    n = int(1e5)

    values = np.random.default_rng(0).uniform(size=n)
    indices = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_simultaneous_sort(
        values, indices, values.shape[0], use_three_way_partition=kind == "3-way"
    )
    dt_ref = perf_counter() - t0

    # worst case pattern (i.e. triggers the quadratic path)
    # for naive 2-way partitioning quicksort:
    values = np.zeros(n)
    indices = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_simultaneous_sort(
        values, indices, values.shape[0], use_three_way_partition=kind == "3-way"
    )
    dt = perf_counter() - t0
    # sorting 100k constant elements should not take more than 10x the time
    # needed to sort 100k random elements:
    assert dt < dt_ref * 10

    # worst case pattern for the better (numpy-style) 2-way partitioning:
    values = np.roll(np.arange(n), -1).astype(np.float32)
    indices = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_simultaneous_sort(
        values, indices, values.shape[0], use_three_way_partition=kind == "3-way"
    )
    dt = perf_counter() - t0
    assert dt < dt_ref * 10

    # worst case pattern for the 3-way partitioning quicksort
    # with median-of-3 pivot:
    k = n // 2
    values = np.array(
        [i if i % 2 == 1 else k + i - 1 for i in range(1, k + 1)]
        + [i for i in range(1, 2 * k + 1) if i % 2 == 0]
    ).astype(np.float64)
    # (very unlikely in real-world non-adversarial data)
    indices = np.arange(n, dtype=np.intp)
    assert values.size == indices.size
    t0 = perf_counter()
    _py_simultaneous_sort(
        values, indices, values.shape[0], use_three_way_partition=kind == "3-way"
    )
    dt = perf_counter() - t0
    assert dt < dt_ref * 10
