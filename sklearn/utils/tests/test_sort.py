# SPDX-License-Identifier: BSD-3-Clause

from time import perf_counter

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._sorting import _py_sort


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
        _py_sort(x_sorted, ind, n, use_three_way_partition=kind == "3-way")
        assert (x_sorted[:-1] <= x_sorted[1:]).all()
        assert_array_equal(x[ind], x_sorted)
        assert_array_equal(np.sort(ind), np.arange(n, dtype=np.intp))


@pytest.mark.parametrize("kind", ["2-way", "3-way"])
def test_simultaneous_sort_not_quadratic(kind):
    n = int(3e5)
    # "killer pattern" (i.e. triggers the quadratic path)
    # for naive 2-way partitioning quicksort:
    dist = np.zeros(n)
    ind = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_sort(dist, ind, dist.shape[0], use_three_way_partition=kind == "3-way")
    dt = perf_counter() - t0
    # sorting 300k elements should take less than 1s unless the sort goes quadratic
    # (it should take something like <10ms)
    assert dt < 1

    # "killer pattern" for the better (numpy-style) 2-way partitioning:
    dist = np.roll(np.arange(n), -1).astype(np.float32)
    ind = np.arange(n, dtype=np.intp)
    t0 = perf_counter()
    _py_sort(dist, ind, dist.shape[0], use_three_way_partition=kind == "3-way")
    dt = perf_counter() - t0
    assert dt < 1

    # The "killer pattern" of the 3-way partitioning quicksort
    # with median-of-3 pivot is:
    # [i if i%2 == 1 else k+i-1 for i in range(1, k+1)]
    # + [i for i in range(1, 2*k + 1) if i%2 == 0]
    # which has 0% chance to exist in real data


def test_sort_log2_build():
    """Non-regression test for gh-30554.

    Using log2 and log in sort correctly sorts feature_values, but the tie breaking is
    different which can results in placing samples in a different order.
    """
    rng = np.random.default_rng(75)
    some = rng.normal(loc=0.0, scale=10.0, size=10).astype(np.float32)
    feature_values = np.concatenate([some] * 5)
    samples = np.arange(50, dtype=np.intp)
    _py_sort(feature_values, samples, 50, use_three_way_partition=True)
    # fmt: off
    # no black reformatting for this specific array
    expected_samples = [
        0, 40, 30, 20, 10, 29, 39, 19, 49,  9, 45, 15, 35,  5, 25, 11, 31,
        41,  1, 21, 22, 12,  2, 42, 32, 23, 13, 43,  3, 33,  6, 36, 46, 16,
        26,  4, 14, 24, 34, 44, 27, 47,  7, 37, 17,  8, 38, 48, 28, 18
    ]
    # fmt: on
    assert_array_equal(samples, expected_samples)
