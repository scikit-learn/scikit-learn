from __future__ import absolute_import

import numpy as np

from skimage._shared import testing
from skimage._shared.testing import assert_array_almost_equal
from skimage.filters import threshold_local, gaussian
from skimage.util.apply_parallel import apply_parallel, dask_available


@testing.skipif(not dask_available, reason="dask not installed")
def test_apply_parallel():
    # data
    a = np.arange(144).reshape(12, 12).astype(float)

    # apply the filter
    expected1 = threshold_local(a, 3)
    result1 = apply_parallel(threshold_local, a, chunks=(6, 6), depth=5,
                             extra_arguments=(3,),
                             extra_keywords={'mode': 'reflect'})

    assert_array_almost_equal(result1, expected1)

    def wrapped_gauss(arr):
        return gaussian(arr, 1, mode='reflect')

    expected2 = gaussian(a, 1, mode='reflect')
    result2 = apply_parallel(wrapped_gauss, a, chunks=(6, 6), depth=5)

    assert_array_almost_equal(result2, expected2)


@testing.skipif(not dask_available, reason="dask not installed")
def test_no_chunks():
    a = np.ones(1 * 4 * 8 * 9).reshape(1, 4, 8, 9)

    def add_42(arr):
        return arr + 42

    expected = add_42(a)
    result = apply_parallel(add_42, a)

    assert_array_almost_equal(result, expected)


@testing.skipif(not dask_available, reason="dask not installed")
def test_apply_parallel_wrap():
    def wrapped(arr):
        return gaussian(arr, 1, mode='wrap')
    a = np.arange(144).reshape(12, 12).astype(float)
    expected = gaussian(a, 1, mode='wrap')
    result = apply_parallel(wrapped, a, chunks=(6, 6), depth=5, mode='wrap')

    assert_array_almost_equal(result, expected)


@testing.skipif(not dask_available, reason="dask not installed")
def test_apply_parallel_nearest():
    def wrapped(arr):
        return gaussian(arr, 1, mode='nearest')
    a = np.arange(144).reshape(12, 12).astype(float)
    expected = gaussian(a, 1, mode='nearest')
    result = apply_parallel(wrapped, a, chunks=(6, 6), depth={0: 5, 1: 5},
                            mode='nearest')

    assert_array_almost_equal(result, expected)
