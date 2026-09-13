import numpy as np
import pytest

from sklearn.tree._partitioner import _py_swap_array_slices


@pytest.mark.parametrize("dtype", [np.float32, np.intp])
def test_py_swap_array_slices_random(dtype, global_random_seed):
    def swap_slices_np(arr, start, end, n):
        """
        Swaps the two slices array[start:start + n] and
        array[start + n:end] while preserving the order in the slices.
        """
        arr = arr.copy()
        arr[start:end] = np.concatenate([arr[start + n : end], arr[start : start + n]])
        return arr

    rng = np.random.default_rng(global_random_seed)

    for _ in range(20):
        size = rng.integers(1, 101)
        arr = rng.permutation(size).astype(dtype)
        n = rng.integers(0, size)
        start = rng.integers(0, size - n)
        end = rng.integers(start + n, size)
        # test the swap of arr[start:start + n] with arr[start + n:end]
        expected = swap_slices_np(arr, start, end, n)

        _py_swap_array_slices(arr, start, end, n)
        np.testing.assert_array_equal(arr, expected)

    # test some edge cases:
    size = 30
    n = 10
    start = rng.integers(0, size - n)
    arr = np.arange(size, dtype=dtype)
    expected = arr.copy()
    # n == end - start should be no-op:
    _py_swap_array_slices(arr, start, start + n, n)
    np.testing.assert_array_equal(arr, expected)
    # n == 0 should be no-op:
    _py_swap_array_slices(arr, start, size, 0)
    np.testing.assert_array_equal(arr, expected)
