import numpy as np
import pytest

from sklearn.tree._utils import _py_swap_array_slices


@pytest.mark.parametrize("dtype", [np.float32, np.intp])
def test_py_swap_array_slices_random(dtype):
    def swap_slices_np(arr, start, end, n):
        """
        Swaps the order of the slices array[start:start + n]
        and array[start + n:end] while preserving the order
        in the slices. Works for any itemsize.
        """
        arr[start:end] = np.concatenate([arr[start + n : end], arr[start : start + n]])

    rng = np.random.default_rng(0)

    for _ in range(10):
        arr = rng.permutation(100).astype(dtype)
        start = rng.integers(40)
        end = rng.integers(60, 100)
        split = rng.integers(end - start)

        expected = arr.copy()
        swap_slices_np(expected, start, end, split)

        _py_swap_array_slices(arr, start, end, split)
        np.testing.assert_array_equal(arr, expected)
