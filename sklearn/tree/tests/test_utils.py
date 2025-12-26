import numpy as np
import pytest

from sklearn.tree._utils import _py_swap_array_slices


@pytest.mark.parametrize("dtype", [np.float32, np.intp])
def test_py_swap_array_slices_random(dtype, global_random_seed):
    def swap_slices_np(arr, start, end, n):
        """
        Swaps the order of the slices array[start:start + n] and
        array[start + n:end] while preserving the order in the slices.
        """
        arr = arr.copy()
        arr[start:end] = np.concatenate([arr[start + n : end], arr[start : start + n]])
        return arr

    rng = np.random.default_rng(global_random_seed)

    for _ in range(20):
        arr = rng.permutation(100).astype(dtype)
        start = rng.integers(40)
        end = rng.integers(60, 100)
        split = rng.integers(end - start)
        # test the swap of arr[start:start + split] with arr[start + split:end]
        expected = swap_slices_np(arr, start, end, split)

        _py_swap_array_slices(arr, start, end, split)
        np.testing.assert_array_equal(arr, expected)
