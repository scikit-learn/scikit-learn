import numpy as np
import pytest

from sklearn.tree._partitioner import _py_swap_array_slices


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
        size = rng.integers(1, 101)
        arr = rng.permutation(size).astype(dtype)
        split = rng.integers(0, size)
        start = rng.integers(0, size - split)
        end = rng.integers(start + split, size)
        # test the swap of arr[start:start + split] with arr[start + split:end]
        expected = swap_slices_np(arr, start, end, split)

        _py_swap_array_slices(arr, start, end, split)
        np.testing.assert_array_equal(arr, expected)
