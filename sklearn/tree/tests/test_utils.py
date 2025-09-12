import numpy as np

from sklearn.tree._utils import _py_swap_array_slices


def test_py_swap_array_slices_random():
    def swap_slices_np(arr, start, end, n):
        """
        Swaps the order of the slices array[start:start + n]
        and array[start + n:end] while preserving the order
        in the slices. Works for any itemsize.
        """
        arr[start:end] = np.concatenate([arr[start + n : end], arr[start : start + n]])

    for _ in range(10):
        arr = np.random.rand(100)
        start = np.random.randint(40)
        end = np.random.randint(60, 100)
        split = np.random.randint(end - start)

        expected = arr.copy()
        swap_slices_np(expected, start, end, split)

        _py_swap_array_slices(arr, start, end, split)
        np.testing.assert_array_equal(arr, expected)
