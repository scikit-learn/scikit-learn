import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from sklearn.utils._sorting import _simultaneous_sort

from sklearn.utils import check_random_state


@pytest.mark.parametrize("kind", ["introsort", "heapsort", "quick_sort"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_simultaneous_sort(kind, dtype, n_pts=201):
    # Sort sanity check
    rng = check_random_state(0)
    values = rng.random_sample(n_pts).astype(dtype, copy=False)
    indices = np.arange(n_pts).astype(np.int64, copy=False)

    values_2 = values.copy()
    indices_2 = indices.copy()

    _simultaneous_sort(values, indices, kind=kind)

    sorted_indices = np.argsort(values_2)
    values_2 = values_2[sorted_indices]
    indices_2 = indices_2[sorted_indices]

    assert_array_almost_equal(values, values_2)
    assert_array_almost_equal(indices, indices_2)
