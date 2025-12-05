import numpy as np

from sklearn.utils._numpy_aligned_allocator import (
    set_aligned_allocation,
    set_numpy_default_allocation,
)


def test_change_numpy_allocator():
    # A no-op
    set_numpy_default_allocation()
    np.arange(3)
    # Change to aligned allocator
    set_aligned_allocation()
    np.arange(3)
    # Change back to default
    set_numpy_default_allocation()
    np.arange(3)
