# Author: Mathis Batoul <batoulmathis@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from sklearn.utils._weight_vector import (
    WeightVector32,
    WeightVector64,
)


def test_type_invariance():
    weights32 = np.random.rand(100).astype(np.float32)
    average_weights32 = np.random.rand(100).astype(np.float32)
    weight_vector32 = WeightVector32(weights32, average_weights32)
    assert weight_vector32.w.dtype is np.dtype(np.float32)
    assert weight_vector32.aw.dtype is np.dtype(np.float32)

    weights64 = np.random.rand(100).astype(np.float64)
    average_weights64 = np.random.rand(100).astype(np.float64)
    weight_vector64 = WeightVector64(weights64, average_weights64)
    assert weight_vector64.w.dtype is np.dtype(np.float64)
    assert weight_vector64.aw.dtype is np.dtype(np.float64)
