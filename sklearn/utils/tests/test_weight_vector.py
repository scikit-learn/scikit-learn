import numpy as np
import pytest
from sklearn.utils._weight_vector import (
    WeightVector32,
    WeightVector64,
)


@pytest.mark.parametrize(
    "dtype, weight_vector_class",
    [
        (np.float32, WeightVector32),
        (np.float64, WeightVector64),
    ],
)
def test_type_invariance(dtype, weight_vector_class):
    weights = np.random.rand(100).astype(dtype)
    average_weights = np.random.rand(100).astype(dtype)
    weight_vector = weight_vector_class(weights, average_weights)
    assert np.asarray(weight_vector.w).dtype is np.dtype(dtype)
    assert np.asarray(weight_vector.aw).dtype is np.dtype(dtype)
