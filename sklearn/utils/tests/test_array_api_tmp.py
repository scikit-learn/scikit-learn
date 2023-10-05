import numpy
import pytest
from numpy.testing import assert_allclose

from sklearn._config import config_context
from sklearn.utils._array_api import (
    _atol_for_type,
    _average,
    _convert_to_numpy,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import (
    _array_api_for_tests,
)


@pytest.mark.parametrize(
    "array_namespace, device, dtype", yield_namespace_device_dtype_combinations()
)
@pytest.mark.parametrize(
    "weights, axis, expected",
    [
        (None, None, 3.5),
        (None, 0, [2.5, 3.5, 4.5]),
        (None, 1, [2, 5]),
        ([0.4, 0.1], 0, [1.6, 2.6, 3.6]),
        ([0.4, 0.2, 0.2], 1, [1.75, 4.75]),
        ([1, 2], 0, [3, 4, 5]),
        ([1, 1, 2], 1, [2.25, 5.25]),
        ([[1, 2, 3], [1, 2, 3]], 0, [2.5, 3.5, 4.5]),
        ([[1, 2, 1], [2, 2, 2]], 1, [2, 5]),
    ],
)
def test_average(array_namespace, device, dtype, weights, axis, expected):
    xp, device, dtype = _array_api_for_tests(array_namespace, device, dtype)
    sample_score = numpy.asarray([[1, 2, 3], [4, 5, 6]], dtype=dtype)
    sample_score = xp.asarray(sample_score, device=device)
    if weights is not None:
        weights = numpy.asarray(weights, dtype=dtype)
        weights = xp.asarray(weights, device=device)

    with config_context(array_api_dispatch=True):
        result = _average(sample_score, axis=axis, weights=weights)

    result = _convert_to_numpy(result, xp)
    assert_allclose(result, expected, atol=_atol_for_type(dtype))
