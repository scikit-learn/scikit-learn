import numpy as np
import pytest

from sklearn.utils._missing import is_scalar_nan


@pytest.mark.parametrize(
    "value, result",
    [
        (float("nan"), True),
        (np.nan, True),
        (float(np.nan), True),
        (np.float32(np.nan), True),
        (np.float64(np.nan), True),
        (0, False),
        (0.0, False),
        (None, False),
        ("", False),
        ("nan", False),
        ([np.nan], False),
        (9867966753463435747313673, False),  # Python int that overflows with C type
    ],
)
def test_is_scalar_nan(value, result):
    assert is_scalar_nan(value) is result
    # make sure that we are returning a Python bool
    assert isinstance(is_scalar_nan(value), bool)
