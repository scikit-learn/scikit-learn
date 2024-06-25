import numpy as np
import pytest

from sklearn.model_selection._search import _get_param_masked_arrays
from sklearn.preprocessing import (
    OneHotEncoder,
    OrdinalEncoder,
)

# Construct these outside the tests so that the same object is used
# for both input and `expected`
one_hot_encoder = OneHotEncoder()
ordinal_encoder = OrdinalEncoder()

# If we construct this directly via `MaskedArray`, the list of tuples
# gets auto-converted to a 2D array.
ma = np.ma.MaskedArray(np.empty(2), mask=True, dtype=object)
ma[0] = (1, 2)
ma[1] = (3, 4)


@pytest.mark.parametrize(
    ("candidate_params", "expected"),
    [
        pytest.param(
            [{"foo": 1}, {"foo": 2}],
            [
                ("param_foo", np.ma.MaskedArray(np.array([1, 2], dtype=np.int64))),
            ],
            id="simple numeric, single param",
        ),
        pytest.param(
            [{"foo": 1, "bar": 3}, {"foo": 2, "bar": 4}, {"foo": 3}],
            [
                ("param_foo", np.ma.MaskedArray(np.array([1, 2, 3]), dtype=np.int64)),
                (
                    "param_bar",
                    np.ma.MaskedArray(
                        np.array([3, 4, 0]), mask=[False, False, True], dtype=np.int64
                    ),
                ),
            ],
            id="simple numeric, one param is missing in one round",
        ),
        pytest.param(
            [{"foo": [[1], [2], [3]]}, {"foo": [[1], [2]]}],
            [
                (
                    "param_foo",
                    np.ma.MaskedArray([[[1], [2], [3]], [[1], [2]]], dtype=object),
                ),
            ],
            id="lists of different lengths",
        ),
        pytest.param(
            [{"foo": (1, 2)}, {"foo": (3, 4)}],
            [
                ("param_foo", ma),
            ],
            id="lists tuples",
        ),
        pytest.param(
            [{"foo": ordinal_encoder}, {"foo": one_hot_encoder}],
            [
                (
                    "param_foo",
                    np.ma.MaskedArray([ordinal_encoder, one_hot_encoder], dtype=object),
                ),
            ],
            id="estimators",
        ),
    ],
)
def test_get_params_masked_array(candidate_params, expected) -> None:
    result = list(_get_param_masked_arrays(candidate_params))
    for (key, value), (expected_key, expected_value) in zip(result, expected):
        assert key == expected_key
        np.testing.assert_array_equal(value, expected_value)
        np.testing.assert_array_equal(value.mask, expected_value.mask)
