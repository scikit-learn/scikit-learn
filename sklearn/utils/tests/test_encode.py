import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._encode import _encode, _encode_labels, _get_counts, _unique


@pytest.mark.parametrize(
    "values, expected",
    [
        (np.array([2, 1, 3, 1, 3], dtype="int64"), np.array([1, 2, 3], dtype="int64")),
        (
            np.array([2, 1, np.nan, 1, np.nan], dtype="float32"),
            np.array([1, 2, np.nan], dtype="float32"),
        ),
        (
            np.array(["b", "a", "c", "a", "c"], dtype=object),
            np.array(["a", "b", "c"], dtype=object),
        ),
        (
            np.array(["b", "a", None, "a", None], dtype=object),
            np.array(["a", "b", None], dtype=object),
        ),
        (np.array(["b", "a", "c", "a", "c"]), np.array(["a", "b", "c"])),
    ],
    ids=["int64", "float32-nan", "object", "object-None", "str"],
)
@pytest.mark.parametrize("encode", [_encode, _encode_labels])
def test_encode_util(values, expected, encode):
    uniques = _unique(values)
    assert_array_equal(uniques, expected)

    result, encoded = _unique(values, return_inverse=True)
    assert_array_equal(result, expected)
    assert_array_equal(encoded, np.array([1, 0, 2, 0, 2]))

    encoded = encode(values, uniques=uniques)
    assert_array_equal(encoded, np.array([1, 0, 2, 0, 2]))

    result, counts = _unique(values, return_counts=True)
    assert_array_equal(result, expected)
    assert_array_equal(counts, np.array([2, 1, 2]))

    result, encoded, counts = _unique(values, return_inverse=True, return_counts=True)
    assert_array_equal(result, expected)
    assert_array_equal(encoded, np.array([1, 0, 2, 0, 2]))
    assert_array_equal(counts, np.array([2, 1, 2]))


def test_encode_unknown_values():
    uniques = np.array([1, 2, 3])
    values = np.array([1, 2, 3, 4])

    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    assert_array_equal(encoded, [0, 1, 2, -1])
    assert_array_equal(diff, [4])

    uniques = np.array(["a", "b", "c"], dtype=object)
    values = np.array(["a", "b", "c", "d"], dtype=object)

    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    assert_array_equal(encoded, [0, 1, 2, -1])
    assert_array_equal(diff, ["d"])


@pytest.mark.parametrize("missing_value", [None, np.nan, float("nan")])
def test_encode_unknown_missing_values(missing_value):
    values = np.array(["d", "c", "a", "b", missing_value], dtype=object)
    uniques = np.array(["c", "a", "b", missing_value], dtype=object)

    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    assert_array_equal(encoded, [-1, 0, 1, 2, 3])
    assert_array_equal(diff, ["d"])

    values = np.array(["d", "c", "a", "b", missing_value], dtype=object)
    uniques = np.array(["c", "a", "b"], dtype=object)

    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    assert_array_equal(encoded, [-1, 0, 1, 2, -1])
    assert_array_equal(diff[:-1], ["d"])
    if missing_value is None:
        assert diff[-1] is None
    else:
        assert np.isnan(diff[-1])

    values = np.array(["a", missing_value], dtype=object)
    uniques = np.array(["a", "b", "z"], dtype=object)

    encoded, diff = _encode(values, uniques=uniques, return_diff=True)
    assert_array_equal(encoded, [0, -1])
    if missing_value is None:
        assert diff[0] is None
    else:
        assert np.isnan(diff[0])


def test_encode_labels_unknown_values():
    with pytest.raises(ValueError, match="y contains previously unseen labels"):
        _encode_labels(np.array([1, 2, 4]), uniques=np.array([1, 2, 3]))


@pytest.mark.parametrize("missing_value", [np.nan, None, float("nan")])
def test_unique_util_missing_values_objects(missing_value):
    # check for _unique and _encode with missing values with object dtypes
    values = np.array(["a", "c", "c", missing_value, "b"], dtype=object)
    expected_uniques = np.array(["a", "b", "c", missing_value], dtype=object)

    uniques = _unique(values)

    if missing_value is None:
        assert_array_equal(uniques, expected_uniques)
    else:  # missing_value == np.nan
        assert_array_equal(uniques[:-1], expected_uniques[:-1])
        assert np.isnan(uniques[-1])

    encoded = _encode(values, uniques=uniques)
    assert_array_equal(encoded, np.array([0, 2, 2, 3, 1]))


def test_unique_util_missing_values_numeric():
    # Check missing values in numerical values
    values = np.array([3, 1, np.nan, 5, 3, np.nan], dtype=float)
    expected_uniques = np.array([1, 3, 5, np.nan], dtype=float)
    expected_inverse = np.array([1, 0, 3, 2, 1, 3])

    uniques = _unique(values)
    assert_array_equal(uniques, expected_uniques)

    uniques, inverse = _unique(values, return_inverse=True)
    assert_array_equal(uniques, expected_uniques)
    assert_array_equal(inverse, expected_inverse)

    encoded = _encode(values, uniques=uniques)
    assert_array_equal(encoded, expected_inverse)


def test_unique_util_with_all_missing_values():
    # test for all types of missing values for object dtype
    values = np.array([np.nan, "a", "c", "c", None, float("nan"), None], dtype=object)

    uniques = _unique(values)
    assert_array_equal(uniques[:-1], ["a", "c", None])
    # last value is nan
    assert np.isnan(uniques[-1])

    expected_inverse = [3, 0, 1, 1, 2, 3, 2]
    _, inverse = _unique(values, return_inverse=True)
    assert_array_equal(inverse, expected_inverse)


def test_encode_with_both_missing_values():
    # test for both types of missing values for object dtype
    values = np.array([np.nan, "a", "c", "c", None, np.nan, None], dtype=object)

    encoded, diff = _encode(
        values, uniques=np.array(["a", "c"], dtype=object), return_diff=True
    )
    assert_array_equal(encoded, [-1, 0, 1, 1, -1, -1, -1])
    assert diff[0] is None
    assert np.isnan(diff[1])


NAN1 = float("nan")
NAN2 = float("nan")


@pytest.mark.parametrize(
    "values, uniques, expected_counts",
    [
        (np.array([1] * 10 + [2] * 4 + [3] * 15), np.array([1, 2, 3]), [10, 4, 15]),
        (
            np.array([1] * 10 + [2] * 4 + [3] * 15),
            np.array([1, 2, 3, 5]),
            [10, 4, 15, 0],
        ),
        (
            np.array([np.nan] * 10 + [2] * 4 + [3] * 15),
            np.array([2, 3, np.nan]),
            [4, 15, 10],
        ),
        (
            np.array(["b"] * 4 + ["a"] * 16 + ["c"] * 20, dtype=object),
            ["a", "b", "c"],
            [16, 4, 20],
        ),
        (
            np.array(["b"] * 4 + ["a"] * 16 + ["c"] * 20, dtype=object),
            ["c", "b", "a"],
            [20, 4, 16],
        ),
        (
            np.array([np.nan] * 4 + ["a"] * 16 + ["c"] * 20, dtype=object),
            ["c", np.nan, "a"],
            [20, 4, 16],
        ),
        (
            np.array(["b"] * 4 + ["a"] * 16 + ["c"] * 20, dtype=object),
            ["a", "b", "c", "e"],
            [16, 4, 20, 0],
        ),
    ],
)
def test_get_counts(values, uniques, expected_counts):
    counts = _get_counts(values, uniques)
    assert_array_equal(counts, expected_counts)


def test_get_counts_multiple_nans():
    """
    When both np.nan and float("nan") are present, they get merged into np.nan.
    """
    values = np.array(
        ["a", np.nan, NAN1, np.nan, NAN2, NAN1, np.nan, "a"],
        dtype=object,
    )
    uniques = np.array(["a", np.nan], dtype=object)
    expected_counts = [2, 6]
    assert_array_equal(
        _get_counts(values, uniques, [np.nan, NAN1, NAN2]), expected_counts
    )

    # Now try it via _unique, to make sure this works end-to-end:
    real_uniques, real_counts = _unique(values, return_counts=True)
    # Comparing two arrays with nan fails cause the nans are not equal to
    # themselves. So compare as Python lists:
    assert list(uniques) == list(real_uniques)
    assert_array_equal(real_counts, expected_counts)
