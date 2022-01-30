import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    concat,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "ops, names",
    [
        ([np.sqrt], ["sqrt"]),
        ([np.abs, np.sqrt], ["absolute", "sqrt"]),
        (np.array([np.sqrt]), ["sqrt"]),
        (np.array([np.abs, np.sqrt]), ["absolute", "sqrt"]),
    ],
)
def test_transform_listlike(string_series, ops, names):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([op(string_series) for op in ops], axis=1)
        expected.columns = names
        result = string_series.transform(ops)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("box", [dict, Series])
def test_transform_dictlike(string_series, box):
    # GH 35964
    with np.errstate(all="ignore"):
        expected = concat([np.sqrt(string_series), np.abs(string_series)], axis=1)
    expected.columns = ["foo", "bar"]
    result = string_series.transform(box({"foo": np.sqrt, "bar": np.abs}))
    tm.assert_frame_equal(result, expected)


def test_transform_dictlike_mixed():
    # GH 40018 - mix of lists and non-lists in values of a dictionary
    df = Series([1, 4])
    result = df.transform({"b": ["sqrt", "abs"], "c": "sqrt"})
    expected = DataFrame(
        [[1.0, 1, 1.0], [2.0, 4, 2.0]],
        columns=MultiIndex([("b", "c"), ("sqrt", "abs")], [(0, 0, 1), (0, 1, 0)]),
    )
    tm.assert_frame_equal(result, expected)
