"""Tests for dataframe detection functions."""

import numpy as np
import pytest

from sklearn.utils._dataframe import is_df_or_series, is_polars_df
from sklearn.utils._testing import _convert_container


@pytest.mark.parametrize("constructor_name", ["pyarrow", "pandas", "polars"])
def test_is_df_or_series(constructor_name):
    df = _convert_container([[1, 4, 2], [3, 3, 6]], constructor_name)

    assert is_df_or_series(df)
    assert not is_df_or_series(np.asarray([1, 2, 3]))


@pytest.mark.parametrize("constructor_name", ["pyarrow", "pandas", "polars"])
def test_is_polars_df_other_libraries(constructor_name):
    df = _convert_container([[1, 4, 2], [3, 3, 6]], constructor_name)
    if constructor_name in ("pyarrow", "pandas"):
        assert not is_polars_df(df)
    else:
        assert is_polars_df(df)


def test_is_polars_df_for_duck_typed_polars_dataframe():
    """Check is_polars_df for object that looks like a polars dataframe"""

    class NotAPolarsDataFrame:
        def __init__(self):
            self.columns = [1, 2, 3]
            self.schema = "my_schema"

    not_a_polars_df = NotAPolarsDataFrame()
    assert not is_polars_df(not_a_polars_df)


def test_is_polars_df():
    """Check that is_polars_df return False for non-dataframe objects."""

    class LooksLikePolars:
        def __init__(self):
            self.columns = ["a", "b"]
            self.schema = ["a", "b"]

    assert not is_polars_df(LooksLikePolars())
