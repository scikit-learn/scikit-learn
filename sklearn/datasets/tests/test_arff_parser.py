import numpy as np
import pytest

from sklearn.datasets._arff_parser import (
    _cast_frame,
    _post_process_frame,
    load_arff_from_gzip_file,
)


def test_cast_frame():
    """Check that the casting is working as expected."""
    pd = pytest.importorskip("pandas")

    X_original = pd.DataFrame(
        {
            "col_int_as_integer": [1, 2, 3],
            "col_int_as_numeric": [1, 2, 3],
            "col_float_as_real": [1.0, 2.0, 3.0],
            "col_float_as_numeric": [1.0, 2.0, 3.0],
            "col_categorical": ["a", "b", "c"],
            "col_string": ["a", "b", "c"],
        }
    )

    columns_info = {
        "col_int_as_integer": {"data_type": "integer"},
        "col_int_as_numeric": {"data_type": "numeric"},
        "col_float_as_real": {"data_type": "real"},
        "col_float_as_numeric": {"data_type": "numeric"},
        "col_categorical": {"data_type": "nominal"},
        "col_string": {"data_type": "string"},
    }
    X_casted = _cast_frame(X_original.copy(), columns_info)
    X_casted_dtypes = X_casted.dtypes

    assert X_casted_dtypes["col_int_as_integer"] == "Int64"
    assert X_casted_dtypes["col_int_as_numeric"].kind == "i"
    assert X_casted_dtypes["col_float_as_real"].kind == "f"
    assert X_casted_dtypes["col_float_as_numeric"].kind == "f"
    assert X_casted_dtypes["col_categorical"].is_dtype("category")
    assert X_casted_dtypes["col_string"].kind == "O"

    X_original = pd.DataFrame(
        {
            "col_int_as_integer": [1, 2, 3, np.nan],
            "col_int_as_numeric": [1, 2, 3, np.nan],
            "col_float_as_real": [1.0, 2.0, 3.0, np.nan],
            "col_float_as_numeric": [1.0, 2.0, 3.0, np.nan],
            "col_categorical": ["a", "b", "c", np.nan],
            "col_string": ["a", "b", "c", np.nan],
        }
    )

    # since we use `np.nan` as missing values sentinel, we expect the integer
    # columns to actually be floats dtypes
    assert X_original.dtypes["col_int_as_integer"].kind == "f"
    assert X_original.dtypes["col_int_as_numeric"].kind == "f"

    X_casted = _cast_frame(X_original.copy(), columns_info)
    X_casted_dtypes = X_casted.dtypes

    assert X_casted_dtypes["col_int_as_integer"] == "Int64"
    # we cannot recover the integer dtype in this case
    assert X_casted_dtypes["col_int_as_numeric"].kind == "f"
    assert X_casted_dtypes["col_float_as_real"].kind == "f"
    assert X_casted_dtypes["col_float_as_numeric"].kind == "f"
    assert X_casted_dtypes["col_categorical"].is_dtype("category")
    # check that `np.nan` is not a category
    pd.testing.assert_index_equal(
        X_casted["col_categorical"].cat.categories, pd.Index(["a", "b", "c"])
    )
    assert X_casted_dtypes["col_string"].kind == "O"


@pytest.mark.parametrize(
    "feature_names, target_names",
    [
        (
            [
                "col_int_as_integer",
                "col_int_as_numeric",
                "col_float_as_real",
                "col_float_as_numeric",
            ],
            ["col_categorical", "col_string"],
        ),
        (
            [
                "col_int_as_integer",
                "col_int_as_numeric",
                "col_float_as_real",
                "col_float_as_numeric",
            ],
            ["col_categorical"],
        ),
        (
            [
                "col_int_as_integer",
                "col_int_as_numeric",
                "col_float_as_real",
                "col_float_as_numeric",
            ],
            [],
        ),
    ],
)
def test_post_process_frame(feature_names, target_names):
    """Check the behaviour of the post-processing function for splitting a dataframe."""
    pd = pytest.importorskip("pandas")

    X_original = pd.DataFrame(
        {
            "col_int_as_integer": [1, 2, 3],
            "col_int_as_numeric": [1, 2, 3],
            "col_float_as_real": [1.0, 2.0, 3.0],
            "col_float_as_numeric": [1.0, 2.0, 3.0],
            "col_categorical": ["a", "b", "c"],
            "col_string": ["a", "b", "c"],
        }
    )

    X, y = _post_process_frame(X_original, feature_names, target_names)
    assert isinstance(X, pd.DataFrame)
    if len(target_names) >= 2:
        assert isinstance(y, pd.DataFrame)
    elif len(target_names) == 1:
        assert isinstance(y, pd.Series)
    else:
        assert y is None


def test_load_arff_from_gzip_file_error_parser():
    """An error will be raised if the parser is not known."""
    # None of the input parameters are required to be accurate since the check
    # of the parser will be carried out first.

    err_msg = "Unknown parser: 'xxx'. Should be 'liac-arff' or 'pandas'"
    with pytest.raises(ValueError, match=err_msg):
        load_arff_from_gzip_file("xxx", "xxx", "xxx", "xxx", "xxx", "xxx")
