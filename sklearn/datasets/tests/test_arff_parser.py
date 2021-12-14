import numpy as np
import pytest

from sklearn.datasets._arff_parser import _cast_frame


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

    assert X_casted_dtypes["col_int_as_integer"].kind == "i"
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

    assert X_casted_dtypes["col_int_as_integer"].kind == "i"
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
