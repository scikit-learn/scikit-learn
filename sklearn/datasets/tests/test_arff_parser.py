import pytest

from sklearn.datasets._arff_parser import (
    _post_process_frame,
    load_arff_from_gzip_file,
)


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
