import pytest

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._feature_names import _get_feature_names


def _construct_array(array_type, column_names):
    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)

    if array_type == "dataframe":
        pd = pytest.importorskip("pandas")
        return pd.DataFrame(X, columns=column_names)
    else:
        xr = pytest.importorskip("xarray")
        return xr.DataArray(X, dims=('index', 'columns'),
                            coords={'columns': column_names})


@pytest.mark.parametrize("array_type", ["dataframe", "dataarray"])
def test_get_feature_names(array_type):
    column_names = [f'col_{i}' for i in range(3)]
    X = _construct_array(array_type, column_names)

    names = _get_feature_names(X)
    assert_array_equal(names, column_names)


@pytest.mark.parametrize("array_type", ["dataframe", "dataarray"])
@pytest.mark.parametrize("column_names", [
    np.array(["one", 2, "tree"], dtype=object),
    np.array([1, 2, 3], dtype=object)
])
def test_get_feature_names_non_str(array_type, column_names):
    X = _construct_array(array_type, column_names)
    msg = "X contains non-string feature names"
    with pytest.raises(ValueError, match=msg):
        _get_feature_names(X)
