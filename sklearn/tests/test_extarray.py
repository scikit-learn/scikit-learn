import pytest

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._extarray import _get_feature_names


def _construct_array(array_type, column_names):
    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)

    if array_type == "dataframe":
        pd = pytest.importorskip("pandas")
        return pd.DataFrame(X, columns=column_names)
    elif array_type == "dataarray":
        xr = pytest.importorskip("xarray")
        return xr.DataArray(X, dims=('index', 'columns'),
                            coords={'columns': column_names})


@pytest.mark.parametrize("array_type", ["dataframe", "dataarray"])
def test_pandas_get_feature_names(array_type):
    column_names = [f'col_{i}' for i in range(3)]
    X = _construct_array(array_type, column_names)
    names = _get_feature_names(X)

    assert_array_equal(names, column_names)
