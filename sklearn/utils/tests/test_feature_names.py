import pytest
import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._feature_names import _get_feature_names


@pytest.mark.parametrize(
    "names",
    [
        list(range(2)),
        range(2),
        [["a", "b"], ["c", "d"]],
    ],
    ids=["list-int", "range", "multi-index"],
)
def test_get_feature_names_pandas_warns(names):
    """Get feature names with pandas dataframes with warnings."""
    pd = pytest.importorskip("pandas")
    X = pd.DataFrame([[1, 2], [4, 5], [5, 6]], columns=names)

    msg = "Feature name support requires all feature names to be strings"
    with pytest.warns(FutureWarning, match=msg):
        names = _get_feature_names(X)
    assert names is None


def test_get_feature_names_pandas():
    """Get feature names with pandas dataframes."""
    pd = pytest.importorskip("pandas")
    columns = [f"col_{i}" for i in range(3)]
    X = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=columns)
    feature_names = _get_feature_names(X)

    assert_array_equal(feature_names, columns)


def test_get_feature_names_numpy():
    """Get feature names return None for numpy arrays."""
    X = np.array([[1, 2, 3], [4, 5, 6]])
    names = _get_feature_names(X)
    assert names is None
