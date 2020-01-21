"""
Testing for the ICE module.
"""

import numpy as np
import pytest

from sklearn.inspection._ice import _grid_from_X
from sklearn.utils._testing import assert_array_equal


def test_grid_from_X():
    # tests for _grid_from_X: sanity check for output, and for shapes.

    # Make sure that the grid is a cartesian product of the input (it will use
    # the unique values instead of the percentiles)
    percentiles = (.05, .95)
    grid_resolution = 100
    X = np.asarray([[1, 2],
                    [3, 4]])
    grid, axes = _grid_from_X(X, percentiles, grid_resolution)
    assert_array_equal(grid, [[1, 2],
                              [1, 4],
                              [3, 2],
                              [3, 4]])
    assert_array_equal(axes, X.T)

    # test shapes of returned objects depending on the number of unique values
    # for a feature.
    rng = np.random.RandomState(0)
    grid_resolution = 15

    # n_unique_values > grid_resolution
    X = rng.normal(size=(20, 2))
    grid, axes = _grid_from_X(X, percentiles, grid_resolution=grid_resolution)
    assert grid.shape == (grid_resolution * grid_resolution, X.shape[1])
    assert np.asarray(axes).shape == (2, grid_resolution)

    # n_unique_values < grid_resolution, will use actual values
    n_unique_values = 12
    X[n_unique_values - 1:, 0] = 12345
    rng.shuffle(X)  # just to make sure the order is irrelevant
    grid, axes = _grid_from_X(X, percentiles, grid_resolution=grid_resolution)
    assert grid.shape == (n_unique_values * grid_resolution, X.shape[1])
    # axes is a list of arrays of different shapes
    assert axes[0].shape == (n_unique_values,)
    assert axes[1].shape == (grid_resolution,)


@pytest.mark.parametrize(
    "grid_resolution, percentiles, err_msg",
    [(2, (0, 0.0001), "percentiles are too close"),
     (100, (1, 2, 3, 4), "'percentiles' must be a sequence of 2 elements"),
     (100, 12345, "'percentiles' must be a sequence of 2 elements"),
     (100, (-1, .95), r"'percentiles' values must be in \[0, 1\]"),
     (100, (.05, 2), r"'percentiles' values must be in \[0, 1\]"),
     (100, (.9, .1), r"percentiles\[0\] must be strictly less than"),
     (1, (0.05, 0.95), "'grid_resolution' must be strictly greater than 1")]
)
def test_grid_from_X_error(grid_resolution, percentiles, err_msg):
    X = np.asarray([[1, 2], [3, 4]])
    with pytest.raises(ValueError, match=err_msg):
        _grid_from_X(
            X, grid_resolution=grid_resolution, percentiles=percentiles
        )
