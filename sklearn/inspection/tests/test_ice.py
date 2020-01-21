"""
Testing for the ICE module.
"""

import numpy as np
import pytest

from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.inspection import individual_conditional_expectation
from sklearn.inspection._ice import _grid_from_X
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.tree import DecisionTreeRegressor
from sklearn.utils._testing import assert_array_equal

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]


# (X, y), n_targets  <-- as expected in the output of partial_dep()
binary_classification_data = (make_classification(n_samples=50,
                                                  random_state=0), 1)
multiclass_classification_data = (make_classification(n_samples=50,
                                                      n_classes=3,
                                                      n_clusters_per_class=1,
                                                      random_state=0), 3)
regression_data = (make_regression(n_samples=50, random_state=0), 1)
multioutput_regression_data = (make_regression(n_samples=50, n_targets=2,
                                               random_state=0), 2)


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


@pytest.mark.parametrize('Estimator, data', [
    (GradientBoostingClassifier, binary_classification_data),
    (GradientBoostingClassifier, multiclass_classification_data),
    (GradientBoostingRegressor, regression_data),
    (DecisionTreeRegressor, regression_data),
    (LinearRegression, regression_data),
    (LinearRegression, multioutput_regression_data),
    (LogisticRegression, binary_classification_data),
    (LogisticRegression, multiclass_classification_data),
    (MultiTaskLasso, multioutput_regression_data),
    ])
@pytest.mark.parametrize('grid_resolution', (5, 10))
@pytest.mark.parametrize('features', ([1], [1, 2]))
def test_output_shape(Estimator, data, grid_resolution, features):
    # Check that individual_conditional_expectation has consistent output
    # shape for different kinds of estimators:
    # - classifiers with binary and multiclass settings
    # - regressors
    # - multi-task regressors

    est = Estimator()

    # n_target corresponds to the number of classes (1 for binary classif) or
    # the number of tasks / outputs in multi task settings. It's equal to 1 for
    # classical regression_data. n_instances corresponds to the number of
    # training samples in X
    (X, y), n_targets = data
    n_instances = X.shape[0]

    est.fit(X, y)
    pdp, axes = individual_conditional_expectation(
        est, X=X, features=features, grid_resolution=grid_resolution
    )

    expected_ice_shape = (n_targets, n_instances,
                          *[grid_resolution for _ in range(len(features))])
    expected_axes_shape = (len(features), grid_resolution)

    assert pdp.shape == expected_ice_shape
    assert axes is not None
    assert np.asarray(axes).shape == expected_axes_shape
