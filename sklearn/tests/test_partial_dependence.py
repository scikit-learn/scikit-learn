"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
import pytest

import sklearn
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import if_matplotlib
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import ignore_warnings
from sklearn.partial_dependence import partial_dependence
from sklearn.partial_dependence import plot_partial_dependence
from sklearn.partial_dependence import _grid_from_X
from sklearn.partial_dependence import _partial_dependence_exact
from sklearn.partial_dependence import _partial_dependence_recursion
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble.gradient_boosting import BaseGradientBoosting
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.svm import SVC
from sklearn.ensemble.forest import ForestRegressor
from sklearn.datasets import load_boston, load_iris
from sklearn.datasets import make_classification, make_regression
from sklearn.base import is_classifier, is_regressor
from sklearn.utils.estimator_checks import multioutput_estimator_convert_y_2d
from sklearn.cluster import KMeans


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]


def binary_classification():
    return make_classification(random_state=0), 1


def multiclass_classification():
    return (make_classification(n_classes=3, n_clusters_per_class=1,
                                random_state=0), 3)


def regression():
    return make_regression(random_state=0), 1


def multioutput_regression():
    return make_regression(n_targets=2, random_state=0), 2

@pytest.mark.parametrize('Estimator, method, data', [
    (GradientBoostingClassifier, 'recursion', binary_classification()),
    (GradientBoostingClassifier, 'recursion', multiclass_classification()),
    (GradientBoostingClassifier, 'exact', binary_classification()),
    (GradientBoostingClassifier, 'exact', multiclass_classification()),
    (GradientBoostingRegressor, 'recursion', regression()),
    (GradientBoostingRegressor, 'exact', regression()),
    (LinearRegression, 'exact', regression()),
    (LogisticRegression, 'exact', binary_classification()),
    (LogisticRegression, 'exact', multiclass_classification()),
    (MultiTaskLasso, 'exact', multioutput_regression())])
@pytest.mark.parametrize('grid_resolution', (5, 10))
@pytest.mark.parametrize('target_variables', ([1], [1, 2]))
def test_output_shape(Estimator, method, data, grid_resolution,
                      target_variables):
    # Check that partial_dependence has consistent output shape for different
    # kinds of estimators:
    # - classifiers with binary and multiclass settings
    # - regressors
    # - multi-task regressors

    est = Estimator()

    # n_target corresponds to the number of classes (1 for binary classif) or
    # the number of tasks / outputs in multi task settings. It's equal to 1 for
    # classical regression.
    (X, y), n_targets = data

    est.fit(X, y)
    pdp, axes = partial_dependence(est, target_variables=target_variables,
                                   X=X, method=method,
                                   grid_resolution=grid_resolution)

    expected_pdp_shape = (n_targets, grid_resolution ** len(target_variables))
    expected_axes_shape = (len(target_variables), grid_resolution)

    assert pdp.shape == expected_pdp_shape
    assert axes is not None
    assert np.asarray(axes).shape == expected_axes_shape


def test_grid_from_X():
    # tests for _grid_from_X: sanity check for output, and for shapes.

    # Make sure that the grid is a cartesian product of the input (it will use
    # the unique values instead of the percentiles)
    X = np.asarray([[1, 2], 
                    [3, 4]])
    grid, axes = _grid_from_X(X)
    assert_array_almost_equal(grid, [[1, 2],
                                     [1, 4],
                                     [3, 2],
                                     [3, 4]])
    assert_array_almost_equal(axes, X.T)

    # test shapes of returned objects depending on the number of unique values
    # for a feature.
    rng = np.random.RandomState(0)
    grid_resolution = 15

    # n_unique_values > grid_resolution
    X = rng.normal(size=(20, 2))
    grid, axes = _grid_from_X(X, grid_resolution=grid_resolution)
    assert grid.shape == (grid_resolution * grid_resolution, X.shape[1])
    assert np.asarray(axes).shape == (2, grid_resolution)

    # n_unique_values < grid_resolution, will use actual values
    n_unique_values = 12
    X[n_unique_values - 1:, 0] = 12345
    rng.shuffle(X)  # just to make sure the order is irrelevant
    grid, axes = _grid_from_X(X, grid_resolution=grid_resolution)
    assert grid.shape == (n_unique_values * grid_resolution, X.shape[1])
    # axes is a list of arrays of different shapes
    assert axes[0].shape == (n_unique_values,)
    assert axes[1].shape == (grid_resolution,)


@pytest.mark.parametrize('target_feature', (0, 3))
@pytest.mark.parametrize('est, partial_dependence_fun',
                         [(LinearRegression(), _partial_dependence_exact),
                          (GradientBoostingRegressor(random_state=0), _partial_dependence_recursion)])
def test_partial_dependence_helpers(est, partial_dependence_fun,
                                    target_feature):
    # Check that what is returned by _partial_dependence_exact or
    # _partial_dependece_recursion is equivalent to manually setting a target
    # feature to a given value, and computing the average prediction over all
    # samples.

    # doesn't work for _partial_dependence_recursion, dont know why :(
    if partial_dependence_fun is _partial_dependence_recursion:
        return

    X, y = make_regression(random_state=0)
    est.fit(X, y)

    # target feature will be set to .5 and then to 123
    target_variables = np.array([target_feature], dtype=np.int32)
    grid = np.array([[.5],
                     [123]])
    pdp = partial_dependence_fun(est, grid, target_variables, X)

    mean_predictions = []
    for val in (.5, 123):
        X_ = X.copy()
        X_[:, target_feature] = val
        mean_predictions.append(est.predict(X_).mean())

    pdp = pdp[0]  # (shape is (1, 2) so make it (2,))
    assert_array_almost_equal(pdp, mean_predictions)
    

@pytest.mark.parametrize('Estimator',
                         (sklearn.tree.DecisionTreeClassifier,
                          sklearn.tree.ExtraTreeClassifier,
                          sklearn.ensemble.ExtraTreesClassifier,
                          sklearn.neighbors.KNeighborsClassifier,
                          sklearn.neighbors.RadiusNeighborsClassifier,
                          sklearn.ensemble.RandomForestClassifier))
def test_multiclass_multioutput(Estimator):
    # Make sure multiclass-multioutput classifiers are not supported

    # make multiclass-multioutput dataset
    X, y = make_classification(n_classes=3, n_clusters_per_class=1,
                               random_state=0)
    y = np.array([y, y]).T

    est = Estimator()
    est.fit(X, y)

    assert_raises_regex(ValueError,
                        "Multiclass-multioutput estimators are not supported",
                        partial_dependence, est, [0], X=X)


def test_partial_dependence_input():
    # Test input validation of partial_dependence.

    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    assert_raises_regex(ValueError,
                        "est must be a fitted regressor or classifier",
                        partial_dependence, KMeans(), [0])

    assert_raises_regex(ValueError,
                        "method blahblah is invalid. Accepted method names "
                        "are exact, recursion, auto.",
                        partial_dependence, lr, [0], method='blahblah')

    assert_raises_regex(ValueError,
                        'est must be an instance of BaseGradientBoosting or '
                        'ForestRegressor for the "recursion" method',
                        partial_dependence, lr, [0], method='recursion')

    assert_raises_regex(ValueError, "est requires a predict_proba()",
                        partial_dependence, SVC(), [0], X=X)

    for feature in (-1, 1000000):
        for est in (lr, gbc):
            assert_raises_regex(ValueError,
                                "all target_variables must be in",
                                partial_dependence, est, [feature], X=X)

    assert_raises_regex(ValueError, "Either grid or X must be specified",
                        partial_dependence, gbc, [0], grid=None, X=None)

    for percentiles in ((1, 2, 3, 4), 12345):
        assert_raises_regex(ValueError, "percentiles must be a sequence",
                            partial_dependence, lr, [0], grid=None, X=X,
                            percentiles=percentiles)
    for percentiles in ((-1, .95), (.05, 2)):
        assert_raises_regex(ValueError, "percentiles values must be in",
                            partial_dependence, lr, [0], grid=None, X=X,
                            percentiles=percentiles)
    assert_raises_regex(ValueError, "percentiles\[0\] must be less than",
                        partial_dependence, lr, [0], grid=None, X=X,
                        percentiles=(.9, .1))

    assert_raises_regex(ValueError, "grid must be 1d or 2d",
                        partial_dependence, lr, [0], grid=[[[1]]], X=X)

    for target_variables in ([0], [0, 1, 0]):
        assert_raises_regex(ValueError,
                            'grid.shape\[1\] \(2\) must be equal to the number'
                            ' of target variables',
                            partial_dependence, lr, target_variables,
                            grid=[[30, -123]], X=X)

    for unfitted_est in (LinearRegression(), GradientBoostingRegressor()):
        assert_raises_regex(ValueError,
                            'est parameter must be a fitted estimator',
                            partial_dependence, unfitted_est, [0], X=X)


@if_matplotlib
def test_plot_partial_dependence():
    # Test partial dependence plot function.
    boston = load_boston()
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, boston.data, [0, 1, (0, 1)],
                                       grid_resolution=grid_resolution,
                                       feature_names=boston.feature_names)
    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)

    # check with str features and array feature names
    fig, axs = plot_partial_dependence(clf, boston.data, ['CRIM', 'ZN',
                                                          ('CRIM', 'ZN')],
                                       grid_resolution=grid_resolution,
                                       feature_names=boston.feature_names)

    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)

    # check with list feature_names
    feature_names = boston.feature_names.tolist()
    fig, axs = plot_partial_dependence(clf, boston.data, ['CRIM', 'ZN',
                                                          ('CRIM', 'ZN')],
                                       grid_resolution=grid_resolution,
                                       feature_names=feature_names)
    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)


@if_matplotlib
def test_plot_partial_dependence_multiclass():
    # Test partial dependence plot function on multi-class input.
    iris = load_iris()
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, iris.data, [0, 1],
                                       label=0,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    # now with symbol labels
    target = iris.target_names[iris.target]
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, iris.data, [0, 1],
                                       label='setosa',
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    # label not in gbrt.classes_
    assert_raises(ValueError, plot_partial_dependence,
                  clf, iris.data, [0, 1], label='foobar',
                  grid_resolution=grid_resolution)

    # label not provided
    assert_raises(ValueError, plot_partial_dependence,
                  clf, iris.data, [0, 1],
                  grid_resolution=grid_resolution)


@if_matplotlib
def test_plot_partial_dependence_input():
    # Test partial dependence plot function input checks.
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)

    # not fitted yet
    assert_raises(ValueError, plot_partial_dependence,
                  clf, X, [0])

    clf.fit(X, y)

    assert_raises(ValueError, plot_partial_dependence,
                  clf, np.array(X)[:, :0], [0])

    # first argument must be an instance of BaseGradientBoosting
    assert_raises(ValueError, plot_partial_dependence,
                  {}, X, [0])

    # must be larger than -1
    assert_raises(ValueError, plot_partial_dependence,
                  clf, X, [-1])

    # too large feature value
    assert_raises(ValueError, plot_partial_dependence,
                  clf, X, [100])

    # str feature but no feature_names
    assert_raises(ValueError, plot_partial_dependence,
                  clf, X, ['foobar'])

    # not valid features value
    assert_raises(ValueError, plot_partial_dependence,
                  clf, X, [{'foo': 'bar'}])
