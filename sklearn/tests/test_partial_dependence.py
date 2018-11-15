"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

import sklearn
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import if_matplotlib
from sklearn.partial_dependence import partial_dependence
from sklearn.partial_dependence import plot_partial_dependence
from sklearn.partial_dependence import _grid_from_X
from sklearn.partial_dependence import _partial_dependence_brute
from sklearn.partial_dependence import _partial_dependence_recursion
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import MultiTaskLasso
from sklearn.svm import SVC
from sklearn.datasets import load_boston, load_iris
from sklearn.datasets import make_classification, make_regression
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
    (GradientBoostingClassifier, 'brute', binary_classification()),
    (GradientBoostingClassifier, 'brute', multiclass_classification()),
    (GradientBoostingRegressor, 'recursion', regression()),
    (GradientBoostingRegressor, 'brute', regression()),
    (LinearRegression, 'brute', regression()),
    (LinearRegression, 'brute', multioutput_regression()),
    (LogisticRegression, 'brute', binary_classification()),
    (LogisticRegression, 'brute', multiclass_classification()),
    (MultiTaskLasso, 'brute', multioutput_regression()),
    ])
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

    assert_raises_regex(ValueError, 'percentiles are too close',
                        _grid_from_X, X, grid_resolution=2,
                        percentiles=(0, 0.0001))

    for percentiles in ((1, 2, 3, 4), 12345):
        assert_raises_regex(ValueError, "percentiles must be a sequence",
                            _grid_from_X, X, percentiles=percentiles)

    for percentiles in ((-1, .95), (.05, 2)):
        assert_raises_regex(ValueError, "percentiles values must be in",
                            _grid_from_X, X, percentiles=percentiles)

    assert_raises_regex(ValueError,
                        r"percentiles\[0\] must be strictly less than",
                        _grid_from_X, X, percentiles=(.9, .1))

    assert_raises_regex(ValueError,
                        'grid_resolution must be strictly greater than 1.',
                        _grid_from_X, X, grid_resolution=1)


@pytest.mark.parametrize('target_feature', (0, 3))
@pytest.mark.parametrize('est, method',
                         [(LinearRegression(), 'brute'),
                          (GradientBoostingRegressor(random_state=0),
                          'recursion')])
def test_partial_dependence_helpers(est, method, target_feature):
    # Check that what is returned by _partial_dependence_brute or
    # _partial_dependece_recursion is equivalent to manually setting a target
    # feature to a given value, and computing the average prediction over all
    # samples.
    # This also checks that the brute method and the recursion give the same
    # output.

    X, y = make_regression(random_state=0)
    # The 'init' estimator for GBDT (here the average prediction) isn't taken
    # into account with the recursion method, for technical reasons. We set
    # the mean to 0 to that this 'bug' doesn't have any effect.
    y = y - y.mean()
    est.fit(X, y)

    # target feature will be set to .5 and then to 123
    target_variables = np.array([target_feature], dtype=np.int32)
    grid = np.array([[.5],
                     [123]])

    if method == 'brute':
        pdp = _partial_dependence_brute(est, grid, target_variables, X)
    else:
        pdp = _partial_dependence_recursion(est, grid, target_variables)

    mean_predictions = []
    for val in (.5, 123):
        X_ = X.copy()
        X_[:, target_feature] = val
        mean_predictions.append(est.predict(X_).mean())

    pdp = pdp[0]  # (shape is (1, 2) so make it (2,))
    assert_array_almost_equal(pdp, mean_predictions, decimal=3)


@pytest.mark.parametrize('Estimator',
                         (sklearn.tree.DecisionTreeClassifier,
                          sklearn.tree.ExtraTreeClassifier,
                          sklearn.ensemble.ExtraTreesClassifier,
                          sklearn.neighbors.KNeighborsClassifier,
                          sklearn.neighbors.RadiusNeighborsClassifier,
                          sklearn.ensemble.RandomForestClassifier))
def test_multiclass_multioutput(Estimator):
    # Make sure error is raised for multiclass-multioutput classifiers

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
                        "are brute, recursion, auto.",
                        partial_dependence, lr, [0], method='blahblah')

    assert_raises_regex(ValueError,
                        'est must be an instance of BaseGradientBoosting '
                        'for the "recursion" method',
                        partial_dependence, lr, [0], method='recursion')

    assert_raises_regex(ValueError, "X is required for brute method",
                        partial_dependence, lr, [0], grid=[[[1]]])

    assert_raises_regex(ValueError, "est requires a predict_proba()",
                        partial_dependence, SVC(), [0], X=X)

    for feature in (-1, 1000000):
        for est in (lr, gbc):
            assert_raises_regex(ValueError,
                                "all target_variables must be in",
                                partial_dependence, est, [feature], X=X)

    assert_raises_regex(ValueError, "Either grid or X must be specified",
                        partial_dependence, gbc, [0], grid=None, X=None)

    assert_raises_regex(ValueError, "grid must be 1d or 2d",
                        partial_dependence, lr, [0], grid=[[[1]]], X=X)

    for target_variables in ([0], [0, 1, 0]):
        assert_raises_regex(ValueError,
                            r'grid.shape\[1\] \(2\) must be equal '
                            r'to the number of target variables',
                            partial_dependence, lr, target_variables,
                            grid=[[30, -123]], X=X)

    for unfitted_est in (LinearRegression(), GradientBoostingRegressor()):
        assert_raises_regex(ValueError,
                            'est parameter must be a fitted estimator',
                            partial_dependence, unfitted_est, [0], X=X)

    # check that array-like objects are accepted
    for est in (lr, gbc):
        partial_dependence(est, [0], grid=[1, 2], X=list(X))
        partial_dependence(est, [0], grid=[[1], [2]], X=list(X))

    partial_dependence(gbc, [0], grid=[1, 2])


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
                                       target=0,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    # now with symbol labels
    target = iris.target_names[iris.target]
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, iris.data, [0, 1],
                                       target='setosa',
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    # label not in gbrt.classes_
    assert_raises(ValueError, plot_partial_dependence,
                  clf, iris.data, [0, 1], target='foobar',
                  grid_resolution=grid_resolution)

    # label not provided
    assert_raises(ValueError, plot_partial_dependence,
                  clf, iris.data, [0, 1],
                  grid_resolution=grid_resolution)


@if_matplotlib
def test_plot_partial_dependence_input():
    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    # check target param for multiclass
    (X_m, y_m), _ = multiclass_classification()
    lr_m = LogisticRegression()
    lr_m.fit(X_m, y_m)
    assert_raises_regex(ValueError,
                        'target must be specified for multi-class',
                        plot_partial_dependence, lr_m, X_m, [0],
                        target=None)
    for target in (-1, 100):
        assert_raises_regex(ValueError,
                            'target not in est.classes_',
                            plot_partial_dependence, lr_m, X_m, [0],
                            target=target)

    # check target param for multioutput
    (X_m, y_m), _ = multioutput_regression()
    lr_m = LinearRegression()
    lr_m.fit(X_m, y_m)
    assert_raises_regex(ValueError,
                        'target must be specified for multi-output',
                        plot_partial_dependence, lr_m, X_m, [0],
                        target=None)
    for target in (-1, 100):
        assert_raises_regex(ValueError,
                            r'target must be in \[0, n_tasks\]',
                            plot_partial_dependence, lr_m, X_m, [0],
                            target=target)

    for feature_names in (None, ['abcd', 'def']):
        assert_raises_regex(ValueError,
                            'Feature foobar not in feature_names',
                            plot_partial_dependence, lr, X,
                            features=['foobar'],
                            feature_names=feature_names)

    for features in([(1, 2, 3)], [1, {}], [tuple()]):
        assert_raises_regex(ValueError,
                            'Each entry in features must be either an int, ',
                            plot_partial_dependence, lr, X,
                            features=features)

    assert_raises_regex(ValueError,
                        'All entries of features must be less than ',
                        plot_partial_dependence, lr, X,
                        features=[123],
                        feature_names=['blah'])
