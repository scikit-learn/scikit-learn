"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import if_matplotlib
from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import ignore_warnings
from sklearn.partial_dependence import partial_dependence
from sklearn.partial_dependence import plot_partial_dependence
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble.gradient_boosting import BaseGradientBoosting
from sklearn.ensemble.forest import ForestRegressor
from sklearn.datasets import load_boston, load_iris
from sklearn.datasets import make_classification, make_regression

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]

# Make some sample data to test output shapes
X_c, y_c = make_classification(n_features=10, n_informative=5, random_state=0)
# Non-negative for MultinomialNB
X_c = X_c + np.abs(X_c.min())
X_r, y_r = make_regression(n_features=10, n_informative=5, random_state=0)

# Load the boston & iris datasets
boston = load_boston()
iris = load_iris()


@ignore_warnings()
def test_output_shape_recursion():
    # Test recursion partial_dependence has same output shape for everything
    for name, Estimator in all_estimators():
        est = Estimator()
        if not (isinstance(est, BaseGradientBoosting) or
                isinstance(est, ForestRegressor)):
            continue
        if est._estimator_type == 'classifier':
            est.fit(X_c, y_c)
            pdp, axes = partial_dependence(est,
                                           target_variables=[1],
                                           X=X_c,
                                           method='recursion',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 10))
            pdp, axes = partial_dependence(est,
                                           target_variables=[1, 2],
                                           X=X_c,
                                           method='recursion',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 100))
        elif est._estimator_type == 'regressor':
            est.fit(X_r, y_r)
            pdp, axes = partial_dependence(est,
                                           target_variables=[1],
                                           X=X_r,
                                           method='recursion',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 10))
            pdp, axes = partial_dependence(est,
                                           target_variables=[1, 2],
                                           X=X_r,
                                           method='recursion',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 100))


@ignore_warnings()
def test_output_shape_exact():
    # Test exact partial_dependence has same output shape for everything
    for name, Estimator in all_estimators():
        est = Estimator()
        if not hasattr(est, '_estimator_type') or 'MultiTask' in name:
            continue
        if est._estimator_type == 'classifier':
            if not hasattr(est, 'predict_proba'):
                continue
            est.fit(X_c, y_c)
            pdp, axes = partial_dependence(est,
                                           target_variables=[1],
                                           X=X_c,
                                           method='exact',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 10))
            pdp, axes = partial_dependence(est,
                                           target_variables=[1, 2],
                                           X=X_c,
                                           method='exact',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 100))
        elif est._estimator_type == 'regressor':
            est.fit(X_r, y_r)
            pdp, axes = partial_dependence(est,
                                           target_variables=[1],
                                           X=X_r,
                                           method='exact',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 10))
            pdp, axes = partial_dependence(est,
                                           target_variables=[1, 2],
                                           X=X_r,
                                           method='exact',
                                           grid_resolution=10)
            assert (pdp.shape == (1, 100))


def test_partial_dependence_classifier():
    # Test partial dependence for classifier
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)

    pdp, axes = partial_dependence(clf, [0], X=X, grid_resolution=5)

    # only 4 grid points instead of 5 because only 4 unique X[:,0] vals
    assert pdp.shape == (1, 4)
    assert axes[0].shape[0] == 4

    # now with our own grid
    X_ = np.asarray(X)
    grid = np.unique(X_[:, 0])
    pdp_2, axes = partial_dependence(clf, [0], grid=grid)

    assert axes is None
    assert_array_equal(pdp, pdp_2)


def test_partial_dependence_multiclass():
    # Test partial dependence for multi-class classifier
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)

    grid_resolution = 25
    n_classes = clf.n_classes_
    pdp, axes = partial_dependence(
        clf, [0], X=iris.data, grid_resolution=grid_resolution)

    assert pdp.shape == (n_classes, grid_resolution)
    assert len(axes) == 1
    assert axes[0].shape[0] == grid_resolution


def test_partial_dependence_regressor():
    # Test partial dependence for regressor
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    pdp, axes = partial_dependence(
        clf, [0], X=boston.data, grid_resolution=grid_resolution)

    assert pdp.shape == (1, grid_resolution)
    assert axes[0].shape[0] == grid_resolution


def test_partial_dependecy_input():
    # Test input validation of partial dependence.
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)

    assert_raises(ValueError, partial_dependence,
                  clf, [0], grid=None, X=None)

    assert_raises(ValueError, partial_dependence,
                  clf, [0], grid=[0, 1], X=X)

    # first argument must be an instance of BaseGradientBoosting
    assert_raises(ValueError, partial_dependence,
                  {}, [0], X=X)

    # Gradient boosting estimator must be fit
    assert_raises(ValueError, partial_dependence,
                  GradientBoostingClassifier(), [0], X=X)

    assert_raises(ValueError, partial_dependence, clf, [-1], X=X)

    assert_raises(ValueError, partial_dependence, clf, [100], X=X)

    # wrong ndim for grid
    grid = np.random.rand(10, 2, 1)
    assert_raises(ValueError, partial_dependence, clf, [0], grid=grid)


@if_matplotlib
def test_plot_partial_dependence():
    # Test partial dependence plot function.
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


@if_matplotlib
def test_plot_partial_dependence_multiclass():
    # Test partial dependence plot function on multi-class input.
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
