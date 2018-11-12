"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
import pytest

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
from sklearn.base import is_classifier, is_regressor

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




@pytest.mark.parametrize('Estimator', all_estimators())
@pytest.mark.parametrize('method', ('recursion', 'exact'))
@pytest.mark.parametrize('multiclass', (True, False))
@pytest.mark.parametrize('grid_resolution', (10,))
@pytest.mark.parametrize('target_variables', ([1], [1, 2]))
def test_output_shape(Estimator, method, multiclass, grid_resolution,
                      target_variables):
    # Check that partial_dependence has consistent output shape

    name, Estimator = Estimator
    est = Estimator()
    if not (is_classifier(est) or is_regressor(est)):
        return
    if method == 'recursion':
        # recursion method only accepts some of the ensemble estimators
        if not (isinstance(est, BaseGradientBoosting) or
                isinstance(est, ForestRegressor)):
            return
    else:
        # why skip multitask? We shouldn't
        if 'MultiTask' in name:
            return
        # classifiers with exact method need predict_proba()
        if is_classifier(est) and not hasattr(est, 'predict_proba'):
            return

    if is_classifier(est):
        if multiclass == True:
            X, y = iris.data, iris.target
            n_classes = 3
        else:
            X, y = X_c, y_c
            n_classes = 1
    else:  # regressor
            X, y = X_r, y_r
            n_classes = 1

    est.fit(X, y)
    pdp, axes = partial_dependence(est,
                                    target_variables=target_variables,
                                    X=X,
                                    method=method,
                                    grid_resolution=grid_resolution)

    expected_pdp_shape = (n_classes,
                            grid_resolution ** len(target_variables))
    expected_axes_shape = (len(target_variables), grid_resolution)

    assert pdp.shape == expected_pdp_shape
    assert axes is not None
    assert np.asarray(axes).shape == expected_axes_shape


def test_grid_from_X():
    pass


def test_partial_dependence_classifier():
    # Test partial dependence for classifier
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)

    pdp, axes = partial_dependence(clf, [0], X=X, grid_resolution=5)

    # only 4 grid points instead of 5 because only 4 unique X[:,0] vals
    assert pdp.shape == (1, 4)
    assert np.asarray(axes).shape == (1, 4)

    # now with our own grid
    X_ = np.asarray(X)
    grid = np.unique(X_[:, 0])
    pdp_2, axes = partial_dependence(clf, [0], grid=grid)

    assert axes is None
    assert_array_almost_equal(pdp, pdp_2)


def test_partial_dependency_input():
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
