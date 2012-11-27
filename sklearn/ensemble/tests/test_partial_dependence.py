"""
Testing for the partial dependence module.
"""

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import if_matplotlib
from sklearn.ensemble.partial_dependence import partial_dependence
from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn import datasets


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

rng = np.random.RandomState(0)
# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_partial_dependence_classifier():
    """Test partial dependence for classifier """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(X, y)

    pdp, axes = partial_dependence(clf, [0], X=X,
                                                     grid_resolution=5)

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
    """Test partial dependence for multi-class classifier """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    clf.fit(iris.data, iris.target)

    grid_resolution = 25
    n_classes = clf.n_classes_
    pdp, axes = partial_dependence(
        clf, [0], X=iris.data, grid_resolution=grid_resolution)

    assert pdp.shape == (n_classes, grid_resolution)
    assert len(axes) == 1
    assert axes[0].shape[0] == grid_resolution


def test_partial_dependence_regressor():
    """Test partial dependence for regressor """
    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    pdp, axes = partial_dependence(
        clf, [0], X=boston.data, grid_resolution=grid_resolution)

    assert pdp.shape == (1, grid_resolution)
    assert axes[0].shape[0] == grid_resolution


def test_partial_dependecy_input():
    """Test input validation of partial dependence. """
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
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
    """Test partial dependence plot function. """
    clf = GradientBoostingRegressor(n_estimators=100, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, boston.data, [0, 1, (0, 1)],
                                       grid_resolution=grid_resolution)
    assert len(axs) == 3
    assert all(ax.has_data for ax in axs)
