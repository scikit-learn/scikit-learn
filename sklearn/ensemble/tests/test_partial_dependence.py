"""
Testing for the partial dependence module.
"""
import pytest

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

from sklearn.utils.testing import assert_raises
from sklearn.ensemble.partial_dependence import partial_dependence
from sklearn.ensemble.partial_dependence import plot_partial_dependence
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn import datasets
from sklearn.utils.testing import ignore_warnings


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
sample_weight = [1, 1, 1, 2, 2, 2]

# also load the boston dataset
boston = datasets.load_boston()

# also load the iris dataset
iris = datasets.load_iris()


@ignore_warnings(category=DeprecationWarning)
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

    # with trivial (no-op) sample weights
    clf.fit(X, y, sample_weight=np.ones(len(y)))

    pdp_w, axes_w = partial_dependence(clf, [0], X=X, grid_resolution=5)

    assert pdp_w.shape == (1, 4)
    assert axes_w[0].shape[0] == 4
    assert_allclose(pdp_w, pdp)

    # with non-trivial sample weights
    clf.fit(X, y, sample_weight=sample_weight)

    pdp_w2, axes_w2 = partial_dependence(clf, [0], X=X, grid_resolution=5)

    assert pdp_w2.shape == (1, 4)
    assert axes_w2[0].shape[0] == 4
    assert np.all(np.abs(pdp_w2 - pdp_w) / np.abs(pdp_w) > 0.1)


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
def test_partial_dependence_regressor():
    # Test partial dependence for regressor
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)

    grid_resolution = 25
    pdp, axes = partial_dependence(
        clf, [0], X=boston.data, grid_resolution=grid_resolution)

    assert pdp.shape == (1, grid_resolution)
    assert axes[0].shape[0] == grid_resolution


@ignore_warnings(category=DeprecationWarning)
def test_partial_dependence_sample_weight():
    # Test near perfect correlation between partial dependence and diagonal
    # when sample weights emphasize y = x predictions
    N = 1000
    rng = np.random.RandomState(123456)
    mask = rng.randint(2, size=N, dtype=bool)

    x = rng.rand(N)
    # set y = x on mask and y = -x outside
    y = x.copy()
    y[~mask] = -y[~mask]
    X = np.c_[mask, x]
    # sample weights to emphasize data points where y = x
    sample_weight = np.ones(N)
    sample_weight[mask] = 1000.

    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(X, y, sample_weight=sample_weight)

    grid = np.arange(0, 1, 0.01)
    pdp = partial_dependence(clf, [1], grid=grid)

    assert np.corrcoef(np.ravel(pdp[0]), grid)[0, 1] > 0.99


@ignore_warnings(category=DeprecationWarning)
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


@ignore_warnings(category=DeprecationWarning)
@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
# matplotlib Python3.7 warning
def test_plot_partial_dependence(pyplot):
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


@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
# matplotlib Python3.7 warning
@ignore_warnings(category=DeprecationWarning)
def test_plot_partial_dependence_input(pyplot):
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


@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
# matplotlib Python3.7 warning
@ignore_warnings(category=DeprecationWarning)
def test_plot_partial_dependence_multiclass(pyplot):
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


@pytest.mark.parametrize(
    "func, params",
    [(partial_dependence, {'target_variables': [0], 'X': boston.data}),
     (plot_partial_dependence, {'X': boston.data, 'features': [0, 1, (0, 1)]})]
)
def test_raise_deprecation_warning(pyplot, func, params):
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)
    grid_resolution = 25

    warn_msg = "The function ensemble.{} has been deprecated".format(
        func.__name__
    )
    with pytest.warns(DeprecationWarning, match=warn_msg):
        func(clf, **params, grid_resolution=grid_resolution)
