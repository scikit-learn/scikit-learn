"""
Testing for the partial dependence module.
"""

import pytest

from sklearn.utils.testing import if_matplotlib
from sklearn.plot import plot_partial_dependence
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import load_boston, load_iris
from sklearn.datasets import make_classification, make_regression


# (X, y), n_targets  <-- as expected in the output of partial_dep()
binary_classification_data = (make_classification(random_state=0), 1)
multiclass_classification_data = (make_classification(n_classes=3,
                                                      n_clusters_per_class=1,
                                                      random_state=0), 3)
regression_data = (make_regression(random_state=0), 1)
multioutput_regression_data = (make_regression(n_targets=2, random_state=0), 2)


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


@if_matplotlib
def test_plot_partial_dependence_multioutput():
    # Test partial dependence plot function on multi-output input.
    (X, y), _ = multioutput_regression_data
    clf = LinearRegression()
    clf.fit(X, y)

    grid_resolution = 25
    fig, axs = plot_partial_dependence(clf, X, [0, 1],
                                       target=0,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)

    fig, axs = plot_partial_dependence(clf, X, [0, 1],
                                       target=1,
                                       grid_resolution=grid_resolution)
    assert len(axs) == 2
    assert all(ax.has_data for ax in axs)


@if_matplotlib
@pytest.mark.filterwarnings('ignore:Default solver will be changed ')  # 0.22
@pytest.mark.filterwarnings('ignore:Default multi_class will be')  # 0.22
def test_plot_partial_dependence_input():
    X, y = make_classification(random_state=0)

    lr = LinearRegression()
    lr.fit(X, y)
    gbc = GradientBoostingClassifier(random_state=0)
    gbc.fit(X, y)

    # check target param for multiclass
    (X_m, y_m), _ = multiclass_classification_data
    lr_m = LogisticRegression()
    lr_m.fit(X_m, y_m)
    with pytest.raises(
            ValueError,
            match='target must be specified for multi-class'):
        plot_partial_dependence(lr_m, X_m, [0], target=None)
    for target in (-1, 100):
        with pytest.raises(
                ValueError,
                match='target not in est.classes_'):
            plot_partial_dependence(lr_m, X_m, [0], target=target)

    # check target param for multioutput
    (X_m, y_m), _ = multioutput_regression_data
    lr_m = LinearRegression()
    lr_m.fit(X_m, y_m)
    with pytest.raises(
            ValueError,
            match='target must be specified for multi-output'):
        plot_partial_dependence(lr_m, X_m, [0], target=None)
    for target in (-1, 100):
        with pytest.raises(
                ValueError,
                match=r'target must be in \[0, n_tasks\]'):
            plot_partial_dependence(lr_m, X_m, [0], target=target)

    for feature_names in (None, ['abcd', 'def']):
        with pytest.raises(
                ValueError,
                match='Feature foobar not in feature_names'):
            plot_partial_dependence(lr, X, features=['foobar'],
                                    feature_names=feature_names)

    for features in([(1, 2, 3)], [1, {}], [tuple()]):
        with pytest.raises(
                ValueError,
                match='Each entry in features must be either an int, '):
            plot_partial_dependence(lr, X, features=features)

    with pytest.raises(
            ValueError,
            match='All entries of features must be less than '):
        plot_partial_dependence(lr, X, features=[123],
                                feature_names=['blah'])

    with pytest.raises(
            ValueError,
            match='feature_names should not contain duplicates'):
        plot_partial_dependence(lr, X, features=[0, 1, 2],
                                feature_names=['a', 'b', 'a'])


@if_matplotlib
def test_plot_partial_dependence_fig():
    # Make sure fig object is correctly used if not None

    import matplotlib.pyplot as plt

    (X, y), _ = regression_data
    clf = LinearRegression()
    clf.fit(X, y)

    fig = plt.figure()
    grid_resolution = 25
    returned_fig, axs = plot_partial_dependence(
        clf, X, [0, 1], target=0, grid_resolution=grid_resolution, fig=fig)

    assert returned_fig is fig
