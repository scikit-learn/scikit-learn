import numpy as np
from scipy.stats.mstats import mquantiles

import pytest
from numpy.testing import assert_allclose

from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LinearRegression
from sklearn.inspection import plot_partial_dependence


@pytest.fixture(scope="module")
def boston():
    return load_boston()


@pytest.fixture(scope="module")
def clf_boston(boston):
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)
    return clf


@pytest.mark.parametrize("grid_resolution", [10, 20])
def test_plot_partial_dependence(grid_resolution, pyplot, clf_boston, boston):
    # Test partial dependence plot function.
    feature_names = boston.feature_names
    disp = plot_partial_dependence(clf_boston, boston.data,
                                   [0, 1, (0, 1)],
                                   grid_resolution=grid_resolution,
                                   feature_names=feature_names,
                                   contour_kw={"cmap": "jet"})
    fig = pyplot.gcf()
    axs = fig.get_axes()
    assert disp.figure_ is fig
    assert len(axs) == 4

    assert disp.bounding_ax_ is not None
    assert disp.axes_.shape == (1, 3)
    assert disp.lines_.shape == (1, 3)
    assert disp.contours_.shape == (1, 3)

    assert disp.lines_[0, 2] is None
    assert disp.contours_[0, 0] is None
    assert disp.contours_[0, 1] is None

    assert disp.features == [(0, ), (1, ), (0, 1)]
    assert np.all(disp.feature_names == feature_names)
    assert len(disp.deciles) == 2
    for i in [0, 1]:
        assert_allclose(disp.deciles[i],
                        mquantiles(boston.data[:, i],
                                   prob=np.arange(0.1, 1.0, 0.1)))

    single_feature_positions = [(0, 0), (0, 1)]

    for i, pos in enumerate(single_feature_positions):
        ax = disp.axes_[pos]
        assert ax.get_xlabel() == boston.feature_names[i]
        assert ax.get_ylabel() == "Partial dependence"
        assert_allclose(ax.get_ylim(), disp.pdp_lim[1])

        line = disp.lines_[pos]

        avg_preds, values = disp.pd_results[i]
        assert avg_preds.shape == (1, grid_resolution)
        target_idx = disp.target_idx

        line_data = line.get_data()
        assert_allclose(line_data[0], values[0])
        assert_allclose(line_data[1], avg_preds[target_idx].ravel())

    # two feature position
    ax = disp.axes_[0, 2]
    coutour = disp.contours_[0, 2]
    expected_levels = np.linspace(*disp.pdp_lim[2], num=8)
    assert_allclose(coutour.levels, expected_levels)
    assert coutour.get_cmap().name == "jet"
    assert ax.get_xlabel() == boston.feature_names[0]
    assert ax.get_ylabel() == boston.feature_names[1]


def test_plot_partial_dependence_str_features(pyplot, clf_boston, boston):
    grid_resolution = 25
    # check with str features and array feature names and single column
    disp = plot_partial_dependence(clf_boston, boston.data,
                                   [('CRIM', 'ZN'), 'ZN'],
                                   grid_resolution=grid_resolution,
                                   feature_names=boston.feature_names,
                                   n_cols=1, line_kw={"alpha": 0.8})
    fig = pyplot.gcf()
    axs = fig.get_axes()
    assert len(axs) == 3

    assert disp.figure_ is fig
    assert disp.axes_.shape == (2, 1)
    assert disp.lines_.shape == (2, 1)
    assert disp.contours_.shape == (2, 1)

    assert disp.lines_[0, 0] is None
    assert disp.contours_[1, 0] is None

    # line
    ax = disp.axes_[1, 0]
    assert ax.get_xlabel() == "ZN"
    assert ax.get_ylabel() == "Partial dependence"

    line = disp.lines_[1, 0]
    avg_preds, values = disp.pd_results[1]
    target_idx = disp.target_idx
    assert line.get_alpha() == 0.8

    line_data = line.get_data()
    assert_allclose(line_data[0], values[0])
    assert_allclose(line_data[1], avg_preds[target_idx].ravel())

    # contour
    ax = disp.axes_[0, 0]
    coutour = disp.contours_[0, 0]
    expect_levels = np.linspace(*disp.pdp_lim[2], num=8)
    assert_allclose(coutour.levels, expect_levels)
    assert ax.get_xlabel() == "CRIM"
    assert ax.get_ylabel() == "ZN"


def test_plot_partial_dependence_custom_axes(pyplot, clf_boston, boston):
    grid_resolution = 25
    fig, (ax1, ax2) = pyplot.subplots(1, 2)
    feature_names = boston.feature_names.tolist()
    disp = plot_partial_dependence(clf_boston, boston.data,
                                   ['CRIM', ('CRIM', 'ZN')],
                                   grid_resolution=grid_resolution,
                                   feature_names=feature_names, ax=[ax1, ax2])
    assert fig is disp.figure_
    assert disp.bounding_ax_ is None
    assert disp.axes_.shape == (2, )
    assert disp.axes_[0] is ax1
    assert disp.axes_[1] is ax2

    ax = disp.axes_[0]
    assert ax.get_xlabel() == "CRIM"
    assert ax.get_ylabel() == "Partial dependence"

    line = disp.lines_[0]
    avg_preds, values = disp.pd_results[0]
    target_idx = disp.target_idx

    line_data = line.get_data()
    assert_allclose(line_data[0], values[0])
    assert_allclose(line_data[1], avg_preds[target_idx].ravel())

    # contour
    ax = disp.axes_[1]
    coutour = disp.contours_[1]
    expect_levels = np.linspace(*disp.pdp_lim[2], num=8)
    assert_allclose(coutour.levels, expect_levels)
    assert ax.get_xlabel() == "CRIM"
    assert ax.get_ylabel() == "ZN"


def test_plot_partial_dependence_passing_numpy_axes(pyplot, clf_boston,
                                                    boston):
    grid_resolution = 25
    feature_names = boston.feature_names.tolist()
    disp1 = plot_partial_dependence(clf_boston, boston.data,
                                    ['CRIM', 'ZN'],
                                    grid_resolution=grid_resolution,
                                    feature_names=feature_names)
    assert disp1.axes_.shape == (1, 2)
    assert len(disp1.axes_[0, 0].get_lines()) == 1
    assert len(disp1.axes_[0, 1].get_lines()) == 1

    lr = LinearRegression()
    lr.fit(boston.data, boston.target)

    disp2 = plot_partial_dependence(lr, boston.data,
                                    ['CRIM', 'ZN'],
                                    grid_resolution=grid_resolution,
                                    feature_names=feature_names,
                                    ax=disp1.axes_)

    assert np.all(disp1.axes_ == disp2.axes_)
    assert len(disp2.axes_[0, 0].get_lines()) == 2
    assert len(disp2.axes_[0, 1].get_lines()) == 2


def test_plot_partial_dependence_incorrent_num_axes(pyplot, clf_boston,
                                                    boston):
    grid_resolution = 25
    fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3)

    msg = r"Expected len\(ax\) == len\(features\), got len\(ax\) = 3"
    with pytest.raises(ValueError, match=msg):
        plot_partial_dependence(clf_boston, boston.data,
                                ['CRIM', ('CRIM', 'ZN')],
                                grid_resolution=grid_resolution,
                                feature_names=boston.feature_names,
                                ax=[ax1, ax2, ax3])

    disp = plot_partial_dependence(clf_boston, boston.data,
                                   ['CRIM', ('CRIM', 'ZN')],
                                   grid_resolution=grid_resolution,
                                   feature_names=boston.feature_names)

    with pytest.raises(ValueError, match=msg):
        disp.plot(ax=[ax1, ax2, ax3])


def test_plot_partial_dependence_with_same_axes(pyplot, clf_boston, boston):
    # The first call to `plot_*` will plot the axes

    grid_resolution = 25
    fig, ax = pyplot.subplots()
    plot_partial_dependence(clf_boston, boston.data, ['CRIM', 'ZN'],
                            grid_resolution=grid_resolution,
                            feature_names=boston.feature_names, ax=ax)

    msg = ("The ax was already used in another plot function, please set "
           "ax=display.axes_ instead")

    with pytest.raises(ValueError, match=msg):
        plot_partial_dependence(clf_boston, boston.data,
                                ['CRIM', 'ZN'],
                                grid_resolution=grid_resolution,
                                feature_names=boston.feature_names, ax=ax)


def test_plot_partial_dependence_multiclass(pyplot):
    grid_resolution = 25
    clf_int = GradientBoostingClassifier(n_estimators=10, random_state=1)
    iris = load_iris()

    # Test partial dependence plot function on multi-class input.
    clf_int.fit(iris.data, iris.target)
    disp_target_0 = plot_partial_dependence(clf_int, iris.data, [0, 1],
                                            target=0,
                                            grid_resolution=grid_resolution)
    assert disp_target_0.figure_ is pyplot.gcf()
    assert disp_target_0.axes_.shape == (1, 2)
    assert disp_target_0.lines_.shape == (1, 2)
    assert disp_target_0.contours_.shape == (1, 2)
    assert all(c is None for c in disp_target_0.contours_.flat)
    assert disp_target_0.target_idx == 0

    # now with symbol labels
    target = iris.target_names[iris.target]
    clf_symbol = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf_symbol.fit(iris.data, target)
    disp_symbol = plot_partial_dependence(clf_symbol, iris.data, [0, 1],
                                          target='setosa',
                                          grid_resolution=grid_resolution)
    assert disp_symbol.figure_ is pyplot.gcf()
    assert disp_symbol.axes_.shape == (1, 2)
    assert disp_symbol.lines_.shape == (1, 2)
    assert disp_symbol.contours_.shape == (1, 2)
    assert all(c is None for c in disp_symbol.contours_.flat)
    assert disp_symbol.target_idx == 0

    for int_result, symbol_result in zip(disp_target_0.pd_results,
                                         disp_symbol.pd_results):
        avg_preds_int, values_int = int_result
        avg_preds_symbol, values_symbol = symbol_result
        assert_allclose(avg_preds_int, avg_preds_symbol)
        assert_allclose(values_int, values_symbol)

    # check that the pd plots are different for another target
    disp_target_1 = plot_partial_dependence(clf_int, iris.data, [0, 1],
                                            target=1,
                                            grid_resolution=grid_resolution)
    target_0_data_y = disp_target_0.lines_[0, 0].get_data()[1]
    target_1_data_y = disp_target_1.lines_[0, 0].get_data()[1]
    assert any(target_0_data_y != target_1_data_y)


multioutput_regression_data = make_regression(n_samples=50, n_targets=2,
                                              random_state=0)


@pytest.mark.parametrize("target", [0, 1])
def test_plot_partial_dependence_multioutput(pyplot, target):
    # Test partial dependence plot function on multi-output input.
    X, y = multioutput_regression_data
    clf = LinearRegression().fit(X, y)

    grid_resolution = 25
    disp = plot_partial_dependence(clf, X, [0, 1], target=target,
                                   grid_resolution=grid_resolution)
    fig = pyplot.gcf()
    axs = fig.get_axes()
    assert len(axs) == 3
    assert disp.target_idx == target
    assert disp.bounding_ax_ is not None

    positions = [(0, 0), (0, 1)]

    for i, pos in enumerate(positions):
        ax = disp.axes_[pos]
        assert ax.get_ylabel() == "Partial dependence"
        assert ax.get_xlabel() == "{}".format(i)


dummy_classification_data = make_classification(random_state=0)


@pytest.mark.parametrize(
    "data, params, err_msg",
    [(multioutput_regression_data, {"target": None, 'features': [0]},
      "target must be specified for multi-output"),
     (multioutput_regression_data, {"target": -1, 'features': [0]},
      r'target must be in \[0, n_tasks\]'),
     (multioutput_regression_data, {"target": 100, 'features': [0]},
      r'target must be in \[0, n_tasks\]'),
     (dummy_classification_data,
     {'features': ['foobar'], 'feature_names': None},
     'Feature foobar not in feature_names'),
     (dummy_classification_data,
     {'features': ['foobar'], 'feature_names': ['abcd', 'def']},
      'Feature foobar not in feature_names'),
     (dummy_classification_data, {'features': [(1, 2, 3)]},
      'Each entry in features must be either an int, '),
     (dummy_classification_data, {'features': [1, {}]},
      'Each entry in features must be either an int, '),
     (dummy_classification_data, {'features': [tuple()]},
      'Each entry in features must be either an int, '),
     (dummy_classification_data,
      {'features': [123], 'feature_names': ['blahblah']},
      'All entries of features must be less than '),
     (dummy_classification_data,
      {'features': [0, 1, 2], 'feature_names': ['a', 'b', 'a']},
      'feature_names should not contain duplicates')]
)
def test_plot_partial_dependence_error(pyplot, data, params, err_msg):
    X, y = data
    estimator = LinearRegression().fit(X, y)

    with pytest.raises(ValueError, match=err_msg):
        plot_partial_dependence(estimator, X, **params)


@pytest.mark.parametrize("params, err_msg", [
    ({'target': 4, 'features': [0]},
     'target not in est.classes_, got 4'),
    ({'target': None, 'features': [0]},
     'target must be specified for multi-class'),
    ({'target': 1, 'features': [4.5]},
     'Each entry in features must be either an int,'),
])
def test_plot_partial_dependence_multiclass_error(pyplot, params, err_msg):
    iris = load_iris()
    clf = GradientBoostingClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)

    with pytest.raises(ValueError, match=err_msg):
        plot_partial_dependence(clf, iris.data, **params)


def test_plot_partial_dependence_fig_deprecated(pyplot):
    # Make sure fig object is correctly used if not None
    X, y = make_regression(n_samples=50, random_state=0)
    clf = LinearRegression()
    clf.fit(X, y)

    fig = pyplot.figure()
    grid_resolution = 25

    msg = ("The fig parameter is deprecated in version 0.22 and will be "
           "removed in version 0.24")
    with pytest.warns(DeprecationWarning, match=msg):
        plot_partial_dependence(
            clf, X, [0, 1], target=0, grid_resolution=grid_resolution, fig=fig)

    assert pyplot.gcf() is fig
