import pytest

from sklearn.datasets import load_boston
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.inspection import (
    plot_individual_conditional_expectation, plot_partial_dependence)


@pytest.fixture(scope="module")
def boston():
    return load_boston()


@pytest.fixture(scope="module")
def clf_boston(boston):
    clf = GradientBoostingRegressor(n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)
    return clf


@pytest.mark.parametrize("plot_method, nrows, ncols", [
    (plot_partial_dependence, 2, 2),
    (plot_individual_conditional_expectation, 3, 1),
])
def test_plot_pdp_ice_incorrect_num_axes(pyplot, clf_boston, boston,
                                         plot_method, nrows, ncols):
    grid_resolution = 5
    fig, axes = pyplot.subplots(nrows, ncols)
    axes_formats = [list(axes.ravel()), tuple(axes.ravel()), axes]

    msg = "Expected ax to have 2 axes, got {}".format(nrows * ncols)

    disp = plot_method(
        clf_boston, boston.data, ['CRIM', 'ZN'],
        grid_resolution=grid_resolution, feature_names=boston.feature_names
    )

    for ax_format in axes_formats:
        with pytest.raises(ValueError, match=msg):
            plot_method(
                clf_boston, boston.data, ['CRIM', 'ZN'],
                grid_resolution=grid_resolution,
                feature_names=boston.feature_names, ax=ax_format
            )

        # with axes object
        with pytest.raises(ValueError, match=msg):
            disp.plot(ax=ax_format)


@pytest.mark.parametrize("plot_method", [
    plot_partial_dependence, plot_individual_conditional_expectation
])
def test_plot_partial_dependence_with_same_axes(pyplot, clf_boston, boston,
                                                plot_method):
    # The first call to plot_pdp or plot_ice method will create two new axes
    # to place in the space of the passed in axes, which results in a total
    # of three axes in the figure.
    # Currently the API does not allow for the second call to plot_pdp or
    # plot_ice to use the same axes again, because it will create two new
    # axes in the space resulting in five axes.
    # To get the expected behavior one needs to pass the generated axes into
    # the second call:
    # disp1 = plot_individual_conditional_expectation(...)
    # disp2 = plot_individual_conditional_expectation(..., ax=disp1.axes_)

    grid_resolution = 25
    fig, ax = pyplot.subplots()
    plot_method(
        clf_boston, boston.data, ['CRIM', 'ZN'],
        grid_resolution=grid_resolution, feature_names=boston.feature_names,
        ax=ax
    )

    msg = ("The ax was already used in another plot function, please set "
           "ax=display.axes_ instead")

    with pytest.raises(ValueError, match=msg):
        plot_method(
            clf_boston, boston.data, ['CRIM', 'ZN'],
            grid_resolution=grid_resolution, feature_names=boston.feature_names,
            ax=ax
        )
