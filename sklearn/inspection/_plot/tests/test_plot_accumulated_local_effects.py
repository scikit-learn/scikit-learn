import numpy as np
import pytest

from sklearn.datasets import make_classification, make_regression
from sklearn.inspection import AccumulatedLocalEffectsDisplay, accumulated_local_effects
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.utils._testing import assert_allclose


@pytest.fixture(scope="module")
def regression_data():
    X, y = make_regression(n_samples=50, n_features=4, random_state=0)
    X[:, 3] = np.random.RandomState(0).randint(0, 4, size=50)
    return X, y, LinearRegression().fit(X, y)


def test_plot_continuous_and_categorical(pyplot, regression_data):
    X, _, reg = regression_data
    display = AccumulatedLocalEffectsDisplay.from_estimator(
        reg,
        X,
        [0, "x3"],
        categorical_features=[3],
        grid_resolution=10,
        line_kw={"color": "red"},
    )

    assert display.axes_.shape == (1, 2)
    line, bars = display.lines_[0, 0], display.bars_[0, 1]
    assert display.bars_[0, 0] is None and display.lines_[0, 1] is None
    assert line.get_color() == "red"

    expected = accumulated_local_effects(reg, X, 0, grid_resolution=10)
    assert_allclose(line.get_xdata(), expected.grid_values[0])
    assert_allclose(line.get_ydata(), expected.effects[0])
    assert len(bars) == 4

    assert display.axes_[0, 0].get_xlabel() == "x0"
    assert display.axes_[0, 0].get_ylabel() == "Accumulated local effect"
    assert display.axes_[0, 1].get_xlabel() == "x3"
    assert display.axes_[0, 1].get_ylabel() == ""


def test_plot_shared_y_limits_and_grid_ticks(pyplot):
    # As in PartialDependenceDisplay, single-feature plots share their y limits and
    # only the first column shows the y label and tick labels.
    rng = np.random.RandomState(0)
    X = rng.uniform(-1, 1, size=(200, 2))
    reg = LinearRegression().fit(X, 10 * X[:, 0] + 0.1 * X[:, 1])
    display = AccumulatedLocalEffectsDisplay.from_estimator(
        reg, X, [0, 1], grid_resolution=6
    )
    first, second = display.axes_[0]
    assert first.get_ylim() == second.get_ylim()
    assert first.get_ylabel() == "Accumulated local effect"
    assert second.get_ylabel() == ""
    assert all(label.get_text() == "" for label in second.get_yticklabels())

    # One tick at each grid point, which are quantiles of the data.
    grid = accumulated_local_effects(reg, X, 0, grid_resolution=6).grid_values[0]
    ticks = first.collections[0].get_segments()
    assert_allclose([segment[0][0] for segment in ticks], grid)


def test_plot_grid_layout(pyplot, regression_data):
    X, _, reg = regression_data
    display = AccumulatedLocalEffectsDisplay.from_estimator(reg, X, [0, 1, 2], n_cols=2)
    assert display.axes_.shape == (2, 2)
    assert not display.axes_[1, 1].axison
    assert display.axes_[1, 0].get_ylabel() == "Accumulated local effect"


def test_plot_with_axes(pyplot, regression_data):
    X, _, reg = regression_data
    _, axes = pyplot.subplots(1, 2)
    display = AccumulatedLocalEffectsDisplay.from_estimator(reg, X, [0, 1], ax=axes)
    assert display.figure_ is axes[0].figure
    assert display.axes_.shape == (1, 2)

    with pytest.raises(ValueError, match="Expected ax to have 3 axes"):
        AccumulatedLocalEffectsDisplay.from_estimator(reg, X, [0, 1, 2], ax=axes)


def test_plot_multiclass_target(pyplot):
    X, y = make_classification(
        n_samples=50, n_classes=3, n_informative=3, random_state=0
    )
    clf = LogisticRegression().fit(X, y)
    display = AccumulatedLocalEffectsDisplay.from_estimator(clf, X, [0], target=2)
    expected = accumulated_local_effects(clf, X, 0)
    assert_allclose(display.lines_[0, 0].get_ydata(), expected.effects[2])


def test_plot_multioutput_target(pyplot):
    X, y = make_regression(n_samples=50, n_features=3, n_targets=2, random_state=0)
    reg = LinearRegression().fit(X, y)
    display = AccumulatedLocalEffectsDisplay.from_estimator(reg, X, [0], target=1)
    expected = accumulated_local_effects(reg, X, 0)
    assert display.target_idx == 1
    assert_allclose(display.lines_[0, 0].get_ydata(), expected.effects[1])


@pytest.mark.parametrize(
    "model, params, err_msg",
    [
        (
            "regression",
            {"features": [(0, 1, 2)]},
            "Each entry in features must be an int, a string, or a pair",
        ),
        ("multiclass", {"features": [0], "target": 5}, "target not in est.classes_"),
        ("multiclass", {"features": [0]}, "target must be specified for multi-class"),
        (
            "multioutput",
            {"features": [0]},
            "target must be specified for multi-output regressors",
        ),
        ("multioutput", {"features": [0], "target": 2}, "target must be in"),
    ],
)
def test_plot_error(pyplot, model, params, err_msg):
    X, y = make_regression(n_samples=50, n_features=3, n_targets=2, random_state=0)
    if model == "regression":
        estimator = LinearRegression().fit(X, y[:, 0])
    elif model == "multiclass":
        estimator = LogisticRegression().fit(X, np.digitize(X[:, 0], [-0.5, 0.5]))
    else:
        estimator = LinearRegression().fit(X, y)
    with pytest.raises(ValueError, match=err_msg):
        AccumulatedLocalEffectsDisplay.from_estimator(estimator, X, **params)


def test_plot_two_way(pyplot):
    rng = np.random.RandomState(0)
    X = rng.uniform(0, 1, size=(300, 3))
    X[:, 1] = X[:, 0] + rng.normal(scale=0.05, size=300)
    reg = LinearRegression().fit(X, X[:, 0] * X[:, 1])

    display = AccumulatedLocalEffectsDisplay.from_estimator(
        reg, X, [0, (0, 1)], grid_resolution=6, heatmap_kw={"cmap": "viridis"}
    )
    assert display.heatmaps_[0, 0] is None
    heatmap = display.heatmaps_[0, 1]
    assert display.lines_[0, 1] is None and display.bars_[0, 1] is None

    expected = accumulated_local_effects(reg, X, [0, 1], grid_resolution=6)
    empty = expected.bin_sizes == 0
    assert empty.any()
    assert np.ma.getmaskarray(heatmap.get_array()).sum() == empty.sum()
    assert heatmap.get_cmap().name == "viridis"
    assert display.axes_[0, 1].get_xlabel() == "x0"
    assert display.axes_[0, 1].get_ylabel() == "x1"
