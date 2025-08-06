import numpy as np
import pytest

from sklearn.linear_model import LogisticRegression
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _deprecate_estimator_name,
    _despine,
    _interval_max_min_ratio,
    _validate_score_name,
    _validate_style_kwargs,
)
from sklearn.utils._response import _get_response_values_binary
from sklearn.utils._testing import assert_allclose


@pytest.mark.parametrize("ax", [None, "Ax"])
@pytest.mark.parametrize(
    "name, expected_name_out", [(None, "TestEstimator"), ("CustomName", "CustomName")]
)
def test_validate_plot_params(pyplot, ax, name, expected_name_out):
    """Check `_validate_plot_params` returns the correct values."""
    display = _BinaryClassifierCurveDisplayMixin()
    display.estimator_name = "TestEstimator"
    if ax:
        _, ax = pyplot.subplots()
    ax_out, _, name_out = display._validate_plot_params(ax=ax, name=name)

    assert name_out == expected_name_out

    if ax:
        assert ax == ax_out


@pytest.mark.parametrize("pos_label", [None, 0])
@pytest.mark.parametrize("name", [None, "CustomName"])
@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
def test_validate_and_get_response_values(pyplot, pos_label, name, response_method):
    """Check `_validate_and_get_response_values` returns the correct values."""
    X = np.array([[0, 0], [1, 1], [2, 2], [3, 3]])
    y = np.array([0, 0, 2, 2])
    estimator = LogisticRegression().fit(X, y)

    y_pred, pos_label, name_out = (
        _BinaryClassifierCurveDisplayMixin._validate_and_get_response_values(
            estimator,
            X,
            y,
            response_method=response_method,
            pos_label=pos_label,
            name=name,
        )
    )

    expected_y_pred, expected_pos_label = _get_response_values_binary(
        estimator, X, response_method=response_method, pos_label=pos_label
    )

    assert_allclose(y_pred, expected_y_pred)
    assert pos_label == expected_pos_label

    # Check name is handled correctly
    expected_name = name if name is not None else "LogisticRegression"
    assert name_out == expected_name


@pytest.mark.parametrize(
    "y_true, error_message",
    [
        (np.array([0, 1, 2]), "The target y is not binary."),
        (np.array([0, 1]), "Found input variables with inconsistent"),
        (np.array([0, 2, 0, 2]), r"y_true takes value in \{0, 2\} and pos_label"),
    ],
)
def test_validate_from_predictions_params_errors(pyplot, y_true, error_message):
    """Check `_validate_from_predictions_params` raises the correct errors."""
    y_pred = np.array([0.1, 0.2, 0.3, 0.4])
    sample_weight = np.ones(4)

    with pytest.raises(ValueError, match=error_message):
        _BinaryClassifierCurveDisplayMixin._validate_from_predictions_params(
            y_true=y_true,
            y_pred=y_pred,
            sample_weight=sample_weight,
            pos_label=None,
        )


@pytest.mark.parametrize("name", [None, "CustomName"])
@pytest.mark.parametrize(
    "pos_label, y_true",
    [
        (None, np.array([0, 1, 0, 1])),
        (2, np.array([0, 2, 0, 2])),
    ],
)
def test_validate_from_predictions_params_returns(pyplot, name, pos_label, y_true):
    """Check `_validate_from_predictions_params` returns the correct values."""
    y_pred = np.array([0.1, 0.2, 0.3, 0.4])
    pos_label_out, name_out = (
        _BinaryClassifierCurveDisplayMixin._validate_from_predictions_params(
            y_true=y_true,
            y_pred=y_pred,
            sample_weight=None,
            pos_label=pos_label,
            name=name,
        )
    )

    # Check name is handled correctly
    expected_name = name if name is not None else "Classifier"
    assert name_out == expected_name

    # Check pos_label is handled correctly
    expected_pos_label = pos_label if pos_label is not None else 1
    assert pos_label_out == expected_pos_label


@pytest.mark.parametrize(
    "params, err_msg",
    [
        (
            {
                # Missing "indices" key
                "cv_results": {"estimator": "dummy"},
                "X": np.array([[1, 2], [3, 4]]),
                "y": np.array([0, 1]),
                "sample_weight": None,
                "pos_label": None,
            },
            "`cv_results` does not contain one of the following",
        ),
        (
            {
                "cv_results": {
                    "estimator": "dummy",
                    "indices": {"test": [[1, 2], [1, 2]], "train": [[3, 4], [3, 4]]},
                },
                # `X` wrong length
                "X": np.array([[1, 2]]),
                "y": np.array([0, 1]),
                "sample_weight": None,
                "pos_label": None,
            },
            "`X` does not contain the correct number of",
        ),
        (
            {
                "cv_results": {
                    "estimator": "dummy",
                    "indices": {"test": [[1, 2], [1, 2]], "train": [[3, 4], [3, 4]]},
                },
                "X": np.array([1, 2, 3, 4]),
                # `y` not binary
                "y": np.array([0, 2, 1, 3]),
                "sample_weight": None,
                "pos_label": None,
            },
            "The target `y` is not binary",
        ),
        (
            {
                "cv_results": {
                    "estimator": "dummy",
                    "indices": {"test": [[1, 2], [1, 2]], "train": [[3, 4], [3, 4]]},
                },
                "X": np.array([1, 2, 3, 4]),
                "y": np.array([0, 1, 0, 1]),
                # `sample_weight` wrong length
                "sample_weight": np.array([0.5]),
                "pos_label": None,
            },
            "Found input variables with inconsistent",
        ),
        (
            {
                "cv_results": {
                    "estimator": "dummy",
                    "indices": {"test": [[1, 2], [1, 2]], "train": [[3, 4], [3, 4]]},
                },
                "X": np.array([1, 2, 3, 4]),
                "y": np.array([2, 3, 2, 3]),
                "sample_weight": None,
                # Not specified when `y` not in {0, 1} or {-1, 1}
                "pos_label": None,
            },
            "y takes value in {2, 3} and pos_label is not specified",
        ),
    ],
)
def test_validate_from_cv_results_params(pyplot, params, err_msg):
    """Check parameter validation is performed correctly."""
    with pytest.raises(ValueError, match=err_msg):
        _BinaryClassifierCurveDisplayMixin()._validate_from_cv_results_params(**params)


@pytest.mark.parametrize(
    "curve_legend_metric, curve_name, expected_label",
    [
        (0.85, None, "AUC = 0.85"),
        (None, "Model A", "Model A"),
        (0.95, "Random Forest", "Random Forest (AUC = 0.95)"),
        (None, None, None),
    ],
)
def test_get_legend_label(curve_legend_metric, curve_name, expected_label):
    """Check `_get_legend_label` returns the correct label."""
    legend_metric_name = "AUC"
    label = _BinaryClassifierCurveDisplayMixin._get_legend_label(
        curve_legend_metric, curve_name, legend_metric_name
    )
    assert label == expected_label


# TODO(1.9) : Remove
@pytest.mark.parametrize("curve_kwargs", [{"alpha": 1.0}, None])
@pytest.mark.parametrize("kwargs", [{}, {"alpha": 1.0}])
def test_validate_curve_kwargs_deprecate_kwargs(curve_kwargs, kwargs):
    """Check `_validate_curve_kwargs` deprecates kwargs correctly."""
    n_curves = 1
    name = None
    legend_metric = {"mean": 0.8, "std": 0.1}
    legend_metric_name = "AUC"

    if curve_kwargs and kwargs:
        with pytest.raises(ValueError, match="Cannot provide both `curve_kwargs`"):
            _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
                n_curves,
                name,
                legend_metric,
                legend_metric_name,
                curve_kwargs,
                **kwargs,
            )
    elif kwargs:
        with pytest.warns(FutureWarning, match=r"`\*\*kwargs` is deprecated and"):
            _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
                n_curves,
                name,
                legend_metric,
                legend_metric_name,
                curve_kwargs,
                **kwargs,
            )
    else:
        # No warning or error should be raised
        _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
            n_curves, name, legend_metric, legend_metric_name, curve_kwargs, **kwargs
        )


def test_validate_curve_kwargs_error():
    """Check `_validate_curve_kwargs` performs parameter validation correctly."""
    n_curves = 3
    legend_metric = {"mean": 0.8, "std": 0.1}
    legend_metric_name = "AUC"
    with pytest.raises(ValueError, match="`curve_kwargs` must be None"):
        _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
            n_curves=n_curves,
            name=None,
            legend_metric=legend_metric,
            legend_metric_name=legend_metric_name,
            curve_kwargs=[{"alpha": 1.0}],
        )
    with pytest.raises(ValueError, match="To avoid labeling individual curves"):
        name = ["one", "two", "three"]
        _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
            n_curves=n_curves,
            name=name,
            legend_metric=legend_metric,
            legend_metric_name=legend_metric_name,
            curve_kwargs=None,
        )
        _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
            n_curves=n_curves,
            name=name,
            legend_metric=legend_metric,
            legend_metric_name=legend_metric_name,
            curve_kwargs={"alpha": 1.0},
        )


@pytest.mark.parametrize("name", [None, "curve_name", ["curve_name"]])
@pytest.mark.parametrize(
    "legend_metric",
    [
        {"mean": 0.8, "std": 0.2},
        {"mean": None, "std": None},
    ],
)
@pytest.mark.parametrize("legend_metric_name", ["AUC", "AP"])
@pytest.mark.parametrize(
    "curve_kwargs",
    [
        None,
        {"color": "red"},
    ],
)
def test_validate_curve_kwargs_single_legend(
    name, legend_metric, legend_metric_name, curve_kwargs
):
    """Check `_validate_curve_kwargs` returns correct kwargs for single legend entry."""
    n_curves = 3
    curve_kwargs_out = _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
        n_curves=n_curves,
        name=name,
        legend_metric=legend_metric,
        legend_metric_name=legend_metric_name,
        curve_kwargs=curve_kwargs,
    )

    assert isinstance(curve_kwargs_out, list)
    assert len(curve_kwargs_out) == n_curves

    expected_label = None
    if isinstance(name, list):
        name = name[0]
    if name is not None:
        expected_label = name
        if legend_metric["mean"] is not None:
            expected_label = expected_label + f" ({legend_metric_name} = 0.80 +/- 0.20)"
    # `name` is None
    elif legend_metric["mean"] is not None:
        expected_label = f"{legend_metric_name} = 0.80 +/- 0.20"

    assert curve_kwargs_out[0]["label"] == expected_label
    # All remaining curves should have None as "label"
    assert curve_kwargs_out[1]["label"] is None
    assert curve_kwargs_out[2]["label"] is None

    # Default multi-curve kwargs
    if curve_kwargs is None:
        assert all(len(kwargs) == 4 for kwargs in curve_kwargs_out)
        assert all(kwargs["alpha"] == 0.5 for kwargs in curve_kwargs_out)
        assert all(kwargs["linestyle"] == "--" for kwargs in curve_kwargs_out)
        assert all(kwargs["color"] == "blue" for kwargs in curve_kwargs_out)
    else:
        assert all(len(kwargs) == 2 for kwargs in curve_kwargs_out)
        assert all(kwargs["color"] == "red" for kwargs in curve_kwargs_out)


@pytest.mark.parametrize("name", [None, "curve_name", ["one", "two", "three"]])
@pytest.mark.parametrize(
    "legend_metric", [{"metric": [1.0, 1.0, 1.0]}, {"metric": [None, None, None]}]
)
@pytest.mark.parametrize("legend_metric_name", ["AUC", "AP"])
def test_validate_curve_kwargs_multi_legend(name, legend_metric, legend_metric_name):
    """Check `_validate_curve_kwargs` returns correct kwargs for multi legend entry."""
    n_curves = 3
    curve_kwargs = [{"color": "red"}, {"color": "yellow"}, {"color": "blue"}]
    curve_kwargs_out = _BinaryClassifierCurveDisplayMixin._validate_curve_kwargs(
        n_curves=n_curves,
        name=name,
        legend_metric=legend_metric,
        legend_metric_name=legend_metric_name,
        curve_kwargs=curve_kwargs,
    )

    assert isinstance(curve_kwargs_out, list)
    assert len(curve_kwargs_out) == n_curves

    expected_labels = [None, None, None]
    if isinstance(name, str):
        expected_labels = "curve_name"
        if legend_metric["metric"][0] is not None:
            expected_labels = expected_labels + f" ({legend_metric_name} = 1.00)"
        expected_labels = [expected_labels] * n_curves
    elif isinstance(name, list) and legend_metric["metric"][0] is None:
        expected_labels = name
    elif isinstance(name, list) and legend_metric["metric"][0] is not None:
        expected_labels = [
            f"{name_single} ({legend_metric_name} = 1.00)" for name_single in name
        ]
    # `name` is None
    elif legend_metric["metric"][0] is not None:
        expected_labels = [f"{legend_metric_name} = 1.00"] * n_curves

    for idx, expected_label in enumerate(expected_labels):
        assert curve_kwargs_out[idx]["label"] == expected_label

    assert all(len(kwargs) == 2 for kwargs in curve_kwargs_out)
    for curve_kwarg, curve_kwarg_out in zip(curve_kwargs, curve_kwargs_out):
        assert curve_kwarg_out["color"] == curve_kwarg["color"]


def metric():
    pass  # pragma: no cover


def neg_metric():
    pass  # pragma: no cover


@pytest.mark.parametrize(
    "score_name, scoring, negate_score, expected_score_name",
    [
        ("accuracy", None, False, "accuracy"),  # do not transform the name
        (None, "accuracy", False, "Accuracy"),  # capitalize the name
        (None, "accuracy", True, "Negative accuracy"),  # add "Negative"
        (None, "neg_mean_absolute_error", False, "Negative mean absolute error"),
        (None, "neg_mean_absolute_error", True, "Mean absolute error"),  # remove "neg_"
        ("MAE", "neg_mean_absolute_error", True, "MAE"),  # keep score_name
        (None, None, False, "Score"),  # default name
        (None, None, True, "Negative score"),  # default name but negated
        ("Some metric", metric, False, "Some metric"),  # do not transform the name
        ("Some metric", metric, True, "Some metric"),  # do not transform the name
        (None, metric, False, "Metric"),  # default name
        (None, metric, True, "Negative metric"),  # default name but negated
        ("Some metric", neg_metric, False, "Some metric"),  # do not transform the name
        ("Some metric", neg_metric, True, "Some metric"),  # do not transform the name
        (None, neg_metric, False, "Negative metric"),  # default name
        (None, neg_metric, True, "Metric"),  # default name but negated
    ],
)
def test_validate_score_name(score_name, scoring, negate_score, expected_score_name):
    """Check that we return the right score name."""
    assert (
        _validate_score_name(score_name, scoring, negate_score) == expected_score_name
    )


# In the following test, we check the value of the max to min ratio
# for parameter value intervals to check that using a decision threshold
# of 5. is a good heuristic to decide between linear and log scales on
# common ranges of parameter values.
@pytest.mark.parametrize(
    "data, lower_bound, upper_bound",
    [
        # Such a range could be clearly displayed with either log scale or linear
        # scale.
        (np.geomspace(0.1, 1, 5), 5, 6),
        # Checking that the ratio is still positive on a negative log scale.
        (-np.geomspace(0.1, 1, 10), 7, 8),
        # Evenly spaced parameter values lead to a ratio of 1.
        (np.linspace(0, 1, 5), 0.9, 1.1),
        # This is not exactly spaced on a log scale but we will benefit from treating
        # it as such for visualization.
        ([1, 2, 5, 10, 20, 50], 20, 40),
    ],
)
def test_inverval_max_min_ratio(data, lower_bound, upper_bound):
    assert lower_bound < _interval_max_min_ratio(data) < upper_bound


@pytest.mark.parametrize(
    "default_kwargs, user_kwargs, expected",
    [
        (
            {"color": "blue", "linewidth": 2},
            {"linestyle": "dashed"},
            {"color": "blue", "linewidth": 2, "linestyle": "dashed"},
        ),
        (
            {"color": "blue", "linestyle": "solid"},
            {"c": "red", "ls": "dashed"},
            {"color": "red", "linestyle": "dashed"},
        ),
        (
            {"label": "xxx", "color": "k", "linestyle": "--"},
            {"ls": "-."},
            {"label": "xxx", "color": "k", "linestyle": "-."},
        ),
        ({}, {}, {}),
        (
            {},
            {
                "ls": "dashed",
                "c": "red",
                "ec": "black",
                "fc": "yellow",
                "lw": 2,
                "mec": "green",
                "mfcalt": "blue",
                "ms": 5,
            },
            {
                "linestyle": "dashed",
                "color": "red",
                "edgecolor": "black",
                "facecolor": "yellow",
                "linewidth": 2,
                "markeredgecolor": "green",
                "markerfacecoloralt": "blue",
                "markersize": 5,
            },
        ),
    ],
)
def test_validate_style_kwargs(default_kwargs, user_kwargs, expected):
    """Check the behaviour of `validate_style_kwargs` with various type of entries."""
    result = _validate_style_kwargs(default_kwargs, user_kwargs)
    assert result == expected, (
        "The validation of style keywords does not provide the expected results: "
        f"Got {result} instead of {expected}."
    )


@pytest.mark.parametrize(
    "default_kwargs, user_kwargs",
    [({}, {"ls": 2, "linestyle": 3}), ({}, {"c": "r", "color": "blue"})],
)
def test_validate_style_kwargs_error(default_kwargs, user_kwargs):
    """Check that `validate_style_kwargs` raises TypeError"""
    with pytest.raises(TypeError):
        _validate_style_kwargs(default_kwargs, user_kwargs)


def test_despine(pyplot):
    ax = pyplot.gca()
    _despine(ax)
    assert ax.spines["top"].get_visible() is False
    assert ax.spines["right"].get_visible() is False
    assert ax.spines["bottom"].get_bounds() == (0, 1)
    assert ax.spines["left"].get_bounds() == (0, 1)


@pytest.mark.parametrize("estimator_name", ["my_est_name", "deprecated"])
@pytest.mark.parametrize("name", [None, "my_name"])
def test_deprecate_estimator_name(estimator_name, name):
    """Check `_deprecate_estimator_name` behaves correctly"""
    version = "1.7"
    version_remove = "1.9"

    if estimator_name == "deprecated":
        name_out = _deprecate_estimator_name(estimator_name, name, version)
        assert name_out == name
    # `estimator_name` is provided and `name` is:
    elif name is None:
        warning_message = (
            f"`estimator_name` is deprecated in {version} and will be removed in "
            f"{version_remove}. Use `name` instead."
        )
        with pytest.warns(FutureWarning, match=warning_message):
            result = _deprecate_estimator_name(estimator_name, name, version)
        assert result == estimator_name
    elif name is not None:
        error_message = (
            f"Cannot provide both `estimator_name` and `name`. `estimator_name` "
            f"is deprecated in {version} and will be removed in {version_remove}. "
        )
        with pytest.raises(ValueError, match=error_message):
            _deprecate_estimator_name(estimator_name, name, version)
