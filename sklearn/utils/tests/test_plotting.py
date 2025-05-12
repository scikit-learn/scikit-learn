import numpy as np
import pytest

from sklearn.linear_model import LogisticRegression
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
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
