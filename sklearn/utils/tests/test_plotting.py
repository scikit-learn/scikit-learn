import numpy as np
import pytest

from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _despine,
    _interval_max_min_ratio,
    _validate_score_name,
    _validate_style_kwargs,
)


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
    print(label)
    assert label == expected_label


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
