import numpy as np
import pytest

from sklearn.utils._plotting import _validate_score_name, _compute_scale_type_ratio


def my_metric():
    pass  # pragma: no cover


@pytest.mark.parametrize(
    "score_name, scoring, negate_score, expected_score_name",
    [
        ("accuracy", None, False, "accuracy"),  # do not transform the name
        (None, "accuracy", False, "Accuracy"),  # capitalize the name
        (None, "neg_mean_absolute_error", True, "Mean absolute error"),  # remove "neg_"
        ("MAE", "neg_mean_absolute_error", True, "MAE"),  # keep score_name
        (None, None, False, "Score"),  # default name
        (None, None, True, "Negative score"),  # default name but negated
        ("my metric", my_metric, False, "my metric"),  # do not transform the name
        (None, my_metric, False, "My metric"),  # default name
        (None, my_metric, True, "My metric"),  # default name but negated
    ],
)
def test_validate_score_name(score_name, scoring, negate_score, expected_score_name):
    """Check that we return the right score name."""
    assert (
        _validate_score_name(score_name, scoring, negate_score) == expected_score_name
    )


@pytest.mark.parametrize(
    "data, lower_bound, upper_bound",
    [
        # Such data generation could either used with log scale or uniform scale.
        (np.logspace(-1, 0, 5), 5, 6),
        # A uniform generation will be close to 1.
        (np.linspace(0, 1, 5), 0.9, 1.1),
        # This is not strictly log generated data but we will benefit from using a
        # log scale.
        ([1, 2, 5, 10, 20, 50], 20, 40),
    ],
)
def test_compute_scale_type_ratio(data, lower_bound, upper_bound):
    print(_compute_scale_type_ratio(data))
    assert lower_bound < _compute_scale_type_ratio(data) < upper_bound
