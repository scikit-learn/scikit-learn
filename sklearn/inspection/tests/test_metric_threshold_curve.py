import pytest
import numpy as np

from functools import partial

from numpy.testing import assert_allclose

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import metric_threshold_curve
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    fbeta_score,
    precision_score,
    recall_score,
)
from sklearn.utils.validation import check_random_state


def test_grid_int_bigger_than_set_then_all():
    """When `threshold_grid` parameter is bigger than the number of unique
    `y_score` then `len(thresholds)` should be equal to `len(set(y_score))`
    and thresholds should be the same from what we get with
    `threshold_grid=None`.
    """

    X, y = make_classification()
    clf = RandomForestClassifier(n_estimators=10, random_state=42).fit(X, y)
    y_score = clf.predict_proba(X)[:, 1]

    _, thresholds_big_int = metric_threshold_curve(
        y, y_score, accuracy_score, threshold_grid=len(set(y_score)) + 1000
    )

    _, thresholds_none = metric_threshold_curve(
        y, y_score, accuracy_score, threshold_grid=None
    )

    assert_allclose(thresholds_big_int, thresholds_none)
    assert len(thresholds_big_int) == len(set(y_score))


def test_binary_clf_curve_multiclass_error():
    rng = check_random_state(404)
    y_true = rng.randint(0, 3, size=10)
    y_pred = rng.rand(10)
    msg = "multiclass format is not supported"
    with pytest.raises(ValueError, match=msg):
        metric_threshold_curve(y_true, y_pred, accuracy_score)


@pytest.mark.parametrize(
    "metric",
    [
        partial(fbeta_score, beta=3),
        partial(fbeta_score, beta=0.5),
        f1_score,
        precision_score,
        recall_score,
        accuracy_score,
    ],
)
def test_metric_threshold_curve_end_points(metric):
    rng = check_random_state(0)
    y_true = np.array([0] * 50 + [1] * 50)
    y_score = rng.normal(3, size=100)
    min_pred, max_score = min(y_score), max(y_score)

    metric_values, _ = metric_threshold_curve(y_true, y_score, metric)

    assert metric_values[0] == metric(y_true, (y_score > min_pred) * 1)
    assert metric_values[-1] == metric(y_true, (y_score > max_score) * 1)


@pytest.mark.parametrize(
    "metric",
    [partial(fbeta_score, beta=3), precision_score, recall_score],
)
def test_zero_sample_weight_equals_excluding(metric):
    rng = check_random_state(0)
    y_true = np.array([0] * 50 + [1] * 50)
    y_score = rng.normal(3, size=100)

    sample_weight = np.array([0] * 20 + [1] * 80)
    metric_values_sw, _ = metric_threshold_curve(
        y_true, y_score, metric, sample_weight=sample_weight
    )

    y_true_exclude = y_true[sample_weight != 0]
    y_score_exclude = y_score[sample_weight != 0]
    metric_values_exclude, _ = metric_threshold_curve(
        y_true_exclude, y_score_exclude, metric
    )

    assert_allclose(metric_values_sw, metric_values_exclude)


def test_len_of_threshold_when_passing_int():
    y = [0] * 500 + [1] * 500
    y_score = list(range(1000))
    _, thresholds = metric_threshold_curve(
        y, y_score, accuracy_score, threshold_grid=13
    )

    assert len(thresholds) == 13


def test_passing_the_grid():
    y = [0] * 500 + [1] * 500
    y_score = list(range(1000))

    grid_sorted = np.array(list(range(200, 300)))
    _, thresholds_sorted = metric_threshold_curve(
        y, y_score, accuracy_score, threshold_grid=grid_sorted
    )

    assert_allclose(grid_sorted, thresholds_sorted)

    grid_not_sorted = grid_sorted[::-1]
    _, thresholds_not_sorted = metric_threshold_curve(
        y, y_score, accuracy_score, threshold_grid=grid_not_sorted
    )

    assert_allclose(grid_sorted, thresholds_not_sorted)
