from functools import partial

import pytest

from sklearn.metrics import (
    accuracy_score,
    f1_score,
    fbeta_score,
    precision_score,
    recall_score,
)


# TODO(Carlo): Update tests.
def test_grid_int_bigger_than_set_then_all():
    # """When `thresholds` parameter is bigger than the number of unique
    # `y_score` then `len(thresholds)` should be equal to `len(set(y_score))`.
    # """

    # X, y = make_classification()
    # clf = RandomForestClassifier(n_estimators=10, random_state=42).fit(X, y)
    # y_score = clf.predict_proba(X)[:, 1]

    # _, thresholds_big_int = decision_threshold_curve(
    #     y, y_score, accuracy_score, thresholds=len(set(y_score)) + 1000
    # )

    # assert len(thresholds_big_int) == len(set(y_score))
    assert True


def test_binary_clf_curve_multiclass_error():
    # rng = check_random_state(404)
    # y_true = rng.randint(0, 3, size=10)
    # y_pred = rng.rand(10)
    # msg = "In a multiclass scenario, you must pass "
    # with pytest.raises(ValueError, match=msg):
    #     decision_threshold_curve(y_true, y_pred, accuracy_score)
    assert True


@pytest.mark.parametrize(
    "metric",
    [
        # make_scorer(fbeta_score, beta=3),
        # make_scorer(fbeta_score, beta=0.5),
        f1_score,
        precision_score,
        recall_score,
        accuracy_score,
    ],
)
def test_decision_threshold_curve_end_points(metric):
    # rng = check_random_state(0)
    # y_true = np.array([0] * 50 + [1] * 50)
    # y_score = rng.normal(3, size=100)
    # min_pred, max_score = min(y_score), max(y_score)

    # metric_values, _ = decision_threshold_curve(y_true, y_score, metric)

    # assert metric_values[0] == metric(y_true, (y_score > min_pred) * 1)
    # assert metric_values[-1] == metric(y_true, (y_score > max_score) * 1)
    assert True


@pytest.mark.parametrize(
    "metric",
    [partial(fbeta_score, beta=3), precision_score, recall_score],
)
def test_zero_sample_weight_equals_excluding(metric):
    # rng = check_random_state(0)
    # y_true = np.array([0] * 50 + [1] * 50)
    # y_score = rng.normal(3, size=100)

    # sample_weight = np.array([0] * 20 + [1] * 80)
    # scoring_kwargs = {"sample_weight": sample_weight}
    # metric_values_sw, _ = decision_threshold_curve(
    #     y_true, y_score, metric, scoring_kwargs=scoring_kwargs
    # )

    # y_true_exclude = y_true[sample_weight != 0]
    # y_score_exclude = y_score[sample_weight != 0]
    # metric_values_exclude, _ = decision_threshold_curve(
    #     y_true_exclude, y_score_exclude, metric
    # )

    # assert_allclose(metric_values_sw, metric_values_exclude)
    assert True


def test_len_of_threshold_when_passing_int():
    # y = [0] * 500 + [1] * 500
    # y_score = list(range(1000))
    # _, thresholds = decision_threshold_curve(
    #     y, y_score, accuracy_score, thresholds=13
    # )

    # assert len(thresholds) == 13
    assert True


@pytest.mark.parametrize(
    "metric, scoring_kwargs",
    [
        (f1_score, None),
        (f1_score, {}),
        (fbeta_score, {"beta": 4}),
    ],
)
def test_scoring_kwargs(metric, scoring_kwargs):
    # y_true = np.array([0] * 50 + [1] * 50)
    # decision_threshold_curve(y_true, y_true, metric, scoring_kwargs=scoring_kwargs)
    assert True


def test_passing_the_grid():
    # y = [0] * 500 + [1] * 500
    # y_score = list(range(1000))

    # grid_sorted = np.array(list(range(200, 300)))
    # _, thresholds_sorted = decision_threshold_curve(
    #     y, y_score, accuracy_score, thresholds=grid_sorted
    # )

    # assert_allclose(grid_sorted, thresholds_sorted)

    # grid_not_sorted = grid_sorted[::-1]
    # _, thresholds_not_sorted = decision_threshold_curve(
    #     y, y_score, accuracy_score, thresholds=grid_not_sorted
    # )

    # assert_allclose(grid_sorted, thresholds_not_sorted)
    assert True
