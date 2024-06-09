import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    recall_score,
)
from sklearn.metrics._classification_threshold import _CurveScorer


def test_curve_scorer():
    """Check the behaviour of the `_CurveScorer` class."""
    X, y = make_classification(random_state=0)
    estimator = LogisticRegression().fit(X, y)
    curve_scorer = _CurveScorer(
        balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={},
    )
    scores, thresholds = curve_scorer(estimator, X, y)

    assert thresholds.shape == scores.shape
    # check that the thresholds are probabilities with extreme values close to 0 and 1.
    # they are not exactly 0 and 1 because they are the extremum of the
    # `estimator.predict_proba(X)` values.
    assert 0 <= thresholds.min() <= 0.01
    assert 0.99 <= thresholds.max() <= 1
    # balanced accuracy should be between 0.5 and 1 when it is not adjusted
    assert 0.5 <= scores.min() <= 1

    # check that passing kwargs to the scorer works
    curve_scorer = _CurveScorer(
        balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={"adjusted": True},
    )
    scores, thresholds = curve_scorer(estimator, X, y)

    # balanced accuracy should be between 0.5 and 1 when it is not adjusted
    assert 0 <= scores.min() <= 0.5

    # check that we can inverse the sign of the score when dealing with `neg_*` scorer
    curve_scorer = _CurveScorer(
        balanced_accuracy_score,
        sign=-1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={"adjusted": True},
    )
    scores, thresholds = curve_scorer(estimator, X, y)

    assert all(scores <= 0)


def test_curve_scorer_pos_label(global_random_seed):
    """Check that we propagate properly the `pos_label` parameter to the scorer."""
    n_samples = 30
    X, y = make_classification(
        n_samples=n_samples, weights=[0.9, 0.1], random_state=global_random_seed
    )
    estimator = LogisticRegression().fit(X, y)

    curve_scorer = _CurveScorer(
        recall_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={"pos_label": 1},
    )
    scores_pos_label_1, thresholds_pos_label_1 = curve_scorer(estimator, X, y)

    curve_scorer = _CurveScorer(
        recall_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={"pos_label": 0},
    )
    scores_pos_label_0, thresholds_pos_label_0 = curve_scorer(estimator, X, y)

    # Since `pos_label` is forwarded to the curve_scorer, the thresholds are not equal.
    assert not (thresholds_pos_label_1 == thresholds_pos_label_0).all()
    # The min-max range for the thresholds is defined by the probabilities of the
    # `pos_label` class (the column of `predict_proba`).
    y_pred = estimator.predict_proba(X)
    assert thresholds_pos_label_0.min() == pytest.approx(y_pred.min(axis=0)[0])
    assert thresholds_pos_label_0.max() == pytest.approx(y_pred.max(axis=0)[0])
    assert thresholds_pos_label_1.min() == pytest.approx(y_pred.min(axis=0)[1])
    assert thresholds_pos_label_1.max() == pytest.approx(y_pred.max(axis=0)[1])

    # The recall cannot be negative and `pos_label=1` should have a higher recall
    # since there is less samples to be considered.
    assert 0.0 < scores_pos_label_0.min() < scores_pos_label_1.min()
    assert scores_pos_label_0.max() == pytest.approx(1.0)
    assert scores_pos_label_1.max() == pytest.approx(1.0)
