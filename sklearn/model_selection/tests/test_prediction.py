import numpy as np
import pytest

from sklearn.datasets import load_breast_cancer, make_classification
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, fbeta_score, f1_score, make_scorer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.utils._testing import assert_allclose, assert_array_equal

from sklearn.model_selection import CutOffClassifier


def test_cutoffclassifier_no_binary():
    """Check that we raise an informative error message for non-binary problem."""
    X, y = make_classification(n_classes=3, n_clusters_per_class=1)
    err_msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=err_msg):
        CutOffClassifier(LogisticRegression()).fit(X, y)


@pytest.mark.parametrize(
    "estimator",
    [LogisticRegression(), SVC(), GradientBoostingClassifier(n_estimators=4)],
)
@pytest.mark.parametrize(
    "response_method", ["predict_proba", "predict_log_proba", "decision_function"]
)
def test_cutoffclassifier_estimator_response_methods(estimator, response_method):
    """Check that `CutOffClassifier` exposes the same response methods as the
    underlying estimator.
    """
    X, y = make_classification(n_samples=100, random_state=0)

    model = CutOffClassifier(estimator)
    assert hasattr(model, response_method) == hasattr(estimator, response_method)

    model.fit(X, y)
    assert hasattr(model, response_method) == hasattr(estimator, response_method)

    if hasattr(model, response_method):
        y_pred_cutoff = getattr(model, response_method)(X)
        y_pred_underlying_estimator = getattr(model.estimator_, response_method)(X)

        assert_allclose(y_pred_cutoff, y_pred_underlying_estimator)


@pytest.mark.parametrize(
    "response_method", ["auto", "decision_function", "predict_proba"]
)
def test_cutoffclassifier_with_objective_value(response_method):
    """Check that `CutOffClassifier` is optimizing a given objective metric."""
    X, y = load_breast_cancer(return_X_y=True)
    # remove feature to degrade performances
    X = X[:, :5]

    # make the problem completely imbalanced such that the balanced accuracy is low
    indices_pos = np.flatnonzero(y == 1)
    indices_pos = indices_pos[: indices_pos.size // 50]
    indices_neg = np.flatnonzero(y == 0)

    X = np.vstack([X[indices_neg], X[indices_pos]])
    y = np.hstack([y[indices_neg], y[indices_pos]])

    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model = CutOffClassifier(
        estimator=lr,
        objective_metric="balanced_accuracy",
        response_method=response_method,
    )
    score_optimized = balanced_accuracy_score(y, model.fit(X, y).predict(X))
    score_baseline = balanced_accuracy_score(y, lr.predict(X))
    assert score_optimized > score_baseline


def test_cutoffclassifier_limit_tpr_tnr():
    """Check that an objective value of 0 give opposite predictions with objective
    metrics `tpr` and `tnr`.
    """
    X, y = load_breast_cancer(return_X_y=True)
    estimator = make_pipeline(StandardScaler(), LogisticRegression())
    clf = CutOffClassifier(
        estimator=estimator, objective_metric="tpr", objective_value=0
    )
    y_pred_tpr = clf.fit(X, y).predict(X)
    clf.set_params(objective_metric="tnr")
    y_pred_tnr = (~clf.fit(X, y).predict(X).astype(bool)).astype(int)
    assert np.mean(y_pred_tnr == y_pred_tpr) > 0.98


def test_cutoffclassifier_metric_with_parameter():
    """Check that we can pass a metric with a parameter in addition check that
    `f_beta with beta=1` is equivalent to `f1`.
    """
    X, y = load_breast_cancer(return_X_y=True)
    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model_fbeta = CutOffClassifier(
        estimator=lr, objective_metric=make_scorer(fbeta_score, beta=1)
    ).fit(X, y)
    model_f1 = CutOffClassifier(
        estimator=lr, objective_metric=make_scorer(f1_score)
    ).fit(X, y)

    assert model_fbeta.decision_threshold_ == pytest.approx(
        model_f1.decision_threshold_
    )


@pytest.mark.parametrize(
    "response_method", ["auto", "decision_function", "predict_proba"]
)
@pytest.mark.parametrize(
    "metric",
    [
        "tpr",
        "tnr",
        make_scorer(balanced_accuracy_score),
        make_scorer(f1_score, pos_label="cancer"),
    ],
)
def test_cutoffclassifier_with_string_targets(response_method, metric):
    """Check that targets represented by str are properly managed.
    Also, check with several metrics to be sure that `pos_label` is properly
    dispatched.
    """
    X, y = load_breast_cancer(return_X_y=True)
    # Encode numeric targets by meaningful strings. We purposely designed the class
    # names such that the `pos_label` is the first alphabetically sorted class and thus
    # encoded as 0.
    classes = np.array(["healthy", "cancer"], dtype=object)
    y = classes[y]
    model = CutOffClassifier(
        estimator=make_pipeline(StandardScaler(), LogisticRegression()),
        objective_metric=metric,
        objective_value=0.5,
        pos_label="cancer",
        response_method=response_method,
    ).fit(X, y)
    assert_array_equal(model.classes_, np.sort(classes))
    y_pred = model.predict(X[[0], :])
    assert y_pred.item(0) in classes

    print(model.decision_threshold_)
