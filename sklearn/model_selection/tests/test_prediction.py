import numpy as np
import pytest

from sklearn.base import BaseEstimator
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import fbeta_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.utils._testing import assert_array_equal

from sklearn.model_selection import CutoffClassifier


class MockNoPredictorClassifier(BaseEstimator):
    """Classifier which does not predict."""
    def fit(self, X, y):
        self.classes_ = np.array([0, 1])
        return self


@pytest.mark.parametrize(
    "Estimator, params, err_type, err_msg",
    [
        (LogisticRegression, {"method": "xxx"}, ValueError,
         "'method' should be one of"),
        (MockNoPredictorClassifier, {"method": "auto"}, TypeError,
         "'base_estimator' must implement one of the"),
        (SVC, {"method": "predict_proba"}, TypeError,
         "'base_estimator' does not implement predict_proba"),
        (LogisticRegression,
         {"objective_metric": "accuracy", "objective_value": 0.5}, ValueError,
         "When 'objective_metric' is a scoring function"),
        (LogisticRegression,
         {"objective_metric": "cost_matrix",
          "objective_value": np.array([[1], [2]])}, ValueError,
         "When 'objective_metric' is a cost matrix"),
    ]
)
def test_cutoffclassifier_valid_params_error(Estimator, params, err_type,
                                             err_msg):
    # check that the proper errors are raised with wrong parameters
    X, y = make_classification(n_samples=200, n_features=6, random_state=42,
                               n_classes=2)
    with pytest.raises(err_type, match=err_msg):
        clf = CutoffClassifier(base_estimator=Estimator(), **params)
        clf.fit(X, y)


def test_cutoffclassifier_error_pos_label():
    # check that we raise when the classes are not in {0, 1} or {-1, 1}
    X, y = load_breast_cancer(return_X_y=True)
    y += 1
    err_msg = "'y_true' takes value in 1, 2 and 'pos_label' is not specified"
    with pytest.raises(ValueError, match=err_msg):
        CutoffClassifier(
            base_estimator=make_pipeline(StandardScaler(),
                                         LogisticRegression())
        ).fit(X, y)


def test_cutoffclassifier_not_binary():
    # check that we only accept binary target
    X, y = load_iris(return_X_y=True)
    with pytest.raises(ValueError, match="Expected target of binary type."):
        CutoffClassifier(
            base_estimator=make_pipeline(StandardScaler(),
                                         LogisticRegression())
        ).fit(X, y)


def test_cutoffclassifier_limit_tpr_tnr():
    # check that an objective value of 0 give opposite predictions in with
    # tpr and tnr
    X, y = load_breast_cancer(return_X_y=True)
    clf = CutoffClassifier(
        base_estimator=make_pipeline(StandardScaler(), LogisticRegression()),
        objective_metric="tpr",
        objective_value=0,
    )
    y_pred_tpr = clf.fit(X, y).predict(X)
    clf.set_params(objective_metric="tnr")
    y_pred_tnr = (~clf.fit(X, y).predict(X).astype(bool)).astype(int)
    assert_array_equal(y_pred_tnr, y_pred_tpr)


@pytest.mark.parametrize(
    "method",
    ["auto", "decision_function", "predict_proba"]
)
def test_cutoffclassifier_with_objective_value(method):
    # check that we can optimize a given metric as a callable
    X, y = load_breast_cancer(return_X_y=True)
    # remove feature to degrade performances
    X = X[:, :5]

    # make the problem completely imbalanced such that the balanced accuracy
    # is low
    indices_pos = np.flatnonzero(y == 1)
    indices_pos = indices_pos[:indices_pos.size // 50]
    indices_neg = np.flatnonzero(y == 0)

    X = np.vstack([X[indices_neg], X[indices_pos]])
    y = np.hstack([y[indices_neg], y[indices_pos]])

    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model = CutoffClassifier(
        base_estimator=lr,
        objective_metric=balanced_accuracy_score,
        method=method,
    )
    score_optimized = balanced_accuracy_score(y, model.fit(X, y).predict(X))
    score_baseline = balanced_accuracy_score(y, lr.predict(X))
    assert score_optimized > score_baseline
    assert_array_equal(model.classes_, [0, 1])


def test_cutoffclassifier_metric_with_parameter():
    # check that we can pass a metric with a parameter
    # in addition check that f_beta with beta=1 is equivalent to f1
    X, y = load_breast_cancer(return_X_y=True)
    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model_fbeta = CutoffClassifier(
        base_estimator=lr, objective_metric=fbeta_score,
        objective_metric_params={"beta": 1}
    ).fit(X, y)
    model_f1 = CutoffClassifier(
        base_estimator=lr, objective_metric=f1_score,
    ).fit(X, y)

    assert (model_fbeta.decision_threshold_ ==
            pytest.approx(model_f1.decision_threshold_))


def test_cutoffclassifier_pretrained_estimator():
    # check that passing a pre-trained estimator is equivalent to training it
    # in the meta-estimator
    X, y = load_breast_cancer(return_X_y=True)
    lr_prefit = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    lr = make_pipeline(StandardScaler(), LogisticRegression())
    model_prefit = CutoffClassifier(base_estimator=lr_prefit).fit(X, y)
    model = CutoffClassifier(base_estimator=lr).fit(X, y)

    assert (model_prefit.decision_threshold_ ==
            pytest.approx(model.decision_threshold_))

    # check that we did not make any clone/copy of the pretrained estimator
    assert model_prefit._estimator is lr_prefit


@pytest.mark.parametrize(
    "method",
    ["auto", "decision_function", "predict_proba"]
)
@pytest.mark.parametrize("metric", [balanced_accuracy_score, f1_score])
@pytest.mark.parametrize("dtype", [None, object])
def test_cutoffclassifier_with_string_targets(method, dtype, metric):
    # check that targets represented by str are properly managed
    # check with several metrics to be sure that `pos_label` is properly
    # dispatched
    X, y = load_breast_cancer(return_X_y=True)
    # replaces y by some strings
    classes = np.array(["healthy", "cancer"])
    if dtype is not None:
        classes = classes.astype(dtype)
    y = classes[y]
    model = CutoffClassifier(
        base_estimator=make_pipeline(StandardScaler(), LogisticRegression()),
        objective_metric=metric,
        pos_label="cancer",
        method=method,
    ).fit(X, y)
    assert_array_equal(np.sort(model.classes_), np.sort(classes))
    y_pred = model.predict(X[[0], :])
    assert y_pred.item(0) in classes


def test_cutoffclassifier_cost_matrix():
    X, y = load_breast_cancer(return_X_y=True)
    cost_matrix = np.array([[0, 0],
                            [0, 1]])
    clf = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model = CutoffClassifier(
        base_estimator=clf,
        objective_metric="cost_matrix",
        objective_value=cost_matrix,
    )
    model.fit(X, y)
    assert clf.predict(X).sum() < model.predict(X).sum()
    assert model.predict(X).sum() > (y.size * 0.9)
