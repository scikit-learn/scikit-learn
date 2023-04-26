import numpy as np
import pytest

from sklearn.base import clone
from sklearn.datasets import load_breast_cancer, load_iris, make_classification
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    fbeta_score,
    f1_score,
    make_scorer,
    precision_recall_curve,
    roc_curve,
)
from sklearn.metrics._scorer import _ContinuousScorer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.utils._mocking import CheckingClassifier
from sklearn.utils._testing import (
    _convert_container,
    assert_allclose,
    assert_array_equal,
)

from sklearn.model_selection import CutOffClassifier
from sklearn.model_selection._prediction import _fit_and_score


@pytest.mark.parametrize(
    "scorer, score_method",
    [
        (
            _ContinuousScorer(
                score_func=balanced_accuracy_score,
                sign=1,
                response_method="predict_proba",
                kwargs={},
            ),
            "balanced_accuracy",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tnr_at_tpr_constraint",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tpr_at_tnr_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_precision_at_recall_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_recall_at_precision_constraint",
        ),
    ],
)
def test_fit_and_score_scorers(scorer, score_method):
    """Check that `_fit_and_score` returns thresholds in ascending order for the
    different accepted scorers."""
    X, y = make_classification(n_samples=100, random_state=0)
    train_idx, val_idx = np.arange(50), np.arange(50, 100)
    classifier = LogisticRegression()

    thresholds, scores = _fit_and_score(
        classifier,
        X,
        y,
        sample_weight=None,
        fit_params={},
        train_idx=train_idx,
        val_idx=val_idx,
        scorer=scorer,
        score_method=score_method,
    )

    if score_method.startswith("max_"):
        assert_array_equal(np.argsort(thresholds), np.arange(len(thresholds)))
        assert isinstance(scores, tuple) and len(scores) == 2
        for sc in scores:
            assert np.logical_and(sc >= 0, sc <= 1).all()
    else:
        assert_array_equal(np.argsort(thresholds), np.arange(len(thresholds)))
        assert isinstance(scores, np.ndarray)
        assert np.logical_and(scores >= 0, scores <= 1).all()


@pytest.mark.parametrize(
    "scorer, score_method, expected_score",
    [
        (
            _ContinuousScorer(
                score_func=balanced_accuracy_score,
                sign=1,
                response_method="predict_proba",
                kwargs={},
            ),
            "balanced_accuracy",
            [0.5, 1.0],
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tnr_at_tpr_constraint",
            [[0.0, 1.0, 1.0], [1.0, 1.0, 0.0]],
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tpr_at_tnr_constraint",
            [[0.0, 1.0, 1.0], [1.0, 1.0, 0.0]],
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_precision_at_recall_constraint",
            [[0.5, 1.0], [1.0, 1.0]],
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_recall_at_precision_constraint",
            [[0.5, 1.0], [1.0, 1.0]],
        ),
    ],
)
def test_fit_and_score_prefit(scorer, score_method, expected_score):
    """Check the behaviour with a prefit classifier."""
    X, y = make_classification(n_samples=100, random_state=0)

    # `train_idx is None` to indicate that the classifier is prefit
    train_idx, val_idx = None, np.arange(50, 100)
    classifier = DecisionTreeClassifier(random_state=0)

    with pytest.raises(NotFittedError):
        _fit_and_score(
            classifier,
            X,
            y,
            sample_weight=None,
            fit_params={},
            train_idx=train_idx,
            val_idx=val_idx,
            scorer=scorer,
            score_method=score_method,
        )

    classifier.fit(X, y)
    # make sure that the classifier memorized the full dataset such that
    # we get perfect predictions and thus match the expected score
    assert classifier.score(X[val_idx], y[val_idx]) == pytest.approx(1.0)

    thresholds, scores = _fit_and_score(
        classifier,
        X,
        y,
        sample_weight=None,
        fit_params={},
        train_idx=train_idx,
        val_idx=val_idx,
        scorer=scorer,
        score_method=score_method,
    )
    assert_array_equal(np.argsort(thresholds), np.arange(len(thresholds)))
    assert_allclose(scores, expected_score)


@pytest.mark.parametrize(
    "scorer, score_method",
    [
        (
            _ContinuousScorer(
                score_func=balanced_accuracy_score,
                sign=1,
                response_method="predict_proba",
                kwargs={},
            ),
            "balanced_accuracy",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tnr_at_tpr_constraint",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tpr_at_tnr_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_precision_at_recall_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_recall_at_precision_constraint",
        ),
    ],
)
def test_fit_and_score_sample_weight(scorer, score_method):
    """Check that we dispatch the sample-weight to fit and score the classifier."""
    X, y = load_iris(return_X_y=True)
    X, y = X[:100], y[:100]  # only 2 classes

    # create a dataset and repeat twice the sample of class #0
    X_repeated, y_repeated = np.vstack([X, X[y == 0]]), np.hstack([y, y[y == 0]])
    # create a sample weight vector that is equivalent to the repeated dataset
    sample_weight = np.ones_like(y)
    sample_weight[:50] *= 2

    classifier = LogisticRegression()
    train_repeated_idx = np.arange(X_repeated.shape[0])
    val_repeated_idx = np.arange(X_repeated.shape[0])
    thresholds_repeated, scores_repeated = _fit_and_score(
        classifier,
        X_repeated,
        y_repeated,
        sample_weight=None,
        fit_params={},
        train_idx=train_repeated_idx,
        val_idx=val_repeated_idx,
        scorer=scorer,
        score_method=score_method,
    )

    train_idx, val_idx = np.arange(X.shape[0]), np.arange(X.shape[0])
    thresholds, scores = _fit_and_score(
        classifier,
        X,
        y,
        sample_weight=sample_weight,
        fit_params={},
        train_idx=train_idx,
        val_idx=val_idx,
        scorer=scorer,
        score_method=score_method,
    )

    assert_allclose(thresholds_repeated, thresholds)
    assert_allclose(scores_repeated, scores)


@pytest.mark.parametrize(
    "scorer, score_method",
    [
        (
            _ContinuousScorer(
                score_func=balanced_accuracy_score,
                sign=1,
                response_method="predict_proba",
                kwargs={},
            ),
            "balanced_accuracy",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tnr_at_tpr_constraint",
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tpr_at_tnr_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_precision_at_recall_constraint",
        ),
        (
            make_scorer(precision_recall_curve, needs_proba=True),
            "max_recall_at_precision_constraint",
        ),
    ],
)
@pytest.mark.parametrize("fit_params_type", ["list", "array"])
def test_fit_and_score_fit_params(scorer, score_method, fit_params_type):
    """Check that we pass `fit_params` to the classifier when calling `fit`."""
    X, y = make_classification(n_samples=100, random_state=0)
    fit_params = {
        "a": _convert_container(y, fit_params_type),
        "b": _convert_container(y, fit_params_type),
    }

    classifier = CheckingClassifier(expected_fit_params=["a", "b"])
    train_idx, val_idx = np.arange(50), np.arange(50, 100)

    _fit_and_score(
        classifier,
        X,
        y,
        sample_weight=None,
        fit_params=fit_params,
        train_idx=train_idx,
        val_idx=val_idx,
        scorer=scorer,
        score_method=score_method,
    )


def test_cutoffclassifier_no_binary():
    """Check that we raise an informative error message for non-binary problem."""
    X, y = make_classification(n_classes=3, n_clusters_per_class=1)
    err_msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=err_msg):
        CutOffClassifier(LogisticRegression()).fit(X, y)


@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        (
            {"cv": "prefit", "refit": True},
            ValueError,
            "When cv='prefit', refit cannot be True.",
        ),
        (
            {"cv": 10, "refit": False},
            ValueError,
            "When cv has several folds, refit cannot be False.",
        ),
        (
            {"cv": "prefit", "refit": False},
            NotFittedError,
            "`estimator` must be fitted.",
        ),
    ],
)
def test_cutoffclassifier_conflict_cv_refit(params, err_type, err_msg):
    """Check that we raise an informative error message when `cv` and `refit`
    cannot be used together.
    """
    X, y = make_classification(n_samples=100, random_state=0)
    with pytest.raises(err_type, match=err_msg):
        CutOffClassifier(LogisticRegression(), **params).fit(X, y)


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
def test_cutoffclassifier_with_constraint_value(response_method):
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


@pytest.mark.parametrize(
    "metrics",
    [
        ("max_tpr_at_tnr_constraint", "max_tnr_at_tpr_constraint"),
        ("max_tnr_at_tpr_constraint", "max_tpr_at_tnr_constraint"),
    ],
)
def test_cutoffclassifier_limit_metric_tradeoff(metrics):
    """Check that an objective value of 0 give opposite predictions with tnr/tpr and
    precision/recall.
    """
    X, y = load_breast_cancer(return_X_y=True)
    estimator = make_pipeline(StandardScaler(), LogisticRegression())
    model = CutOffClassifier(
        estimator=estimator,
        objective_metric=metrics[0],
        constraint_value=0,
    )
    y_pred_1 = model.fit(X, y).predict(X)
    model.set_params(objective_metric=metrics[1])
    y_pred_2 = (~model.fit(X, y).predict(X).astype(bool)).astype(int)
    assert np.mean(y_pred_1 == y_pred_2) == pytest.approx(1.0)


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
        "max_tnr_at_tpr_constraint",
        "max_tpr_at_tnr_constraint",
        "max_precision_at_recall_constraint",
        "max_recall_at_precision_constraint",
        make_scorer(balanced_accuracy_score),
        make_scorer(f1_score, pos_label="cancer"),
        {"tp": 5, "tn": 1, "fp": -1, "fn": -1},
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
    classes = np.array(["cancer", "healthy"], dtype=object)
    y = classes[y]
    model = CutOffClassifier(
        estimator=make_pipeline(StandardScaler(), LogisticRegression()),
        objective_metric=metric,
        constraint_value=0.9,
        pos_label="cancer",
        response_method=response_method,
        n_thresholds=100,
    ).fit(X, y)
    assert_array_equal(model.classes_, np.sort(classes))
    y_pred = model.predict(X)
    assert_array_equal(np.sort(np.unique(y_pred)), np.sort(classes))


@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_cutoffclassifier_refit(with_sample_weight, global_random_seed):
    """Check the behaviour of the `refit` parameter."""
    rng = np.random.RandomState(global_random_seed)
    X, y = make_classification(n_samples=100, random_state=0)
    if with_sample_weight:
        sample_weight = rng.randn(X.shape[0])
        sample_weight = np.abs(sample_weight, out=sample_weight)
    else:
        sample_weight = None

    # check that `estimator_` if fitted on the full dataset when `refit=True`
    estimator = LogisticRegression()
    model = CutOffClassifier(estimator, refit=True).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is not estimator
    estimator.fit(X, y, sample_weight=sample_weight)
    assert_allclose(model.estimator_.coef_, estimator.coef_)
    assert_allclose(model.estimator_.intercept_, estimator.intercept_)

    # check that `estimator_` was not altered when `refit=False` and `cv="prefit"`
    estimator = LogisticRegression().fit(X, y, sample_weight=sample_weight)
    coef = estimator.coef_.copy()
    model = CutOffClassifier(estimator, cv="prefit", refit=False).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is estimator
    assert_allclose(model.estimator_.coef_, coef)

    # check that we train `estimator_` on the training split of a given cross-validation
    estimator = LogisticRegression()
    cv = [
        (np.arange(50), np.arange(50, 100)),
    ]  # single split
    model = CutOffClassifier(estimator, cv=cv, refit=False).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is not estimator
    if with_sample_weight:
        sw_train = sample_weight[cv[0][0]]
    else:
        sw_train = None
    estimator.fit(X[cv[0][0]], y[cv[0][0]], sample_weight=sw_train)
    assert_allclose(model.estimator_.coef_, estimator.coef_)


@pytest.mark.parametrize(
    "objective_metric",
    [
        "max_tnr_at_tpr_constraint",
        "max_tpr_at_tnr_constraint",
        "max_precision_at_recall_constraint",
        "max_recall_at_precision_constraint",
        "balanced_accuracy",
    ],
)
@pytest.mark.parametrize("fit_params_type", ["list", "array"])
def test_cutoffclassifier_fit_params(objective_metric, fit_params_type):
    """Check that we pass `fit_params` to the classifier when calling `fit`."""
    X, y = make_classification(n_samples=100, random_state=0)
    fit_params = {
        "a": _convert_container(y, fit_params_type),
        "b": _convert_container(y, fit_params_type),
    }

    classifier = CheckingClassifier(expected_fit_params=["a", "b"])
    model = CutOffClassifier(
        classifier, objective_metric=objective_metric, constraint_value=0.5
    )
    model.fit(X, y, **fit_params)


@pytest.mark.parametrize(
    "objective_metric, constraint_value",
    [
        ("max_tnr_at_tpr_constraint", 0.5),
        ("max_tpr_at_tnr_constraint", 0.5),
        ("max_precision_at_recall_constraint", 0.5),
        ("max_recall_at_precision_constraint", 0.5),
    ],
)
@pytest.mark.parametrize(
    "response_method", ["auto", "decision_function", "predict_proba"]
)
def test_cutoffclassifier_response_method_scorer_with_constraint_metric(
    objective_metric, constraint_value, response_method, global_random_seed
):
    """Check that we use the proper scorer and forwarding the requested response method
    for TNR/TPR and precision/recall metrics.
    """
    X, y = make_classification(n_samples=100, random_state=global_random_seed)
    classifier = LogisticRegression()

    model = CutOffClassifier(
        classifier,
        objective_metric=objective_metric,
        constraint_value=constraint_value,
        response_method=response_method,
    )
    model.fit(X, y)

    if response_method in ("auto", "predict_proba"):
        # "auto" will fall back  in priority on `predict_proba` if `estimator`
        # supports it.
        # we expect the decision threshold to be in [0, 1]
        if objective_metric in (
            "max_tnr_at_tpr_constraint",
            "max_precision_at_recall_constraint",
        ):
            assert 0.5 <= model.decision_threshold_ <= 1
        else:  # "max_tpr_at_tnr_constraint" or "max_recall_at_precision_constraint"
            assert 0 <= model.decision_threshold_ <= 0.5
    else:  # "decision_function"
        # we expect the decision function to be centered in 0.0 and to be larger than
        # -1 and 1.
        if objective_metric in (
            "max_tnr_at_tpr_constraint",
            "max_precision_at_recall_constraint",
        ):
            assert 0 < model.decision_threshold_ < 20
        else:  # "max_tpr_at_tnr_constraint" or "max_recall_at_precision_constraint"
            assert -20 < model.decision_threshold_ < 0


def test_cutoffclassifier_objective_metric_dict(global_random_seed):
    """Check that we can pass a custom objective metric."""
    X, y = make_classification(n_samples=500, random_state=global_random_seed)
    classifier = LogisticRegression()

    # we need to set a small number of thresholds to avoid ties and picking a too low
    # threshold.
    n_thresholds = 5

    # affect a high gain to true negative and force the classifier to mainly
    # predict the negative class.
    costs_and_again = {"tp": 0, "tn": 10, "fp": 0, "fn": 0}
    model = CutOffClassifier(
        classifier, objective_metric=costs_and_again, n_thresholds=n_thresholds
    )
    model.fit(X, y)

    assert model.decision_threshold_ > 0.99
    assert np.mean(model.predict(X) == 0) > 0.9

    # use the true positive now
    costs_and_again = {"tp": 10, "tn": 0, "fp": 0, "fn": 0}
    model = CutOffClassifier(
        classifier, objective_metric=costs_and_again, n_thresholds=n_thresholds
    )
    model.fit(X, y)

    assert model.decision_threshold_ < 0.01
    assert np.mean(model.predict(X) == 1) > 0.9

    # flipping the `pos_label` to zero should force the classifier to always predict 0
    # and thus have a low threshold
    pos_label = 0
    model = CutOffClassifier(
        classifier,
        objective_metric=costs_and_again,
        n_thresholds=n_thresholds,
        pos_label=pos_label,
    )
    model.fit(X, y)

    assert model.decision_threshold_ < 0.01
    assert np.mean(model.predict(X) == 0) > 0.9


def test_cutoffclassifier_sample_weight_costs_and_again():
    """Check that we dispatch the `sample_weight` to the scorer when computing the
    confusion matrix."""
    X, y = load_iris(return_X_y=True)
    X, y = X[:100], y[:100]  # only 2 classes

    # create a dataset and repeat twice the sample of class #0
    X_repeated, y_repeated = np.vstack([X, X[y == 0]]), np.hstack([y, y[y == 0]])
    # create a sample weight vector that is equivalent to the repeated dataset
    sample_weight = np.ones_like(y)
    sample_weight[:50] *= 2

    # we use a prefit classifier to simplify the test
    cv = "prefit"
    estimator = LogisticRegression().fit(X, y)
    costs_and_again = {"tp": 1, "tn": 1, "fp": -1, "fn": -1}

    model_repeat = CutOffClassifier(estimator, cv=cv, objective_metric=costs_and_again)
    model_repeat.fit(X_repeated, y_repeated, sample_weight=None)

    model_sw = CutOffClassifier(estimator, cv=cv, objective_metric=costs_and_again)
    model_sw.fit(X, y, sample_weight=sample_weight)

    assert model_repeat.objective_score_ == pytest.approx(model_sw.objective_score_)


def test_cutoffclassifier_cv_zeros_sample_weights_equivalence():
    """Check that passing removing some sample from the dataset `X` is
    equivalent to passing a `sample_weight` with a factor 0."""
    X, y = load_iris(return_X_y=True)
    # Scale the data to avoid any convergence issue
    X = StandardScaler().fit_transform(X)
    # Only use 2 classes and select samples such that 2-fold cross-validation
    # split will lead to an equivalence with a `sample_weight` of 0
    X = np.vstack((X[:40], X[50:90]))
    y = np.hstack((y[:40], y[50:90]))
    sample_weight = np.zeros_like(y)
    sample_weight[::2] = 1

    estimator = LogisticRegression()
    model_without_weights = CutOffClassifier(estimator, cv=2)
    model_with_weights = clone(model_without_weights)

    model_with_weights.fit(X, y, sample_weight=sample_weight)
    model_without_weights.fit(X[::2], y[::2])

    assert_allclose(
        model_with_weights.estimator_.coef_, model_without_weights.estimator_.coef_
    )

    y_pred_with_weights = model_with_weights.predict_proba(X)
    y_pred_without_weights = model_without_weights.predict_proba(X)
    assert_allclose(y_pred_with_weights, y_pred_without_weights)


def test_cutoffclassifier_error_constant_learner():
    """Check that we raise an error message when providing an estimator that predicts
    only a single class."""
    X, y = make_classification(random_state=0)
    estimator = DummyClassifier(strategy="constant", constant=1)
    err_msg = "The provided estimator makes constant predictions."
    with pytest.raises(ValueError, match=err_msg):
        CutOffClassifier(estimator).fit(X, y)


# TODO write non-regression test when pos_label corresponds to idx #0
# Before we did not think to potentially remap the the output of the comparison
# to the original labels.
