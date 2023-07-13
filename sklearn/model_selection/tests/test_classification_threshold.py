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
    confusion_matrix,
    f1_score,
    fbeta_score,
    make_scorer,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_curve,
)
from sklearn.metrics._scorer import _ContinuousScorer
from sklearn.model_selection import TunedThresholdClassifier
from sklearn.model_selection._classification_threshold import _fit_and_score
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
            [[0.0, 1.0], [1.0, 1.0]],
        ),
        (
            make_scorer(roc_curve, needs_proba=True),
            "max_tpr_at_tnr_constraint",
            [[0.0, 1.0], [1.0, 1.0]],
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

    classifier = CheckingClassifier(expected_fit_params=["a", "b"], random_state=0)
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


def test_tunedthresholdclassifier_no_binary():
    """Check that we raise an informative error message for non-binary problem."""
    X, y = make_classification(n_classes=3, n_clusters_per_class=1)
    err_msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=err_msg):
        TunedThresholdClassifier(LogisticRegression()).fit(X, y)


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
@pytest.mark.parametrize("strategy", ["optimum", "constant"])
def test_tunedthresholdclassifier_conflict_cv_refit(
    strategy, params, err_type, err_msg
):
    """Check that we raise an informative error message when `cv` and `refit`
    cannot be used together.
    """
    X, y = make_classification(n_samples=100, random_state=0)
    with pytest.raises(err_type, match=err_msg):
        TunedThresholdClassifier(LogisticRegression(), strategy=strategy, **params).fit(
            X, y
        )


@pytest.mark.parametrize(
    "estimator",
    [LogisticRegression(), SVC(), GradientBoostingClassifier(n_estimators=4)],
)
@pytest.mark.parametrize(
    "response_method", ["predict_proba", "predict_log_proba", "decision_function"]
)
@pytest.mark.parametrize("strategy", ["optimum", "constant"])
def test_tunedthresholdclassifier_estimator_response_methods(
    estimator, strategy, response_method
):
    """Check that `TunedThresholdClassifier` exposes the same response methods as the
    underlying estimator.
    """
    X, y = make_classification(n_samples=100, random_state=0)

    model = TunedThresholdClassifier(estimator, strategy=strategy)
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
def test_tunedthresholdclassifier_with_constraint_value(response_method):
    """Check that `TunedThresholdClassifier` is optimizing a given objective metric."""
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
    n_thresholds = 100
    model = TunedThresholdClassifier(
        estimator=lr,
        objective_metric="balanced_accuracy",
        response_method=response_method,
        n_thresholds=n_thresholds,
    )
    score_optimized = balanced_accuracy_score(y, model.fit(X, y).predict(X))
    score_baseline = balanced_accuracy_score(y, lr.predict(X))
    assert score_optimized > score_baseline
    assert model.decision_thresholds_.shape == (n_thresholds,)
    assert model.objective_scores_.shape == (n_thresholds,)


@pytest.mark.parametrize(
    "metrics",
    [
        ("max_tpr_at_tnr_constraint", "max_tnr_at_tpr_constraint"),
        ("max_tnr_at_tpr_constraint", "max_tpr_at_tnr_constraint"),
    ],
)
def test_tunedthresholdclassifier_limit_metric_tradeoff(metrics):
    """Check that an objective value of 0 give opposite predictions with tnr/tpr and
    precision/recall.
    """
    X, y = load_breast_cancer(return_X_y=True)
    estimator = make_pipeline(StandardScaler(), LogisticRegression())
    model = TunedThresholdClassifier(
        estimator=estimator,
        objective_metric=metrics[0],
        constraint_value=0,
    )
    y_pred_1 = model.fit(X, y).predict(X)
    model.set_params(objective_metric=metrics[1])
    y_pred_2 = (~model.fit(X, y).predict(X).astype(bool)).astype(int)
    assert np.mean(y_pred_1 == y_pred_2) > 0.98


def test_tunedthresholdclassifier_metric_with_parameter():
    """Check that we can pass a metric with a parameter in addition check that
    `f_beta with beta=1` is equivalent to `f1`.
    """
    X, y = load_breast_cancer(return_X_y=True)
    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model_fbeta = TunedThresholdClassifier(
        estimator=lr, objective_metric=make_scorer(fbeta_score, beta=1)
    ).fit(X, y)
    model_f1 = TunedThresholdClassifier(
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
def test_tunedthresholdclassifier_with_string_targets(response_method, metric):
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
    model = TunedThresholdClassifier(
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


@pytest.mark.parametrize("strategy", ["optimum", "constant"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_tunedthresholdclassifier_refit(
    strategy, with_sample_weight, global_random_seed
):
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
    model = TunedThresholdClassifier(estimator, strategy=strategy, refit=True).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is not estimator
    estimator.fit(X, y, sample_weight=sample_weight)
    assert_allclose(model.estimator_.coef_, estimator.coef_)
    assert_allclose(model.estimator_.intercept_, estimator.intercept_)

    # check that `estimator_` was not altered when `refit=False` and `cv="prefit"`
    estimator = LogisticRegression().fit(X, y, sample_weight=sample_weight)
    coef = estimator.coef_.copy()
    model = TunedThresholdClassifier(
        estimator, strategy=strategy, cv="prefit", refit=False
    ).fit(X, y, sample_weight=sample_weight)

    assert model.estimator_ is estimator
    assert_allclose(model.estimator_.coef_, coef)

    # check that we train `estimator_` on the training split of a given cross-validation
    estimator = LogisticRegression()
    cv = [
        (np.arange(50), np.arange(50, 100)),
    ]  # single split
    model = TunedThresholdClassifier(
        estimator, strategy=strategy, cv=cv, refit=False
    ).fit(X, y, sample_weight=sample_weight)

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
def test_tunedthresholdclassifier_fit_params(objective_metric, fit_params_type):
    """Check that we pass `fit_params` to the classifier when calling `fit`."""
    X, y = make_classification(n_samples=100, random_state=0)
    fit_params = {
        "a": _convert_container(y, fit_params_type),
        "b": _convert_container(y, fit_params_type),
    }

    classifier = CheckingClassifier(expected_fit_params=["a", "b"], random_state=0)
    model = TunedThresholdClassifier(
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
def test_tunedthresholdclassifier_response_method_scorer_with_constraint_metric(
    objective_metric, constraint_value, response_method, global_random_seed
):
    """Check that we use the proper scorer and forwarding the requested response method
    for TNR/TPR and precision/recall metrics.
    """
    X, y = make_classification(n_samples=100, random_state=global_random_seed)
    classifier = LogisticRegression()

    n_thresholds = 100
    model = TunedThresholdClassifier(
        classifier,
        objective_metric=objective_metric,
        constraint_value=constraint_value,
        response_method=response_method,
        n_thresholds=n_thresholds,
    )
    model.fit(X, y)
    assert model.decision_thresholds_.shape == (n_thresholds,)
    assert all(score.shape == (n_thresholds,) for score in model.objective_scores_)

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


def test_tunedthresholdclassifier_objective_metric_dict(global_random_seed):
    """Check that we can pass a custom objective metric."""
    X, y = make_classification(n_samples=500, random_state=global_random_seed)
    classifier = LogisticRegression()

    # we need to set a small number of thresholds to avoid ties and picking a too low
    # threshold.
    n_thresholds = 5

    # affect a high gain to true negative and force the classifier to mainly
    # predict the negative class.
    costs_and_again = {"tp": 0, "tn": 10, "fp": 0, "fn": 0}
    model = TunedThresholdClassifier(
        classifier, objective_metric=costs_and_again, n_thresholds=n_thresholds
    )
    model.fit(X, y)

    assert model.decision_thresholds_.shape == (n_thresholds,)
    assert model.objective_scores_.shape == (n_thresholds,)

    assert model.decision_threshold_ > 0.99
    assert np.mean(model.predict(X) == 0) > 0.9

    # use the true positive now
    costs_and_again = {"tp": 10, "tn": 0, "fp": 0, "fn": 0}
    model = TunedThresholdClassifier(
        classifier, objective_metric=costs_and_again, n_thresholds=n_thresholds
    )
    model.fit(X, y)

    assert model.decision_thresholds_.shape == (n_thresholds,)
    assert model.objective_scores_.shape == (n_thresholds,)

    assert model.decision_threshold_ < 0.01
    assert np.mean(model.predict(X) == 1) > 0.9

    # flipping the `pos_label` to zero should force the classifier to always predict 0
    # and thus have a low threshold
    pos_label = 0
    model = TunedThresholdClassifier(
        classifier,
        objective_metric=costs_and_again,
        n_thresholds=n_thresholds,
        pos_label=pos_label,
    )
    model.fit(X, y)

    assert model.decision_thresholds_.shape == (n_thresholds,)
    assert model.objective_scores_.shape == (n_thresholds,)

    assert model.decision_threshold_ < 0.01
    assert np.mean(model.predict(X) == 0) > 0.9


def test_tunedthresholdclassifier_sample_weight_costs_and_gain():
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

    model_repeat = TunedThresholdClassifier(
        estimator, cv=cv, objective_metric=costs_and_again
    )
    model_repeat.fit(X_repeated, y_repeated, sample_weight=None)

    model_sw = TunedThresholdClassifier(
        estimator, cv=cv, objective_metric=costs_and_again
    )
    model_sw.fit(X, y, sample_weight=sample_weight)

    assert model_repeat.objective_score_ == pytest.approx(model_sw.objective_score_)


def test_tunedthresholdclassifier_cv_zeros_sample_weights_equivalence():
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
    model_without_weights = TunedThresholdClassifier(estimator, cv=2)
    model_with_weights = clone(model_without_weights)

    model_with_weights.fit(X, y, sample_weight=sample_weight)
    model_without_weights.fit(X[::2], y[::2])

    assert_allclose(
        model_with_weights.estimator_.coef_, model_without_weights.estimator_.coef_
    )

    y_pred_with_weights = model_with_weights.predict_proba(X)
    y_pred_without_weights = model_without_weights.predict_proba(X)
    assert_allclose(y_pred_with_weights, y_pred_without_weights)


def test_tunedthresholdclassifier_error_constant_learner():
    """Check that we raise an error message when providing an estimator that predicts
    only a single class."""
    X, y = make_classification(random_state=0)
    estimator = DummyClassifier(strategy="constant", constant=1)
    err_msg = "The provided estimator makes constant predictions."
    with pytest.raises(ValueError, match=err_msg):
        TunedThresholdClassifier(estimator).fit(X, y)


@pytest.mark.parametrize(
    "objective_metric",
    ["max_precision_at_recall_constraint", "max_recall_at_precision_constraint"],
)
@pytest.mark.parametrize("pos_label", [0, 1])
def test_tunedthresholdclassifier_pos_label_precision_recall(
    objective_metric, pos_label
):
    """Check that `pos_label` is dispatched correctly by checking the precision and
    recall score found during the optimization and the one found at `predict` time."""
    X, y = make_classification(n_samples=5_000, weights=[0.6, 0.4], random_state=42)

    # prefit the estimator to avoid variability due to the cross-validation
    estimator = LogisticRegression().fit(X, y)

    constraint_value = 0.7
    model = TunedThresholdClassifier(
        estimator,
        objective_metric=objective_metric,
        constraint_value=constraint_value,
        cv="prefit",
        pos_label=pos_label,
    ).fit(X, y)

    precision = precision_score(y, model.predict(X), pos_label=pos_label)
    recall = recall_score(y, model.predict(X), pos_label=pos_label)

    # due to internal interpolation, the scores will vary slightly
    if objective_metric == "max_precision_at_recall_constraint":
        assert recall == pytest.approx(model.objective_score_[0], abs=1e-3)
        assert precision == pytest.approx(model.objective_score_[1], abs=1e-3)
    else:
        assert precision == pytest.approx(model.objective_score_[0], abs=1e-3)
        assert recall == pytest.approx(model.objective_score_[1], abs=1e-3)


@pytest.mark.parametrize(
    "objective_metric", ["max_tnr_at_tpr_constraint", "max_tpr_at_tnr_constraint"]
)
@pytest.mark.parametrize("pos_label", [0, 1])
def test_tunedthresholdclassifier_pos_label_tnr_tpr(objective_metric, pos_label):
    """Check that `pos_label` is dispatched correctly by checking the TNR and TPR
    score found during the optimization and the one found at `predict` time."""
    X, y = make_classification(n_samples=5_000, weights=[0.6, 0.4], random_state=42)

    # prefit the estimator to avoid variability due to the cross-validation
    estimator = LogisticRegression().fit(X, y)

    constraint_value = 0.7
    model = TunedThresholdClassifier(
        estimator,
        objective_metric=objective_metric,
        constraint_value=constraint_value,
        cv="prefit",
        pos_label=pos_label,
    ).fit(X, y)

    def tnr_tpr_score(y_true, y_pred, pos_label=pos_label):
        cm = confusion_matrix(y_true, y_pred)
        if pos_label == 0:
            cm = cm[::-1, ::-1]
        tn, fp, fn, tp = cm.ravel()
        tnr = tn / (tn + fp)
        tpr = tp / (tp + fn)
        return tnr, tpr

    tnr, tpr = tnr_tpr_score(y, model.predict(X), pos_label=pos_label)
    # due to internal interpolation, the scores will vary slightly
    if objective_metric == "max_tnr_at_tpr_constraint":
        assert tpr == pytest.approx(model.objective_score_[0], abs=0.05)
        assert tnr == pytest.approx(model.objective_score_[1], abs=0.05)
    else:
        assert tnr == pytest.approx(model.objective_score_[0], abs=0.05)
        assert tpr == pytest.approx(model.objective_score_[1], abs=0.05)


@pytest.mark.parametrize(
    "metric_type",
    ["string", "scorer_without_pos_label", "scorer_with_pos_label"],
)
@pytest.mark.parametrize("pos_label", [0, 1])
def test_tunedthresholdclassifier_pos_label_single_metric(pos_label, metric_type):
    """Check that `pos_label` is dispatched correctly when getting a scorer linked to
    a known metric. By default, the scorer in scikit-learn only have a default value
    for `pos_label` which is 1.
    """
    X, y = make_classification(n_samples=100, weights=[0.6, 0.4], random_state=42)

    # prefit the estimator to avoid variability due to the cross-validation
    estimator = LogisticRegression().fit(X, y)

    if metric_type == "string":
        objective_metric = "precision"
    elif metric_type == "scorer_without_pos_label":
        objective_metric = make_scorer(precision_score)
    else:  # metric_type == "scorer_with_pos_label"
        objective_metric = make_scorer(precision_score, pos_label=pos_label)

    model = TunedThresholdClassifier(
        estimator,
        objective_metric=objective_metric,
        cv="prefit",
        pos_label=pos_label,
        n_thresholds=500,
    ).fit(X, y)

    precision = precision_score(y, model.predict(X), pos_label=pos_label)
    assert precision == pytest.approx(model.objective_score_, abs=1e-3)


@pytest.mark.parametrize(
    "predict_method",
    ["predict", "predict_proba", "decision_function", "predict_log_proba"],
)
def test_tunedthresholdclassifier_constant_strategy(predict_method):
    """Check the behavior when `strategy='contant'."""
    X, y = make_classification(n_samples=100, weights=[0.6, 0.4], random_state=42)

    # With a constant strategy and a threshold at 0.5, we should get the same than the
    # original model
    estimator = LogisticRegression().fit(X, y)
    constant_threshold = 0.5
    tuned_model = TunedThresholdClassifier(
        estimator, strategy="constant", constant_threshold=constant_threshold
    ).fit(X, y)
    assert tuned_model.decision_threshold_ == pytest.approx(constant_threshold)
    for attribute in ("decision_thresholds_", "objective_score_", "objective_scores_"):
        assert getattr(tuned_model, attribute) is None

    assert_allclose(
        getattr(tuned_model, predict_method)(X), getattr(estimator, predict_method)(X)
    )
