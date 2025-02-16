import numpy as np
import pytest

from sklearn import config_context
from sklearn.base import clone
from sklearn.datasets import (
    load_breast_cancer,
    load_iris,
    make_classification,
    make_multilabel_classification,
)
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    fbeta_score,
    make_scorer,
)
from sklearn.metrics._scorer import _CurveScorer
from sklearn.model_selection import (
    FixedThresholdClassifier,
    StratifiedShuffleSplit,
    TunedThresholdClassifierCV,
)
from sklearn.model_selection._classification_threshold import (
    _fit_and_score_over_thresholds,
)
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


def test_fit_and_score_over_thresholds_curve_scorers():
    """Check that `_fit_and_score_over_thresholds` returns thresholds in ascending order
    for the different accepted curve scorers."""
    X, y = make_classification(n_samples=100, random_state=0)
    train_idx, val_idx = np.arange(50), np.arange(50, 100)
    classifier = LogisticRegression()

    curve_scorer = _CurveScorer(
        score_func=balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={},
    )
    scores, thresholds = _fit_and_score_over_thresholds(
        classifier,
        X,
        y,
        fit_params={},
        train_idx=train_idx,
        val_idx=val_idx,
        curve_scorer=curve_scorer,
        score_params={},
    )

    assert np.all(thresholds[:-1] <= thresholds[1:])
    assert isinstance(scores, np.ndarray)
    assert np.logical_and(scores >= 0, scores <= 1).all()


def test_fit_and_score_over_thresholds_prefit():
    """Check the behaviour with a prefit classifier."""
    X, y = make_classification(n_samples=100, random_state=0)

    # `train_idx is None` to indicate that the classifier is prefit
    train_idx, val_idx = None, np.arange(50, 100)
    classifier = DecisionTreeClassifier(random_state=0).fit(X, y)
    # make sure that the classifier memorized the full dataset such that
    # we get perfect predictions and thus match the expected score
    assert classifier.score(X[val_idx], y[val_idx]) == pytest.approx(1.0)

    curve_scorer = _CurveScorer(
        score_func=balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=2,
        kwargs={},
    )
    scores, thresholds = _fit_and_score_over_thresholds(
        classifier,
        X,
        y,
        fit_params={},
        train_idx=train_idx,
        val_idx=val_idx,
        curve_scorer=curve_scorer,
        score_params={},
    )
    assert np.all(thresholds[:-1] <= thresholds[1:])
    assert_allclose(scores, [0.5, 1.0])


@config_context(enable_metadata_routing=True)
def test_fit_and_score_over_thresholds_sample_weight():
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
    curve_scorer = _CurveScorer(
        score_func=balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={},
    )
    scores_repeated, thresholds_repeated = _fit_and_score_over_thresholds(
        classifier,
        X_repeated,
        y_repeated,
        fit_params={},
        train_idx=train_repeated_idx,
        val_idx=val_repeated_idx,
        curve_scorer=curve_scorer,
        score_params={},
    )

    train_idx, val_idx = np.arange(X.shape[0]), np.arange(X.shape[0])
    scores, thresholds = _fit_and_score_over_thresholds(
        classifier.set_fit_request(sample_weight=True),
        X,
        y,
        fit_params={"sample_weight": sample_weight},
        train_idx=train_idx,
        val_idx=val_idx,
        curve_scorer=curve_scorer.set_score_request(sample_weight=True),
        score_params={"sample_weight": sample_weight},
    )

    assert_allclose(thresholds_repeated, thresholds)
    assert_allclose(scores_repeated, scores)


@pytest.mark.parametrize("fit_params_type", ["list", "array"])
@config_context(enable_metadata_routing=True)
def test_fit_and_score_over_thresholds_fit_params(fit_params_type):
    """Check that we pass `fit_params` to the classifier when calling `fit`."""
    X, y = make_classification(n_samples=100, random_state=0)
    fit_params = {
        "a": _convert_container(y, fit_params_type),
        "b": _convert_container(y, fit_params_type),
    }

    classifier = CheckingClassifier(expected_fit_params=["a", "b"], random_state=0)
    classifier.set_fit_request(a=True, b=True)
    train_idx, val_idx = np.arange(50), np.arange(50, 100)

    curve_scorer = _CurveScorer(
        score_func=balanced_accuracy_score,
        sign=1,
        response_method="predict_proba",
        thresholds=10,
        kwargs={},
    )
    _fit_and_score_over_thresholds(
        classifier,
        X,
        y,
        fit_params=fit_params,
        train_idx=train_idx,
        val_idx=val_idx,
        curve_scorer=curve_scorer,
        score_params={},
    )


@pytest.mark.parametrize(
    "data",
    [
        make_classification(n_classes=3, n_clusters_per_class=1, random_state=0),
        make_multilabel_classification(random_state=0),
    ],
)
def test_tuned_threshold_classifier_no_binary(data):
    """Check that we raise an informative error message for non-binary problem."""
    err_msg = "Only binary classification is supported."
    with pytest.raises(ValueError, match=err_msg):
        TunedThresholdClassifierCV(LogisticRegression()).fit(*data)


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
def test_tuned_threshold_classifier_conflict_cv_refit(params, err_type, err_msg):
    """Check that we raise an informative error message when `cv` and `refit`
    cannot be used together.
    """
    X, y = make_classification(n_samples=100, random_state=0)
    with pytest.raises(err_type, match=err_msg):
        TunedThresholdClassifierCV(LogisticRegression(), **params).fit(X, y)


@pytest.mark.parametrize(
    "estimator",
    [LogisticRegression(), SVC(), GradientBoostingClassifier(n_estimators=4)],
)
@pytest.mark.parametrize(
    "response_method", ["predict_proba", "predict_log_proba", "decision_function"]
)
@pytest.mark.parametrize(
    "ThresholdClassifier", [FixedThresholdClassifier, TunedThresholdClassifierCV]
)
def test_threshold_classifier_estimator_response_methods(
    ThresholdClassifier, estimator, response_method
):
    """Check that `TunedThresholdClassifierCV` exposes the same response methods as the
    underlying estimator.
    """
    X, y = make_classification(n_samples=100, random_state=0)

    model = ThresholdClassifier(estimator=estimator)
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
def test_tuned_threshold_classifier_without_constraint_value(response_method):
    """Check that `TunedThresholdClassifierCV` is optimizing a given objective
    metric."""
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
    thresholds = 100
    model = TunedThresholdClassifierCV(
        estimator=lr,
        scoring="balanced_accuracy",
        response_method=response_method,
        thresholds=thresholds,
        store_cv_results=True,
    )
    score_optimized = balanced_accuracy_score(y, model.fit(X, y).predict(X))
    score_baseline = balanced_accuracy_score(y, lr.predict(X))
    assert score_optimized > score_baseline
    assert model.cv_results_["thresholds"].shape == (thresholds,)
    assert model.cv_results_["scores"].shape == (thresholds,)


def test_tuned_threshold_classifier_metric_with_parameter():
    """Check that we can pass a metric with a parameter in addition check that
    `f_beta` with `beta=1` is equivalent to `f1` and different from `f_beta` with
    `beta=2`.
    """
    X, y = load_breast_cancer(return_X_y=True)
    lr = make_pipeline(StandardScaler(), LogisticRegression()).fit(X, y)
    model_fbeta_1 = TunedThresholdClassifierCV(
        estimator=lr, scoring=make_scorer(fbeta_score, beta=1)
    ).fit(X, y)
    model_fbeta_2 = TunedThresholdClassifierCV(
        estimator=lr, scoring=make_scorer(fbeta_score, beta=2)
    ).fit(X, y)
    model_f1 = TunedThresholdClassifierCV(
        estimator=lr, scoring=make_scorer(f1_score)
    ).fit(X, y)

    assert model_fbeta_1.best_threshold_ == pytest.approx(model_f1.best_threshold_)
    assert model_fbeta_1.best_threshold_ != pytest.approx(model_fbeta_2.best_threshold_)


@pytest.mark.parametrize(
    "response_method", ["auto", "decision_function", "predict_proba"]
)
@pytest.mark.parametrize(
    "metric",
    [
        make_scorer(balanced_accuracy_score),
        make_scorer(f1_score, pos_label="cancer"),
    ],
)
def test_tuned_threshold_classifier_with_string_targets(response_method, metric):
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
    model = TunedThresholdClassifierCV(
        estimator=make_pipeline(StandardScaler(), LogisticRegression()),
        scoring=metric,
        response_method=response_method,
        thresholds=100,
    ).fit(X, y)
    assert_array_equal(model.classes_, np.sort(classes))
    y_pred = model.predict(X)
    assert_array_equal(np.unique(y_pred), np.sort(classes))


@pytest.mark.parametrize("with_sample_weight", [True, False])
@config_context(enable_metadata_routing=True)
def test_tuned_threshold_classifier_refit(with_sample_weight, global_random_seed):
    """Check the behaviour of the `refit` parameter."""
    rng = np.random.RandomState(global_random_seed)
    X, y = make_classification(n_samples=100, random_state=0)
    if with_sample_weight:
        sample_weight = rng.randn(X.shape[0])
        sample_weight = np.abs(sample_weight, out=sample_weight)
    else:
        sample_weight = None

    # check that `estimator_` if fitted on the full dataset when `refit=True`
    estimator = LogisticRegression().set_fit_request(sample_weight=True)
    model = TunedThresholdClassifierCV(estimator, refit=True).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is not estimator
    estimator.fit(X, y, sample_weight=sample_weight)
    assert_allclose(model.estimator_.coef_, estimator.coef_)
    assert_allclose(model.estimator_.intercept_, estimator.intercept_)

    # check that `estimator_` was not altered when `refit=False` and `cv="prefit"`
    estimator = LogisticRegression().set_fit_request(sample_weight=True)
    estimator.fit(X, y, sample_weight=sample_weight)
    coef = estimator.coef_.copy()
    model = TunedThresholdClassifierCV(estimator, cv="prefit", refit=False).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is estimator
    assert_allclose(model.estimator_.coef_, coef)

    # check that we train `estimator_` on the training split of a given cross-validation
    estimator = LogisticRegression().set_fit_request(sample_weight=True)
    cv = [
        (np.arange(50), np.arange(50, 100)),
    ]  # single split
    model = TunedThresholdClassifierCV(estimator, cv=cv, refit=False).fit(
        X, y, sample_weight=sample_weight
    )

    assert model.estimator_ is not estimator
    if with_sample_weight:
        sw_train = sample_weight[cv[0][0]]
    else:
        sw_train = None
    estimator.fit(X[cv[0][0]], y[cv[0][0]], sample_weight=sw_train)
    assert_allclose(model.estimator_.coef_, estimator.coef_)


@pytest.mark.parametrize("fit_params_type", ["list", "array"])
@config_context(enable_metadata_routing=True)
def test_tuned_threshold_classifier_fit_params(fit_params_type):
    """Check that we pass `fit_params` to the classifier when calling `fit`."""
    X, y = make_classification(n_samples=100, random_state=0)
    fit_params = {
        "a": _convert_container(y, fit_params_type),
        "b": _convert_container(y, fit_params_type),
    }

    classifier = CheckingClassifier(expected_fit_params=["a", "b"], random_state=0)
    classifier.set_fit_request(a=True, b=True)
    model = TunedThresholdClassifierCV(classifier)
    model.fit(X, y, **fit_params)


@config_context(enable_metadata_routing=True)
def test_tuned_threshold_classifier_cv_zeros_sample_weights_equivalence():
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

    estimator = LogisticRegression().set_fit_request(sample_weight=True)
    model_without_weights = TunedThresholdClassifierCV(estimator, cv=2)
    model_with_weights = clone(model_without_weights)

    model_with_weights.fit(X, y, sample_weight=sample_weight)
    model_without_weights.fit(X[::2], y[::2])

    assert_allclose(
        model_with_weights.estimator_.coef_, model_without_weights.estimator_.coef_
    )

    y_pred_with_weights = model_with_weights.predict_proba(X)
    y_pred_without_weights = model_without_weights.predict_proba(X)
    assert_allclose(y_pred_with_weights, y_pred_without_weights)


def test_tuned_threshold_classifier_thresholds_array():
    """Check that we can pass an array to `thresholds` and it is used as candidate
    threshold internally."""
    X, y = make_classification(random_state=0)
    estimator = LogisticRegression()
    thresholds = np.linspace(0, 1, 11)
    tuned_model = TunedThresholdClassifierCV(
        estimator,
        thresholds=thresholds,
        response_method="predict_proba",
        store_cv_results=True,
    ).fit(X, y)
    assert_allclose(tuned_model.cv_results_["thresholds"], thresholds)


@pytest.mark.parametrize("store_cv_results", [True, False])
def test_tuned_threshold_classifier_store_cv_results(store_cv_results):
    """Check that if `cv_results_` exists depending on `store_cv_results`."""
    X, y = make_classification(random_state=0)
    estimator = LogisticRegression()
    tuned_model = TunedThresholdClassifierCV(
        estimator, store_cv_results=store_cv_results
    ).fit(X, y)
    if store_cv_results:
        assert hasattr(tuned_model, "cv_results_")
    else:
        assert not hasattr(tuned_model, "cv_results_")


def test_tuned_threshold_classifier_cv_float():
    """Check the behaviour when `cv` is set to a float."""
    X, y = make_classification(random_state=0)

    # case where `refit=False` and cv is a float: the underlying estimator will be fit
    # on the training set given by a ShuffleSplit. We check that we get the same model
    # coefficients.
    test_size = 0.3
    estimator = LogisticRegression()
    tuned_model = TunedThresholdClassifierCV(
        estimator, cv=test_size, refit=False, random_state=0
    ).fit(X, y)
    tuned_model.fit(X, y)

    cv = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=0)
    train_idx, val_idx = next(cv.split(X, y))
    cloned_estimator = clone(estimator).fit(X[train_idx], y[train_idx])

    assert_allclose(tuned_model.estimator_.coef_, cloned_estimator.coef_)

    # case where `refit=True`, then the underlying estimator is fitted on the full
    # dataset.
    tuned_model.set_params(refit=True).fit(X, y)
    cloned_estimator = clone(estimator).fit(X, y)

    assert_allclose(tuned_model.estimator_.coef_, cloned_estimator.coef_)


def test_tuned_threshold_classifier_error_constant_predictor():
    """Check that we raise a ValueError if the underlying classifier returns constant
    probabilities such that we cannot find any threshold.
    """
    X, y = make_classification(random_state=0)
    estimator = DummyClassifier(strategy="constant", constant=1)
    tuned_model = TunedThresholdClassifierCV(estimator, response_method="predict_proba")
    err_msg = "The provided estimator makes constant predictions"
    with pytest.raises(ValueError, match=err_msg):
        tuned_model.fit(X, y)


@pytest.mark.parametrize(
    "response_method", ["auto", "predict_proba", "decision_function"]
)
def test_fixed_threshold_classifier_equivalence_default(response_method):
    """Check that `FixedThresholdClassifier` has the same behaviour as the vanilla
    classifier.
    """
    X, y = make_classification(random_state=0)
    classifier = LogisticRegression().fit(X, y)
    classifier_default_threshold = FixedThresholdClassifier(
        estimator=clone(classifier), response_method=response_method
    )
    classifier_default_threshold.fit(X, y)

    # emulate the response method that should take into account the `pos_label`
    if response_method in ("auto", "predict_proba"):
        y_score = classifier_default_threshold.predict_proba(X)[:, 1]
        threshold = 0.5
    else:  # response_method == "decision_function"
        y_score = classifier_default_threshold.decision_function(X)
        threshold = 0.0

    y_pred_lr = (y_score >= threshold).astype(int)
    assert_allclose(classifier_default_threshold.predict(X), y_pred_lr)


@pytest.mark.parametrize(
    "response_method, threshold", [("predict_proba", 0.7), ("decision_function", 2.0)]
)
@pytest.mark.parametrize("pos_label", [0, 1])
def test_fixed_threshold_classifier(response_method, threshold, pos_label):
    """Check that applying `predict` lead to the same prediction as applying the
    threshold to the output of the response method.
    """
    X, y = make_classification(n_samples=50, random_state=0)
    logistic_regression = LogisticRegression().fit(X, y)
    model = FixedThresholdClassifier(
        estimator=clone(logistic_regression),
        threshold=threshold,
        response_method=response_method,
        pos_label=pos_label,
    ).fit(X, y)

    # check that the underlying estimator is the same
    assert_allclose(model.estimator_.coef_, logistic_regression.coef_)

    # emulate the response method that should take into account the `pos_label`
    if response_method == "predict_proba":
        y_score = model.predict_proba(X)[:, pos_label]
    else:  # response_method == "decision_function"
        y_score = model.decision_function(X)
        y_score = y_score if pos_label == 1 else -y_score

    # create a mapping from boolean values to class labels
    map_to_label = np.array([0, 1]) if pos_label == 1 else np.array([1, 0])
    y_pred_lr = map_to_label[(y_score >= threshold).astype(int)]
    assert_allclose(model.predict(X), y_pred_lr)

    for method in ("predict_proba", "predict_log_proba", "decision_function"):
        assert_allclose(
            getattr(model, method)(X), getattr(logistic_regression, method)(X)
        )
        assert_allclose(
            getattr(model.estimator_, method)(X),
            getattr(logistic_regression, method)(X),
        )


@config_context(enable_metadata_routing=True)
def test_fixed_threshold_classifier_metadata_routing():
    """Check that everything works with metadata routing."""
    X, y = make_classification(random_state=0)
    sample_weight = np.ones_like(y)
    sample_weight[::2] = 2
    classifier = LogisticRegression().set_fit_request(sample_weight=True)
    classifier.fit(X, y, sample_weight=sample_weight)
    classifier_default_threshold = FixedThresholdClassifier(estimator=clone(classifier))
    classifier_default_threshold.fit(X, y, sample_weight=sample_weight)
    assert_allclose(classifier_default_threshold.estimator_.coef_, classifier.coef_)


@pytest.mark.parametrize(
    "method", ["predict_proba", "decision_function", "predict", "predict_log_proba"]
)
def test_fixed_threshold_classifier_fitted_estimator(method):
    """Check that if the underlying estimator is already fitted, no fit is required."""
    X, y = make_classification(random_state=0)
    classifier = LogisticRegression().fit(X, y)
    fixed_threshold_classifier = FixedThresholdClassifier(estimator=classifier)
    # This should not raise an error
    getattr(fixed_threshold_classifier, method)(X)


def test_fixed_threshold_classifier_classes_():
    """Check that the classes_ attribute is properly set."""
    X, y = make_classification(random_state=0)
    with pytest.raises(
        AttributeError, match="The underlying estimator is not fitted yet."
    ):
        FixedThresholdClassifier(estimator=LogisticRegression()).classes_

    classifier = LogisticRegression().fit(X, y)
    fixed_threshold_classifier = FixedThresholdClassifier(estimator=classifier)
    assert_array_equal(fixed_threshold_classifier.classes_, classifier.classes_)
