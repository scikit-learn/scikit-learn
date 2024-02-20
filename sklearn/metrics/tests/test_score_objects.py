import numbers
import os
import pickle
import shutil
import tempfile
from copy import deepcopy
from functools import partial
from unittest.mock import Mock

import joblib
import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import config_context
from sklearn.base import BaseEstimator
from sklearn.cluster import KMeans
from sklearn.datasets import (
    load_diabetes,
    make_blobs,
    make_classification,
    make_multilabel_classification,
    make_regression,
)
from sklearn.linear_model import LogisticRegression, Perceptron, Ridge
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    check_scoring,
    f1_score,
    fbeta_score,
    get_scorer,
    get_scorer_names,
    jaccard_score,
    log_loss,
    make_scorer,
    matthews_corrcoef,
    precision_score,
    r2_score,
    recall_score,
    roc_auc_score,
    top_k_accuracy_score,
)
from sklearn.metrics import cluster as cluster_module
from sklearn.metrics._scorer import (
    _check_multimetric_scoring,
    _MultimetricScorer,
    _PassthroughScorer,
    _Scorer,
)
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.svm import LinearSVC
from sklearn.tests.metadata_routing_common import (
    assert_request_equal,
    assert_request_is_empty,
)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils._testing import (
    assert_almost_equal,
    assert_array_equal,
    ignore_warnings,
)
from sklearn.utils.metadata_routing import MetadataRouter

REGRESSION_SCORERS = [
    "explained_variance",
    "r2",
    "neg_mean_absolute_error",
    "neg_mean_squared_error",
    "neg_mean_absolute_percentage_error",
    "neg_mean_squared_log_error",
    "neg_median_absolute_error",
    "neg_root_mean_squared_error",
    "neg_root_mean_squared_log_error",
    "mean_absolute_error",
    "mean_absolute_percentage_error",
    "mean_squared_error",
    "median_absolute_error",
    "max_error",
    "neg_mean_poisson_deviance",
    "neg_mean_gamma_deviance",
]

CLF_SCORERS = [
    "accuracy",
    "balanced_accuracy",
    "top_k_accuracy",
    "f1",
    "f1_weighted",
    "f1_macro",
    "f1_micro",
    "roc_auc",
    "average_precision",
    "precision",
    "precision_weighted",
    "precision_macro",
    "precision_micro",
    "recall",
    "recall_weighted",
    "recall_macro",
    "recall_micro",
    "neg_log_loss",
    "neg_brier_score",
    "jaccard",
    "jaccard_weighted",
    "jaccard_macro",
    "jaccard_micro",
    "roc_auc_ovr",
    "roc_auc_ovo",
    "roc_auc_ovr_weighted",
    "roc_auc_ovo_weighted",
    "matthews_corrcoef",
    "positive_likelihood_ratio",
    "neg_negative_likelihood_ratio",
]

# All supervised cluster scorers (They behave like classification metric)
CLUSTER_SCORERS = [
    "adjusted_rand_score",
    "rand_score",
    "homogeneity_score",
    "completeness_score",
    "v_measure_score",
    "mutual_info_score",
    "adjusted_mutual_info_score",
    "normalized_mutual_info_score",
    "fowlkes_mallows_score",
]

MULTILABEL_ONLY_SCORERS = [
    "precision_samples",
    "recall_samples",
    "f1_samples",
    "jaccard_samples",
]

REQUIRE_POSITIVE_Y_SCORERS = ["neg_mean_poisson_deviance", "neg_mean_gamma_deviance"]


def _require_positive_y(y):
    """Make targets strictly positive"""
    offset = abs(y.min()) + 1
    y = y + offset
    return y


def _make_estimators(X_train, y_train, y_ml_train):
    # Make estimators that make sense to test various scoring methods
    sensible_regr = DecisionTreeRegressor(random_state=0)
    # some of the regressions scorers require strictly positive input.
    sensible_regr.fit(X_train, _require_positive_y(y_train))
    sensible_clf = DecisionTreeClassifier(random_state=0)
    sensible_clf.fit(X_train, y_train)
    sensible_ml_clf = DecisionTreeClassifier(random_state=0)
    sensible_ml_clf.fit(X_train, y_ml_train)
    return dict(
        [(name, sensible_regr) for name in REGRESSION_SCORERS]
        + [(name, sensible_clf) for name in CLF_SCORERS]
        + [(name, sensible_clf) for name in CLUSTER_SCORERS]
        + [(name, sensible_ml_clf) for name in MULTILABEL_ONLY_SCORERS]
    )


X_mm, y_mm, y_ml_mm = None, None, None
ESTIMATORS = None
TEMP_FOLDER = None


def setup_module():
    # Create some memory mapped data
    global X_mm, y_mm, y_ml_mm, TEMP_FOLDER, ESTIMATORS
    TEMP_FOLDER = tempfile.mkdtemp(prefix="sklearn_test_score_objects_")
    X, y = make_classification(n_samples=30, n_features=5, random_state=0)
    _, y_ml = make_multilabel_classification(n_samples=X.shape[0], random_state=0)
    filename = os.path.join(TEMP_FOLDER, "test_data.pkl")
    joblib.dump((X, y, y_ml), filename)
    X_mm, y_mm, y_ml_mm = joblib.load(filename, mmap_mode="r")
    ESTIMATORS = _make_estimators(X_mm, y_mm, y_ml_mm)


def teardown_module():
    global X_mm, y_mm, y_ml_mm, TEMP_FOLDER, ESTIMATORS
    # GC closes the mmap file descriptors
    X_mm, y_mm, y_ml_mm, ESTIMATORS = None, None, None, None
    shutil.rmtree(TEMP_FOLDER)


class EstimatorWithFit(BaseEstimator):
    """Dummy estimator to test scoring validators"""

    def fit(self, X, y):
        return self


class EstimatorWithFitAndScore:
    """Dummy estimator to test scoring validators"""

    def fit(self, X, y):
        return self

    def score(self, X, y):
        return 1.0


class EstimatorWithFitAndPredict:
    """Dummy estimator to test scoring validators"""

    def fit(self, X, y):
        self.y = y
        return self

    def predict(self, X):
        return self.y


class DummyScorer:
    """Dummy scorer that always returns 1."""

    def __call__(self, est, X, y):
        return 1


def test_all_scorers_repr():
    # Test that all scorers have a working repr
    for name in get_scorer_names():
        repr(get_scorer(name))


def check_scoring_validator_for_single_metric_usecases(scoring_validator):
    # Test all branches of single metric usecases
    estimator = EstimatorWithFitAndScore()
    estimator.fit([[1]], [1])
    scorer = scoring_validator(estimator)
    assert isinstance(scorer, _PassthroughScorer)
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFitAndPredict()
    estimator.fit([[1]], [1])
    pattern = (
        r"If no scoring is specified, the estimator passed should have"
        r" a 'score' method\. The estimator .* does not\."
    )
    with pytest.raises(TypeError, match=pattern):
        scoring_validator(estimator)

    scorer = scoring_validator(estimator, scoring="accuracy")
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFit()
    scorer = scoring_validator(estimator, scoring="accuracy")
    assert isinstance(scorer, _Scorer)
    assert scorer._response_method == "predict"

    # Test the allow_none parameter for check_scoring alone
    if scoring_validator is check_scoring:
        estimator = EstimatorWithFit()
        scorer = scoring_validator(estimator, allow_none=True)
        assert scorer is None


@pytest.mark.parametrize(
    "scoring",
    (
        ("accuracy",),
        ["precision"],
        {"acc": "accuracy", "precision": "precision"},
        ("accuracy", "precision"),
        ["precision", "accuracy"],
        {
            "accuracy": make_scorer(accuracy_score),
            "precision": make_scorer(precision_score),
        },
    ),
    ids=[
        "single_tuple",
        "single_list",
        "dict_str",
        "multi_tuple",
        "multi_list",
        "dict_callable",
    ],
)
def test_check_scoring_and_check_multimetric_scoring(scoring):
    check_scoring_validator_for_single_metric_usecases(check_scoring)
    # To make sure the check_scoring is correctly applied to the constituent
    # scorers

    estimator = LinearSVC(dual="auto", random_state=0)
    estimator.fit([[1], [2], [3]], [1, 1, 0])

    scorers = _check_multimetric_scoring(estimator, scoring)
    assert isinstance(scorers, dict)
    assert sorted(scorers.keys()) == sorted(list(scoring))
    assert all([isinstance(scorer, _Scorer) for scorer in list(scorers.values())])
    assert all(scorer._response_method == "predict" for scorer in scorers.values())

    if "acc" in scoring:
        assert_almost_equal(
            scorers["acc"](estimator, [[1], [2], [3]], [1, 0, 0]), 2.0 / 3.0
        )
    if "accuracy" in scoring:
        assert_almost_equal(
            scorers["accuracy"](estimator, [[1], [2], [3]], [1, 0, 0]), 2.0 / 3.0
        )
    if "precision" in scoring:
        assert_almost_equal(
            scorers["precision"](estimator, [[1], [2], [3]], [1, 0, 0]), 0.5
        )


@pytest.mark.parametrize(
    "scoring, msg",
    [
        (
            (make_scorer(precision_score), make_scorer(accuracy_score)),
            "One or more of the elements were callables",
        ),
        ([5], "Non-string types were found"),
        ((make_scorer(precision_score),), "One or more of the elements were callables"),
        ((), "Empty list was given"),
        (("f1", "f1"), "Duplicate elements were found"),
        ({4: "accuracy"}, "Non-string types were found in the keys"),
        ({}, "An empty dict was passed"),
    ],
    ids=[
        "tuple of callables",
        "list of int",
        "tuple of one callable",
        "empty tuple",
        "non-unique str",
        "non-string key dict",
        "empty dict",
    ],
)
def test_check_scoring_and_check_multimetric_scoring_errors(scoring, msg):
    # Make sure it raises errors when scoring parameter is not valid.
    # More weird corner cases are tested at test_validation.py
    estimator = EstimatorWithFitAndPredict()
    estimator.fit([[1]], [1])

    with pytest.raises(ValueError, match=msg):
        _check_multimetric_scoring(estimator, scoring=scoring)


def test_check_scoring_gridsearchcv():
    # test that check_scoring works on GridSearchCV and pipeline.
    # slightly redundant non-regression test.

    grid = GridSearchCV(LinearSVC(dual="auto"), param_grid={"C": [0.1, 1]}, cv=3)
    scorer = check_scoring(grid, scoring="f1")
    assert isinstance(scorer, _Scorer)
    assert scorer._response_method == "predict"

    pipe = make_pipeline(LinearSVC(dual="auto"))
    scorer = check_scoring(pipe, scoring="f1")
    assert isinstance(scorer, _Scorer)
    assert scorer._response_method == "predict"

    # check that cross_val_score definitely calls the scorer
    # and doesn't make any assumptions about the estimator apart from having a
    # fit.
    scores = cross_val_score(
        EstimatorWithFit(), [[1], [2], [3]], [1, 0, 1], scoring=DummyScorer(), cv=3
    )
    assert_array_equal(scores, 1)


@pytest.mark.parametrize(
    "scorer_name, metric",
    [
        ("f1", f1_score),
        ("f1_weighted", partial(f1_score, average="weighted")),
        ("f1_macro", partial(f1_score, average="macro")),
        ("f1_micro", partial(f1_score, average="micro")),
        ("precision", precision_score),
        ("precision_weighted", partial(precision_score, average="weighted")),
        ("precision_macro", partial(precision_score, average="macro")),
        ("precision_micro", partial(precision_score, average="micro")),
        ("recall", recall_score),
        ("recall_weighted", partial(recall_score, average="weighted")),
        ("recall_macro", partial(recall_score, average="macro")),
        ("recall_micro", partial(recall_score, average="micro")),
        ("jaccard", jaccard_score),
        ("jaccard_weighted", partial(jaccard_score, average="weighted")),
        ("jaccard_macro", partial(jaccard_score, average="macro")),
        ("jaccard_micro", partial(jaccard_score, average="micro")),
        ("top_k_accuracy", top_k_accuracy_score),
        ("matthews_corrcoef", matthews_corrcoef),
    ],
)
def test_classification_binary_scores(scorer_name, metric):
    # check consistency between score and scorer for scores supporting
    # binary classification.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(dual="auto", random_state=0)
    clf.fit(X_train, y_train)

    score = get_scorer(scorer_name)(clf, X_test, y_test)
    expected_score = metric(y_test, clf.predict(X_test))
    assert_almost_equal(score, expected_score)


@pytest.mark.parametrize(
    "scorer_name, metric",
    [
        ("accuracy", accuracy_score),
        ("balanced_accuracy", balanced_accuracy_score),
        ("f1_weighted", partial(f1_score, average="weighted")),
        ("f1_macro", partial(f1_score, average="macro")),
        ("f1_micro", partial(f1_score, average="micro")),
        ("precision_weighted", partial(precision_score, average="weighted")),
        ("precision_macro", partial(precision_score, average="macro")),
        ("precision_micro", partial(precision_score, average="micro")),
        ("recall_weighted", partial(recall_score, average="weighted")),
        ("recall_macro", partial(recall_score, average="macro")),
        ("recall_micro", partial(recall_score, average="micro")),
        ("jaccard_weighted", partial(jaccard_score, average="weighted")),
        ("jaccard_macro", partial(jaccard_score, average="macro")),
        ("jaccard_micro", partial(jaccard_score, average="micro")),
    ],
)
def test_classification_multiclass_scores(scorer_name, metric):
    # check consistency between score and scorer for scores supporting
    # multiclass classification.
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=30, random_state=0
    )

    # use `stratify` = y to ensure train and test sets capture all classes
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=0, stratify=y
    )

    clf = DecisionTreeClassifier(random_state=0)
    clf.fit(X_train, y_train)
    score = get_scorer(scorer_name)(clf, X_test, y_test)
    expected_score = metric(y_test, clf.predict(X_test))
    assert score == pytest.approx(expected_score)


def test_custom_scorer_pickling():
    # test that custom scorer can be pickled
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(dual="auto", random_state=0)
    clf.fit(X_train, y_train)

    scorer = make_scorer(fbeta_score, beta=2)
    score1 = scorer(clf, X_test, y_test)
    unpickled_scorer = pickle.loads(pickle.dumps(scorer))
    score2 = unpickled_scorer(clf, X_test, y_test)
    assert score1 == pytest.approx(score2)

    # smoke test the repr:
    repr(fbeta_score)


def test_regression_scorers():
    # Test regression scorers.
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = Ridge()
    clf.fit(X_train, y_train)
    score1 = get_scorer("r2")(clf, X_test, y_test)
    score2 = r2_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)


def test_thresholded_scorers():
    # Test scorers that take thresholds.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)
    score1 = get_scorer("roc_auc")(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    score3 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)
    assert_almost_equal(score1, score3)

    logscore = get_scorer("neg_log_loss")(clf, X_test, y_test)
    logloss = log_loss(y_test, clf.predict_proba(X_test))
    assert_almost_equal(-logscore, logloss)

    # same for an estimator without decision_function
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    score1 = get_scorer("roc_auc")(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)

    # test with a regressor (no decision_function)
    reg = DecisionTreeRegressor()
    reg.fit(X_train, y_train)
    err_msg = "DecisionTreeRegressor has none of the following attributes"
    with pytest.raises(AttributeError, match=err_msg):
        get_scorer("roc_auc")(reg, X_test, y_test)

    # Test that an exception is raised on more than two classes
    X, y = make_blobs(random_state=0, centers=3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf.fit(X_train, y_train)
    with pytest.raises(ValueError, match="multi_class must be in \\('ovo', 'ovr'\\)"):
        get_scorer("roc_auc")(clf, X_test, y_test)

    # test error is raised with a single class present in model
    # (predict_proba shape is not suitable for binary auc)
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = DecisionTreeClassifier()
    clf.fit(X_train, np.zeros_like(y_train))
    with pytest.raises(ValueError, match="need classifier with two classes"):
        get_scorer("roc_auc")(clf, X_test, y_test)

    # for proba scorers
    with pytest.raises(ValueError, match="need classifier with two classes"):
        get_scorer("neg_log_loss")(clf, X_test, y_test)


def test_thresholded_scorers_multilabel_indicator_data():
    # Test that the scorer work with multilabel-indicator format
    # for multilabel and multi-output multi-class classifier
    X, y = make_multilabel_classification(allow_unlabeled=False, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # Multi-output multi-class predict_proba
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    y_proba = clf.predict_proba(X_test)
    score1 = get_scorer("roc_auc")(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, np.vstack([p[:, -1] for p in y_proba]).T)
    assert_almost_equal(score1, score2)

    # Multilabel predict_proba
    clf = OneVsRestClassifier(DecisionTreeClassifier())
    clf.fit(X_train, y_train)
    score1 = get_scorer("roc_auc")(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test))
    assert_almost_equal(score1, score2)

    # Multilabel decision function
    clf = OneVsRestClassifier(LinearSVC(dual="auto", random_state=0))
    clf.fit(X_train, y_train)
    score1 = get_scorer("roc_auc")(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    assert_almost_equal(score1, score2)


def test_supervised_cluster_scorers():
    # Test clustering scorers against gold standard labeling.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    km = KMeans(n_clusters=3, n_init="auto")
    km.fit(X_train)
    for name in CLUSTER_SCORERS:
        score1 = get_scorer(name)(km, X_test, y_test)
        score2 = getattr(cluster_module, name)(y_test, km.predict(X_test))
        assert_almost_equal(score1, score2)


@ignore_warnings
def test_raises_on_score_list():
    # Test that when a list of scores is returned, we raise proper errors.
    X, y = make_blobs(random_state=0)
    f1_scorer_no_average = make_scorer(f1_score, average=None)
    clf = DecisionTreeClassifier()
    with pytest.raises(ValueError):
        cross_val_score(clf, X, y, scoring=f1_scorer_no_average)
    grid_search = GridSearchCV(
        clf, scoring=f1_scorer_no_average, param_grid={"max_depth": [1, 2]}
    )
    with pytest.raises(ValueError):
        grid_search.fit(X, y)


@ignore_warnings
def test_classification_scorer_sample_weight():
    # Test that classification scorers support sample_weight or raise sensible
    # errors

    # Unlike the metrics invariance test, in the scorer case it's harder
    # to ensure that, on the classifier output, weighted and unweighted
    # scores really should be unequal.
    X, y = make_classification(random_state=0)
    _, y_ml = make_multilabel_classification(n_samples=X.shape[0], random_state=0)
    split = train_test_split(X, y, y_ml, random_state=0)
    X_train, X_test, y_train, y_test, y_ml_train, y_ml_test = split

    sample_weight = np.ones_like(y_test)
    sample_weight[:10] = 0

    # get sensible estimators for each metric
    estimator = _make_estimators(X_train, y_train, y_ml_train)

    for name in get_scorer_names():
        scorer = get_scorer(name)
        if name in REGRESSION_SCORERS:
            # skip the regression scores
            continue
        if name == "top_k_accuracy":
            # in the binary case k > 1 will always lead to a perfect score
            scorer._kwargs = {"k": 1}
        if name in MULTILABEL_ONLY_SCORERS:
            target = y_ml_test
        else:
            target = y_test
        try:
            weighted = scorer(
                estimator[name], X_test, target, sample_weight=sample_weight
            )
            ignored = scorer(estimator[name], X_test[10:], target[10:])
            unweighted = scorer(estimator[name], X_test, target)
            # this should not raise. sample_weight should be ignored if None.
            _ = scorer(estimator[name], X_test[:10], target[:10], sample_weight=None)
            assert weighted != unweighted, (
                f"scorer {name} behaves identically when called with "
                f"sample weights: {weighted} vs {unweighted}"
            )
            assert_almost_equal(
                weighted,
                ignored,
                err_msg=(
                    f"scorer {name} behaves differently "
                    "when ignoring samples and setting "
                    f"sample_weight to 0: {weighted} vs {ignored}"
                ),
            )

        except TypeError as e:
            assert "sample_weight" in str(e), (
                f"scorer {name} raises unhelpful exception when called "
                f"with sample weights: {str(e)}"
            )


@ignore_warnings
def test_regression_scorer_sample_weight():
    # Test that regression scorers support sample_weight or raise sensible
    # errors

    # Odd number of test samples req for neg_median_absolute_error
    X, y = make_regression(n_samples=101, n_features=20, random_state=0)
    y = _require_positive_y(y)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    sample_weight = np.ones_like(y_test)
    # Odd number req for neg_median_absolute_error
    sample_weight[:11] = 0

    reg = DecisionTreeRegressor(random_state=0)
    reg.fit(X_train, y_train)

    for name in get_scorer_names():
        scorer = get_scorer(name)
        if name not in REGRESSION_SCORERS:
            # skip classification scorers
            continue
        try:
            weighted = scorer(reg, X_test, y_test, sample_weight=sample_weight)
            ignored = scorer(reg, X_test[11:], y_test[11:])
            unweighted = scorer(reg, X_test, y_test)
            assert weighted != unweighted, (
                f"scorer {name} behaves identically when called with "
                f"sample weights: {weighted} vs {unweighted}"
            )
            assert_almost_equal(
                weighted,
                ignored,
                err_msg=(
                    f"scorer {name} behaves differently "
                    "when ignoring samples and setting "
                    f"sample_weight to 0: {weighted} vs {ignored}"
                ),
            )

        except TypeError as e:
            assert "sample_weight" in str(e), (
                f"scorer {name} raises unhelpful exception when called "
                f"with sample weights: {str(e)}"
            )


@pytest.mark.parametrize("name", get_scorer_names())
def test_scorer_memmap_input(name):
    # Non-regression test for #6147: some score functions would
    # return singleton memmap when computed on memmap data instead of scalar
    # float values.

    if name in REQUIRE_POSITIVE_Y_SCORERS:
        y_mm_1 = _require_positive_y(y_mm)
        y_ml_mm_1 = _require_positive_y(y_ml_mm)
    else:
        y_mm_1, y_ml_mm_1 = y_mm, y_ml_mm

    # UndefinedMetricWarning for P / R scores
    with ignore_warnings():
        scorer, estimator = get_scorer(name), ESTIMATORS[name]
        if name in MULTILABEL_ONLY_SCORERS:
            score = scorer(estimator, X_mm, y_ml_mm_1)
        else:
            score = scorer(estimator, X_mm, y_mm_1)
        assert isinstance(score, numbers.Number), name


def test_scoring_is_not_metric():
    with pytest.raises(ValueError, match="make_scorer"):
        check_scoring(LogisticRegression(), scoring=f1_score)
    with pytest.raises(ValueError, match="make_scorer"):
        check_scoring(LogisticRegression(), scoring=roc_auc_score)
    with pytest.raises(ValueError, match="make_scorer"):
        check_scoring(Ridge(), scoring=r2_score)
    with pytest.raises(ValueError, match="make_scorer"):
        check_scoring(KMeans(), scoring=cluster_module.adjusted_rand_score)
    with pytest.raises(ValueError, match="make_scorer"):
        check_scoring(KMeans(), scoring=cluster_module.rand_score)


@pytest.mark.parametrize(
    (
        "scorers,expected_predict_count,"
        "expected_predict_proba_count,expected_decision_func_count"
    ),
    [
        (
            {
                "a1": "accuracy",
                "a2": "accuracy",
                "ll1": "neg_log_loss",
                "ll2": "neg_log_loss",
                "ra1": "roc_auc",
                "ra2": "roc_auc",
            },
            1,
            1,
            1,
        ),
        (["roc_auc", "accuracy"], 1, 0, 1),
        (["neg_log_loss", "accuracy"], 1, 1, 0),
    ],
)
def test_multimetric_scorer_calls_method_once(
    scorers,
    expected_predict_count,
    expected_predict_proba_count,
    expected_decision_func_count,
):
    X, y = np.array([[1], [1], [0], [0], [0]]), np.array([0, 1, 1, 1, 0])

    mock_est = Mock()
    mock_est._estimator_type = "classifier"
    fit_func = Mock(return_value=mock_est, name="fit")
    fit_func.__name__ = "fit"
    predict_func = Mock(return_value=y, name="predict")
    predict_func.__name__ = "predict"

    pos_proba = np.random.rand(X.shape[0])
    proba = np.c_[1 - pos_proba, pos_proba]
    predict_proba_func = Mock(return_value=proba, name="predict_proba")
    predict_proba_func.__name__ = "predict_proba"
    decision_function_func = Mock(return_value=pos_proba, name="decision_function")
    decision_function_func.__name__ = "decision_function"

    mock_est.fit = fit_func
    mock_est.predict = predict_func
    mock_est.predict_proba = predict_proba_func
    mock_est.decision_function = decision_function_func
    # add the classes that would be found during fit
    mock_est.classes_ = np.array([0, 1])

    scorer_dict = _check_multimetric_scoring(LogisticRegression(), scorers)
    multi_scorer = _MultimetricScorer(scorers=scorer_dict)
    results = multi_scorer(mock_est, X, y)

    assert set(scorers) == set(results)  # compare dict keys

    assert predict_func.call_count == expected_predict_count
    assert predict_proba_func.call_count == expected_predict_proba_count
    assert decision_function_func.call_count == expected_decision_func_count


@pytest.mark.parametrize(
    "scorers",
    [
        (["roc_auc", "neg_log_loss"]),
        (
            {
                "roc_auc": make_scorer(
                    roc_auc_score,
                    response_method=["predict_proba", "decision_function"],
                ),
                "neg_log_loss": make_scorer(log_loss, response_method="predict_proba"),
            }
        ),
    ],
)
def test_multimetric_scorer_calls_method_once_classifier_no_decision(scorers):
    predict_proba_call_cnt = 0

    class MockKNeighborsClassifier(KNeighborsClassifier):
        def predict_proba(self, X):
            nonlocal predict_proba_call_cnt
            predict_proba_call_cnt += 1
            return super().predict_proba(X)

    X, y = np.array([[1], [1], [0], [0], [0]]), np.array([0, 1, 1, 1, 0])

    # no decision function
    clf = MockKNeighborsClassifier(n_neighbors=1)
    clf.fit(X, y)

    scorer_dict = _check_multimetric_scoring(clf, scorers)
    scorer = _MultimetricScorer(scorers=scorer_dict)
    scorer(clf, X, y)

    assert predict_proba_call_cnt == 1


def test_multimetric_scorer_calls_method_once_regressor_threshold():
    predict_called_cnt = 0

    class MockDecisionTreeRegressor(DecisionTreeRegressor):
        def predict(self, X):
            nonlocal predict_called_cnt
            predict_called_cnt += 1
            return super().predict(X)

    X, y = np.array([[1], [1], [0], [0], [0]]), np.array([0, 1, 1, 1, 0])

    # no decision function
    clf = MockDecisionTreeRegressor()
    clf.fit(X, y)

    scorers = {"neg_mse": "neg_mean_squared_error", "r2": "r2"}
    scorer_dict = _check_multimetric_scoring(clf, scorers)
    scorer = _MultimetricScorer(scorers=scorer_dict)
    scorer(clf, X, y)

    assert predict_called_cnt == 1


def test_multimetric_scorer_sanity_check():
    # scoring dictionary returned is the same as calling each scorer separately
    scorers = {
        "a1": "accuracy",
        "a2": "accuracy",
        "ll1": "neg_log_loss",
        "ll2": "neg_log_loss",
        "ra1": "roc_auc",
        "ra2": "roc_auc",
    }

    X, y = make_classification(random_state=0)

    clf = DecisionTreeClassifier()
    clf.fit(X, y)

    scorer_dict = _check_multimetric_scoring(clf, scorers)
    multi_scorer = _MultimetricScorer(scorers=scorer_dict)

    result = multi_scorer(clf, X, y)

    separate_scores = {
        name: get_scorer(name)(clf, X, y)
        for name in ["accuracy", "neg_log_loss", "roc_auc"]
    }

    for key, value in result.items():
        score_name = scorers[key]
        assert_allclose(value, separate_scores[score_name])


@pytest.mark.parametrize("raise_exc", [True, False])
def test_multimetric_scorer_exception_handling(raise_exc):
    """Check that the calling of the `_MultimetricScorer` returns
    exception messages in the result dict for the failing scorers
    in case of `raise_exc` is `False` and if `raise_exc` is `True`,
    then the proper exception is raised.
    """
    scorers = {
        "failing_1": "neg_mean_squared_log_error",
        "non_failing": "neg_median_absolute_error",
        "failing_2": "neg_mean_squared_log_error",
    }

    X, y = make_classification(
        n_samples=50, n_features=2, n_redundant=0, random_state=0
    )
    y *= -1  # neg_mean_squared_log_error fails if y contains negative values

    clf = DecisionTreeClassifier().fit(X, y)

    scorer_dict = _check_multimetric_scoring(clf, scorers)
    multi_scorer = _MultimetricScorer(scorers=scorer_dict, raise_exc=raise_exc)

    error_msg = (
        "Mean Squared Logarithmic Error cannot be used when targets contain"
        " negative values."
    )

    if raise_exc:
        with pytest.raises(ValueError, match=error_msg):
            multi_scorer(clf, X, y)
    else:
        result = multi_scorer(clf, X, y)

        exception_message_1 = result["failing_1"]
        score = result["non_failing"]
        exception_message_2 = result["failing_2"]

        assert isinstance(exception_message_1, str) and error_msg in exception_message_1
        assert isinstance(score, float)
        assert isinstance(exception_message_2, str) and error_msg in exception_message_2


@pytest.mark.parametrize(
    "scorer_name, metric",
    [
        ("roc_auc_ovr", partial(roc_auc_score, multi_class="ovr")),
        ("roc_auc_ovo", partial(roc_auc_score, multi_class="ovo")),
        (
            "roc_auc_ovr_weighted",
            partial(roc_auc_score, multi_class="ovr", average="weighted"),
        ),
        (
            "roc_auc_ovo_weighted",
            partial(roc_auc_score, multi_class="ovo", average="weighted"),
        ),
    ],
)
def test_multiclass_roc_proba_scorer(scorer_name, metric):
    scorer = get_scorer(scorer_name)
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=20, random_state=0
    )
    lr = LogisticRegression(multi_class="multinomial").fit(X, y)
    y_proba = lr.predict_proba(X)
    expected_score = metric(y, y_proba)

    assert scorer(lr, X, y) == pytest.approx(expected_score)


def test_multiclass_roc_proba_scorer_label():
    scorer = make_scorer(
        roc_auc_score,
        multi_class="ovo",
        labels=[0, 1, 2],
        response_method="predict_proba",
    )
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=20, random_state=0
    )
    lr = LogisticRegression(multi_class="multinomial").fit(X, y)
    y_proba = lr.predict_proba(X)

    y_binary = y == 0
    expected_score = roc_auc_score(
        y_binary, y_proba, multi_class="ovo", labels=[0, 1, 2]
    )

    assert scorer(lr, X, y_binary) == pytest.approx(expected_score)


@pytest.mark.parametrize(
    "scorer_name",
    ["roc_auc_ovr", "roc_auc_ovo", "roc_auc_ovr_weighted", "roc_auc_ovo_weighted"],
)
def test_multiclass_roc_no_proba_scorer_errors(scorer_name):
    # Perceptron has no predict_proba
    scorer = get_scorer(scorer_name)
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=20, random_state=0
    )
    lr = Perceptron().fit(X, y)
    msg = "Perceptron has none of the following attributes: predict_proba."
    with pytest.raises(AttributeError, match=msg):
        scorer(lr, X, y)


@pytest.fixture
def string_labeled_classification_problem():
    """Train a classifier on binary problem with string target.

    The classifier is trained on a binary classification problem where the
    minority class of interest has a string label that is intentionally not the
    greatest class label using the lexicographic order. In this case, "cancer"
    is the positive label, and `classifier.classes_` is
    `["cancer", "not cancer"]`.

    In addition, the dataset is imbalanced to better identify problems when
    using non-symmetric performance metrics such as f1-score, average precision
    and so on.

    Returns
    -------
    classifier : estimator object
        Trained classifier on the binary problem.
    X_test : ndarray of shape (n_samples, n_features)
        Data to be used as testing set in tests.
    y_test : ndarray of shape (n_samples,), dtype=object
        Binary target where labels are strings.
    y_pred : ndarray of shape (n_samples,), dtype=object
        Prediction of `classifier` when predicting for `X_test`.
    y_pred_proba : ndarray of shape (n_samples, 2), dtype=np.float64
        Probabilities of `classifier` when predicting for `X_test`.
    y_pred_decision : ndarray of shape (n_samples,), dtype=np.float64
        Decision function values of `classifier` when predicting on `X_test`.
    """
    from sklearn.datasets import load_breast_cancer
    from sklearn.utils import shuffle

    X, y = load_breast_cancer(return_X_y=True)
    # create an highly imbalanced classification task
    idx_positive = np.flatnonzero(y == 1)
    idx_negative = np.flatnonzero(y == 0)
    idx_selected = np.hstack([idx_negative, idx_positive[:25]])
    X, y = X[idx_selected], y[idx_selected]
    X, y = shuffle(X, y, random_state=42)
    # only use 2 features to make the problem even harder
    X = X[:, :2]
    y = np.array(["cancer" if c == 1 else "not cancer" for c in y], dtype=object)
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        stratify=y,
        random_state=0,
    )
    classifier = LogisticRegression().fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    y_pred_proba = classifier.predict_proba(X_test)
    y_pred_decision = classifier.decision_function(X_test)

    return classifier, X_test, y_test, y_pred, y_pred_proba, y_pred_decision


def test_average_precision_pos_label(string_labeled_classification_problem):
    # check that _Scorer will lead to the right score when passing
    # `pos_label`. Currently, only `average_precision_score` is defined to
    # be such a scorer.
    (
        clf,
        X_test,
        y_test,
        _,
        y_pred_proba,
        y_pred_decision,
    ) = string_labeled_classification_problem

    pos_label = "cancer"
    # we need to select the positive column or reverse the decision values
    y_pred_proba = y_pred_proba[:, 0]
    y_pred_decision = y_pred_decision * -1
    assert clf.classes_[0] == pos_label

    # check that when calling the scoring function, probability estimates and
    # decision values lead to the same results
    ap_proba = average_precision_score(y_test, y_pred_proba, pos_label=pos_label)
    ap_decision_function = average_precision_score(
        y_test, y_pred_decision, pos_label=pos_label
    )
    assert ap_proba == pytest.approx(ap_decision_function)

    # create a scorer which would require to pass a `pos_label`
    # check that it fails if `pos_label` is not provided
    average_precision_scorer = make_scorer(
        average_precision_score,
        response_method=("decision_function", "predict_proba"),
    )
    err_msg = "pos_label=1 is not a valid label. It should be one of "
    with pytest.raises(ValueError, match=err_msg):
        average_precision_scorer(clf, X_test, y_test)

    # otherwise, the scorer should give the same results than calling the
    # scoring function
    average_precision_scorer = make_scorer(
        average_precision_score,
        response_method=("decision_function", "predict_proba"),
        pos_label=pos_label,
    )
    ap_scorer = average_precision_scorer(clf, X_test, y_test)

    assert ap_scorer == pytest.approx(ap_proba)

    # The above scorer call is using `clf.decision_function`. We will force
    # it to use `clf.predict_proba`.
    clf_without_predict_proba = deepcopy(clf)

    def _predict_proba(self, X):
        raise NotImplementedError

    clf_without_predict_proba.predict_proba = partial(
        _predict_proba, clf_without_predict_proba
    )
    # sanity check
    with pytest.raises(NotImplementedError):
        clf_without_predict_proba.predict_proba(X_test)

    ap_scorer = average_precision_scorer(clf_without_predict_proba, X_test, y_test)
    assert ap_scorer == pytest.approx(ap_proba)


def test_brier_score_loss_pos_label(string_labeled_classification_problem):
    # check that _Scorer leads to the right score when `pos_label` is
    # provided. Currently only the `brier_score_loss` is defined to be such
    # a scorer.
    clf, X_test, y_test, _, y_pred_proba, _ = string_labeled_classification_problem

    pos_label = "cancer"
    assert clf.classes_[0] == pos_label

    # brier score loss is symmetric
    brier_pos_cancer = brier_score_loss(y_test, y_pred_proba[:, 0], pos_label="cancer")
    brier_pos_not_cancer = brier_score_loss(
        y_test, y_pred_proba[:, 1], pos_label="not cancer"
    )
    assert brier_pos_cancer == pytest.approx(brier_pos_not_cancer)

    brier_scorer = make_scorer(
        brier_score_loss,
        response_method="predict_proba",
        pos_label=pos_label,
    )
    assert brier_scorer(clf, X_test, y_test) == pytest.approx(brier_pos_cancer)


@pytest.mark.parametrize(
    "score_func", [f1_score, precision_score, recall_score, jaccard_score]
)
def test_non_symmetric_metric_pos_label(
    score_func, string_labeled_classification_problem
):
    # check that _Scorer leads to the right score when `pos_label` is
    # provided. We check for all possible metric supported.
    # Note: At some point we may end up having "scorer tags".
    clf, X_test, y_test, y_pred, _, _ = string_labeled_classification_problem

    pos_label = "cancer"
    assert clf.classes_[0] == pos_label

    score_pos_cancer = score_func(y_test, y_pred, pos_label="cancer")
    score_pos_not_cancer = score_func(y_test, y_pred, pos_label="not cancer")

    assert score_pos_cancer != pytest.approx(score_pos_not_cancer)

    scorer = make_scorer(score_func, pos_label=pos_label)
    assert scorer(clf, X_test, y_test) == pytest.approx(score_pos_cancer)


@pytest.mark.parametrize(
    "scorer",
    [
        make_scorer(
            average_precision_score,
            response_method=("decision_function", "predict_proba"),
            pos_label="xxx",
        ),
        make_scorer(brier_score_loss, response_method="predict_proba", pos_label="xxx"),
        make_scorer(f1_score, pos_label="xxx"),
    ],
    ids=["non-thresholded scorer", "probability scorer", "thresholded scorer"],
)
def test_scorer_select_proba_error(scorer):
    # check that we raise the proper error when passing an unknown
    # pos_label
    X, y = make_classification(
        n_classes=2, n_informative=3, n_samples=20, random_state=0
    )
    lr = LogisticRegression().fit(X, y)
    assert scorer._kwargs["pos_label"] not in np.unique(y).tolist()

    err_msg = "is not a valid label"
    with pytest.raises(ValueError, match=err_msg):
        scorer(lr, X, y)


def test_get_scorer_return_copy():
    # test that get_scorer returns a copy
    assert get_scorer("roc_auc") is not get_scorer("roc_auc")


def test_scorer_no_op_multiclass_select_proba():
    # check that calling a _Scorer on a multiclass problem do not raise
    # even if `y_true` would be binary during the scoring.
    # `_select_proba_binary` should not be called in this case.
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=20, random_state=0
    )
    lr = LogisticRegression().fit(X, y)

    mask_last_class = y == lr.classes_[-1]
    X_test, y_test = X[~mask_last_class], y[~mask_last_class]
    assert_array_equal(np.unique(y_test), lr.classes_[:-1])

    scorer = make_scorer(
        roc_auc_score,
        response_method="predict_proba",
        multi_class="ovo",
        labels=lr.classes_,
    )
    scorer(lr, X_test, y_test)


@pytest.mark.parametrize("name", get_scorer_names())
def test_scorer_set_score_request_raises(name):
    """Test that set_score_request is only available when feature flag is on."""
    # Make sure they expose the routing methods.
    scorer = get_scorer(name)
    with pytest.raises(RuntimeError, match="This method is only available"):
        scorer.set_score_request()


@pytest.mark.usefixtures("enable_slep006")
@pytest.mark.parametrize("name", get_scorer_names(), ids=get_scorer_names())
def test_scorer_metadata_request(name):
    """Testing metadata requests for scorers.

    This test checks many small things in a large test, to reduce the
    boilerplate required for each section.
    """
    # Make sure they expose the routing methods.
    scorer = get_scorer(name)
    assert hasattr(scorer, "set_score_request")
    assert hasattr(scorer, "get_metadata_routing")

    # Check that by default no metadata is requested.
    assert_request_is_empty(scorer.get_metadata_routing())

    weighted_scorer = scorer.set_score_request(sample_weight=True)
    # set_score_request should mutate the instance, rather than returning a
    # new instance
    assert weighted_scorer is scorer

    # make sure the scorer doesn't request anything on methods other than
    # `score`, and that the requested value on `score` is correct.
    assert_request_is_empty(weighted_scorer.get_metadata_routing(), exclude="score")
    assert (
        weighted_scorer.get_metadata_routing().score.requests["sample_weight"] is True
    )

    # make sure putting the scorer in a router doesn't request anything by
    # default
    router = MetadataRouter(owner="test").add(
        method_mapping="score", scorer=get_scorer(name)
    )
    # make sure `sample_weight` is refused if passed.
    with pytest.raises(TypeError, match="got unexpected argument"):
        router.validate_metadata(params={"sample_weight": 1}, method="score")
    # make sure `sample_weight` is not routed even if passed.
    routed_params = router.route_params(params={"sample_weight": 1}, caller="score")
    assert not routed_params.scorer.score

    # make sure putting weighted_scorer in a router requests sample_weight
    router = MetadataRouter(owner="test").add(
        scorer=weighted_scorer, method_mapping="score"
    )
    router.validate_metadata(params={"sample_weight": 1}, method="score")
    routed_params = router.route_params(params={"sample_weight": 1}, caller="score")
    assert list(routed_params.scorer.score.keys()) == ["sample_weight"]


@pytest.mark.usefixtures("enable_slep006")
def test_metadata_kwarg_conflict():
    """This test makes sure the right warning is raised if the user passes
    some metadata both as a constructor to make_scorer, and during __call__.
    """
    X, y = make_classification(
        n_classes=3, n_informative=3, n_samples=20, random_state=0
    )
    lr = LogisticRegression().fit(X, y)

    scorer = make_scorer(
        roc_auc_score,
        response_method="predict_proba",
        multi_class="ovo",
        labels=lr.classes_,
    )
    with pytest.warns(UserWarning, match="already set as kwargs"):
        scorer.set_score_request(labels=True)

    with pytest.warns(UserWarning, match="There is an overlap"):
        scorer(lr, X, y, labels=lr.classes_)


@pytest.mark.usefixtures("enable_slep006")
def test_PassthroughScorer_metadata_request():
    """Test that _PassthroughScorer properly routes metadata.

    _PassthroughScorer should behave like a consumer, mirroring whatever is the
    underlying score method.
    """
    scorer = _PassthroughScorer(
        estimator=LinearSVC()
        .set_score_request(sample_weight="alias")
        .set_fit_request(sample_weight=True)
    )
    # Test that _PassthroughScorer doesn't change estimator's routing.
    assert_request_equal(
        scorer.get_metadata_routing(),
        {"fit": {"sample_weight": True}, "score": {"sample_weight": "alias"}},
    )


@pytest.mark.usefixtures("enable_slep006")
def test_multimetric_scoring_metadata_routing():
    # Test that _MultimetricScorer properly routes metadata.
    def score1(y_true, y_pred):
        return 1

    def score2(y_true, y_pred, sample_weight="test"):
        # make sure sample_weight is not passed
        assert sample_weight == "test"
        return 1

    def score3(y_true, y_pred, sample_weight=None):
        # make sure sample_weight is passed
        assert sample_weight is not None
        return 1

    scorers = {
        "score1": make_scorer(score1),
        "score2": make_scorer(score2).set_score_request(sample_weight=False),
        "score3": make_scorer(score3).set_score_request(sample_weight=True),
    }

    X, y = make_classification(
        n_samples=50, n_features=2, n_redundant=0, random_state=0
    )

    clf = DecisionTreeClassifier().fit(X, y)

    scorer_dict = _check_multimetric_scoring(clf, scorers)
    multi_scorer = _MultimetricScorer(scorers=scorer_dict)
    # this should fail, because metadata routing is not enabled and w/o it we
    # don't support different metadata for different scorers.
    # TODO: remove when enable_metadata_routing is deprecated
    with config_context(enable_metadata_routing=False):
        with pytest.raises(TypeError, match="got an unexpected keyword argument"):
            multi_scorer(clf, X, y, sample_weight=1)

    # This passes since routing is done.
    multi_scorer(clf, X, y, sample_weight=1)


def test_kwargs_without_metadata_routing_error():
    # Test that kwargs are not supported in scorers if metadata routing is not
    # enabled.
    # TODO: remove when enable_metadata_routing is deprecated
    def score(y_true, y_pred, param=None):
        return 1  # pragma: no cover

    X, y = make_classification(
        n_samples=50, n_features=2, n_redundant=0, random_state=0
    )

    clf = DecisionTreeClassifier().fit(X, y)
    scorer = make_scorer(score)
    with config_context(enable_metadata_routing=False):
        with pytest.raises(
            ValueError, match="is only supported if enable_metadata_routing=True"
        ):
            scorer(clf, X, y, param="blah")


def test_get_scorer_multilabel_indicator():
    """Check that our scorer deal with multi-label indicator matrices.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/26817
    """
    X, Y = make_multilabel_classification(n_samples=72, n_classes=3, random_state=0)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, random_state=0)

    estimator = KNeighborsClassifier().fit(X_train, Y_train)

    score = get_scorer("average_precision")(estimator, X_test, Y_test)
    assert score > 0.8


@pytest.mark.parametrize(
    "scorer, expected_repr",
    [
        (
            get_scorer("accuracy"),
            "make_scorer(accuracy_score, response_method='predict')",
        ),
        (
            get_scorer("neg_log_loss"),
            (
                "make_scorer(log_loss, greater_is_better=False,"
                " response_method='predict_proba')"
            ),
        ),
        (
            get_scorer("roc_auc"),
            (
                "make_scorer(roc_auc_score, response_method="
                "('decision_function', 'predict_proba'))"
            ),
        ),
        (
            make_scorer(fbeta_score, beta=2),
            "make_scorer(fbeta_score, response_method='predict', beta=2)",
        ),
    ],
)
def test_make_scorer_repr(scorer, expected_repr):
    """Check the representation of the scorer."""
    assert repr(scorer) == expected_repr


# TODO(1.6): rework this test after the deprecation of `needs_proba` and
# `needs_threshold`
@pytest.mark.filterwarnings("ignore:.*needs_proba.*:FutureWarning")
@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        # response_method should not be set if needs_* are set
        (
            {"response_method": "predict_proba", "needs_proba": True},
            ValueError,
            "You cannot set both `response_method`",
        ),
        (
            {"response_method": "predict_proba", "needs_threshold": True},
            ValueError,
            "You cannot set both `response_method`",
        ),
        # cannot set both needs_proba and needs_threshold
        (
            {"needs_proba": True, "needs_threshold": True},
            ValueError,
            "You cannot set both `needs_proba` and `needs_threshold`",
        ),
    ],
)
def test_make_scorer_error(params, err_type, err_msg):
    """Check that `make_scorer` raises errors if the parameter used."""
    with pytest.raises(err_type, match=err_msg):
        make_scorer(lambda y_true, y_pred: 1, **params)


# TODO(1.6): remove the following test
@pytest.mark.parametrize(
    "deprecated_params, new_params, warn_msg",
    [
        (
            {"needs_proba": True},
            {"response_method": "predict_proba"},
            "The `needs_threshold` and `needs_proba` parameter are deprecated",
        ),
        (
            {"needs_proba": True, "needs_threshold": False},
            {"response_method": "predict_proba"},
            "The `needs_threshold` and `needs_proba` parameter are deprecated",
        ),
        (
            {"needs_threshold": True},
            {"response_method": ("decision_function", "predict_proba")},
            "The `needs_threshold` and `needs_proba` parameter are deprecated",
        ),
        (
            {"needs_threshold": True, "needs_proba": False},
            {"response_method": ("decision_function", "predict_proba")},
            "The `needs_threshold` and `needs_proba` parameter are deprecated",
        ),
        (
            {"needs_threshold": False, "needs_proba": False},
            {"response_method": "predict"},
            "The `needs_threshold` and `needs_proba` parameter are deprecated",
        ),
    ],
)
def test_make_scorer_deprecation(deprecated_params, new_params, warn_msg):
    """Check that we raise a deprecation warning when using `needs_proba` or
    `needs_threshold`."""
    X, y = make_classification(n_samples=150, n_features=10, random_state=0)
    classifier = LogisticRegression().fit(X, y)

    # check deprecation of needs_proba
    with pytest.warns(FutureWarning, match=warn_msg):
        deprecated_roc_auc_scorer = make_scorer(roc_auc_score, **deprecated_params)
    roc_auc_scorer = make_scorer(roc_auc_score, **new_params)

    assert deprecated_roc_auc_scorer(classifier, X, y) == pytest.approx(
        roc_auc_scorer(classifier, X, y)
    )


@pytest.mark.parametrize("pass_estimator", [True, False])
def test_get_scorer_multimetric(pass_estimator):
    """Check that check_scoring is compatible with multi-metric configurations."""
    X, y = make_classification(n_samples=150, n_features=10, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)

    if pass_estimator:
        check_scoring_ = check_scoring
    else:
        check_scoring_ = partial(check_scoring, clf)

    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)

    expected_results = {
        "r2": r2_score(y_test, y_pred),
        "roc_auc": roc_auc_score(y_test, y_proba[:, 1]),
        "accuracy": accuracy_score(y_test, y_pred),
    }

    for container in [set, list, tuple]:
        scoring = check_scoring_(scoring=container(["r2", "roc_auc", "accuracy"]))
        result = scoring(clf, X_test, y_test)

        assert result.keys() == expected_results.keys()
        for name in result:
            assert result[name] == pytest.approx(expected_results[name])

    def double_accuracy(y_true, y_pred):
        return 2 * accuracy_score(y_true, y_pred)

    custom_scorer = make_scorer(double_accuracy, response_method="predict")

    # dict with different names
    dict_scoring = check_scoring_(
        scoring={
            "my_r2": "r2",
            "my_roc_auc": "roc_auc",
            "double_accuracy": custom_scorer,
        }
    )
    dict_result = dict_scoring(clf, X_test, y_test)
    assert len(dict_result) == 3
    assert dict_result["my_r2"] == pytest.approx(expected_results["r2"])
    assert dict_result["my_roc_auc"] == pytest.approx(expected_results["roc_auc"])
    assert dict_result["double_accuracy"] == pytest.approx(
        2 * expected_results["accuracy"]
    )


def test_multimetric_scorer_repr():
    """Check repr for multimetric scorer"""
    multi_metric_scorer = check_scoring(scoring=["accuracy", "r2"])

    assert str(multi_metric_scorer) == 'MultiMetricScorer("accuracy", "r2")'


@pytest.mark.parametrize("enable_metadata_routing", [True, False])
def test_metadata_routing_multimetric_metadata_routing(enable_metadata_routing):
    """Test multimetric scorer works with and without metadata routing enabled when
    there is no actual metadata to pass.

    Non-regression test for https://github.com/scikit-learn/scikit-learn/issues/28256
    """
    X, y = make_classification(n_samples=50, n_features=10, random_state=0)
    estimator = EstimatorWithFitAndPredict().fit(X, y)

    multimetric_scorer = _MultimetricScorer(scorers={"acc": get_scorer("accuracy")})
    with config_context(enable_metadata_routing=enable_metadata_routing):
        multimetric_scorer(estimator, X, y)
