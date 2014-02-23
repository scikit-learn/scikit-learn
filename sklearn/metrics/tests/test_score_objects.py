import pickle

import numpy as np

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import ignore_warnings

from sklearn.metrics import (accuracy_score, f1_score, r2_score, roc_auc_score,
                             fbeta_score, log_loss, mean_squared_error,
                             average_precision_score)
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.scorer import check_scoring, _evaluate_scorers
from sklearn.metrics import make_scorer, SCORERS
from sklearn.svm import LinearSVC, SVC
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.datasets import make_blobs
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification
from sklearn.datasets import load_diabetes
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn.multiclass import OneVsRestClassifier

# FIXME: temporary, to demonstrate ranking with several relevance levels.
def dcg_score(y_true, y_score, k=10, gains="exponential"):
    order = np.argsort(y_score)[::-1]
    y_true = np.take(y_true, order[:k])

    if gains == "exponential":
        gains = 2 ** y_true - 1
    elif gains == "linear":
        gains = y_true
    else:
        raise ValueError("Invalid gains option.")

    # highest rank is 1 so +2 instead of +1
    discounts = np.log2(np.arange(len(y_true)) + 2)
    return np.sum(gains / discounts)

dcg_scorer = make_scorer(dcg_score, needs_threshold=True)


class EstimatorWithoutFit(object):
    """Dummy estimator to test check_scoring"""
    pass


class EstimatorWithFit(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        return self


class EstimatorWithFitAndScore(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        return self

    def score(self, X, y):
        return 1.0


class EstimatorWithFitAndPredict(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        self.y = y
        return self

    def predict(self, X):
        return self.y


def test_check_scoring():
    """Test all branches of check_scoring"""
    estimator = EstimatorWithoutFit()
    assert_raises(TypeError, check_scoring, estimator)

    estimator = EstimatorWithFitAndScore()
    estimator.fit([[1]], [1])
    scorer = check_scoring(estimator)
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFitAndPredict()
    estimator.fit([[1]], [1])
    assert_raises(TypeError, check_scoring, estimator)

    scorer = check_scoring(estimator, "accuracy")
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFit()
    assert_raises(TypeError, check_scoring, estimator)


class EstimatorWithoutFit(object):
    """Dummy estimator to test check_scoring"""
    pass


class EstimatorWithFit(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        return self


class EstimatorWithFitAndScore(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        return self

    def score(self, X, y):
        return 1.0


class EstimatorWithFitAndPredict(object):
    """Dummy estimator to test check_scoring"""
    def fit(self, X, y):
        self.y = y
        return self

    def predict(self, X):
        return self.y


def test_check_scoring():
    """Test all branches of check_scoring"""
    estimator = EstimatorWithoutFit()
    pattern = (r"estimator should a be an estimator implementing 'fit' method,"
               r" .* was passed")
    assert_raises_regexp(TypeError, pattern, check_scoring, estimator)

    estimator = EstimatorWithFitAndScore()
    estimator.fit([[1]], [1])
    scorer = check_scoring(estimator)
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFitAndPredict()
    estimator.fit([[1]], [1])
    pattern = (r"If no scoring is specified, the estimator passed should have"
               r" a 'score' method\. The estimator .* does not\.")
    assert_raises_regexp(TypeError, pattern, check_scoring, estimator)

    scorer = check_scoring(estimator, "accuracy")
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFit()
    pattern = (r"The estimator passed should have a 'score'"
               r" or a 'predict' method\. The estimator .* does not\.")
    assert_raises_regexp(TypeError, pattern, check_scoring, estimator,
                         "accuracy")

    estimator = EstimatorWithFit()
    scorer = check_scoring(estimator, allow_none=True)
    assert_true(scorer is None)


def test_make_scorer():
    """Sanity check on the make_scorer factory function."""
    f = lambda *args: 0
    assert_raises(ValueError, make_scorer, f, needs_threshold=True,
                  needs_proba=True)


def test_classification_scores():
    """Test classification scorers."""
    X, y = make_blobs(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score1 = SCORERS['f1'](clf, X_test, y_test)
    score2 = f1_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)

    # test fbeta score that takes an argument
    scorer = make_scorer(fbeta_score, beta=2)
    score1 = scorer(clf, X_test, y_test)
    score2 = fbeta_score(y_test, clf.predict(X_test), beta=2)
    assert_almost_equal(score1, score2)

    # test that custom scorer can be pickled
    unpickled_scorer = pickle.loads(pickle.dumps(scorer))
    score3 = unpickled_scorer(clf, X_test, y_test)
    assert_almost_equal(score1, score3)

    # smoke test the repr:
    repr(fbeta_score)


def test_regression_scorers():
    """Test regression scorers."""
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = Ridge()
    clf.fit(X_train, y_train)
    score1 = SCORERS['r2'](clf, X_test, y_test)
    score2 = r2_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)


def test_thresholded_scorers():
    """Test scorers that take thresholds."""
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    score3 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)
    assert_almost_equal(score1, score3)

    logscore = SCORERS['log_loss'](clf, X_test, y_test)
    logloss = log_loss(y_test, clf.predict_proba(X_test))
    assert_almost_equal(-logscore, logloss)

    # same for an estimator without decision_function
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)

    # Test that an exception is raised on more than two classes
    X, y = make_blobs(random_state=0, centers=3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf.fit(X_train, y_train)
    assert_raises(ValueError, SCORERS['roc_auc'], clf, X_test, y_test)


def test_thresholded_scorers_multilabel_indicator_data():
    """Test that the scorer work with multilabel-indicator format
    for multilabel and multi-output multi-class classifier
    """
    X, y = make_multilabel_classification(return_indicator=True,
                                          allow_unlabeled=False,
                                          random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # Multilabel predict_proba
    clf = OneVsRestClassifier(DecisionTreeClassifier())
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test))
    assert_almost_equal(score1, score2)

    # Multilabel decision function
    clf = OneVsRestClassifier(LinearSVC(random_state=0))
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    assert_almost_equal(score1, score2)


def test_unsupervised_scorers():
    """Test clustering scorers against gold standard labeling."""
    # We don't have any real unsupervised Scorers yet.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    km = KMeans(n_clusters=3)
    km.fit(X_train)
    score1 = SCORERS['adjusted_rand_score'](km, X_test, y_test)
    score2 = adjusted_rand_score(y_test, km.predict(X_test))
    assert_almost_equal(score1, score2)


def test_evaluate_scorers_binary():
    X, y = make_classification(n_classes=2, random_state=0)

    # Test a classifier with decision_function.
    for clf in (SVC(), LinearSVC()):
        clf.fit(X, y)

        s1, s2 = _evaluate_scorers(clf, X, y, [SCORERS["f1"],
                                              SCORERS["roc_auc"]])
        df = clf.decision_function(X)
        y_pred = clf.predict(X)

        assert_almost_equal(s1, f1_score(y, y_pred))
        assert_almost_equal(s2, roc_auc_score(y, df))

    # Test a classifier with predict_proba.
    clf = LogisticRegression()
    clf.fit(X, y)

    s1, s2 = _evaluate_scorers(clf, X, y, [SCORERS["f1"],
                                          SCORERS["roc_auc"]])
    y_proba = clf.predict_proba(X)[:, 1]
    y_pred = clf.predict(X)

    assert_almost_equal(s1, f1_score(y, y_pred))
    assert_almost_equal(s2, roc_auc_score(y, y_proba))


def test_evaluate_scorers_multiclass():
    iris = load_iris()
    X, y = iris.data, iris.target

    # Test a classifier with decision_function.
    clf = LinearSVC()
    clf.fit(X, y)

    s1, s2 = _evaluate_scorers(clf, X, y, [SCORERS["f1"],
                                          SCORERS["accuracy"]])
    y_pred = clf.predict(X)

    assert_almost_equal(s1, f1_score(y, y_pred))
    assert_almost_equal(s2, accuracy_score(y, y_pred))

    # Test a classifier with predict_proba.
    clf = LogisticRegression()
    clf.fit(X, y)

    s1, s2, s3 = _evaluate_scorers(clf, X, y, [SCORERS["f1"],
                                              SCORERS["accuracy"],
                                              SCORERS["log_loss"]])
    y_proba = clf.predict_proba(X)
    y_pred = clf.predict(X)

    assert_almost_equal(s1, f1_score(y, y_pred))
    assert_almost_equal(s2, accuracy_score(y, y_pred))
    assert_almost_equal(s3, -log_loss(y, y_proba))


def test_evaluate_scorers_regression():
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target

    reg = Ridge()
    reg.fit(X, y)

    s1, s2 = _evaluate_scorers(reg, X, y, [SCORERS["r2"],
                                          SCORERS["mean_squared_error"]])
    y_pred = reg.predict(X)

    assert_almost_equal(s1, r2_score(y, y_pred))
    assert_almost_equal(s2, -mean_squared_error(y, y_pred))


def test_evaluate_scorers_ranking_by_regression():
    X, y = make_classification(n_classes=2, random_state=0)

    reg = DecisionTreeRegressor()
    reg.fit(X, y)

    s1, s2 = _evaluate_scorers(reg, X, y, [SCORERS["roc_auc"],
                                          SCORERS["average_precision"]])
    y_pred = reg.predict(X)

    assert_almost_equal(s1, roc_auc_score(y, y_pred))
    assert_almost_equal(s2, average_precision_score(y, y_pred))

    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target

    reg.fit(X, y)

    s1, s2 = _evaluate_scorers(reg, X, y, [SCORERS["r2"],
                                          dcg_scorer])
    y_pred = reg.predict(X)

    assert_almost_equal(s1, r2_score(y, y_pred))
    assert_almost_equal(s2, dcg_score(y, y_pred))


def test_evaluate_scorers_exceptions():
    clf = LinearSVC()
    # log_loss needs probabilities but LinearSVC does not have predict_proba.
    assert_raises(ValueError, _evaluate_scorers, clf, [], [],
                  [SCORERS["log_loss"]])


@ignore_warnings
def test_raises_on_score_list():
    """Test that when a list of scores is returned, we raise proper errors."""
    X, y = make_blobs(random_state=0)
    f1_scorer_no_average = make_scorer(f1_score, average=None)
    clf = DecisionTreeClassifier()
    assert_raises(ValueError, cross_val_score, clf, X, y,
                  scoring=f1_scorer_no_average)
    grid_search = GridSearchCV(clf, scoring=f1_scorer_no_average,
                               param_grid={'max_depth': [1, 2]})
    assert_raises(ValueError, grid_search.fit, X, y)
