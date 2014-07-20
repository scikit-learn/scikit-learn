import pickle

import numpy as np

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_not_equal

from sklearn.metrics import (f1_score, r2_score, roc_auc_score, fbeta_score,
                             log_loss)
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.scorer import check_scoring
from sklearn.metrics import make_scorer, SCORERS
from sklearn.svm import LinearSVC
from sklearn.cluster import KMeans
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import make_blobs
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification
from sklearn.datasets import load_diabetes
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn.multiclass import OneVsRestClassifier


REGRESSION_SCORERS = ['r2', 'mean_absolute_error', 'mean_squared_error']
CLF_SCORERS = ['accuracy', 'f1', 'roc_auc', 'average_precision', 'precision',
               'recall', 'log_loss',
               'adjusted_rand_score'  # not really, but works
               ]


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

    # Multi-output multi-class predict_proba
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    y_proba = clf.predict_proba(X_test)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, np.vstack(p[:, -1] for p in y_proba).T)
    assert_almost_equal(score1, score2)

    # Multi-output multi-class decision_function
    # TODO Is there any yet?
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    clf._predict_proba = clf.predict_proba
    clf.predict_proba = None
    clf.decision_function = lambda X: [p[:, 1] for p in clf._predict_proba(X)]

    y_proba = clf.decision_function(X_test)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = roc_auc_score(y_test, np.vstack(p for p in y_proba).T)
    assert_almost_equal(score1, score2)

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


def test_scorer_sample_weight():
    """Test that scorers support sample_weight or raise sensible errors"""

    # Unlike the metrics invariance test, in the scorer case it's harder
    # to ensure that, on the classifier output, weighted and unweighted
    # scores really should be unequal.
    X, y = make_classification(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    sample_weight = np.ones_like(y_test)
    sample_weight[:10] = 0

    # get sensible estimators for each metric
    sensible_regr = DummyRegressor(strategy='median')
    sensible_regr.fit(X_train, y_train)
    sensible_clf = DecisionTreeClassifier()
    sensible_clf.fit(X_train, y_train)
    estimator = dict([(name, sensible_regr)
                      for name in REGRESSION_SCORERS] +
                     [(name, sensible_clf)
                      for name in CLF_SCORERS])

    for name, scorer in SCORERS.items():
        try:
            weighted = scorer(estimator[name], X_test, y_test,
                              sample_weight=sample_weight)
            ignored = scorer(estimator[name], X_test[10:], y_test[10:])
            unweighted = scorer(estimator[name], X_test, y_test)
            assert_not_equal(weighted, unweighted,
                             msg="scorer {0} behaves identically when "
                             "called with sample weights: {1} vs "
                             "{2}".format(name, weighted, unweighted))
            assert_almost_equal(weighted, ignored,
                         err_msg="scorer {0} behaves differently when "
                         "ignoring samples and setting sample_weight to 0: "
                         "{1} vs {2}".format(name, weighted, ignored))

        except TypeError as e:
            assert_true("sample_weight" in str(e),
                        "scorer {0} raises unhelpful exception when called "
                        "with sample weights: {1}".format(name, str(e)))
