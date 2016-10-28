import pickle
import tempfile
import shutil
import os
import numbers

import numpy as np

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_warns_message

from sklearn.base import BaseEstimator
from sklearn.metrics import (f1_score, r2_score, roc_auc_score, fbeta_score,
                             log_loss, precision_score, recall_score)
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.scorer import (check_scoring, _PredictScorer,
                                    _passthrough_scorer)
from sklearn.metrics import make_scorer, get_scorer, SCORERS
from sklearn.svm import LinearSVC
from sklearn.pipeline import make_pipeline
from sklearn.cluster import KMeans
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.datasets import make_blobs
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.multiclass import OneVsRestClassifier
from sklearn.externals import joblib


REGRESSION_SCORERS = ['r2', 'neg_mean_absolute_error',
                      'neg_mean_squared_error', 'neg_median_absolute_error',
                      'mean_absolute_error',
                      'mean_squared_error', 'median_absolute_error']

CLF_SCORERS = ['accuracy', 'f1', 'f1_weighted', 'f1_macro', 'f1_micro',
               'roc_auc', 'average_precision', 'precision',
               'precision_weighted', 'precision_macro', 'precision_micro',
               'recall', 'recall_weighted', 'recall_macro', 'recall_micro',
               'neg_log_loss', 'log_loss',
               'adjusted_rand_score'  # not really, but works
               ]

MULTILABEL_ONLY_SCORERS = ['precision_samples', 'recall_samples', 'f1_samples']


def _make_estimators(X_train, y_train, y_ml_train):
    # Make estimators that make sense to test various scoring methods
    sensible_regr = DummyRegressor(strategy='median')
    sensible_regr.fit(X_train, y_train)
    sensible_clf = DecisionTreeClassifier(random_state=0)
    sensible_clf.fit(X_train, y_train)
    sensible_ml_clf = DecisionTreeClassifier(random_state=0)
    sensible_ml_clf.fit(X_train, y_ml_train)
    return dict(
        [(name, sensible_regr) for name in REGRESSION_SCORERS] +
        [(name, sensible_clf) for name in CLF_SCORERS] +
        [(name, sensible_ml_clf) for name in MULTILABEL_ONLY_SCORERS]
    )


X_mm, y_mm, y_ml_mm = None, None, None
ESTIMATORS = None
TEMP_FOLDER = None


def setup_module():
    # Create some memory mapped data
    global X_mm, y_mm, y_ml_mm, TEMP_FOLDER, ESTIMATORS
    TEMP_FOLDER = tempfile.mkdtemp(prefix='sklearn_test_score_objects_')
    X, y = make_classification(n_samples=30, n_features=5, random_state=0)
    _, y_ml = make_multilabel_classification(n_samples=X.shape[0],
                                             random_state=0)
    filename = os.path.join(TEMP_FOLDER, 'test_data.pkl')
    joblib.dump((X, y, y_ml), filename)
    X_mm, y_mm, y_ml_mm = joblib.load(filename, mmap_mode='r')
    ESTIMATORS = _make_estimators(X_mm, y_mm, y_ml_mm)


def teardown_module():
    global X_mm, y_mm, y_ml_mm, TEMP_FOLDER, ESTIMATORS
    # GC closes the mmap file descriptors
    X_mm, y_mm, y_ml_mm, ESTIMATORS = None, None, None, None
    shutil.rmtree(TEMP_FOLDER)


class EstimatorWithoutFit(object):
    """Dummy estimator to test check_scoring"""
    pass


class EstimatorWithFit(BaseEstimator):
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


class DummyScorer(object):
    """Dummy scorer that always returns 1."""
    def __call__(self, est, X, y):
        return 1


def test_all_scorers_repr():
    # Test that all scorers have a working repr
    for name, scorer in SCORERS.items():
        repr(scorer)


def test_check_scoring():
    # Test all branches of check_scoring
    estimator = EstimatorWithoutFit()
    pattern = (r"estimator should be an estimator implementing 'fit' method,"
               r" .* was passed")
    assert_raises_regexp(TypeError, pattern, check_scoring, estimator)

    estimator = EstimatorWithFitAndScore()
    estimator.fit([[1]], [1])
    scorer = check_scoring(estimator)
    assert_true(scorer is _passthrough_scorer)
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFitAndPredict()
    estimator.fit([[1]], [1])
    pattern = (r"If no scoring is specified, the estimator passed should have"
               r" a 'score' method\. The estimator .* does not\.")
    assert_raises_regexp(TypeError, pattern, check_scoring, estimator)

    scorer = check_scoring(estimator, "accuracy")
    assert_almost_equal(scorer(estimator, [[1]], [1]), 1.0)

    estimator = EstimatorWithFit()
    scorer = check_scoring(estimator, "accuracy")
    assert_true(isinstance(scorer, _PredictScorer))

    estimator = EstimatorWithFit()
    scorer = check_scoring(estimator, allow_none=True)
    assert_true(scorer is None)


def test_check_scoring_gridsearchcv():
    # test that check_scoring works on GridSearchCV and pipeline.
    # slightly redundant non-regression test.

    grid = GridSearchCV(LinearSVC(), param_grid={'C': [.1, 1]})
    scorer = check_scoring(grid, "f1")
    assert_true(isinstance(scorer, _PredictScorer))

    pipe = make_pipeline(LinearSVC())
    scorer = check_scoring(pipe, "f1")
    assert_true(isinstance(scorer, _PredictScorer))

    # check that cross_val_score definitely calls the scorer
    # and doesn't make any assumptions about the estimator apart from having a
    # fit.
    scores = cross_val_score(EstimatorWithFit(), [[1], [2], [3]], [1, 0, 1],
                             scoring=DummyScorer())
    assert_array_equal(scores, 1)


def test_make_scorer():
    # Sanity check on the make_scorer factory function.
    f = lambda *args: 0
    assert_raises(ValueError, make_scorer, f, needs_threshold=True,
                  needs_proba=True)


def test_classification_scores():
    # Test classification scorers.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)

    for prefix, metric in [('f1', f1_score), ('precision', precision_score),
                           ('recall', recall_score)]:

        score1 = get_scorer('%s_weighted' % prefix)(clf, X_test, y_test)
        score2 = metric(y_test, clf.predict(X_test), pos_label=None,
                        average='weighted')
        assert_almost_equal(score1, score2)

        score1 = get_scorer('%s_macro' % prefix)(clf, X_test, y_test)
        score2 = metric(y_test, clf.predict(X_test), pos_label=None,
                        average='macro')
        assert_almost_equal(score1, score2)

        score1 = get_scorer('%s_micro' % prefix)(clf, X_test, y_test)
        score2 = metric(y_test, clf.predict(X_test), pos_label=None,
                        average='micro')
        assert_almost_equal(score1, score2)

        score1 = get_scorer('%s' % prefix)(clf, X_test, y_test)
        score2 = metric(y_test, clf.predict(X_test), pos_label=1)
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
    # Test regression scorers.
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = Ridge()
    clf.fit(X_train, y_train)
    score1 = get_scorer('r2')(clf, X_test, y_test)
    score2 = r2_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)


def test_thresholded_scorers():
    # Test scorers that take thresholds.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    score3 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)
    assert_almost_equal(score1, score3)

    logscore = get_scorer('neg_log_loss')(clf, X_test, y_test)
    logloss = log_loss(y_test, clf.predict_proba(X_test))
    assert_almost_equal(-logscore, logloss)

    # same for an estimator without decision_function
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)

    # test with a regressor (no decision_function)
    reg = DecisionTreeRegressor()
    reg.fit(X_train, y_train)
    score1 = get_scorer('roc_auc')(reg, X_test, y_test)
    score2 = roc_auc_score(y_test, reg.predict(X_test))
    assert_almost_equal(score1, score2)

    # Test that an exception is raised on more than two classes
    X, y = make_blobs(random_state=0, centers=3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf.fit(X_train, y_train)
    assert_raises(ValueError, get_scorer('roc_auc'), clf, X_test, y_test)


def test_thresholded_scorers_multilabel_indicator_data():
    # Test that the scorer work with multilabel-indicator format
    # for multilabel and multi-output multi-class classifier
    X, y = make_multilabel_classification(allow_unlabeled=False,
                                          random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    # Multi-output multi-class predict_proba
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    y_proba = clf.predict_proba(X_test)
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
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
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, np.vstack(p for p in y_proba).T)
    assert_almost_equal(score1, score2)

    # Multilabel predict_proba
    clf = OneVsRestClassifier(DecisionTreeClassifier())
    clf.fit(X_train, y_train)
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.predict_proba(X_test))
    assert_almost_equal(score1, score2)

    # Multilabel decision function
    clf = OneVsRestClassifier(LinearSVC(random_state=0))
    clf.fit(X_train, y_train)
    score1 = get_scorer('roc_auc')(clf, X_test, y_test)
    score2 = roc_auc_score(y_test, clf.decision_function(X_test))
    assert_almost_equal(score1, score2)


def test_unsupervised_scorers():
    # Test clustering scorers against gold standard labeling.
    # We don't have any real unsupervised Scorers yet.
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    km = KMeans(n_clusters=3)
    km.fit(X_train)
    score1 = get_scorer('adjusted_rand_score')(km, X_test, y_test)
    score2 = adjusted_rand_score(y_test, km.predict(X_test))
    assert_almost_equal(score1, score2)


@ignore_warnings
def test_raises_on_score_list():
    # Test that when a list of scores is returned, we raise proper errors.
    X, y = make_blobs(random_state=0)
    f1_scorer_no_average = make_scorer(f1_score, average=None)
    clf = DecisionTreeClassifier()
    assert_raises(ValueError, cross_val_score, clf, X, y,
                  scoring=f1_scorer_no_average)
    grid_search = GridSearchCV(clf, scoring=f1_scorer_no_average,
                               param_grid={'max_depth': [1, 2]})
    assert_raises(ValueError, grid_search.fit, X, y)


@ignore_warnings
def test_scorer_sample_weight():
    # Test that scorers support sample_weight or raise sensible errors

    # Unlike the metrics invariance test, in the scorer case it's harder
    # to ensure that, on the classifier output, weighted and unweighted
    # scores really should be unequal.
    X, y = make_classification(random_state=0)
    _, y_ml = make_multilabel_classification(n_samples=X.shape[0],
                                             random_state=0)
    split = train_test_split(X, y, y_ml, random_state=0)
    X_train, X_test, y_train, y_test, y_ml_train, y_ml_test = split

    sample_weight = np.ones_like(y_test)
    sample_weight[:10] = 0

    # get sensible estimators for each metric
    estimator = _make_estimators(X_train, y_train, y_ml_train)

    for name, scorer in SCORERS.items():
        if name in MULTILABEL_ONLY_SCORERS:
            target = y_ml_test
        else:
            target = y_test
        try:
            weighted = scorer(estimator[name], X_test, target,
                              sample_weight=sample_weight)
            ignored = scorer(estimator[name], X_test[10:], target[10:])
            unweighted = scorer(estimator[name], X_test, target)
            assert_not_equal(weighted, unweighted,
                             msg="scorer {0} behaves identically when "
                             "called with sample weights: {1} vs "
                             "{2}".format(name, weighted, unweighted))
            assert_almost_equal(weighted, ignored,
                                err_msg="scorer {0} behaves differently when "
                                "ignoring samples and setting sample_weight to"
                                " 0: {1} vs {2}".format(name, weighted,
                                                        ignored))

        except TypeError as e:
            assert_true("sample_weight" in str(e),
                        "scorer {0} raises unhelpful exception when called "
                        "with sample weights: {1}".format(name, str(e)))


@ignore_warnings  # UndefinedMetricWarning for P / R scores
def check_scorer_memmap(scorer_name):
    scorer, estimator = SCORERS[scorer_name], ESTIMATORS[scorer_name]
    if scorer_name in MULTILABEL_ONLY_SCORERS:
        score = scorer(estimator, X_mm, y_ml_mm)
    else:
        score = scorer(estimator, X_mm, y_mm)
    assert isinstance(score, numbers.Number), scorer_name


def test_scorer_memmap_input():
    # Non-regression test for #6147: some score functions would
    # return singleton memmap when computed on memmap data instead of scalar
    # float values.
    for name in SCORERS.keys():
        yield check_scorer_memmap, name


def test_deprecated_names():
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)

    for name in ('mean_absolute_error', 'mean_squared_error',
                 'median_absolute_error', 'log_loss'):
        warning_msg = "Scoring method %s was renamed to" % name
        for scorer in (get_scorer(name), SCORERS[name]):
            assert_warns_message(DeprecationWarning,
                                 warning_msg,
                                 scorer, clf, X, y)

        assert_warns_message(DeprecationWarning,
                             warning_msg,
                             cross_val_score, clf, X, y, scoring=name)


def test_scoring_is_not_metric():
    assert_raises_regexp(ValueError, 'make_scorer', check_scoring,
                         LogisticRegression(), f1_score)
    assert_raises_regexp(ValueError, 'make_scorer', check_scoring,
                         LogisticRegression(), roc_auc_score)
    assert_raises_regexp(ValueError, 'make_scorer', check_scoring,
                         Ridge(), r2_score)
    assert_raises_regexp(ValueError, 'make_scorer', check_scoring,
                         KMeans(), adjusted_rand_score)
