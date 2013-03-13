import pickle

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true

from sklearn.metrics import f1_score, r2_score, auc_score, fbeta_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import SCORERS, Scorer, PRFScorer, WrapScorer, EstimatorScorer
from sklearn.svm import LinearSVC
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import make_blobs, load_diabetes
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.grid_search import GridSearchCV

# TODO: test scorers without ground truth


def test_classification_scores():
    X, y = make_blobs(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score1 = SCORERS['f1'](clf, X_test, y_test)
    score2 = f1_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)

    # test fbeta score that takes an argument
    scorer = Scorer(fbeta_score, beta=2)
    score1 = scorer(clf, X_test, y_test)
    score2 = fbeta_score(y_test, clf.predict(X_test), beta=2)
    assert_almost_equal(score1, score2)

    # test that custom scorer can be pickled
    unpickled_scorer = pickle.loads(pickle.dumps(scorer))
    score3 = unpickled_scorer(clf, X_test, y_test)
    assert_almost_equal(score1, score3)

    # smoke test the repr:
    repr(fbeta_score)


def test_regression_scores():
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = Ridge()
    clf.fit(X_train, y_train)
    score1 = SCORERS['r2'](clf, X_test, y_test)
    score2 = r2_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)


def test_thresholded_scores():
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = auc_score(y_test, clf.decision_function(X_test))
    score3 = auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)
    assert_almost_equal(score1, score3)

    # same for an estimator without decision_function
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    score1 = SCORERS['roc_auc'](clf, X_test, y_test)
    score2 = auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)

    # Test that an exception is raised on more than two classes
    X, y = make_blobs(random_state=0, centers=3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf.fit(X_train, y_train)
    assert_raises(ValueError, SCORERS['roc_auc'], clf, X_test, y_test)


def test_unsupervised_scores():
    # test clustering where there is some true y.
    # We don't have any real unsupervised SCORERS yet
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    km = KMeans(n_clusters=3)
    km.fit(X_train)
    score1 = SCORERS['ari'](km, X_test, y_test)
    score2 = adjusted_rand_score(y_test, km.predict(X_test))
    assert_almost_equal(score1, score2)


def test_raises_on_score_list():
    # test that when a list of scores is returned, we raise proper errors.
    X, y = make_blobs(random_state=0)
    f1_scorer_no_average = Scorer(f1_score, average=None)
    clf = DecisionTreeClassifier()
    assert_raises(ValueError, cross_val_score, clf, X, y,
                  scoring=f1_scorer_no_average)
    grid_search = GridSearchCV(clf, scoring=f1_scorer_no_average,
                               param_grid={'max_depth': [1, 2]})
    assert_raises(ValueError, grid_search.fit, X, y)


def test_calc_scores():
    """Test that the score returned by __call__ is named 'score' by calc_scores"""
    scorer = SCORERS['roc_auc']
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score = scorer(clf, X_test, y_test)
    scores = dict(scorer.calc_scores(clf, X_test, y_test))
    assert_true('score' in scores)
    assert_equal(score, scores['score'])


def test_prf_scorer():
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)

    f1_scorer = PRFScorer()
    f1_score = f1_scorer(clf, X_test, y_test)
    f1_scores = dict(f1_scorer.calc_scores(clf, X_test, y_test))

    f2_scorer = PRFScorer(beta=2.)
    f2_score = f2_scorer(clf, X_test, y_test)
    f2_scores = dict(f2_scorer.calc_scores(clf, X_test, y_test))

    def F(p, r, beta):
        return (1 + beta * beta) * p * r / (beta * beta * p + r)

    assert_equal(f1_score, f1_scores['score'])
    assert_equal(f2_score, f2_scores['score'])
    assert_almost_equal(f1_score, F(f1_scores['precision'], f1_scores['recall'], 1.))
    assert_almost_equal(f2_score, F(f2_scores['precision'], f2_scores['recall'], 2.))


def test_estimator_scorer():
    scorer = EstimatorScorer()
    X, y = make_blobs(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score = scorer(clf, X_test, y_test)
    assert_equal(clf.score(X_test, y_test), score)
    assert_equal(score, dict(scorer.calc_scores(clf, X_test, y_test))['score'])


def test_wrap_scorer():
    scorer = WrapScorer(lambda clf, X, y: clf.score(X, y) * 100)
    X, y = make_blobs(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score = scorer(clf, X_test, y_test)
    assert_equal(clf.score(X_test, y_test) * 100, score)
    assert_equal(score, dict(scorer.calc_scores(clf, X_test, y_test))['score'])
