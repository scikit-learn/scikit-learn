from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.metrics import f1_score, r2_score, auc_score, fbeta_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.score_objects import scorers, AsScorer
from sklearn.svm import LinearSVC
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import make_blobs, load_diabetes
from sklearn.cross_validation import train_test_split


def test_classification_scores():
    X, y = make_blobs(random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LinearSVC(random_state=0)
    clf.fit(X_train, y_train)
    score1 = scorers['f1'](clf, X_test, y_test)
    score2 = f1_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)

    # test fbeta score that takes an argument
    scorer = AsScorer(fbeta_score, beta=2)
    score1 = scorer(clf, X_test, y_test)
    score2 = fbeta_score(y_test, clf.predict(X_test), beta=2)
    assert_almost_equal(score1, score2)


def test_regression_scores():
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = Ridge()
    clf.fit(X_train, y_train)
    score1 = scorers['r2'](clf, X_test, y_test)
    score2 = r2_score(y_test, clf.predict(X_test))
    assert_almost_equal(score1, score2)


def test_thresholded_scores():
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf = LogisticRegression(random_state=0)
    clf.fit(X_train, y_train)
    score1 = scorers['roc_auc'](clf, X_test, y_test)
    score2 = auc_score(y_test, clf.decision_function(X_test))
    score3 = auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)
    assert_almost_equal(score1, score3)

    # same for an estimator without decision_function
    clf = DecisionTreeClassifier()
    clf.fit(X_train, y_train)
    score1 = scorers['roc_auc'](clf, X_test, y_test)
    score2 = auc_score(y_test, clf.predict_proba(X_test)[:, 1])
    assert_almost_equal(score1, score2)

    # Test that an exception is raised on more than two classes
    X, y = make_blobs(random_state=0, centers=3)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    clf.fit(X_train, y_train)
    assert_raises(ValueError, scorers['roc_auc'], clf, X_test, y_test)


def test_unsupervised_scores():
    # test clustering where there is some true y.
    # We don't have any real unsupervised scorers yet
    X, y = make_blobs(random_state=0, centers=2)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    km = KMeans(n_clusters=3)
    km.fit(X_train)
    score1 = scorers['ari'](km, X_test, y_test)
    score2 = adjusted_rand_score(y_test, km.predict(X_test))
    assert_almost_equal(score1, score2)
