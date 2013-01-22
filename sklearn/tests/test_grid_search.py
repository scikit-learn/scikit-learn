"""
Testing for grid search module (sklearn.grid_search)

"""

from cStringIO import StringIO
import sys

import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal

from sklearn.base import BaseEstimator
from sklearn.grid_search import GridSearchCV
from sklearn.datasets.samples_generator import make_classification, make_blobs
from sklearn.svm import LinearSVC, SVC
from sklearn.cluster import KMeans, MeanShift
from sklearn.metrics import f1_score, precision_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.cross_validation import KFold


class MockClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation"""
    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y):
        assert_true(len(X) == len(Y))
        return self

    def predict(self, T):
        return T.shape[0]

    def score(self, X=None, Y=None):
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score


class MockListClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation.

    Checks that GridSearchCV didn't convert X to array.
    """
    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y):
        assert_true(len(X) == len(Y))
        assert_true(isinstance(X, list))
        return self

    def predict(self, T):
        return T.shape[0]

    def score(self, X=None, Y=None):
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score


X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
y = np.array([1, 1, 2, 2])


def test_grid_search():
    """Test that the best estimator contains the right value for foo_param"""
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, verbose=3)
    # make sure it selects the smallest parameter in case of ties
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    grid_search.fit(X, y)
    sys.stdout = old_stdout
    assert_equal(grid_search.best_estimator_.foo_param, 2)

    for i, foo_i in enumerate([1, 2, 3]):
        assert_true(grid_search.grid_scores_[i][0] == {'foo_param': foo_i})
    # Smoke test the score:
    grid_search.score(X, y)


def test_no_refit():
    """Test that grid search can be used for model selection only"""
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, refit=False)
    grid_search.fit(X, y)
    assert_true(hasattr(grid_search, "best_params_"))


def test_grid_search_error():
    """Test that grid search will capture errors on data with different
    length"""
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, X_[:180], y_)


def test_grid_search_one_grid_point():
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)
    param_dict = {"C": [1.0], "kernel": ["rbf"], "gamma": [0.1]}

    clf = SVC()
    cv = GridSearchCV(clf, param_dict)
    cv.fit(X_, y_)

    clf = SVC(C=1.0, kernel="rbf", gamma=0.1)
    clf.fit(X_, y_)

    assert_array_equal(clf.dual_coef_, cv.best_estimator_.dual_coef_)


def test_grid_search_bad_param_grid():
    param_dict = {"C": 1.0}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)

    param_dict = {"C": []}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)

    param_dict = {"C": np.ones(6).reshape(3, 2)}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)


def test_grid_search_sparse():
    """Test that grid search works with both dense and sparse matrices"""
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(X_[:180], y_[:180])
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator_.C

    X_ = sp.csr_matrix(X_)
    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(X_[:180].tocoo(), y_[:180])
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator_.C

    assert_true(np.mean(y_pred == y_pred2) >= .9)
    assert_equal(C, C2)


def test_grid_search_sparse_score_func():
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, score_func=f1_score)
    cv.fit(X_[:180], y_[:180])
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator_.C

    X_ = sp.csr_matrix(X_)
    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, score_func=f1_score)
    cv.fit(X_[:180], y_[:180])
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator_.C

    assert_array_equal(y_pred, y_pred2)
    assert_equal(C, C2)
    # Smoke test the score
    #np.testing.assert_allclose(f1_score(cv.predict(X_[:180]), y[:180]),
    #                        cv.score(X_[:180], y[:180]))

    # test loss_func
    def f1_loss(y_true_, y_pred_):
        return -f1_score(y_true_, y_pred_)
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, loss_func=f1_loss)
    cv.fit(X_[:180], y_[:180])
    y_pred3 = cv.predict(X_[180:])
    C3 = cv.best_estimator_.C

    assert_array_equal(y_pred, y_pred3)
    assert_equal(C, C3)


def test_grid_search_precomputed_kernel():
    """Test that grid search works when the input features are given in the
    form of a precomputed kernel matrix """
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    # compute the training kernel matrix corresponding to the linear kernel
    K_train = np.dot(X_[:180], X_[:180].T)
    y_train = y_[:180]

    clf = SVC(kernel='precomputed')
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(K_train, y_train)

    assert_true(cv.best_score_ >= 0)

    # compute the test kernel matrix
    K_test = np.dot(X_[180:], X_[:180].T)
    y_test = y_[180:]

    y_pred = cv.predict(K_test)

    assert_true(np.mean(y_pred == y_test) >= 0)

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    assert_raises(ValueError, cv.fit, K_train.tolist(), y_train)


def test_grid_search_precomputed_kernel_error_nonsquare():
    """Test that grid search returns an error with a non-square precomputed
    training kernel matrix"""
    K_train = np.zeros((10, 20))
    y_train = np.ones((10, ))
    clf = SVC(kernel='precomputed')
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, K_train, y_train)


def test_grid_search_precomputed_kernel_error_kernel_function():
    """Test that grid search returns an error when using a kernel_function"""
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)
    kernel_function = lambda x1, x2: np.dot(x1, x2.T)
    clf = SVC(kernel=kernel_function)
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, X_, y_)


class BrokenClassifier(BaseEstimator):
    """Broken classifier that cannot be fit twice"""

    def __init__(self, parameter=None):
        self.parameter = parameter

    def fit(self, X, y):
        assert_true(not hasattr(self, 'has_been_fit_'))
        self.has_been_fit_ = True

    def predict(self, X):
        return np.zeros(X.shape[0])


def test_refit():
    """Regression test for bug in refitting

    Simulates re-fitting a broken estimator; this used to break with
    sparse SVMs.
    """
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = GridSearchCV(BrokenClassifier(), [{'parameter': [0, 1]}],
                       score_func=precision_score, refit=True)
    clf.fit(X, y)


def test_X_as_list():
    """Pass X as list in GridSearchCV
    """
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = MockListClassifier()
    cv = KFold(n=len(X), n_folds=3)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=cv)
    grid_search.fit(X.tolist(), y).score(X, y)


def test_unsupervised_grid_search():
    # test grid-search with unsupervised estimator
    X, y = make_blobs(random_state=0)
    km = KMeans(random_state=0)
    grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]),
                               score_func=adjusted_rand_score)
    grid_search.fit(X)
    # most number of clusters should be best
    assert_equal(grid_search.best_params_["n_clusters"], 4)


def test_bad_estimator():
    # test grid-search with unsupervised estimator
    ms = MeanShift()
    assert_raises(TypeError, GridSearchCV, ms,
                  param_grid=dict(gamma=[.1, 1, 10]),
                  score_func=adjusted_rand_score)
