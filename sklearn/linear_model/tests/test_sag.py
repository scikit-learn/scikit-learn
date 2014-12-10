import math

import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import SAGRegressor, SAGClassifier
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises_regexp
from sklearn.datasets import make_blobs


# this is used for sag classification
def log_dloss(p, y):
    z = p * y
    # approximately equal and saves the computation of the log
    if z > 18.0:
        return math.exp(-z) * -y
    if z < -18.0:
        return -y
    return -y / (math.exp(z) + 1.0)


# this is used for sag regression
def squared_dloss(p, y):
    return p - y


class SparseSAGClassifier(SAGClassifier):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGClassifier.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.decision_function(self, X)

    def predict_proba(self, X):
        X = sp.csr_matrix(X)
        return SAGClassifier.predict_proba(self, X)


class SparseSAGRegressor(SAGRegressor):

    def fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.fit(self, X, y, *args, **kw)

    def partial_fit(self, X, y, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.partial_fit(self, X, y, *args, **kw)

    def decision_function(self, X, *args, **kw):
        X = sp.csr_matrix(X)
        return SAGRegressor.decision_function(self, X, *args, **kw)


def sag(X, y, eta, alpha, intercept_init=0.0,
        n_iter=1, dloss=None):
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(X.shape[1])
    sum_gradient = np.zeros(X.shape[1])
    gradient_memory = np.zeros((n_samples, n_features))
    intercept = intercept_init
    rng = np.random.RandomState(77)
    decay = 1.0
    seen = set()

    # sparse data has a fixed decay of .01
    # if (sp.issparse(X)):
    #     decay = .01

    for epoch in range(n_iter):
        for k in range(n_samples):
            # idx = int(rng.rand(1) * n_samples)
            idx = k
            entry = X[idx]
            seen.add(idx)
            p = np.dot(entry, weights) + intercept
            gradient = dloss(p, y[idx])
            update = entry * gradient + alpha * weights
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            intercept -= eta * gradient * decay
            weights -= eta * sum_gradient / len(seen)

    return weights, intercept


def sag_sparse(X, y, eta, alpha, intercept_init=0.0, n_iter=1,
               dloss=None, class_weight=None):
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(n_features)
    sum_gradient = np.zeros(n_features)
    last_updated = np.zeros(n_features, dtype=np.int)
    gradient_memory = np.zeros((n_samples, n_features))
    rng = np.random.RandomState(77)
    intercept = intercept_init
    wscale = 1.0
    decay = 1.0
    seen = set()

    c_sum = np.zeros(n_iter * n_samples + 1)

    # sparse data has a fixed decay of .01
    if (sp.issparse(X)):
        decay = .01

    counter = 0
    for epoch in range(n_iter):
        for k in range(n_samples):
            # idx = int(rng.rand(1) * n_samples)
            idx = k
            entry = X[idx]
            seen.add(idx)

            for j in range(n_features):
                weights[j] -= ((c_sum[counter] - c_sum[last_updated[j]]) *
                               sum_gradient[j])
                last_updated[j] = counter

            p = (wscale * np.dot(entry, weights)) + intercept
            gradient = dloss(p, y[idx])

            if class_weight:
                gradient *= class_weight[y[idx]]

            update = entry * gradient
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            wscale *= (1.0 - alpha * eta)
            c_sum[counter + 1] = c_sum[counter] + eta / (wscale * len(seen))

            intercept -= eta * gradient * decay
            counter += 1

    for k in range(n_features):
        weights[k] -= (c_sum[counter] -
                       c_sum[last_updated[k]]) * sum_gradient[k]
    weights *= wscale

    return weights, intercept


def test_sag_computed_correctly():
    """tests the if the sag regressor computed correctly"""
    eta = .001
    alpha = .1
    n_features = 20
    n_samples = 100
    max_iter = 1000
    tol = .001
    n_iter = 51
    clf1 = SAGRegressor(eta0=eta, alpha=alpha,
                        max_iter=max_iter, tol=tol, random_state=77)
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)

    clf1.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_regressor_warm_start():
    """tests the regressor warmstart"""
    eta = .001
    alpha = .1
    n_features = 20
    n_samples = 100
    n_iter = 53
    max_iter = 1000
    tol = .001
    clf1 = SAGRegressor(eta0=eta, alpha=alpha, warm_start=True,
                        max_iter=max_iter, tol=tol, random_state=77)
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)

    clf1.fit(X, y)
    clf1.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_auto_eta():
    """tests the auto eta computed correctly"""
    X = np.array([[1, 2, 3], [2, 3, 4], [2, 3, 2]])
    y = [.5, .6, .7]
    alpha = 1.0
    n_iter = 10
    tol = .01
    max_iter = 1000
    # sum the squares of the second sample because that's the largest
    eta = 4 + 9 + 16
    eta = 1.0 / (eta + alpha)

    clf1 = SAGRegressor(eta0='auto', alpha=alpha,
                        max_iter=max_iter, tol=tol, random_state=77)
    clf1.fit(X, y)
    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_sag_regressor():
    """tests the if the sag regressor performs well"""
    xmin, xmax = -5, 5
    n_samples = 100
    tol = .001
    max_iter = 20
    rng = np.random.RandomState(0)
    X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

    # simple linear function without noise
    y = 0.5 * X.ravel()

    clf = SAGRegressor(eta0=.001, alpha=0.1, max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    score = clf.score(X, y)
    assert_greater(score, 0.99)

    # simple linear function with noise
    y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

    clf = SAGRegressor(eta0=.001, alpha=0.1, max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    score = clf.score(X, y)
    assert_greater(score, 0.5)


def test_sag_computed_correctly():
    """tests the binary classifier computed correctly"""
    eta = .001
    alpha = .1
    n_samples = 50
    n_iter = 59
    tol = .01
    max_iter = 1000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77)
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=1)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_sag_multiclass_computed_correctly():
    """tests the multiclass classifier is computed correctly"""
    eta = .001
    alpha = .1
    n_samples = 50
    tol = .05
    max_iter = 1000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77)
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    classes = np.unique(y)
    itrs = [16, 12, 16]

    clf1.fit(X, y)

    coef = []
    intercept = []
    for cl, it in zip(classes, itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha,
                                            dloss=log_dloss,
                                            n_iter=it)
        coef.append(spweights)
        intercept.append(spintercept)

    coef = np.vstack(coef)
    intercept = np.array(intercept)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)


def test_classifier_results():
    """tests classifier results match target"""
    eta = .2
    alpha = .1
    n_features = 20
    n_samples = 10
    tol = .01
    max_iter = 2000
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)
    y = np.sign(y)
    clf = SAGClassifier(eta0=eta, alpha=alpha,  max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    pred = clf.predict(X)
    assert_almost_equal(pred, y, decimal=12)


def test_binary_classifier_warm_start():
    """tests binary classifier with a warm start"""
    eta = .001
    alpha = .1
    n_samples = 50
    n_iter = 59
    tol = .01
    max_iter = 2000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77,
                         warm_start=True)
    rng = np.random.RandomState(0)
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1.fit(X, y)
    clf1.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_multiclass_classifier_warm_start():
    """tests multiclass classifier with a warm start"""
    eta = .001
    alpha = .1
    n_samples = 20
    tol = .1
    max_iter = 3000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77,
                         warm_start=True)
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    classes = np.unique(y)

    clf1.fit(X, y)
    clf1.fit(X, y)

    itrs = [13, 11, 12]

    coef = []
    intercept = []
    for cl, it in zip(classes, itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha,
                                            n_iter=it,
                                            dloss=log_dloss)
        coef.append(spweights)
        intercept.append(spintercept)

    coef = np.vstack(coef)
    intercept = np.array(intercept)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)


def test_binary_classifier_class_weight():
    """tests binary classifier with classweights for each class"""
    eta = .001
    alpha = .1
    n_samples = 50
    n_iter = 9
    tol = .1
    max_iter = 1000
    rng = np.random.RandomState(0)
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)
    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    class_weight = {1: .45, -1: .55}

    clf1 = SAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter, tol=tol,
                         random_state=77, warm_start=True,
                         class_weight=class_weight)

    clf1.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss,
                                        class_weight=class_weight)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)


def test_multiclass_classifier_class_weight():
    """tests multiclass with classweights for each class"""

    eta = .001
    alpha = .1
    n_samples = 20
    tol = .1
    max_iter = 1000
    class_weight = {0: .45, 1: .55, 2: .75}

    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    classes = np.unique(y)

    clf1 = SAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter, tol=tol,
                         random_state=77, warm_start=True,
                         class_weight=class_weight)

    clf1.fit(X, y)

    itrs = [9, 9, 10]

    coef = []
    intercept = []
    for cl, it in zip(classes, itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1
        cl_weight = {-1: 1.0, 1: 1.0}
        cl_weight[1] = class_weight[cl]

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha, n_iter=it,
                                            dloss=log_dloss,
                                            class_weight=cl_weight)
        coef.append(spweights)
        intercept.append(spintercept)

    coef = np.vstack(coef)
    intercept = np.array(intercept)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)


def test_classifier_single_class():
    """tests value error thrown with only one class"""
    X = [[1, 2], [3, 4]]
    Y = [1, 1]

    assert_raises_regexp(ValueError,
                         "The number of class labels must be "
                         "greater than one.",
                         SAGClassifier().fit,
                         X, Y)
