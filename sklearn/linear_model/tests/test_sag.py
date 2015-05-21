import math

import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import SAGRegressor, SAGClassifier
from sklearn.linear_model import LogisticRegression, Ridge
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


def log_loss(p, y):
    return np.mean(np.log(1. + np.exp(-y * p)))


# this is used for sag regression
def squared_dloss(p, y):
    return p - y


def squared_loss(p, y):
    return np.mean(0.5 * (p - y) * (p - y))


# function for measuring the log loss
def get_pobj(w, alpha, myX, myy, loss):
    w = w.ravel()
    pred = np.dot(myX, w)
    p = loss(pred, myy)
    p += alpha * w.dot(w) / 2.
    return p


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


def sag(X, y, eta, alpha, n_iter=1, dloss=None, sparse=False,
        class_weight=None, fit_intercept=True):
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(X.shape[1])
    sum_gradient = np.zeros(X.shape[1])
    gradient_memory = np.zeros((n_samples, n_features))

    intercept = 0.0
    intercept_sum_gradient = 0.0
    intercept_gradient_memory = np.zeros(n_samples)

    rng = np.random.RandomState(77)
    decay = 1.0
    seen = set()

    # sparse data has a fixed decay of .01
    if sparse:
        decay = .01

    for epoch in range(n_iter):
        for k in range(n_samples):
            idx = int(rng.rand(1) * n_samples)
            # idx = k
            entry = X[idx]
            seen.add(idx)
            p = np.dot(entry, weights) + intercept
            gradient = dloss(p, y[idx])
            if class_weight:
                gradient *= class_weight[y[idx]]
            update = entry * gradient + alpha * weights
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            if fit_intercept:
                intercept_sum_gradient += (gradient -
                                           intercept_gradient_memory[idx])
                intercept_gradient_memory[idx] = gradient
                intercept -= eta * intercept_sum_gradient / len(seen) * decay

            weights -= eta * sum_gradient / len(seen)

    return weights, intercept


def sag_sparse(X, y, eta, alpha, n_iter=1,
               dloss=None, class_weight=None, sparse=False,
               fit_intercept=True):
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(n_features)
    sum_gradient = np.zeros(n_features)
    last_updated = np.zeros(n_features, dtype=np.int)
    gradient_memory = np.zeros(n_samples)
    rng = np.random.RandomState(77)
    intercept = 0.0
    intercept_sum_gradient = 0.0
    wscale = 1.0
    decay = 1.0
    seen = set()

    c_sum = np.zeros(n_iter * n_samples)

    # sparse data has a fixed decay of .01
    if sparse:
        decay = .01

    counter = 0
    for epoch in range(n_iter):
        for k in range(n_samples):
            # idx = k
            idx = int(rng.rand(1) * n_samples)
            entry = X[idx]
            seen.add(idx)

            if counter >= 1:
                for j in range(n_features):
                    if last_updated[j] == 0:
                        weights[j] -= c_sum[counter - 1] * sum_gradient[j]
                    else:
                        weights[j] -= ((c_sum[counter - 1] -
                                        c_sum[last_updated[j] - 1]) *
                                       sum_gradient[j])
                    last_updated[j] = counter

            p = (wscale * np.dot(entry, weights)) + intercept
            gradient = dloss(p, y[idx])

            if class_weight:
                gradient *= class_weight[y[idx]]

            update = entry * gradient
            sum_gradient += update - (gradient_memory[idx] * entry)

            if fit_intercept:
                intercept_sum_gradient += gradient - gradient_memory[idx]
                intercept -= eta * (intercept_sum_gradient / len(seen)) * decay

            gradient_memory[idx] = gradient

            wscale *= (1.0 - alpha * eta)
            if counter == 0:
                c_sum[0] = eta / (wscale * len(seen))
            else:
                c_sum[counter] = (c_sum[counter - 1] +
                                  eta / (wscale * len(seen)))

            if counter >= 1 and wscale < 1e-9:
                for j in range(n_features):
                    if last_updated[j] == 0:
                        weights[j] -= c_sum[counter] * sum_gradient[j]
                    else:
                        weights[j] -= ((c_sum[counter] -
                                        c_sum[last_updated[j] - 1]) *
                                       sum_gradient[j])
                    last_updated[j] = counter + 1
                c_sum[counter] = 0
                weights *= wscale
                wscale = 1.0

            counter += 1

    for j in range(n_features):
        if last_updated[j] == 0:
            weights[j] -= c_sum[counter - 1] * sum_gradient[j]
        else:
            weights[j] -= ((c_sum[counter - 1] -
                            c_sum[last_updated[j] - 1]) *
                           sum_gradient[j])
    weights *= wscale
    return weights, intercept


def get_eta(X, alpha, fit_intercept, classification=True):
    if classification:
        return (4.0 / (np.max(np.sum(X * X, axis=1))
                + fit_intercept + 4.0 * alpha))
    else:
        return 1.0 / (np.max(np.sum(X * X, axis=1)) + fit_intercept + alpha)


def test_classifier_matching():
    n_samples = 40
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)
    y[y == 0] = -1
    alpha = 1.1
    n_iter = 100
    fit_intercept = True
    eta = get_eta(X, alpha, fit_intercept)
    clf = SAGClassifier(fit_intercept=fit_intercept, tol=.00000000001,
                        alpha=alpha, max_iter=n_iter, random_state=10)
    clf.fit(X, y)

    weights, intercept = sag_sparse(X, y, eta, alpha, n_iter=n_iter,
                                    dloss=log_dloss,
                                    fit_intercept=fit_intercept)
    weights2, intercept2 = sag(X, y, eta, alpha, n_iter=n_iter,
                               dloss=log_dloss,
                               fit_intercept=fit_intercept)
    weights = np.atleast_2d(weights)
    intercept = np.atleast_1d(intercept)
    weights2 = np.atleast_2d(weights2)
    intercept2 = np.atleast_1d(intercept2)

    assert_array_almost_equal(weights, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept, clf.intercept_, decimal=10)
    assert_array_almost_equal(weights2, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept2, clf.intercept_, decimal=10)


def test_regressor_matching():
    n_samples = 30
    n_features = 10

    np.random.seed(10)
    X = np.random.random((n_samples, n_features))
    true_w = np.random.random(n_features)
    y = X.dot(true_w)

    alpha = 1.1
    n_iter = 100
    fit_intercept = True

    eta = get_eta(X, alpha, fit_intercept, classification=False)
    clf = SAGRegressor(fit_intercept=fit_intercept, tol=.00000000001,
                       alpha=alpha, max_iter=n_iter, random_state=10)
    clf.fit(X, y)

    weights, intercept = sag_sparse(X, y, eta, alpha, n_iter=n_iter,
                                    dloss=squared_dloss,
                                    fit_intercept=fit_intercept)
    weights2, intercept2 = sag(X, y, eta, alpha, n_iter=n_iter,
                               dloss=squared_dloss,
                               fit_intercept=fit_intercept)
    weights = np.atleast_2d(weights)
    intercept = np.atleast_1d(intercept)
    weights2 = np.atleast_2d(weights2)
    intercept2 = np.atleast_1d(intercept2)

    assert_array_almost_equal(weights, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept, clf.intercept_, decimal=10)
    assert_array_almost_equal(weights2, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept2, clf.intercept_, decimal=10)


def test_sag_pobj_matches_logistic_regression():
    """tests if the sag pobj matches log reg"""
    n_samples = 500
    alpha = 1.0
    n_iter = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    clf1 = SAGClassifier(fit_intercept=False, tol=.0000001, alpha=alpha,
                         max_iter=n_iter, random_state=10)
    clf2 = SparseSAGClassifier(fit_intercept=False, tol=.0000001, alpha=alpha,
                               max_iter=n_iter, random_state=10)
    clf3 = LogisticRegression(fit_intercept=False, tol=.0000001,
                              C=1.0 / (alpha * n_samples), max_iter=n_iter,
                              random_state=10)

    clf1.fit(X, y)
    clf2.fit(X, y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, log_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, log_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, log_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj2, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj1, decimal=4)


def test_sag_pobj_matches_ridge_regression():
    """tests if the sag pobj matches ridge reg"""
    n_samples = 500
    n_features = 10
    alpha = 1.0
    n_iter = 100
    np.random.seed(10)
    X = np.random.random((n_samples, n_features))
    true_w = np.random.random(n_features)
    y = X.dot(true_w)

    clf1 = SAGRegressor(fit_intercept=False, tol=.00001, alpha=alpha,
                        max_iter=n_iter, random_state=10)
    clf2 = SparseSAGRegressor(fit_intercept=False, tol=.00001, alpha=alpha,
                              max_iter=n_iter, random_state=10)
    clf3 = Ridge(fit_intercept=False, tol=.00001,
                 alpha=alpha * n_samples, max_iter=n_iter, solver="lsqr")

    clf1.fit(X, y)
    clf2.fit(X, y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, log_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, log_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, log_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj2, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj1, decimal=4)


def test_sag_regressor_computed_correctly():
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
    clf2 = SparseSAGRegressor(eta0=eta, alpha=alpha,
                              max_iter=max_iter, tol=tol, random_state=77)
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)

    clf1.fit(X, y)
    clf2.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)

    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=n_iter,
                                          dloss=squared_dloss, sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


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
    clf2 = SparseSAGRegressor(eta0=eta, alpha=alpha, warm_start=True,
                              max_iter=max_iter, tol=tol, random_state=77)
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)

    clf1.fit(X, y)
    clf1.fit(X, y)

    clf2.fit(X, y)
    clf2.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)
    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=n_iter,
                                          dloss=squared_dloss, sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept, decimal=1)


def test_auto_eta():
    """tests the auto eta computed correctly"""
    X = np.array([[1, 2, 3], [2, 3, 4], [2, 3, 2]])
    y = [.5, .6, .7]
    alpha = 1.0
    n_iter = 8
    sp_n_iter = 18
    tol = .01
    max_iter = 1000
    # sum the squares of the second sample because that's the largest
    eta = 4 + 9 + 16
    eta = 1.0 / (eta + alpha)

    clf1 = SAGRegressor(eta0='auto', alpha=alpha,
                        max_iter=max_iter, tol=tol, random_state=77)
    clf2 = SparseSAGRegressor(eta0='auto', alpha=alpha,
                              max_iter=max_iter, tol=.001,
                              random_state=77)
    clf1.fit(X, y)
    clf2.fit(X, y)
    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=squared_dloss)

    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=sp_n_iter,
                                          dloss=squared_dloss, sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


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
    clf2 = SparseSAGRegressor(eta0=.001, alpha=0.1, max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    clf2.fit(X, y)
    score = clf.score(X, y)
    score2 = clf2.score(X, y)
    assert_greater(score, 0.99)
    assert_greater(score2, 0.99)

    # simple linear function with noise
    y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

    clf = SAGRegressor(eta0=.001, alpha=0.1, max_iter=max_iter, tol=tol)
    clf2 = SparseSAGRegressor(eta0=.001, alpha=0.1, max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    clf2.fit(X, y)
    score = clf.score(X, y)
    score2 = clf2.score(X, y)
    score2 = clf2.score(X, y)
    assert_greater(score, 0.5)
    assert_greater(score2, 0.5)


def test_sag_classifier_computed_correctly():
    """tests the binary classifier computed correctly"""
    eta = .001
    alpha = .1
    n_samples = 50
    n_iter = 59
    tol = .01
    max_iter = 1000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77)
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha,
                               max_iter=max_iter, tol=tol, random_state=77)
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1.fit(X, y)
    clf2.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss)
    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=n_iter,
                                          dloss=log_dloss, sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=1)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=1)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


def test_sag_multiclass_computed_correctly():
    """tests the multiclass classifier is computed correctly"""
    eta = .001
    alpha = .1
    n_samples = 50
    tol = .01
    max_iter = 1000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77)
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha,
                               max_iter=max_iter, tol=tol, random_state=77)
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    y[y == 0] = -1
    classes = np.unique(y)
    itrs = [84, 53, 49]
    sp_itrs = [48, 55, 49]

    clf1.fit(X, y)
    clf2.fit(X, y)

    coef = []
    intercept = []
    coef2 = []
    intercept2 = []
    for cl, it, sp_it in zip(classes, itrs, sp_itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha,
                                            dloss=log_dloss,
                                            n_iter=it)
        spweights2, spintercept2 = sag_sparse(X, y_encoded, eta, alpha,
                                              dloss=log_dloss,
                                              n_iter=sp_it, sparse=True)
        coef.append(spweights)
        intercept.append(spintercept)

        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef = np.vstack(coef)
    intercept = np.array(intercept)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

        assert_array_almost_equal(clf2.coef_[i].ravel(),
                                  coef2[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf2.intercept_[i], intercept2[i], decimal=1)


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
    clf = SAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter, tol=tol)
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha,
                               max_iter=max_iter, tol=tol)
    clf.fit(X, y)
    clf2.fit(X, y)
    pred = clf.predict(X)
    pred2 = clf2.predict(X)
    assert_almost_equal(pred, y, decimal=12)
    assert_almost_equal(pred2, y, decimal=12)


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
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter,
                               tol=tol, random_state=77, warm_start=True)
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1.fit(X, y)
    clf1.fit(X, y)
    clf2.fit(X, y)
    clf2.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss)
    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=61,
                                          dloss=log_dloss, sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


def test_multiclass_classifier_warm_start():
    """tests multiclass classifier with a warm start"""
    eta = .001
    alpha = .1
    n_samples = 20
    tol = .01
    max_iter = 3000
    clf1 = SAGClassifier(eta0=eta, alpha=alpha,
                         max_iter=max_iter, tol=tol, random_state=77,
                         warm_start=True)
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha,
                               max_iter=max_iter, tol=tol, random_state=77,
                               warm_start=True)
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    classes = np.unique(y)

    clf1.fit(X, y)
    clf1.fit(X, y)
    clf2.fit(X, y)
    clf2.fit(X, y)

    itrs = [65, 63, 66]
    sp_itrs = [39, 63, 66]

    coef = []
    intercept = []
    coef2 = []
    intercept2 = []
    for cl, it, sp_it in zip(classes, itrs, sp_itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha,
                                            n_iter=it,
                                            dloss=log_dloss)
        spweights2, spintercept2 = sag_sparse(X, y_encoded, eta, alpha,
                                              n_iter=sp_it,
                                              dloss=log_dloss, sparse=True)
        coef.append(spweights)
        intercept.append(spintercept)
        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef = np.vstack(coef)
    intercept = np.array(intercept)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

        assert_array_almost_equal(clf2.coef_[i].ravel(),
                                  coef2[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf2.intercept_[i], intercept2[i], decimal=1)


def test_binary_classifier_class_weight():
    """tests binary classifier with classweights for each class"""
    eta = .001
    alpha = .1
    n_samples = 50
    n_iter = 9
    tol = .1
    max_iter = 1000
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
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter,
                               tol=tol, random_state=77, warm_start=True,
                               class_weight=class_weight)
    clf1.fit(X, y)
    clf2.fit(X, y)

    spweights, spintercept = sag_sparse(X, y, eta, alpha,
                                        n_iter=n_iter,
                                        dloss=log_dloss,
                                        class_weight=class_weight)
    spweights2, spintercept2 = sag_sparse(X, y, eta, alpha,
                                          n_iter=n_iter,
                                          dloss=log_dloss,
                                          class_weight=class_weight,
                                          sparse=True)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


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
    clf2 = SparseSAGClassifier(eta0=eta, alpha=alpha, max_iter=max_iter,
                               tol=tol, random_state=77, warm_start=True,
                               class_weight=class_weight)
    clf1.fit(X, y)
    clf2.fit(X, y)

    itrs = [9, 9, 10]

    coef = []
    intercept = []
    coef2 = []
    intercept2 = []
    for cl, it in zip(classes, itrs):
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1
        cl_weight = {-1: 1.0, 1: 1.0}
        cl_weight[1] = class_weight[cl]

        spweights, spintercept = sag_sparse(X, y_encoded, eta, alpha,
                                            n_iter=it, dloss=log_dloss,
                                            class_weight=cl_weight)
        spweights2, spintercept2 = sag_sparse(X, y_encoded, eta, alpha,
                                              n_iter=it, dloss=log_dloss,
                                              class_weight=cl_weight,
                                              sparse=True)
        coef.append(spweights)
        intercept.append(spintercept)
        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef = np.vstack(coef)
    intercept = np.array(intercept)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept[i], decimal=1)

        assert_array_almost_equal(clf2.coef_[i].ravel(),
                                  coef2[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf2.intercept_[i], intercept2[i], decimal=1)


def test_classifier_single_class():
    """tests value error thrown with only one class"""
    X = [[1, 2], [3, 4]]
    Y = [1, 1]

    assert_raises_regexp(ValueError,
                         "The number of class labels must be "
                         "greater than one.",
                         SAGClassifier().fit,
                         X, Y)
