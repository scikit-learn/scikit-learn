# Authors: Danny Sullivan <dbsullivan23@gmail.com>
#          Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#
# Licence: BSD 3 clause

import math
import numpy as np
import scipy.sparse as sp

from sklearn.linear_model.sag import get_auto_step_size
from sklearn.linear_model.sag_fast import _multinomial_grad_loss_all_samples
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.linear_model.base import make_dataset
from sklearn.linear_model.logistic import _multinomial_loss_grad

from sklearn.utils.extmath import logsumexp
from sklearn.utils.extmath import row_norms
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils import compute_class_weight
from sklearn.utils import check_random_state
from sklearn.preprocessing import LabelEncoder, LabelBinarizer
from sklearn.datasets import make_blobs, load_iris
from sklearn.base import clone

iris = load_iris()


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


def sag(X, y, step_size, alpha, n_iter=1, dloss=None, sparse=False,
        sample_weight=None, fit_intercept=True):
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
            if sample_weight is not None:
                gradient *= sample_weight[idx]
            update = entry * gradient + alpha * weights
            sum_gradient += update - gradient_memory[idx]
            gradient_memory[idx] = update

            if fit_intercept:
                intercept_sum_gradient += (gradient -
                                           intercept_gradient_memory[idx])
                intercept_gradient_memory[idx] = gradient
                intercept -= (step_size * intercept_sum_gradient
                              / len(seen) * decay)

            weights -= step_size * sum_gradient / len(seen)

    return weights, intercept


def sag_sparse(X, y, step_size, alpha, n_iter=1,
               dloss=None, sample_weight=None, sparse=False,
               fit_intercept=True):
    if step_size * alpha == 1.:
        raise ZeroDivisionError("Sparse sag does not handle the case "
                                "step_size * alpha == 1")
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

            if sample_weight is not None:
                gradient *= sample_weight[idx]

            update = entry * gradient
            sum_gradient += update - (gradient_memory[idx] * entry)

            if fit_intercept:
                intercept_sum_gradient += gradient - gradient_memory[idx]
                intercept -= (step_size * intercept_sum_gradient
                              / len(seen) * decay)

            gradient_memory[idx] = gradient

            wscale *= (1.0 - alpha * step_size)
            if counter == 0:
                c_sum[0] = step_size / (wscale * len(seen))
            else:
                c_sum[counter] = (c_sum[counter - 1] +
                                  step_size / (wscale * len(seen)))

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


def get_step_size(X, alpha, fit_intercept, classification=True):
    if classification:
        return (4.0 / (np.max(np.sum(X * X, axis=1))
                + fit_intercept + 4.0 * alpha))
    else:
        return 1.0 / (np.max(np.sum(X * X, axis=1)) + fit_intercept + alpha)


@ignore_warnings
def test_classifier_matching():
    n_samples = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)
    y[y == 0] = -1
    alpha = 1.1
    n_iter = 80
    fit_intercept = True
    step_size = get_step_size(X, alpha, fit_intercept)
    clf = LogisticRegression(solver="sag", fit_intercept=fit_intercept,
                             tol=1e-11, C=1. / alpha / n_samples,
                             max_iter=n_iter, random_state=10)
    clf.fit(X, y)

    weights, intercept = sag_sparse(X, y, step_size, alpha, n_iter=n_iter,
                                    dloss=log_dloss,
                                    fit_intercept=fit_intercept)
    weights2, intercept2 = sag(X, y, step_size, alpha, n_iter=n_iter,
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


@ignore_warnings
def test_regressor_matching():
    n_samples = 10
    n_features = 5

    rng = np.random.RandomState(10)
    X = rng.normal(size=(n_samples, n_features))
    true_w = rng.normal(size=n_features)
    y = X.dot(true_w)

    alpha = 1.
    n_iter = 100
    fit_intercept = True

    step_size = get_step_size(X, alpha, fit_intercept, classification=False)
    clf = Ridge(fit_intercept=fit_intercept, tol=.00000000001, solver='sag',
                alpha=alpha * n_samples, max_iter=n_iter)
    clf.fit(X, y)

    weights1, intercept1 = sag_sparse(X, y, step_size, alpha, n_iter=n_iter,
                                      dloss=squared_dloss,
                                      fit_intercept=fit_intercept)
    weights2, intercept2 = sag(X, y, step_size, alpha, n_iter=n_iter,
                               dloss=squared_dloss,
                               fit_intercept=fit_intercept)

    assert_array_almost_equal(weights1, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept1, clf.intercept_, decimal=10)
    assert_array_almost_equal(weights2, clf.coef_, decimal=10)
    assert_array_almost_equal(intercept2, clf.intercept_, decimal=10)


@ignore_warnings
def test_sag_pobj_matches_logistic_regression():
    """tests if the sag pobj matches log reg"""
    n_samples = 100
    alpha = 1.0
    max_iter = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)

    clf1 = LogisticRegression(solver='sag', fit_intercept=False, tol=.0000001,
                              C=1. / alpha / n_samples, max_iter=max_iter,
                              random_state=10)
    clf2 = clone(clf1)
    clf3 = LogisticRegression(fit_intercept=False, tol=.0000001,
                              C=1. / alpha / n_samples, max_iter=max_iter,
                              random_state=10)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, log_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, log_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, log_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj2, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj1, decimal=4)


@ignore_warnings
def test_sag_pobj_matches_ridge_regression():
    """tests if the sag pobj matches ridge reg"""
    n_samples = 100
    n_features = 10
    alpha = 1.0
    n_iter = 100
    fit_intercept = False
    rng = np.random.RandomState(10)
    X = rng.normal(size=(n_samples, n_features))
    true_w = rng.normal(size=n_features)
    y = X.dot(true_w)

    clf1 = Ridge(fit_intercept=fit_intercept, tol=.00000000001, solver='sag',
                 alpha=alpha, max_iter=n_iter, random_state=42)
    clf2 = clone(clf1)
    clf3 = Ridge(fit_intercept=fit_intercept, tol=.00001, solver='lsqr',
                 alpha=alpha, max_iter=n_iter, random_state=42)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, squared_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, squared_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, squared_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj1, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj2, decimal=4)


@ignore_warnings
def test_sag_regressor_computed_correctly():
    """tests if the sag regressor is computed correctly"""
    alpha = .1
    n_features = 10
    n_samples = 40
    max_iter = 50
    tol = .000001
    fit_intercept = True
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w) + 2.
    step_size = get_step_size(X, alpha, fit_intercept, classification=False)

    clf1 = Ridge(fit_intercept=fit_intercept, tol=tol, solver='sag',
                 alpha=alpha * n_samples, max_iter=max_iter)
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)

    spweights1, spintercept1 = sag_sparse(X, y, step_size, alpha,
                                          n_iter=max_iter,
                                          dloss=squared_dloss,
                                          fit_intercept=fit_intercept)

    spweights2, spintercept2 = sag_sparse(X, y, step_size, alpha,
                                          n_iter=max_iter,
                                          dloss=squared_dloss, sparse=True,
                                          fit_intercept=fit_intercept)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights1.ravel(),
                              decimal=3)
    assert_almost_equal(clf1.intercept_, spintercept1, decimal=1)

    # TODO: uncomment when sparse Ridge with intercept will be fixed (#4710)
    #assert_array_almost_equal(clf2.coef_.ravel(),
    #                          spweights2.ravel(),
    #                          decimal=3)
    #assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)'''


@ignore_warnings
def test_get_auto_step_size():
    X = np.array([[1, 2, 3], [2, 3, 4], [2, 3, 2]], dtype=np.float64)
    alpha = 1.2
    fit_intercept = False
    # sum the squares of the second sample because that's the largest
    max_squared_sum = 4 + 9 + 16
    max_squared_sum_ = row_norms(X, squared=True).max()
    assert_almost_equal(max_squared_sum, max_squared_sum_, decimal=4)

    for fit_intercept in (True, False):
        step_size_sqr = 1.0 / (max_squared_sum + alpha + int(fit_intercept))
        step_size_log = 4.0 / (max_squared_sum + 4.0 * alpha +
                               int(fit_intercept))

        step_size_sqr_ = get_auto_step_size(max_squared_sum_, alpha, "squared",
                                            fit_intercept)
        step_size_log_ = get_auto_step_size(max_squared_sum_, alpha, "log",
                                            fit_intercept)

        assert_almost_equal(step_size_sqr, step_size_sqr_, decimal=4)
        assert_almost_equal(step_size_log, step_size_log_, decimal=4)

    msg = 'Unknown loss function for SAG solver, got wrong instead of'
    assert_raise_message(ValueError, msg, get_auto_step_size,
                         max_squared_sum_, alpha, "wrong", fit_intercept)


@ignore_warnings
def test_sag_regressor():
    """tests if the sag regressor performs well"""
    xmin, xmax = -5, 5
    n_samples = 20
    tol = .001
    max_iter = 20
    alpha = 0.1
    rng = np.random.RandomState(0)
    X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

    # simple linear function without noise
    y = 0.5 * X.ravel()

    clf1 = Ridge(tol=tol, solver='sag', max_iter=max_iter,
                 alpha=alpha * n_samples)
    clf2 = clone(clf1)
    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)
    score1 = clf1.score(X, y)
    score2 = clf2.score(X, y)
    assert_greater(score1, 0.99)
    assert_greater(score2, 0.99)

    # simple linear function with noise
    y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

    clf1 = Ridge(tol=tol, solver='sag', max_iter=max_iter,
                 alpha=alpha * n_samples)
    clf2 = clone(clf1)
    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)
    score1 = clf1.score(X, y)
    score2 = clf2.score(X, y)
    score2 = clf2.score(X, y)
    assert_greater(score1, 0.5)
    assert_greater(score2, 0.5)


@ignore_warnings
def test_sag_classifier_computed_correctly():
    """tests if the binary classifier is computed correctly"""
    alpha = .1
    n_samples = 50
    n_iter = 50
    tol = .00001
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0,
                      cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1 = LogisticRegression(solver='sag', C=1. / alpha / n_samples,
                              max_iter=n_iter, tol=tol, random_state=77,
                              fit_intercept=fit_intercept)
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)

    spweights, spintercept = sag_sparse(X, y, step_size, alpha, n_iter=n_iter,
                                        dloss=log_dloss,
                                        fit_intercept=fit_intercept)
    spweights2, spintercept2 = sag_sparse(X, y, step_size, alpha,
                                          n_iter=n_iter,
                                          dloss=log_dloss, sparse=True,
                                          fit_intercept=fit_intercept)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


@ignore_warnings
def test_sag_multiclass_computed_correctly():
    """tests if the multiclass classifier is computed correctly"""
    alpha = .1
    n_samples = 20
    tol = .00001
    max_iter = 40
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)

    clf1 = LogisticRegression(solver='sag', C=1. / alpha / n_samples,
                              max_iter=max_iter, tol=tol, random_state=77,
                              fit_intercept=fit_intercept)
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)

    coef1 = []
    intercept1 = []
    coef2 = []
    intercept2 = []
    for cl in classes:
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights1, spintercept1 = sag_sparse(X, y_encoded, step_size, alpha,
                                              dloss=log_dloss, n_iter=max_iter,
                                              fit_intercept=fit_intercept)
        spweights2, spintercept2 = sag_sparse(X, y_encoded, step_size, alpha,
                                              dloss=log_dloss, n_iter=max_iter,
                                              sparse=True,
                                              fit_intercept=fit_intercept)
        coef1.append(spweights1)
        intercept1.append(spintercept1)

        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef1 = np.vstack(coef1)
    intercept1 = np.array(intercept1)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef1[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept1[i], decimal=1)

        assert_array_almost_equal(clf2.coef_[i].ravel(),
                                  coef2[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf2.intercept_[i], intercept2[i], decimal=1)


@ignore_warnings
def test_classifier_results():
    """tests if classifier results match target"""
    alpha = .1
    n_features = 20
    n_samples = 10
    tol = .01
    max_iter = 200
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)
    y = np.sign(y)
    clf1 = LogisticRegression(solver='sag', C=1. / alpha / n_samples,
                              max_iter=max_iter, tol=tol, random_state=77)
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)
    pred1 = clf1.predict(X)
    pred2 = clf2.predict(X)
    assert_almost_equal(pred1, y, decimal=12)
    assert_almost_equal(pred2, y, decimal=12)


@ignore_warnings
def test_binary_classifier_class_weight():
    """tests binary classifier with classweights for each class"""
    alpha = .1
    n_samples = 50
    n_iter = 20
    tol = .00001
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=10,
                      cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    class_weight = {1: .45, -1: .55}
    clf1 = LogisticRegression(solver='sag', C=1. / alpha / n_samples,
                              max_iter=n_iter, tol=tol, random_state=77,
                              fit_intercept=fit_intercept,
                              class_weight=class_weight)
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)

    le = LabelEncoder()
    class_weight_ = compute_class_weight(class_weight, np.unique(y), y)
    sample_weight = class_weight_[le.fit_transform(y)]
    spweights, spintercept = sag_sparse(X, y, step_size, alpha, n_iter=n_iter,
                                        dloss=log_dloss,
                                        sample_weight=sample_weight,
                                        fit_intercept=fit_intercept)
    spweights2, spintercept2 = sag_sparse(X, y, step_size, alpha,
                                          n_iter=n_iter,
                                          dloss=log_dloss, sparse=True,
                                          sample_weight=sample_weight,
                                          fit_intercept=fit_intercept)

    assert_array_almost_equal(clf1.coef_.ravel(),
                              spweights.ravel(),
                              decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(),
                              spweights2.ravel(),
                              decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


@ignore_warnings
def test_multiclass_classifier_class_weight():
    """tests multiclass with classweights for each class"""
    alpha = .1
    n_samples = 20
    tol = .00001
    max_iter = 50
    class_weight = {0: .45, 1: .55, 2: .75}
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0,
                      cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)

    clf1 = LogisticRegression(solver='sag', C=1. / alpha / n_samples,
                              max_iter=max_iter, tol=tol, random_state=77,
                              fit_intercept=fit_intercept,
                              class_weight=class_weight)
    clf2 = clone(clf1)
    clf1.fit(X, y)
    clf2.fit(sp.csr_matrix(X), y)

    le = LabelEncoder()
    class_weight_ = compute_class_weight(class_weight, np.unique(y), y)
    sample_weight = class_weight_[le.fit_transform(y)]

    coef1 = []
    intercept1 = []
    coef2 = []
    intercept2 = []
    for cl in classes:
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights1, spintercept1 = sag_sparse(X, y_encoded, step_size, alpha,
                                              n_iter=max_iter, dloss=log_dloss,
                                              sample_weight=sample_weight)
        spweights2, spintercept2 = sag_sparse(X, y_encoded, step_size, alpha,
                                              n_iter=max_iter, dloss=log_dloss,
                                              sample_weight=sample_weight,
                                              sparse=True)
        coef1.append(spweights1)
        intercept1.append(spintercept1)
        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef1 = np.vstack(coef1)
    intercept1 = np.array(intercept1)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_array_almost_equal(clf1.coef_[i].ravel(),
                                  coef1[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf1.intercept_[i], intercept1[i], decimal=1)

        assert_array_almost_equal(clf2.coef_[i].ravel(),
                                  coef2[i].ravel(),
                                  decimal=2)
        assert_almost_equal(clf2.intercept_[i], intercept2[i], decimal=1)


def test_classifier_single_class():
    """tests if ValueError is thrown with only one class"""
    X = [[1, 2], [3, 4]]
    y = [1, 1]

    assert_raise_message(ValueError,
                         "This solver needs samples of at least 2 classes "
                         "in the data",
                         LogisticRegression(solver='sag').fit,
                         X, y)


def test_step_size_alpha_error():
    X = [[0, 0], [0, 0]]
    y = [1, -1]
    fit_intercept = False
    alpha = 1.
    msg = ("Current sag implementation does not handle the case"
           " step_size * alpha_scaled == 1")

    clf1 = LogisticRegression(solver='sag', C=1. / alpha,
                              fit_intercept=fit_intercept)
    assert_raise_message(ZeroDivisionError, msg, clf1.fit, X, y)

    clf2 = Ridge(fit_intercept=fit_intercept, solver='sag', alpha=alpha)
    assert_raise_message(ZeroDivisionError, msg, clf2.fit, X, y)


def test_multinomial_loss():
    # test if the multinomial loss and gradient computations are consistent
    X, y = iris.data, iris.target.astype(np.float64)
    n_samples, n_features = X.shape
    n_classes = len(np.unique(y))

    rng = check_random_state(42)
    weights = rng.randn(n_features, n_classes)
    intercept = rng.randn(n_classes)
    sample_weights = rng.randn(n_samples)
    np.abs(sample_weights, sample_weights)

    # compute loss and gradient like in multinomial SAG
    dataset, _ = make_dataset(X, y, sample_weights, random_state=42)
    loss_1, grad_1 = _multinomial_grad_loss_all_samples(dataset, weights,
                                                        intercept, n_samples,
                                                        n_features, n_classes)
    # compute loss and gradient like in multinomial LogisticRegression
    lbin = LabelBinarizer()
    Y_bin = lbin.fit_transform(y)
    weights_intercept = np.vstack((weights, intercept)).T.ravel()
    loss_2, grad_2, _ = _multinomial_loss_grad(weights_intercept, X, Y_bin,
                                               0.0, sample_weights)
    grad_2 = grad_2.reshape(n_classes, -1)
    grad_2 = grad_2[:, :-1].T

    # comparison
    assert_array_almost_equal(grad_1, grad_2)
    assert_almost_equal(loss_1, loss_2)


def test_multinomial_loss_ground_truth():
    # n_samples, n_features, n_classes = 4, 2, 3
    n_classes = 3
    X = np.array([[1.1, 2.2], [2.2, -4.4], [3.3, -2.2], [1.1, 1.1]])
    y = np.array([0, 1, 2, 0])
    lbin = LabelBinarizer()
    Y_bin = lbin.fit_transform(y)

    weights = np.array([[0.1, 0.2, 0.3], [1.1, 1.2, -1.3]])
    intercept = np.array([1., 0, -.2])
    sample_weights = np.array([0.8, 1, 1, 0.8])

    prediction = np.dot(X, weights) + intercept
    logsumexp_prediction = logsumexp(prediction, axis=1)
    p = prediction - logsumexp_prediction[:, np.newaxis]
    loss_1 = -(sample_weights[:, np.newaxis] * p * Y_bin).sum()
    diff = sample_weights[:, np.newaxis] * (np.exp(p) - Y_bin)
    grad_1 = np.dot(X.T, diff)

    weights_intercept = np.vstack((weights, intercept)).T.ravel()
    loss_2, grad_2, _ = _multinomial_loss_grad(weights_intercept, X, Y_bin,
                                               0.0, sample_weights)
    grad_2 = grad_2.reshape(n_classes, -1)
    grad_2 = grad_2[:, :-1].T

    assert_almost_equal(loss_1, loss_2)
    assert_array_almost_equal(grad_1, grad_2)

    # ground truth
    loss_gt = 11.680360354325961
    grad_gt = np.array([[-0.557487, -1.619151, +2.176638],
                        [-0.903942, +5.258745, -4.354803]])
    assert_almost_equal(loss_1, loss_gt)
    assert_array_almost_equal(grad_1, grad_gt)
