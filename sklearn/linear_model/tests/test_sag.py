# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import math
import re
import warnings

import numpy as np
import pytest

from sklearn.base import clone
from sklearn.datasets import load_iris, make_blobs, make_classification
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.linear_model._sag import _build_sag_alias_tables, get_auto_step_size
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import check_random_state, compute_class_weight
from sklearn.utils._testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
)
from sklearn.utils.extmath import row_norms
from sklearn.utils.fixes import CSR_CONTAINERS

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
    return np.mean(np.log(1.0 + np.exp(-y * p)))


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
    p += alpha * w.dot(w) / 2.0
    return p


def sag(
    X,
    y,
    step_size,
    alpha,
    n_iter=1,
    dloss=None,
    sparse=False,
    sample_weight=None,
    fit_intercept=True,
    saga=False,
):
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
        decay = 0.01

    for epoch in range(n_iter):
        for k in range(n_samples):
            idx = int(rng.rand() * n_samples)
            # idx = k
            entry = X[idx]
            seen.add(idx)
            p = np.dot(entry, weights) + intercept
            gradient = dloss(p, y[idx])
            if sample_weight is not None:
                gradient *= sample_weight[idx]
            update = entry * gradient + alpha * weights
            gradient_correction = update - gradient_memory[idx]
            sum_gradient += gradient_correction
            gradient_memory[idx] = update
            if saga:
                weights -= gradient_correction * step_size * (1 - 1.0 / len(seen))

            if fit_intercept:
                gradient_correction = gradient - intercept_gradient_memory[idx]
                intercept_gradient_memory[idx] = gradient
                intercept_sum_gradient += gradient_correction
                gradient_correction *= step_size * (1.0 - 1.0 / len(seen))
                if saga:
                    intercept -= (
                        step_size * intercept_sum_gradient / len(seen) * decay
                    ) + gradient_correction
                else:
                    intercept -= step_size * intercept_sum_gradient / len(seen) * decay

            weights -= step_size * sum_gradient / len(seen)

    return weights, intercept


def sag_sparse(
    X,
    y,
    step_size,
    alpha,
    n_iter=1,
    dloss=None,
    sample_weight=None,
    sparse=False,
    fit_intercept=True,
    saga=False,
    random_state=0,
):
    if step_size * alpha == 1.0:
        raise ZeroDivisionError(
            "Sparse sag does not handle the case step_size * alpha == 1"
        )
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(n_features)
    sum_gradient = np.zeros(n_features)
    last_updated = np.zeros(n_features, dtype=int)
    gradient_memory = np.zeros(n_samples)
    rng = check_random_state(random_state)
    intercept = 0.0
    intercept_sum_gradient = 0.0
    wscale = 1.0
    decay = 1.0
    seen = set()

    c_sum = np.zeros(n_iter * n_samples)

    # sparse data has a fixed decay of .01
    if sparse:
        decay = 0.01

    counter = 0
    for epoch in range(n_iter):
        for k in range(n_samples):
            # idx = k
            idx = int(rng.rand() * n_samples)
            entry = X[idx]
            seen.add(idx)

            if counter >= 1:
                for j in range(n_features):
                    if last_updated[j] == 0:
                        weights[j] -= c_sum[counter - 1] * sum_gradient[j]
                    else:
                        weights[j] -= (
                            c_sum[counter - 1] - c_sum[last_updated[j] - 1]
                        ) * sum_gradient[j]
                    last_updated[j] = counter

            p = (wscale * np.dot(entry, weights)) + intercept
            gradient = dloss(p, y[idx])

            if sample_weight is not None:
                gradient *= sample_weight[idx]

            update = entry * gradient
            gradient_correction = update - (gradient_memory[idx] * entry)
            sum_gradient += gradient_correction
            if saga:
                for j in range(n_features):
                    weights[j] -= (
                        gradient_correction[j]
                        * step_size
                        * (1 - 1.0 / len(seen))
                        / wscale
                    )

            if fit_intercept:
                gradient_correction = gradient - gradient_memory[idx]
                intercept_sum_gradient += gradient_correction
                gradient_correction *= step_size * (1.0 - 1.0 / len(seen))
                if saga:
                    intercept -= (
                        step_size * intercept_sum_gradient / len(seen) * decay
                    ) + gradient_correction
                else:
                    intercept -= step_size * intercept_sum_gradient / len(seen) * decay

            gradient_memory[idx] = gradient

            wscale *= 1.0 - alpha * step_size
            if counter == 0:
                c_sum[0] = step_size / (wscale * len(seen))
            else:
                c_sum[counter] = c_sum[counter - 1] + step_size / (wscale * len(seen))

            if counter >= 1 and wscale < 1e-9:
                for j in range(n_features):
                    if last_updated[j] == 0:
                        weights[j] -= c_sum[counter] * sum_gradient[j]
                    else:
                        weights[j] -= (
                            c_sum[counter] - c_sum[last_updated[j] - 1]
                        ) * sum_gradient[j]
                    last_updated[j] = counter + 1
                c_sum[counter] = 0
                weights *= wscale
                wscale = 1.0

            counter += 1

    for j in range(n_features):
        if last_updated[j] == 0:
            weights[j] -= c_sum[counter - 1] * sum_gradient[j]
        else:
            weights[j] -= (
                c_sum[counter - 1] - c_sum[last_updated[j] - 1]
            ) * sum_gradient[j]
    weights *= wscale
    return weights, intercept


def get_step_size(X, alpha, fit_intercept, classification=True):
    if classification:
        return 4.0 / (np.max(np.sum(X * X, axis=1)) + fit_intercept + 4.0 * alpha)
    else:
        return 1.0 / (np.max(np.sum(X * X, axis=1)) + fit_intercept + alpha)


def test_classifier_matching():
    n_samples = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0, cluster_std=0.1)
    # y must be 0 or 1
    alpha = 1.1
    fit_intercept = True
    step_size = get_step_size(X, alpha, fit_intercept)
    for solver in ["sag", "saga"]:
        if solver == "sag":
            n_iter = 80
        else:
            # SAGA variance w.r.t. stream order is higher
            n_iter = 300
        clf = LogisticRegression(
            solver=solver,
            fit_intercept=fit_intercept,
            tol=1e-11,
            C=1.0 / alpha / n_samples,
            max_iter=n_iter,
            random_state=10,
        )
        clf.fit(X, y)

        weights, intercept = sag_sparse(
            X,
            2 * y - 1,  # y must be -1 or +1
            step_size,
            alpha,
            n_iter=n_iter,
            dloss=log_dloss,
            fit_intercept=fit_intercept,
            saga=solver == "saga",
        )
        weights2, intercept2 = sag(
            X,
            2 * y - 1,  # y must be -1 or +1
            step_size,
            alpha,
            n_iter=n_iter,
            dloss=log_dloss,
            fit_intercept=fit_intercept,
            saga=solver == "saga",
        )
        weights = np.atleast_2d(weights)
        intercept = np.atleast_1d(intercept)
        weights2 = np.atleast_2d(weights2)
        intercept2 = np.atleast_1d(intercept2)

        assert_array_almost_equal(weights, clf.coef_, decimal=9)
        assert_array_almost_equal(intercept, clf.intercept_, decimal=9)
        assert_array_almost_equal(weights2, clf.coef_, decimal=9)
        assert_array_almost_equal(intercept2, clf.intercept_, decimal=9)


def test_regressor_matching():
    n_samples = 10
    n_features = 5

    rng = np.random.RandomState(10)
    X = rng.normal(size=(n_samples, n_features))
    true_w = rng.normal(size=n_features)
    y = X.dot(true_w)

    alpha = 1.0
    n_iter = 100
    fit_intercept = True

    step_size = get_step_size(X, alpha, fit_intercept, classification=False)
    clf = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00000000001,
        solver="sag",
        alpha=alpha * n_samples,
        max_iter=n_iter,
    )
    clf.fit(X, y)

    weights1, intercept1 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
    )
    weights2, intercept2 = sag(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
    )

    assert_allclose(weights1, clf.coef_)
    assert_allclose(intercept1, clf.intercept_)
    assert_allclose(weights2, clf.coef_)
    assert_allclose(intercept2, clf.intercept_)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_pobj_matches_logistic_regression(csr_container):
    """tests if the sag pobj matches log reg"""
    n_samples = 100
    alpha = 1.0
    max_iter = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0, cluster_std=0.1)

    clf1 = LogisticRegression(
        solver="sag",
        fit_intercept=False,
        tol=0.0000001,
        C=1.0 / alpha / n_samples,
        max_iter=max_iter,
        random_state=10,
    )
    clf2 = clone(clf1)
    clf3 = LogisticRegression(
        fit_intercept=False,
        tol=0.0000001,
        C=1.0 / alpha / n_samples,
        max_iter=max_iter,
        random_state=10,
    )

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, log_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, log_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, log_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj2, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj1, decimal=4)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_pobj_matches_ridge_regression(csr_container):
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

    clf1 = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00000000001,
        solver="sag",
        alpha=alpha,
        max_iter=n_iter,
        random_state=42,
    )
    clf2 = clone(clf1)
    clf3 = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00001,
        solver="lsqr",
        alpha=alpha,
        max_iter=n_iter,
        random_state=42,
    )

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    clf3.fit(X, y)

    pobj1 = get_pobj(clf1.coef_, alpha, X, y, squared_loss)
    pobj2 = get_pobj(clf2.coef_, alpha, X, y, squared_loss)
    pobj3 = get_pobj(clf3.coef_, alpha, X, y, squared_loss)

    assert_array_almost_equal(pobj1, pobj2, decimal=4)
    assert_array_almost_equal(pobj1, pobj3, decimal=4)
    assert_array_almost_equal(pobj3, pobj2, decimal=4)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_regressor_computed_correctly(csr_container):
    """tests if the sag regressor is computed correctly"""
    alpha = 0.1
    n_features = 10
    n_samples = 40
    max_iter = 100
    tol = 0.000001
    fit_intercept = True
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w) + 2.0
    step_size = get_step_size(X, alpha, fit_intercept, classification=False)

    clf1 = Ridge(
        fit_intercept=fit_intercept,
        tol=tol,
        solver="sag",
        alpha=alpha * n_samples,
        max_iter=max_iter,
        random_state=rng,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)

    spweights1, spintercept1 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=max_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
        random_state=rng,
    )

    spweights2, spintercept2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=max_iter,
        dloss=squared_dloss,
        sparse=True,
        fit_intercept=fit_intercept,
        random_state=rng,
    )

    assert_array_almost_equal(clf1.coef_.ravel(), spweights1.ravel(), decimal=3)
    assert_almost_equal(clf1.intercept_, spintercept1, decimal=1)

    # TODO: uncomment when sparse Ridge with intercept will be fixed (#4710)
    # assert_array_almost_equal(clf2.coef_.ravel(),
    #                          spweights2.ravel(),
    #                          decimal=3)
    # assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)'''


def test_get_auto_step_size():
    X = np.array([[1, 2, 3], [2, 3, 4], [2, 3, 2]], dtype=np.float64)
    alpha = 1.2
    fit_intercept = False
    # sum the squares of the second sample because that's the largest
    max_squared_sum = 4 + 9 + 16
    max_squared_sum_ = row_norms(X, squared=True).max()
    n_samples = X.shape[0]
    assert_almost_equal(max_squared_sum, max_squared_sum_, decimal=4)

    for saga in [True, False]:
        for fit_intercept in (True, False):
            if saga:
                L_sqr = max_squared_sum + alpha + int(fit_intercept)
                L_log = (max_squared_sum + 4.0 * alpha + int(fit_intercept)) / 4.0
                mun_sqr = min(2 * n_samples * alpha, L_sqr)
                mun_log = min(2 * n_samples * alpha, L_log)
                step_size_sqr = 1 / (2 * L_sqr + mun_sqr)
                step_size_log = 1 / (2 * L_log + mun_log)
            else:
                step_size_sqr = 1.0 / (max_squared_sum + alpha + int(fit_intercept))
                step_size_log = 4.0 / (
                    max_squared_sum + 4.0 * alpha + int(fit_intercept)
                )

            step_size_sqr_ = get_auto_step_size(
                max_squared_sum_,
                alpha,
                "squared",
                fit_intercept,
                n_samples=n_samples,
                is_saga=saga,
            )
            step_size_log_ = get_auto_step_size(
                max_squared_sum_,
                alpha,
                "log",
                fit_intercept,
                n_samples=n_samples,
                is_saga=saga,
            )

            assert_almost_equal(step_size_sqr, step_size_sqr_, decimal=4)
            assert_almost_equal(step_size_log, step_size_log_, decimal=4)

    msg = "Unknown loss function for SAG solver, got wrong instead of"
    with pytest.raises(ValueError, match=msg):
        get_auto_step_size(max_squared_sum_, alpha, "wrong", fit_intercept)


@pytest.mark.parametrize("seed", range(3))  # locally tested with 1000 seeds
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_regressor(seed, csr_container):
    """tests if the sag regressor performs well"""
    xmin, xmax = -5, 5
    n_samples = 300
    tol = 0.001
    max_iter = 100
    alpha = 0.1
    rng = np.random.RandomState(seed)
    X = np.linspace(xmin, xmax, n_samples).reshape(n_samples, 1)

    # simple linear function without noise
    y = 0.5 * X.ravel()

    clf1 = Ridge(
        tol=tol,
        solver="sag",
        max_iter=max_iter,
        alpha=alpha * n_samples,
        random_state=rng,
    )
    clf2 = clone(clf1)
    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    score1 = clf1.score(X, y)
    score2 = clf2.score(X, y)
    assert score1 > 0.98
    assert score2 > 0.98

    # simple linear function with noise
    y = 0.5 * X.ravel() + rng.randn(n_samples, 1).ravel()

    clf1 = Ridge(
        tol=tol,
        solver="sag",
        max_iter=max_iter,
        alpha=alpha * n_samples,
        random_state=rng,
    )
    clf2 = clone(clf1)
    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    score1 = clf1.score(X, y)
    score2 = clf2.score(X, y)
    assert score1 > 0.45
    assert score2 > 0.45


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_classifier_computed_correctly(csr_container):
    """tests if the binary classifier is computed correctly"""
    alpha = 0.1
    n_samples = 50
    n_iter = 50
    tol = 0.00001
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0, cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    clf1 = LogisticRegression(
        solver="sag",
        C=1.0 / alpha / n_samples,
        max_iter=n_iter,
        tol=tol,
        random_state=77,
        fit_intercept=fit_intercept,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)

    spweights, spintercept = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=log_dloss,
        fit_intercept=fit_intercept,
    )
    spweights2, spintercept2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=log_dloss,
        sparse=True,
        fit_intercept=fit_intercept,
    )

    assert_array_almost_equal(clf1.coef_.ravel(), spweights.ravel(), decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(), spweights2.ravel(), decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sag_multiclass_computed_correctly(csr_container):
    """tests if the multiclass classifier is computed correctly"""
    alpha = 0.1
    n_samples = 20
    tol = 1e-5
    max_iter = 70
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=3, random_state=0, cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)

    clf1 = OneVsRestClassifier(
        LogisticRegression(
            solver="sag",
            C=1.0 / alpha / n_samples,
            max_iter=max_iter,
            tol=tol,
            random_state=77,
            fit_intercept=fit_intercept,
        )
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)

    coef1 = []
    intercept1 = []
    coef2 = []
    intercept2 = []
    for cl in classes:
        y_encoded = np.ones(n_samples)
        y_encoded[y != cl] = -1

        spweights1, spintercept1 = sag_sparse(
            X,
            y_encoded,
            step_size,
            alpha,
            dloss=log_dloss,
            n_iter=max_iter,
            fit_intercept=fit_intercept,
        )
        spweights2, spintercept2 = sag_sparse(
            X,
            y_encoded,
            step_size,
            alpha,
            dloss=log_dloss,
            n_iter=max_iter,
            sparse=True,
            fit_intercept=fit_intercept,
        )
        coef1.append(spweights1)
        intercept1.append(spintercept1)

        coef2.append(spweights2)
        intercept2.append(spintercept2)

    coef1 = np.vstack(coef1)
    intercept1 = np.array(intercept1)
    coef2 = np.vstack(coef2)
    intercept2 = np.array(intercept2)

    for i, cl in enumerate(classes):
        assert_allclose(clf1.estimators_[i].coef_.ravel(), coef1[i], rtol=1e-2)
        assert_allclose(clf1.estimators_[i].intercept_, intercept1[i], rtol=1e-1)

        assert_allclose(clf2.estimators_[i].coef_.ravel(), coef2[i], rtol=1e-2)
        # Note the very crude accuracy, i.e. high rtol.
        assert_allclose(clf2.estimators_[i].intercept_, intercept2[i], rtol=5e-1)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_classifier_results(csr_container):
    """tests if classifier results match target"""
    alpha = 0.1
    n_features = 20
    n_samples = 10
    tol = 0.01
    max_iter = 200
    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    w = rng.normal(size=n_features)
    y = np.dot(X, w)
    y = np.sign(y)
    clf1 = LogisticRegression(
        solver="sag",
        C=1.0 / alpha / n_samples,
        max_iter=max_iter,
        tol=tol,
        random_state=77,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)
    pred1 = clf1.predict(X)
    pred2 = clf2.predict(X)
    assert_almost_equal(pred1, y, decimal=12)
    assert_almost_equal(pred2, y, decimal=12)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_binary_classifier_class_weight(csr_container):
    """tests binary classifier with classweights for each class"""
    alpha = 0.1
    n_samples = 50
    n_iter = 20
    tol = 0.00001
    fit_intercept = True
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=10, cluster_std=0.1)
    step_size = get_step_size(X, alpha, fit_intercept, classification=True)
    classes = np.unique(y)
    y_tmp = np.ones(n_samples)
    y_tmp[y != classes[1]] = -1
    y = y_tmp

    class_weight = {1: 0.45, -1: 0.55}
    clf1 = LogisticRegression(
        solver="sag",
        C=1.0 / alpha / n_samples,
        max_iter=n_iter,
        tol=tol,
        random_state=77,
        fit_intercept=fit_intercept,
        class_weight=class_weight,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)

    le = LabelEncoder()
    class_weight_ = compute_class_weight(class_weight, classes=np.unique(y), y=y)
    sample_weight = class_weight_[le.fit_transform(y)]
    spweights, spintercept = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=log_dloss,
        sample_weight=sample_weight,
        fit_intercept=fit_intercept,
    )
    spweights2, spintercept2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        n_iter=n_iter,
        dloss=log_dloss,
        sparse=True,
        sample_weight=sample_weight,
        fit_intercept=fit_intercept,
    )

    assert_array_almost_equal(clf1.coef_.ravel(), spweights.ravel(), decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(), spweights2.ravel(), decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


def test_classifier_single_class():
    """tests if ValueError is thrown with only one class"""
    X = [[1, 2], [3, 4]]
    y = [1, 1]

    msg = "This solver needs samples of at least 2 classes in the data"
    with pytest.raises(ValueError, match=msg):
        LogisticRegression(solver="sag").fit(X, y)


def test_step_size_alpha_error():
    X = [[0, 0], [0, 0]]
    y = [1, -1]
    fit_intercept = False
    alpha = 1.0
    msg = re.escape(
        "Current sag implementation does not handle the case"
        " step_size * alpha_scaled == 1"
    )

    clf1 = LogisticRegression(solver="sag", C=1.0 / alpha, fit_intercept=fit_intercept)
    with pytest.raises(ZeroDivisionError, match=msg):
        clf1.fit(X, y)

    clf2 = Ridge(fit_intercept=fit_intercept, solver="sag", alpha=alpha)
    with pytest.raises(ZeroDivisionError, match=msg):
        clf2.fit(X, y)


@pytest.mark.parametrize("solver", ["sag", "saga"])
def test_sag_classifier_raises_error(solver):
    # Following #13316, the error handling behavior changed in cython sag. This
    # is simply a non-regression test to make sure numerical errors are
    # properly raised.

    # Train a classifier on a simple problem
    rng = np.random.RandomState(42)
    X, y = make_classification(random_state=rng)
    clf = LogisticRegression(solver=solver, random_state=rng, warm_start=True)
    clf.fit(X, y)

    # Trigger a numerical error by:
    # - corrupting the fitted coefficients of the classifier
    # - fit it again starting from its current state thanks to warm_start
    clf.coef_[:] = np.nan

    with pytest.raises(ValueError, match="Floating-point under-/overflow"):
        clf.fit(X, y)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_build_sag_alias_tables(dtype):
    """The alias tables must draw indices proportionally to sample_weight
    and provide a matching inverse-probability importance correction.

    Non-regression test for gh-21305: the SAG/SAGA solver used to draw
    samples uniformly regardless of sample_weight.
    """
    rng = np.random.RandomState(0)
    n_samples = 500
    sample_weight = rng.uniform(low=0.1, high=100, size=n_samples)

    alias_index, prob_threshold, importance_weight = _build_sag_alias_tables(
        sample_weight, dtype
    )

    assert alias_index.dtype == np.int32
    assert prob_threshold.dtype == dtype
    assert importance_weight.dtype == dtype

    # Empirically draw many samples using the alias method and check the
    # resulting frequencies match sample_weight / sample_weight.sum().
    draws = rng.randint(0, n_samples, size=200_000)
    u = rng.uniform(size=200_000)
    picked = np.where(u < prob_threshold[draws], draws, alias_index[draws])
    empirical_freq = np.bincount(picked, minlength=n_samples) / picked.shape[0]
    expected_freq = sample_weight / sample_weight.sum()
    assert_allclose(empirical_freq, expected_freq, atol=5e-3)

    # The importance weight must exactly cancel sample_weight on average,
    # i.e. E_i[sample_weight[i] * importance_weight[i]] == mean(sample_weight),
    # which is what keeps the SAG/SAGA gradient estimator unbiased.
    expected_importance = sample_weight.sum() / (n_samples * sample_weight)
    assert_allclose(importance_weight, expected_importance, rtol=1e-5)


def test_sag_sample_weight_convergence():
    """LogisticRegression(solver='saga') with a very unbalanced
    sample_weight must converge to (essentially) the same predictions as
    fitting on the expanded (duplicated) dataset, without raising
    ConvergenceWarning.

    Non-regression test for gh-21305: SAGA used to draw samples uniformly
    at random regardless of sample_weight, which is fine when weights are
    close to uniform but makes the gradient estimate very high-variance
    (and convergence correspondingly poor or absent) when a few samples
    carry most of the weight.

    Note this is specific to L2 (or no) penalty: plain SAG has no
    mechanism to stay unbiased when a sample is drawn rarely (see
    _sag.py:sag_solver's use_weighted_sampling), so skewing the draw
    towards heavy samples can slow it down instead of helping; and L1's
    proximal-operator thresholding interacts with highly skewed
    sample_weight in a way that can still need many more epochs than the
    default max_iter to satisfy the strict epoch-level tol criterion (it
    reaches a much better solution than before this fix in the same
    number of epochs, just not always within tol). Both are distinct,
    pre-existing rough edges this fix does not address.
    """
    rng = np.random.RandomState(0)
    n_unique, n_features = 2000, 10

    X_unique = rng.randn(n_unique, n_features)
    coef_true = rng.randn(n_features)
    y_unique = (X_unique @ coef_true + 0.1 * rng.randn(n_unique) > 0).astype(int)

    weights = rng.randint(1, 6, size=n_unique)
    # A subset of rows with a much larger multiplicity than the rest,
    # similar in spirit to the reporter's data (min=1, max=37681, mean=3.61
    # over 39341 rows in gh-21305), while keeping the weight mass spread
    # over enough rows to stay well clear of the near-separable-data /
    # unregularized-coefficient-growth regime, which is a distinct,
    # unrelated slow-convergence phenomenon.
    heavy_idx = rng.choice(n_unique, size=20, replace=False)
    weights[heavy_idx] = rng.randint(50, 400, size=20)

    X_dup = np.repeat(X_unique, weights, axis=0)
    y_dup = np.repeat(y_unique, weights)

    # A moderate, default-scale C: the "duplicated rows" and "unique rows +
    # sample_weight" formulations aren't *quite* the same optimization
    # problem regardless of C (the penalty is scaled by 1 / n_samples, and
    # n_samples differs between the two representations), but an extremely
    # large C makes both formulations approach the unregularized MLE, which
    # on randomly generated (imperfectly separable, but not by a lot)
    # data can need many more epochs to satisfy the tol-based stopping
    # criterion regardless of sampling -- an unrelated slow-convergence
    # phenomenon this fix does not address. See the "extreme weights"
    # scenario in this same PR's manual testing for a case where this
    # dominates even more starkly.
    common_params = dict(
        solver="saga",
        penalty="l2",
        C=1.0,
        max_iter=2000,
        tol=1e-4,
        random_state=0,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        clf_dup = LogisticRegression(**common_params).fit(X_dup, y_dup)
        clf_weighted = LogisticRegression(**common_params).fit(
            X_unique, y_unique, sample_weight=weights.astype(float)
        )

    # Raw coef_ can wobble more than this at the margins (e.g. near-zero /
    # L1-zeroed coordinates) since SAG/SAGA only guarantee a finite-tol
    # stochastic solution, not an exact one; predicted-probability
    # agreement is the more meaningful invariant here.
    proba_dup = clf_dup.predict_proba(X_unique)[:, 1]
    proba_weighted = clf_weighted.predict_proba(X_unique)[:, 1]
    assert_allclose(proba_dup, proba_weighted, atol=0.05)


def test_sag_l1_sample_weight_quality_improves():
    """With L1 penalty and highly-skewed sample_weight, SAGA may still need
    more than the default max_iter to satisfy the strict epoch-level tol
    criterion (a separate rough edge in how the proximal operator
    interacts with skewed weights, see test_sag_sample_weight_convergence's
    docstring), but drawing samples proportionally to their weight must
    substantially reduce the (weighted) log loss reached within a fixed,
    identical epoch budget, compared to the previous ("draw uniformly,
    ignore sample_weight for sampling") behavior.

    Non-regression test for gh-21305.
    """
    # A fixed seed (rather than global_random_seed): the size of the
    # improvement is sensitive to how pathological this particular draw of
    # sample_weight happens to be, so this uses one hand-verified to show a
    # large (>20x), robustly-reproducible gap rather than parametrizing
    # over arbitrary seeds where a milder weight skew wouldn't reliably do so.
    rng = np.random.RandomState(0)
    n_samples, n_features = 2000, 10
    X = rng.randn(n_samples, n_features)
    coef_true = rng.randn(n_features)
    y = (X @ coef_true + 0.1 * rng.randn(n_samples) > 0).astype(np.float64)

    weights = rng.randint(1, 6, size=n_samples).astype(float)
    heavy_idx = rng.choice(n_samples, size=20, replace=False)
    weights[heavy_idx] = rng.randint(50, 400, size=20)

    alpha_scaled = 1.0 / n_samples
    beta_scaled = 1.0 / n_samples
    max_squared_sum = row_norms(X, squared=True).max()
    step_size = get_auto_step_size(
        max_squared_sum, alpha_scaled, "log", True, n_samples=n_samples, is_saga=True
    )

    def weighted_log_loss(coef, intercept):
        z = X @ coef + intercept
        p = np.clip(1 / (1 + np.exp(-z)), 1e-12, 1 - 1e-12)
        return np.average(-(y * np.log(p) + (1 - y) * np.log(1 - p)), weights=weights)

    def run(use_weighted_sampling, max_iter=2000):
        from sklearn.linear_model._base import make_dataset
        from sklearn.linear_model._sag_fast import sag64

        dataset, intercept_decay = make_dataset(X, y, weights, 0)
        coef = np.zeros((n_features, 1), dtype=np.float64, order="C")
        intercept = np.zeros(1, dtype=np.float64)
        if use_weighted_sampling:
            alias_index, prob_threshold, importance_weight = _build_sag_alias_tables(
                weights, np.float64
            )
        else:
            alias_index = np.zeros(1, dtype=np.int32)
            prob_threshold = np.zeros(1, dtype=np.float64)
            importance_weight = np.zeros(1, dtype=np.float64)
        sag64(
            dataset,
            coef,
            intercept,
            n_samples,
            n_features,
            1,
            1e-4,
            max_iter,
            "log",
            step_size,
            alpha_scaled,
            beta_scaled,
            np.zeros((n_features, 1), dtype=np.float64, order="C"),
            np.zeros((n_samples, 1), dtype=np.float64, order="C"),
            np.zeros(n_samples, dtype=np.int32, order="C"),
            0,
            True,
            np.zeros(1, dtype=np.float64),
            intercept_decay,
            True,
            use_weighted_sampling,
            alias_index,
            prob_threshold,
            importance_weight,
            12345,
            0,
        )
        return weighted_log_loss(coef.ravel(), intercept[0])

    loss_before_fix = run(use_weighted_sampling=False)
    loss_after_fix = run(use_weighted_sampling=True)

    assert loss_after_fix < 0.5 * loss_before_fix


def test_sag_weighted_sampling_noop_for_uniform_weights():
    """When sample_weight is uniform (including the default, None), the SAG
    solver must take the exact original (unweighted) code path: fitting
    with sample_weight=None and with an explicit all-ones sample_weight
    must produce bit-identical coefficients for a fixed random_state.
    """
    X, y = make_classification(n_samples=200, random_state=0)

    clf_none = LogisticRegression(solver="saga", random_state=0, max_iter=200).fit(X, y)
    clf_ones = LogisticRegression(solver="saga", random_state=0, max_iter=200).fit(
        X, y, sample_weight=np.ones(X.shape[0])
    )

    assert_allclose(clf_none.coef_, clf_ones.coef_, rtol=0, atol=0)
    assert clf_none.n_iter_ == clf_ones.n_iter_
