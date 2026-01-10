# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import math
import re

import numpy as np
import pytest

from sklearn.base import clone
from sklearn.datasets import load_iris, make_blobs, make_classification
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.linear_model._sag import get_auto_step_size, sag_solver
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import check_random_state, compute_class_weight
from sklearn.utils._testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
)
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
    sample_weight=None,
    max_iter=1,
    dloss=None,
    decay=1.0,
    fit_intercept=True,
    saga=False,
    tol=0,
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

    for epoch in range(max_iter):
        previous_weights = weights.copy()
        for k in range(n_samples):
            idx = int(rng.rand() * n_samples)
            # idx = k
            entry = X[idx]
            seen.add(idx)
            S_seen = len(seen)
            if sample_weight is not None:
                S_seen = sample_weight[list(seen)].sum()
            if S_seen == 0:
                continue
            p = np.dot(entry, weights) + intercept
            gradient = dloss(p, y[idx])
            update = entry * gradient + alpha * weights
            gradient_correction = update - gradient_memory[idx]
            if sample_weight is not None:
                gradient_correction *= sample_weight[idx]
            sum_gradient += gradient_correction
            gradient_memory[idx] = update
            weights -= step_size * sum_gradient / S_seen
            if saga:
                weights -= gradient_correction * step_size * (1 - 1.0 / S_seen)

            if fit_intercept:
                update = gradient
                gradient_correction = update - intercept_gradient_memory[idx]
                if sample_weight is not None:
                    gradient_correction *= sample_weight[idx]
                intercept_sum_gradient += gradient_correction
                intercept_gradient_memory[idx] = update
                intercept -= step_size * intercept_sum_gradient / S_seen * decay
                if saga:
                    intercept -= gradient_correction * step_size * (1 - 1.0 / S_seen)

        # stopping criteria
        max_weight = np.abs(weights).max()
        max_change = np.abs(weights - previous_weights).max()
        if (max_weight != 0 and max_change / max_weight <= tol) or (
            max_weight == 0 and max_change == 0
        ):
            break

    n_iter = epoch + 1
    return weights, intercept, n_iter


def sag_sparse(
    X,
    y,
    step_size,
    alpha,
    sample_weight=None,
    max_iter=1,
    dloss=None,
    decay=1.0,
    fit_intercept=True,
    saga=False,
    tol=0,
    random_state=0,
):
    if step_size * alpha == 1.0:
        raise ZeroDivisionError(
            "Sparse sag does not handle the case step_size * alpha == 1"
        )
    n_samples, n_features = X.shape[0], X.shape[1]

    weights = np.zeros(n_features)
    actual_weights = np.zeros(n_features)
    sum_gradient = np.zeros(n_features)
    last_updated = np.zeros(n_features, dtype=int)
    gradient_memory = np.zeros(n_samples)
    rng = check_random_state(random_state)
    intercept = 0.0
    intercept_sum_gradient = 0.0
    wscale = 1.0
    decay = 1.0
    seen = set()

    c_sum = np.zeros(max_iter * n_samples)

    counter = 0
    for epoch in range(max_iter):
        previous_weights = actual_weights
        for k in range(n_samples):
            # idx = k
            idx = int(rng.rand() * n_samples)
            entry = X[idx]
            seen.add(idx)
            S_seen = len(seen)
            if sample_weight is not None:
                S_seen = sample_weight[list(seen)].sum()
            if S_seen == 0:
                continue

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
                        gradient_correction[j] * step_size * (1 - 1.0 / S_seen) / wscale
                    )

            if fit_intercept:
                gradient_correction = gradient - gradient_memory[idx]
                intercept_sum_gradient += gradient_correction
                gradient_correction *= step_size * (1.0 - 1.0 / S_seen)
                if saga:
                    intercept -= (
                        step_size * intercept_sum_gradient / S_seen * decay
                    ) + gradient_correction
                else:
                    intercept -= step_size * intercept_sum_gradient / S_seen * decay

            gradient_memory[idx] = gradient

            wscale *= 1.0 - alpha * step_size
            if counter == 0:
                c_sum[0] = step_size / (wscale * S_seen)
            else:
                c_sum[counter] = c_sum[counter - 1] + step_size / (wscale * S_seen)

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
        # Actual weights is wscale times the just-in-time updates for all features
        actual_weights = weights.copy()
        for j in range(n_features):
            if last_updated[j] == 0:
                actual_weights[j] -= c_sum[counter - 1] * sum_gradient[j]
            else:
                actual_weights[j] -= (
                    c_sum[counter - 1] - c_sum[last_updated[j] - 1]
                ) * sum_gradient[j]
        actual_weights *= wscale

        # stopping criteria
        max_weight = np.abs(actual_weights).max()
        max_change = np.abs(actual_weights - previous_weights).max()
        if (max_weight != 0 and max_change / max_weight <= tol) or (
            max_weight == 0 and max_change == 0
        ):
            break
    # Actual weights is wscale times the just-in-time updates for all features
    for j in range(n_features):
        if last_updated[j] == 0:
            weights[j] -= c_sum[counter - 1] * sum_gradient[j]
        else:
            weights[j] -= (
                c_sum[counter - 1] - c_sum[last_updated[j] - 1]
            ) * sum_gradient[j]
    weights *= wscale

    n_iter = epoch + 1
    return weights, intercept, n_iter


def get_step_size(
    X, alpha, fit_intercept, classification=True, sample_weight=None, is_saga=False
):
    # Lipschitz smoothness constant for f_i(w) = s_i (loss_i(w)) + alpha ||w||^2):
    # L_i = s_i ( kappa * (||x_i||^2 + fit_intercept) + alpha )
    # where kappa = 1/4 for classification and 1 for regression
    kappa = 0.25 if classification else 1.0
    L = kappa * (np.sum(X * X, axis=1) + fit_intercept) + alpha
    if sample_weight is not None:
        L *= sample_weight
    L = L.max()
    # step_size = 1/L for SAG https://arxiv.org/abs/1309.2388
    # step_size = 1/3L for SAGA https://arxiv.org/abs/1407.0202
    step_size = 1 / (3 * L) if is_saga else 1 / L
    return step_size


def test_classifier_matching():
    n_samples = 20
    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0, cluster_std=0.1)
    # y must be 0 or 1
    alpha = 1.1
    fit_intercept = True
    step_size = get_step_size(X, alpha, fit_intercept)
    for solver in ["sag", "saga"]:
        if solver == "sag":
            max_iter = 80
        else:
            # SAGA variance w.r.t. stream order is higher
            max_iter = 300
        clf = LogisticRegression(
            solver=solver,
            fit_intercept=fit_intercept,
            tol=1e-11,
            C=1.0 / alpha / n_samples,
            max_iter=max_iter,
            random_state=10,
        )
        clf.fit(X, y)

        weights, intercept, n_iter = sag_sparse(
            X,
            2 * y - 1,  # y must be -1 or +1
            step_size,
            alpha,
            max_iter=max_iter,
            dloss=log_dloss,
            fit_intercept=fit_intercept,
            saga=solver == "saga",
        )
        weights2, intercept2, n_iter2 = sag(
            X,
            2 * y - 1,  # y must be -1 or +1
            step_size,
            alpha,
            max_iter=max_iter,
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
    max_iter = 100
    fit_intercept = True

    step_size = get_step_size(X, alpha, fit_intercept, classification=False)
    clf = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00000000001,
        solver="sag",
        alpha=alpha * n_samples,
        max_iter=max_iter,
    )
    clf.fit(X, y)

    weights1, intercept1, n_iter1 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
    )
    weights2, intercept2, n_iter2 = sag(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
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
    max_iter = 100
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
        max_iter=max_iter,
        random_state=42,
    )
    clf2 = clone(clf1)
    clf3 = Ridge(
        fit_intercept=fit_intercept,
        tol=0.00001,
        solver="lsqr",
        alpha=alpha,
        max_iter=max_iter,
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

    spweights1, spintercept1, sp_n_iter1 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=squared_dloss,
        fit_intercept=fit_intercept,
        random_state=rng,
    )

    spweights2, spintercept2, sp_n_iter2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=squared_dloss,
        decay=0.01,
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


@pytest.mark.parametrize("saga", [True, False])
@pytest.mark.parametrize("fit_intercept", [True, False])
@pytest.mark.parametrize("classification", [True, False])
@pytest.mark.parametrize("sample_weight", [None, "random"])
def test_get_auto_step_size(saga, fit_intercept, classification, sample_weight):
    alpha = 1.2
    n_samples = 100
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, 10)
    squared_sum = np.sum(X * X, axis=1)
    if sample_weight == "random":
        sample_weight = rng.rand(n_samples)

    step_size = get_step_size(
        X, alpha, fit_intercept, classification, sample_weight, saga
    )
    loss = "log" if classification else "squared"
    auto_step_size = get_auto_step_size(
        squared_sum, alpha, loss, fit_intercept, sample_weight, saga
    )
    assert_allclose(step_size, auto_step_size)

    msg = "Unknown loss function for SAG solver, got wrong instead of"
    with pytest.raises(ValueError, match=msg):
        get_auto_step_size(squared_sum, alpha, "wrong", fit_intercept)


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
@pytest.mark.xfail()
def test_sag_classifier_computed_correctly(csr_container):
    """tests if the binary classifier is computed correctly"""
    alpha = 0.1
    n_samples = 100
    max_iter = 50
    tol = 1e-10
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
        max_iter=max_iter,
        tol=tol,
        random_state=77,
        fit_intercept=fit_intercept,
    )
    clf2 = clone(clf1)

    clf1.fit(X, y)
    clf2.fit(csr_container(X), y)

    spweights, spintercept, sp_n_iter = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=log_dloss,
        fit_intercept=fit_intercept,
    )
    spweights2, spintercept2, sp_n_iter2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=log_dloss,
        decay=0.05,
        fit_intercept=fit_intercept,
    )

    assert_array_almost_equal(clf1.coef_.ravel(), spweights.ravel(), decimal=2)
    assert_almost_equal(clf1.intercept_, spintercept, decimal=1)

    assert_array_almost_equal(clf2.coef_.ravel(), spweights2.ravel(), decimal=2)
    assert_almost_equal(clf2.intercept_, spintercept2, decimal=1)


@pytest.mark.filterwarnings("ignore:The max_iter was reached")
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
@pytest.mark.xfail()
def test_sag_multiclass_computed_correctly(csr_container):
    """tests if the multiclass classifier is computed correctly"""
    alpha = 0.1
    n_samples = 100
    tol = 1e-10
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

        spweights1, spintercept1, sp_n_iter1 = sag_sparse(
            X,
            y_encoded,
            step_size,
            alpha,
            dloss=log_dloss,
            max_iter=max_iter,
            fit_intercept=fit_intercept,
        )
        spweights2, spintercept2, sp_n_iter2 = sag_sparse(
            X,
            y_encoded,
            step_size,
            alpha,
            dloss=log_dloss,
            max_iter=max_iter,
            decay=0.01,
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
@pytest.mark.xfail()
def test_binary_classifier_class_weight(csr_container):
    """tests binary classifier with classweights for each class"""
    alpha = 0.1
    n_samples = 50
    max_iter = 20
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
        max_iter=max_iter,
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
    spweights, spintercept, sp_n_iter = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=log_dloss,
        sample_weight=sample_weight,
        fit_intercept=fit_intercept,
    )
    spweights2, spintercept2, sp_n_iter2 = sag_sparse(
        X,
        y,
        step_size,
        alpha,
        max_iter=max_iter,
        dloss=log_dloss,
        decay=0.01,
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


@pytest.mark.parametrize("solver", [sag, sag_sparse, sag_solver])
@pytest.mark.parametrize("decay", [1.0, 0.05])
@pytest.mark.parametrize("saga", [True, False])
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_sag_weighted_classification_convergence(solver, decay, saga, fit_intercept):
    # FIXME
    if solver == sag_solver:
        pytest.xfail("sag_solver fail convergence test")
    n_samples = 100
    max_iter = 1000
    tol = 1e-10
    alpha = 1.1

    X, y = make_blobs(n_samples=n_samples, centers=2, random_state=0, cluster_std=0.1)
    n_features = X.shape[1]
    y = 1 * (y >= 1)
    sample_weights = np.random.randint(0, 3, size=n_samples)

    est = LogisticRegression(
        max_iter=max_iter,
        tol=tol,
        fit_intercept=fit_intercept,
        solver="lbfgs",
        penalty="l2",
        C=1 / (sample_weights.sum() * alpha),
    )
    est.fit(X, y, sample_weight=sample_weights)
    true_weights = est.coef_[0]
    true_intercept = est.intercept_

    y = 2 * y - 1

    if solver != sag_solver:
        sag_kwargs = dict(
            dloss=log_dloss,
            max_iter=max_iter,
            decay=decay,
            tol=tol,
            fit_intercept=fit_intercept,
            saga=saga,
        )

        step_size = get_step_size(
            X,
            alpha,
            fit_intercept,
            classification=True,
            sample_weight=sample_weights,
            is_saga=saga,
        )

        weights, intercept, n_iter = solver(
            X, y, step_size, alpha, sample_weight=sample_weights, **sag_kwargs
        )
    else:
        sag_kwargs = dict(
            loss="log",
            max_iter=max_iter,
            tol=tol,
            warm_start_mem={
                "coef": np.zeros((n_features + fit_intercept, 1), dtype=X.dtype)
            },
            is_saga=saga,
        )
        weights, n_iter, warm_start_mem = solver(
            X, y, sample_weights, alpha=alpha, **sag_kwargs
        )
        if fit_intercept:
            intercept = weights[-1]
            weights = weights[:-1]
    assert weights.shape == (n_features,)
    assert_allclose(weights, true_weights)
    assert_allclose(intercept, true_intercept)


@pytest.mark.parametrize("solver", [sag, sag_sparse, sag_solver])
@pytest.mark.parametrize("decay", [1.0, 0.05])
@pytest.mark.parametrize("saga", [True, False])
@pytest.mark.parametrize("fit_intercept", [True, False])
def test_sag_weighted_regression_convergence(solver, decay, saga, fit_intercept):
    # FIXME
    if saga and fit_intercept:
        pytest.xfail("saga + fit_intercept fail convergence test")
    if solver == sag_solver:
        pytest.xfail("sag_solver fail convergence test")
    n_samples = 15
    max_iter = 1000
    tol = 1e-11
    alpha = 1.1
    rng = np.random.RandomState(42)
    X = rng.rand(n_samples, n_samples * 2)
    n_features = X.shape[1]
    y = rng.randint(0, 3, size=n_samples)
    # Use random integers (including zero) as weights.
    sample_weights = rng.randint(0, 5, size=n_samples)
    # FIXME: zero sample_weight lead to division by zero in sag_sparse
    if solver == sag_sparse:
        sample_weights += 1
    est = Ridge(
        max_iter=max_iter,
        tol=tol,
        fit_intercept=fit_intercept,
        solver="auto",
        alpha=(sample_weights.sum() * alpha),
    )
    est.fit(X, y, sample_weight=sample_weights)
    true_weights = est.coef_.ravel()
    true_intercept = est.intercept_

    if solver != sag_solver:
        sag_kwargs = dict(
            dloss=squared_dloss,
            max_iter=max_iter,
            decay=decay,
            tol=tol,
            fit_intercept=fit_intercept,
            saga=saga,
        )
        step_size = get_step_size(
            X,
            alpha,
            fit_intercept,
            classification=False,
            sample_weight=sample_weights,
            is_saga=saga,
        )
        weights, intercept, n_iter = solver(
            X, y, step_size, alpha, sample_weight=sample_weights, **sag_kwargs
        )
    else:
        sag_kwargs = dict(
            loss="squared",
            max_iter=max_iter,
            tol=tol,
            warm_start_mem={
                "coef": np.zeros((n_features + fit_intercept, 1), dtype=X.dtype)
            },
            is_saga=saga,
        )
        weights, n_iter, warm_start_mem = solver(
            X, y, sample_weights, alpha=alpha, **sag_kwargs
        )
        if fit_intercept:
            intercept = weights[-1]
            weights = weights[:-1]
    assert weights.shape == (n_features,)
    assert_allclose(weights, true_weights, atol=1e-10)
    assert_allclose(intercept, true_intercept)
