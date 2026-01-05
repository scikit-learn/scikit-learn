# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# Most of the solvers are tested elsewhere. Here, we only test certrain low level
# aspects.
import numpy as np
import pytest
from scipy import sparse

from sklearn._loss import HalfMultinomialLoss, HalfSquaredError
from sklearn.linear_model._cd_fast import enet_coordinate_descent_gram
from sklearn.linear_model._glm._newton_solver import (
    enet_coordinate_descent_multinomial_py,
    min_norm_subgradient,
)
from sklearn.linear_model._linear_loss import (
    LinearModelLoss,
    Multinomial_LDL_Decomposition,
)
from sklearn.utils._testing import assert_allclose
from sklearn.utils.extmath import softmax


def test_min_norm_subgradient():
    coef = np.array([-2.0, 1, 0, 1])  # coef[-1] might be the intercept

    # The test is mostly independent on the choice of base_loss, with the exception of
    # shapes.
    linear_loss = LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True)
    alpha = 8.0  # l1_reg_strength
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [2 - 8, 2 + 8, 0, 2])

    grad = np.full_like(coef, -9.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [-9 - 8, -9 + 8, -1 * (9 - 8), -9])

    linear_loss = LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=False)
    alpha = 1.0
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [2 - 1, 2 + 1, +1 * (2 - 1), 2 + 1])

    mnsg = min_norm_subgradient(alpha, -grad, coef, linear_loss)
    assert_allclose(mnsg, [-2 - 1, -2 + 1, -1 * (2 - 1), -2 + 1])

    grad = np.full_like(coef, -9.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [-9 - 1, -9 + 1, -1 * (9 - 1), -9 + 1])

    n_classes = 3
    coef = np.repeat(coef, n_classes)
    linear_loss = LinearModelLoss(base_loss=HalfMultinomialLoss(), fit_intercept=True)
    alpha = 1.0
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    expected = np.repeat([2 - 1, 2 + 1, +1 * (2 - 1), 2], 3)
    assert_allclose(mnsg, expected)

    # Same as above but with 2-dim coef, coef.shape = (n_classes, n_features)
    coef = coef.reshape(n_classes, -1, order="F")
    grad = grad.reshape(n_classes, -1, order="F")
    expected = expected.reshape(n_classes, -1, order="F")
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, expected)

    # alpha = 0
    mnsg = min_norm_subgradient(0, grad, coef, linear_loss)
    assert mnsg is grad


@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("cd_multinomial", [enet_coordinate_descent_multinomial_py])
@pytest.mark.parametrize("sparse_X", [False, True])
@pytest.mark.parametrize("alpha", [0, 1])
@pytest.mark.parametrize("max_iter", [1, 5])
def test_cd_multinomial(cd_multinomial, sparse_X, alpha, max_iter):
    """Test that enet_coordinate_descent_multinomial is equivalent to
    enet_coordinate_descent_gram."""
    n_samples, n_features, n_classes = 5, 2, 3
    beta = 1
    W = np.zeros((n_classes, n_features), dtype=np.float64, order="F")
    y = np.arange(n_samples, dtype=np.float64) % n_classes
    X = np.zeros((n_samples, n_features), order="F")

    W.flat = np.arange(n_classes * n_features)
    W -= np.mean(W, axis=0)[None, :]
    X[:, 0] = np.arange(n_samples)
    X[:, 1] = 1
    raw = np.asfortranarray(X @ W.T)
    proba = np.asfortranarray(softmax(raw))
    lin_loss = LinearModelLoss(
        HalfMultinomialLoss(n_classes=n_classes), fit_intercept=False
    )
    grad_pointwise = lin_loss.base_loss.gradient(y_true=y, raw_prediction=raw)
    grad_pointwise /= n_samples
    LDL = Multinomial_LDL_Decomposition(proba=proba)
    sqrt_D_Lt_raw = LDL.sqrt_D_Lt_matmul(raw.copy()) / np.sqrt(n_samples)
    b = sqrt_D_Lt_raw - LDL.inverse_L_sqrt_D_matmul(grad_pointwise.copy()) * np.sqrt(
        n_samples
    )

    # As in NewtonCDGramSolver, but correct y=b.
    wflat = W.flatten(order="F")  # returns copy
    G, H, _ = lin_loss.gradient_hessian(coef=wflat, X=X, y=y, l2_reg_strength=0)
    # Note that we need to change the order of features and classes to mimic
    # the CD loop order in enet_coordinate_descent_multinomial: first classes, then
    # innermost features. Therefore we need to rearrange/reshape w, q and Q from F-
    # to C-order.
    w_gram, gap_gram, tol, n_iter = enet_coordinate_descent_gram(
        w=wflat.reshape(n_classes, -1, order="F").ravel(order="C"),
        alpha=alpha,
        beta=beta,
        Q=H.reshape(n_classes, n_features, n_classes, n_features, order="F").reshape(
            n_classes * n_features, n_classes * n_features, order="C"
        ),
        q=(H @ wflat - G).reshape(n_classes, -1, order="F").ravel(order="C"),
        y=b.ravel(),  # used in dual gap
        max_iter=max_iter,
        tol=1e-20,
        rng=np.random.RandomState(0),
        random=False,
        positive=False,
        early_stopping=False,
    )
    # As in NewtonCDSolver
    w_cd, gap_cd, tol, n_iter = cd_multinomial(
        W=W.copy(order="F"),
        alpha=alpha,
        beta=beta,
        X=sparse.csc_array(X) if sparse_X else X,
        sample_weight=None,
        raw_prediction=raw,
        grad_pointwise=grad_pointwise,
        proba=proba,
        fit_intercept=False,
        max_iter=max_iter,
        tol=1e-20,
        early_stopping=False,
    )
    assert_allclose(w_cd, w_gram.reshape(n_classes, -1, order="C"))
    assert gap_cd == pytest.approx(gap_gram, abs=5e-9)

    # Now the same, but WITH INTERCEPT, such that we have the same raw_prediction and
    # same proba as above.
    X = np.asfortranarray(X[:, 0:1])
    assert_allclose(X @ W[:, :-1].T + W[:, -1], raw)
    # A0 = sqrt(D) L' 1
    A0 = np.zeros((n_samples, n_classes, n_classes))
    for k in range(n_classes):
        A0[:, k, k] = LDL.sqrt_d[:, k]
        for l in range(k + 1, n_classes):
            A0[:, k, l] = -LDL.sqrt_d[:, k] * LDL.p[:, l] * LDL.q_inv[:, k]
    A0 /= np.sqrt(n_samples)
    A0tA0 = np.einsum("ijk,ijl->kl", A0, A0)
    # is the same as
    # A0tA0 = np.tensordot(A0, A0, axes=([0, 1], [0, 1]))
    lin_loss = LinearModelLoss(
        HalfMultinomialLoss(n_classes=n_classes), fit_intercept=True
    )
    G, H, _ = lin_loss.gradient_hessian(
        coef=W.ravel(order="F"), X=X, y=y, l2_reg_strength=0
    )
    Hcoef = H @ W.ravel(order="F")
    H0_coef = Hcoef[-n_classes:]
    H0 = H[-n_classes:, :-n_classes]  # shape (n_classes, n_classes * n_features)
    H00 = H[-n_classes:, -n_classes:]
    H00_inv = np.linalg.pinv(H00)
    H = H[:-n_classes, :-n_classes]  # only the part without intercepts
    q = Hcoef[:-n_classes] - G[:-n_classes]
    q0 = Hcoef[-n_classes:] - G[-n_classes:]
    H_centered = H - H0.T @ H00_inv @ H0
    q_centered = q - H0.T @ H00_inv @ q0

    assert_allclose(A0tA0, H00)

    H00_pinv_H0 = H00_inv @ H0
    raw_centered = raw - W[:, -1]
    raw_centered -= (H00_pinv_H0 @ W[:, :-1].ravel(order="F"))[None, :]
    sqrt_D_Lt_raw = LDL.sqrt_D_Lt_matmul(raw_centered) / np.sqrt(n_samples)
    t = H00_inv @ grad_pointwise.sum(axis=0)  # H00^(-1) 1' g
    g_pointwise_centered = grad_pointwise.copy()
    for k in range(n_classes):
        h = proba[:, k] * (1 - proba[:, k]) / n_samples
        g_pointwise_centered[:, k] -= h * t[k]
        for l in range(k + 1, n_classes):
            h = -proba[:, k] * proba[:, l] / n_samples
            g_pointwise_centered[:, k] -= h * t[l]
            g_pointwise_centered[:, l] -= h * t[k]
    b = sqrt_D_Lt_raw - LDL.inverse_L_sqrt_D_matmul(g_pointwise_centered) * np.sqrt(
        n_samples
    )

    # residual
    R = b - sqrt_D_Lt_raw
    # rotated residual
    LD_R = LDL.L_sqrt_D_matmul(R / np.sqrt(n_samples))
    # With fit_intercept, if should sum to zero over samples.
    assert_allclose(np.sum(LD_R, axis=0), 0, atol=1e-14)

    w_gram, gap_gram, tol, n_iter = enet_coordinate_descent_gram(
        w=W[:, :-1].ravel(order="F").copy(),
        alpha=alpha,
        beta=beta,
        Q=H_centered,  # only 1 feature, so no reshaping
        q=q_centered,
        y=b.ravel(),  # used in dual gap
        max_iter=max_iter,
        tol=1e-20,
        rng=np.random.RandomState(0),
        random=False,
        positive=False,
        early_stopping=False,
    )
    w_cd, gap_cd, tol, n_iter = cd_multinomial(
        W=W.copy(order="F"),
        alpha=alpha,
        beta=beta,
        X=sparse.csc_array(X) if sparse_X else X,
        sample_weight=None,
        raw_prediction=raw,
        grad_pointwise=grad_pointwise,
        proba=proba,
        fit_intercept=True,
        max_iter=max_iter,
        tol=1e-20,
        early_stopping=False,
    )
    assert_allclose(w_cd[:, :-1], w_gram.reshape(n_classes, -1, order="F"))
    assert gap_cd == pytest.approx(gap_gram, abs=5e-9)
