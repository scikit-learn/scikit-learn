# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# Most of the solvers are tested elsewhere. Here, we only test certain low level
# aspects.
import warnings

import numpy as np
import pytest
from scipy import sparse

from sklearn import config_context
from sklearn._loss import HalfMultinomialLoss, HalfPoissonLoss, HalfSquaredError
from sklearn._loss.loss import HalfPoissonLossArrayAPI
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model._cd_fast import (
    enet_coordinate_descent_gram,
    enet_coordinate_descent_multinomial,
)
from sklearn.linear_model._glm._newton_solver import (
    NewtonCDGramSolver,
    NewtonCholeskySolver,
    min_norm_subgradient,
)
from sklearn.linear_model._linear_loss import (
    LinearModelLoss,
    Multinomial_LDL_Decomposition,
)
from sklearn.utils._array_api import _convert_to_numpy
from sklearn.utils._testing import (
    _array_api_for_tests,
    assert_allclose,
    assert_array_equal,
)
from sklearn.utils.extmath import softmax
from sklearn.utils.fixes import parse_version, sp_version


def test_newton_cholesky_array_api_zero_curvature():
    """Refresh the full gradient, including the intercept, at zero curvature."""
    xp, device = _array_api_for_tests("array_api_strict", "device1", "float64")
    X = xp.asarray([[1.0], [2.0]], device=device)
    y = xp.asarray([1.0, 2.0], device=device)
    coef = xp.asarray([-1000.0, 0.0], device=device)
    linear_loss = LinearModelLoss(
        base_loss=HalfPoissonLossArrayAPI(xp=xp, device=device), fit_intercept=True
    )
    solver = NewtonCholeskySolver(coef=coef, linear_loss=linear_loss)
    # The predicted rates underflow to zero, while the gradient is finite.
    solver.raw_prediction = X @ coef[:-1] + coef[-1]

    with config_context(array_api_dispatch=True):
        solver.update_gradient_hessian(X, y, sample_weight=None)

    assert solver.hessian_warning
    assert_allclose(_convert_to_numpy(solver.gradient, xp), [-2.5, -1.5])


def test_min_norm_subgradient_raises():
    """Test that min_norm_subgradient raises on wrong input."""
    alpha = 1
    linear_loss = LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True)
    msg = "Shapes of 'gradient' and 'coef' must match"

    gradient = np.arange(4)
    coef = gradient[None, :]
    with pytest.raises(ValueError, match=msg):
        min_norm_subgradient(alpha, gradient, coef, linear_loss)

    coef = np.arange(5)
    with pytest.raises(ValueError, match=msg):
        min_norm_subgradient(alpha, gradient, coef, linear_loss)


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


def _newton_solver_before_line_search(
    coef_newton,
    coef=None,
    y=(1.0, 2.0),
    tol=1e-4,
    fit_intercept=False,
    l1_reg_strength=0,
):
    """Return a tiny Poisson problem and its solver, ready to call line_search.

    Without intercept and penalty, the loss is 1/2 * sum(exp(coef) - y * coef) with
    minimum at coef = log(y).
    """
    y, X = np.asarray(y, dtype=float), np.eye(2)
    linear_loss = LinearModelLoss(
        base_loss=HalfPoissonLoss(), fit_intercept=fit_intercept
    )
    if coef is None:
        coef = linear_loss.init_zero_coef(X)
    # line_search is the one of the base class NewtonSolver.
    sol = NewtonCDGramSolver(
        coef=coef, linear_loss=linear_loss, l1_reg_strength=l1_reg_strength, tol=tol
    )
    sol.setup(X=X, y=y, sample_weight=None)
    sol.iteration = 1
    sol.use_fallback_lbfgs_solve = False
    # For coef = 0, y = [1, 2] and without intercept:
    # loss = 0.5, gradient = [0, -0.5], hessian = [[0.5, 0], [0, 0.5]]
    sol.update_gradient_hessian(X=X, y=y, sample_weight=None)
    sol.coef_newton = np.asarray(coef_newton, dtype=float)
    sol.gradient_times_newton = sol.gradient @ sol.coef_newton
    return X, y, sol


def test_newton_solver_line_search_huge_step():
    """Test that line search copes with non-finite and huge loss values.

    The Newton step is a descent direction, but way too large. For step size
    alpha=1, the loss is exp(1000) = inf, for alpha=1/2 it is exp(500) ~ 1e217.
    The minimum along this direction is at alpha = log(2) / 1000 ~ 0.0007.
    Note: This can happens in the first Newton step with a log-link.
    """
    coef_newton = np.array([0.0, 1000])
    X, y, sol = _newton_solver_before_line_search(coef_newton=coef_newton)
    loss_old = sol.loss_value
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        # Non-finite pointwise gradients may trigger "invalid value in matmul".
        warnings.simplefilter("ignore", RuntimeWarning)
        sol.line_search(X=X, y=y, sample_weight=None)

    assert not sol.use_fallback_lbfgs_solve
    assert sol.loss_value < loss_old
    # The largest step size 1/2**k that satisfies the Armijo condition is k=10.
    assert_allclose(sol.coef, coef_newton / 2**10)


def test_newton_solver_line_search_l1_excludes_intercept():
    """Test that the L1 term of the Armijo condition does not include the intercept.

    We take a step only for the intercept. It is a descent direction but too large:
    the loss is exp(alpha) - 1.5 * alpha, i.e. 1 for alpha=0, 1.22 for alpha=1 and
    0.90 for alpha=1/2.
    If the change of the absolute value of the intercept is added to the Armijo
    term, weighted with the large L1 penalty strength, then the Armijo term gets
    positive and the full step is accepted despite the increased loss.
    """
    coef_newton = np.array([0.0, 0.0, 1.0])
    X, y, sol = _newton_solver_before_line_search(
        coef_newton=coef_newton, fit_intercept=True, l1_reg_strength=1000
    )
    loss_old = sol.loss_value
    assert sol.gradient_times_newton == pytest.approx(-0.5)
    sol.line_search(X=X, y=y, sample_weight=None)

    assert sol.loss_value < loss_old
    assert_allclose(sol.coef, coef_newton / 2)


@pytest.mark.parametrize("l1_reg_strength", [0, 0.25])
def test_newton_solver_line_search_inexact_step_at_machine_precision(l1_reg_strength):
    """Test line search for an inexact Newton step very close to the minimum.

    Close to the minimum, the change of the loss can be smaller than what floating
    point numbers can resolve. On top, an inexact Newton step, e.g. from an inner
    coordinate descent solver that stopped early, might decrease the loss without
    decreasing the norm of the gradient. The line search must then rely on the
    directional derivative and accept the step.
    """
    # We choose y such that the minimum is at coef = log([2, 4]), also with L1
    # penalty. There, the hessian is diag([1, 2]). We start at a point where the
    # gradient, or the minimum norm subgradient in case of L1, is g = [eps, 0].
    eps = 1e-9
    y = np.array([2.0, 4.0]) + 2 * l1_reg_strength
    coef = np.log([2.0, 4.0]) + np.array([eps, 0])
    # The exact Newton step is [-eps, 0]. Our inexact step has
    #   g @ coef_newton = -0.5 * eps**2
    #   coef_newton @ hessian @ coef_newton = 0.53125 * eps**2
    # The quadratic approximation of the loss, which is very good for such tiny steps,
    # decreases by (-0.5 + 0.53125 / 2) * eps**2 ~ -0.23e-18. This is far below the
    # resolution of the loss value of about -0.5.
    # The gradient at coef + alpha * coef_newton is eps * [1 - alpha/2, 3/4 * alpha].
    # Its 1-norm, eps * (1 + alpha / 4), increases for all step sizes alpha > 0.
    coef_newton = np.array([-eps / 2, 0.375 * eps])
    # We set a tiny tol to avoid that the line search stops because of
    # max|gradient| <= tol.
    X, y, sol = _newton_solver_before_line_search(
        coef_newton=coef_newton,
        coef=coef,
        y=y,
        tol=1e-12,
        l1_reg_strength=l1_reg_strength,
    )
    g = min_norm_subgradient(l1_reg_strength, sol.gradient, sol.coef, sol.linear_loss)
    assert_allclose(g, [eps, 0], rtol=1e-6, atol=1e-15)
    coef_old = sol.coef.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        sol.line_search(X=X, y=y, sample_weight=None)

    assert not sol.use_fallback_lbfgs_solve
    # The full step is accepted.
    assert_array_equal(sol.coef, coef_old + coef_newton)


@pytest.mark.parametrize(
    "coef_newton",
    [
        [1, -1],  # not a descent direction
        [0, np.nan],
        [0, np.inf],
    ],
)
def test_newton_solver_line_search_failure_restores_state(coef_newton):
    """Test that a failed line search restores coef, loss and gradient.

    The fallback solver (lbfgs) should not start from a failed, possibly non-finite,
    trial point of the line search.
    """
    X, y, sol = _newton_solver_before_line_search(coef_newton=coef_newton)
    coef_old = sol.coef.copy()
    loss_old = sol.loss_value
    gradient_old = sol.gradient.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with pytest.warns(ConvergenceWarning, match="Line search of Newton solver"):
            sol.line_search(X=X, y=y, sample_weight=None)

    assert sol.use_fallback_lbfgs_solve
    assert_array_equal(sol.coef, coef_old)
    assert sol.loss_value == loss_old
    assert_array_equal(sol.gradient, gradient_old)

    # The fallback solver then finds the solution exp(coef) = y.
    sol.fallback_lbfgs_solve(X=X, y=y, sample_weight=None)
    assert sol.converged
    assert_allclose(sol.coef, np.log(y), atol=1e-3)


@pytest.mark.skipif(
    sp_version < parse_version("1.17"),
    reason="scipy version 1.17 required for sparse slicing",
)
@pytest.mark.filterwarnings("ignore::sklearn.exceptions.ConvergenceWarning")
@pytest.mark.parametrize("sparse_X", [False, True])
@pytest.mark.parametrize("alpha", [0, 1])
@pytest.mark.parametrize("max_iter", [1, 5])
def test_cd_multinomial(sparse_X, alpha, max_iter):
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
    w_cd, gap_cd, tol, n_iter = enet_coordinate_descent_multinomial(
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
    raw_centered = raw - W[:, -1] - (H00_pinv_H0 @ W[:, :-1].ravel(order="F"))[None, :]
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
    w_cd, gap_cd, tol, n_iter = enet_coordinate_descent_multinomial(
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
