# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
Newton solver for Generalized Linear Models
"""

import warnings
from abc import ABC, abstractmethod

import numpy as np
import scipy.linalg
import scipy.optimize

from sklearn._loss.loss import HalfSquaredError
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model._cd_fast import (
    enet_coordinate_descent_gram,
)
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.utils.fixes import _get_additional_lbfgs_options_dict
from sklearn.utils.optimize import _check_optimize_result


def min_norm_subgradient(alpha, gradient, coef, linear_loss):
    """Component-wise minimum norm subgradient of the loss function.

    The subgradient of f(x) + alpha ||x||_1 with the minimum L2 norm is defined
    component-wise, i.e. for component j:

        mn_sg(f(x) + alpha ||x||_1)_j =
            f'(x)_j + alpha                        if x_j > 0
            f'(x)_j - alpha                        if x_j < 0
            sign(f'(x)_j) max(|f'(x)| - alpha, 0)  if x_j = 0

    for positive alpha. See between equations (26) and (27) in
    Yuan, Ho, Lin (2011).

    Parameters
    ----------
    alpha : float
        The L1 penalty strength.

    gradient : ndarray of shape coef.shape
        The gradient of f(coef), i.e. f'(coef)

    coef : ndarray of shape (n_dof,) or (n_classes * n_dof,) or (n_classes, n_dof)

    linear_loss: LinearModelLoss
        Contains some meta information.

    Returns
    -------
    ndarray of shape gradient.shape

    References
    ----------
    - Yuan, G., Ho, C., & Lin, C. (2011). "An improved GLMNET for l1-regularized
      logistic regression." Journal of machine learning research.
      https://doi.org/10.1145/2020408.2020421
    """
    if alpha == 0:
        return gradient

    if gradient.shape != coef.shape or coef.ndim > 2:
        raise ValueError(
            "Shapes of 'gradient' and 'coef' must match, dimension must be smaller 3; "
            f"got {gradient.shape=} and {coef.shape=}."
        )
    if gradient.ndim == 1:
        result = gradient.copy()
    else:
        result = gradient.flatten(order="F")  # flatten returns a copy
        coef = coef.ravel(order="F")
    # mn_sgrad = minimum norm subgradient
    if linear_loss.fit_intercept:
        # For the intercept we just leave the gradient untouched.
        if linear_loss.base_loss.is_multiclass:
            n_classes = linear_loss.base_loss.n_classes
            mn_sgrad = result[:-n_classes]  # mn_sgrad is a view
        else:
            mn_sgrad = result[:-1]  # mn_sgrad is a view
    else:
        mn_sgrad = result
    weights, intercept = linear_loss.weight_intercept(coef)
    # For multiclass, weights.shape = (n_classes, n_features), flatten again.
    if linear_loss.base_loss.is_multiclass:
        weights = weights.ravel(order="F")
    mn_sgrad[weights > 0] += alpha
    mn_sgrad[weights < 0] -= alpha
    mask = weights == 0
    mn_sgrad[mask] = np.sign(mn_sgrad[mask]) * np.fmax(
        np.abs(mn_sgrad[mask]) - alpha, 0
    )
    if gradient.ndim == 1:
        return result
    else:
        return result.reshape(gradient.shape, order="F")


class NewtonSolver(ABC):
    """Newton solver for GLMs.

    This class implements Newton/2nd-order optimization routines for GLMs. Each Newton
    iteration aims at finding the Newton step which is done by the inner solver. With
    Hessian H, gradient g and coefficients coef, one step solves:

        H @ coef_newton = -g

    For our GLM / LinearModelLoss, we have gradient g and Hessian H:

        g = X.T @ loss.gradient + l2_reg_strength * coef
        H = X.T @ diag(loss.hessian) @ X + l2_reg_strength * identity

    Backtracking line search updates coef = coef_old + t * coef_newton for some t in
    (0, 1].

    This is a base class, actual implementations (child classes) may deviate from the
    above pattern and use structure specific tricks.

    Usage pattern:
        - initialize solver: sol = NewtonSolver(...)
        - solve the problem: sol.solve(X, y, sample_weight)

    References
    ----------
    - Jorge Nocedal, Stephen J. Wright. (2006) "Numerical Optimization"
      2nd edition
      https://doi.org/10.1007/978-0-387-40065-5

    - Stephen P. Boyd, Lieven Vandenberghe. (2004) "Convex Optimization."
      Cambridge University Press, 2004.
      https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf

    - Lee, J., Sun, Y., & Saunders, M.A. (2012). "Proximal Newton-Type Methods for
      Minimizing Composite Functions." SIAM J. Optim., 24, 1420-1443.
      https://doi.org/10.1137/130921428
      https://arxiv.org/abs/1206.1623

    - Yuan, G., Ho, C., & Lin, C. (2011). "An improved GLMNET for l1-regularized
      logistic regression." Journal of machine learning research.
      https://doi.org/10.1145/2020408.2020421

    Parameters
    ----------
    coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
        Initial coefficients of a linear model.
        If shape (n_classes * n_dof,), the classes of one feature are contiguous,
        i.e. one reconstructs the 2d-array via
        coef.reshape((n_classes, -1), order="F").

    linear_loss : LinearModelLoss
        The loss to be minimized.

    l1_reg_strength : float, default=0.0
        L1 regularization strength. This needs special care as it is non-smooth.

    l2_reg_strength : float, default=0.0
        L2 regularization strength.

    tol : float, default=1e-4
        The optimization problem is solved when each of the following condition is
        fulfilled:
        1. maximum |gradient| <= tol
        2. Newton decrement d: 1/2 * d^2 <= tol

    max_iter : int, default=100
        Maximum number of Newton steps allowed.

    n_threads : int, default=1
        Number of OpenMP threads to use for the computation of the Hessian and gradient
        of the loss function.

    Attributes
    ----------
    coef_old : ndarray of shape coef.shape
        Coefficient of previous iteration.

    coef_newton : ndarray of shape coef.shape
        Newton step.

    gradient : ndarray of shape coef.shape
        Gradient of the loss w.r.t. the coefficients.

    gradient_old : ndarray of shape coef.shape
        Gradient of previous iteration.

    loss_value : float
        Value of objective function = loss + penalty.

    loss_value_old : float
        Value of objective function of previous itertion.

    raw_prediction : ndarray of shape (n_samples,) or (n_samples, n_classes)

    converged : bool
        Indicator for convergence of the solver.

    iteration : int
        Number of Newton steps, i.e. calls to inner_solve

    use_fallback_lbfgs_solve : bool
        If set to True, the solver will resort to call LBFGS to finish the optimisation
        procedure in case of convergence issues.

    gradient_times_newton : float
        gradient @ coef_newton, set in inner_solve and used by line_search. If the
        Newton step is a descent direction, this is negative.
    """

    def __init__(
        self,
        *,
        coef,
        linear_loss=LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True),
        l1_reg_strength=0.0,
        l2_reg_strength=0.0,
        tol=1e-4,
        max_iter=100,
        n_threads=1,
        verbose=0,
    ):
        self.coef = coef
        self.linear_loss = linear_loss
        self.l1_reg_strength = l1_reg_strength
        self.l2_reg_strength = l2_reg_strength
        self.tol = tol
        self.max_iter = max_iter
        self.n_threads = n_threads
        self.verbose = verbose

    def setup(self, X, y, sample_weight):
        """Precomputations

        If None, initializes:
            - self.coef
        Sets:
            - self.raw_prediction
            - self.loss_value
        """
        _, _, self.raw_prediction = self.linear_loss.weight_intercept_raw(self.coef, X)
        self.loss_value = self.linear_loss.loss(
            coef=self.coef,
            X=X,
            y=y,
            sample_weight=sample_weight,
            l1_reg_strength=self.l1_reg_strength,
            l2_reg_strength=self.l2_reg_strength,
            n_threads=self.n_threads,
            raw_prediction=self.raw_prediction,
        )

    @abstractmethod
    def update_gradient_hessian(self, X, y, sample_weight):
        """Update gradient and Hessian."""

    @abstractmethod
    def inner_solve(self, X, y, sample_weight):
        """Compute Newton step.

        Sets:
            - self.coef_newton
            - self.gradient_times_newton
        """

    @abstractmethod
    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""

    def fallback_lbfgs_solve(self, X, y, sample_weight):
        """Fallback solver in case of emergency.

        If a solver detects convergence problems, it may fall back to this methods in
        the hope to exit with success instead of raising an error.

        Sets:
            - self.coef
            - self.converged
        """
        coef_shape = self.coef.shape
        self.coef = self.coef.ravel(order="F")  # scipy minimize expects 1d arrays
        max_iter = self.max_iter - self.iteration
        opt_res = scipy.optimize.minimize(
            self.linear_loss.loss_gradient,
            self.coef,
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": max_iter,
                "maxls": 50,  # default is 20
                "gtol": self.tol,
                "ftol": 64 * np.finfo(np.float64).eps,
                **_get_additional_lbfgs_options_dict("iprint", self.verbose - 1),
            },
            args=(X, y, sample_weight, 0, self.l2_reg_strength, self.n_threads),
        )
        self.iteration += _check_optimize_result("lbfgs", opt_res, max_iter=max_iter)
        self.coef = opt_res.x
        self.converged = opt_res.status == 0
        if len(coef_shape) > 1:
            self.coef = self.coef.reshape(coef_shape, order="F")

    def line_search(self, X, y, sample_weight):
        """Backtracking line search with Armijo condition.

        Backtracking line search with sufficient decrease condition, i.e. Armijo
        condition, is enough, no curvature condition needed, see Algorithm 3.1 of
        Nocedal & Wright 2nd ed.

        Sets:
            - self.coef_old
            - self.coef
            - self.loss_value_old
            - self.loss_value
            - self.gradient_old
            - self.gradient
            - self.raw_prediction
        """
        # line search parameters
        beta, sigma = 0.5, 0.00048828125  # 1/2, 1/2**11
        # Remember: dtype follows X, also the one of self.loss_value. For Array API
        # support, self.loss_value might be float instead of np.floatXX.
        eps = 16 * np.finfo(X.dtype).eps
        t = 1  # step size

        self.coef_old = self.coef
        self.loss_value_old = self.loss_value
        self.gradient_old = self.gradient

        # gradient_times_newton = self.gradient @ self.coef_newton
        # was computed in inner_solve.
        armijo_term = sigma * self.gradient_times_newton
        if self.l1_reg_strength > 0:
            # add ||coef||_1 - ||coef_old||_1, see Eq. 2.19 and 2.14 of
            # Lee, Sun, Saunder (2012) "Proximal Newton-Type Methods"
            L1 = self.linear_loss.l1_penalty(
                self.coef + self.coef_newton, self.l1_reg_strength
            )
            L1_old = self.linear_loss.l1_penalty(self.coef_old, self.l1_reg_strength)
            armijo_term += sigma * (L1 - L1_old)

        _, _, raw_prediction_newton = self.linear_loss.weight_intercept_raw(
            self.coef_newton, X
        )

        # np.sum(np.abs(self.gradient_old))
        sum_abs_grad_old = None
        g_max_abs_old = None

        is_verbose = self.verbose >= 2
        if is_verbose:
            print("  Backtracking Line Search")
            print(f"    eps=16 * finfo.eps={eps}")

        for i in range(21):  # until and including t = beta**20 ~ 1e-6
            self.coef = self.coef_old + t * self.coef_newton
            raw = self.raw_prediction + t * raw_prediction_newton
            self.loss_value, self.gradient = self.linear_loss.loss_gradient(
                coef=self.coef,
                X=X,
                y=y,
                sample_weight=sample_weight,
                l1_reg_strength=self.l1_reg_strength,
                l2_reg_strength=self.l2_reg_strength,
                n_threads=self.n_threads,
                raw_prediction=raw,
            )
            # Note: If coef_newton is too large, loss_gradient may produce inf values,
            # potentially accompanied by a RuntimeWarning.
            # This case will be captured by the Armijo condition.

            # 1. Check Armijo / sufficient decrease condition.
            # The smaller (more negative) the better.
            loss_improvement = self.loss_value - self.loss_value_old
            check = loss_improvement <= t * armijo_term
            if is_verbose:
                print(
                    f"    line search iteration={i + 1}, step size={t}\n"
                    f"      check loss improvement <= armijo term: {loss_improvement} "
                    f"<= {t * armijo_term} {check}"
                )
            if check:
                break
            # 2. Deal with too small loss improvements (Armijo not satisfied) that
            #    nevertheless result in a good reduction of the gradient norm which is
            #    our final stopping criterion.
            #    The author could not find justification in literature. It might be one
            #    of those improvements only practical implementations incorporate.
            #    The value of 1/8 = 0.125 is pure expert judgement, no derivation.
            if g_max_abs_old is None:
                mnsg_old = min_norm_subgradient(
                    self.l1_reg_strength,
                    self.gradient_old,
                    self.coef_old,
                    self.linear_loss,
                )
                sum_abs_grad_old = np.linalg.norm(mnsg_old.ravel(), ord=1)
                g_max_abs_old = np.linalg.norm(mnsg_old.ravel(), ord=np.inf)
            mnsg = min_norm_subgradient(
                self.l1_reg_strength,
                self.gradient,
                self.coef,
                self.linear_loss,
            )
            sum_abs_grad = np.linalg.norm(mnsg.ravel(), ord=1)
            g_max_abs = np.linalg.norm(mnsg.ravel(), ord=np.inf)
            check = (loss_improvement < 0) and (g_max_abs < 0.125 * g_max_abs_old)
            if is_verbose:
                print(
                    f"      check loss improvement < 0: {loss_improvement} < 0; and "
                    "max |gradient| <= 1/8 * max |gradient_old|: "
                    f"{g_max_abs} <= {0.125 * g_max_abs_old} {check}"
                )
            if check:
                break
            # 3. Deal with relative loss differences around machine precision.
            tiny_loss = np.abs(self.loss_value_old * eps)
            check = np.abs(loss_improvement) <= tiny_loss
            if is_verbose:
                print(
                    "      check loss |improvement| <= eps * |loss_old|:"
                    f" {np.abs(loss_improvement)} <= {tiny_loss} {check}"
                )
            if check:
                # 3.1 Check sum of absolute gradients as alternative condition.
                check = sum_abs_grad < sum_abs_grad_old
                if is_verbose:
                    print(
                        "      check sum(|gradient|) < sum(|gradient_old|): "
                        f"{sum_abs_grad} < {sum_abs_grad_old} {check}"
                    )
                if check:
                    break

            t *= beta
        else:
            warnings.warn(
                (
                    f"Line search of Newton solver {self.__class__.__name__} at"
                    f" iteration #{self.iteration} did not converge after 21 line "
                    "search refinement iterations. It will now resort to lbfgs instead."
                ),
                ConvergenceWarning,
            )
            if self.verbose:
                print("  Line search did not converge and resorts to lbfgs instead.")
            self.use_fallback_lbfgs_solve = True
            return

        self.raw_prediction = raw
        if is_verbose:
            print(
                f"    line search successful after {i + 1} iterations with "
                f"loss={self.loss_value}."
            )

    def check_convergence(self, X, y, sample_weight):
        """Check for convergence.

        Sets self.converged.
        """
        if self.verbose:
            print("  Check Convergence")
        # Note: Checking maximum relative change of coefficient <= tol is a bad
        # convergence criterion because even a large step could have brought us close
        # to the true minimum.
        # coef_step = self.coef - self.coef_old
        # change = np.max(np.abs(coef_step) / np.maximum(1, np.abs(self.coef_old)))
        # check = change <= tol

        # 1. Criterion: maximum |gradient| <= tol
        #    With L1 penalty: maximum |minimum-norm subgradient| <= tol
        #    Note: The gradient was already updated in line_search()
        #    Note: Yuan, Ho, Lin (2011) argue for the L1-norm instead of the
        #          infinity-norm. We do not find their argument that much convincing
        #          and, as above, use the infinity norm.
        mnsg = min_norm_subgradient(
            self.l1_reg_strength, self.gradient, self.coef, self.linear_loss
        )  # equals gradient if l1_reg_strength = 0.
        g_max_abs = np.linalg.norm(mnsg.ravel(), ord=np.inf)
        check = g_max_abs <= self.tol
        if self.verbose:
            if self.l1_reg_strength == 0:
                print(f"    1. max |gradient| {g_max_abs} <= {self.tol} {check}")
            else:
                print(
                    "    1. max |minimum L2-norm subgradient|"
                    f" {g_max_abs} <= {self.tol} {check}"
                )
        if not check:
            return

        # 2. Criterion: For Newton decrement d, check 1/2 * d^2 <= tol
        #       d = sqrt(grad @ hessian^-1 @ grad)
        #         = sqrt(coef_newton @ hessian @ coef_newton)
        #    See Boyd, Vanderberghe (2009) "Convex Optimization" Chapter 9.5.1. and
        #    Eq. 2.14, 2.15 of Lee, Sun, Saunder (2012) "Proximal Newton-Type Methods".
        d2 = self.compute_d2(X, sample_weight=sample_weight)
        check = 0.5 * d2 <= self.tol
        if self.verbose:
            print(f"    2. Newton decrement {0.5 * d2} <= {self.tol} {check}")
        if not check:
            return

        if self.verbose:
            loss_value = self.linear_loss.loss(
                coef=self.coef,
                X=X,
                y=y,
                sample_weight=sample_weight,
                l1_reg_strength=self.l1_reg_strength,
                l2_reg_strength=self.l2_reg_strength,
                n_threads=self.n_threads,
            )
            print(f"  Solver did converge at loss = {loss_value}.")
        self.converged = True

    def finalize(self, X, y, sample_weight):
        """Finalize the solvers results.

        Some solvers may need this, others not.
        """
        pass

    def solve(self, X, y, sample_weight):
        """Solve the optimization problem.

        This is the main routine.

        Order of calls:
            self.setup()
            while iteration:
                self.update_gradient_hessian()
                self.inner_solve()
                self.line_search()
                self.check_convergence()
            self.finalize()

        Returns
        -------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Solution of the optimization problem.
        """
        if self.verbose:
            print(self.__class__.__name__)
        # setup usually:
        #   - initializes self.coef if needed
        #   - initializes and calculates self.raw_predictions, self.loss_value
        self.setup(X=X, y=y, sample_weight=sample_weight)

        self.iteration = 1
        self.converged = False
        self.use_fallback_lbfgs_solve = False

        while self.iteration <= self.max_iter and not self.converged:
            if self.verbose:
                print(f"Newton iter={self.iteration}")

            self.use_fallback_lbfgs_solve = False  # Fallback solver.

            # 1. Update Hessian and gradient
            self.update_gradient_hessian(X=X, y=y, sample_weight=sample_weight)

            # TODO:
            # if iteration == 1:
            # We might stop early, e.g. we already are close to the optimum,
            # usually detected by zero gradients at this stage.

            # 2. Inner solver
            #    Calculate Newton step/direction
            #    This usually sets self.coef_newton and self.gradient_times_newton.
            self.inner_solve(X=X, y=y, sample_weight=sample_weight)
            if self.use_fallback_lbfgs_solve:
                break

            # 3. Backtracking line search
            #    This usually sets self.coef_old, self.coef, self.loss_value_old
            #    self.loss_value, self.gradient_old, self.gradient,
            #    self.raw_prediction.
            self.line_search(X=X, y=y, sample_weight=sample_weight)
            if self.use_fallback_lbfgs_solve:
                break

            # 4. Check convergence
            #    Sets self.converged.
            self.check_convergence(X=X, y=y, sample_weight=sample_weight)

            # 5. Next iteration
            self.iteration += 1

        if not self.converged:
            if self.use_fallback_lbfgs_solve:
                # Note: The fallback solver circumvents check_convergence and relies on
                # the convergence checks of lbfgs instead. Enough warnings have been
                # raised on the way.
                self.fallback_lbfgs_solve(X=X, y=y, sample_weight=sample_weight)
            else:
                warnings.warn(
                    (
                        f"Newton solver did not converge after {self.iteration - 1} "
                        "iterations."
                    ),
                    ConvergenceWarning,
                )

        self.iteration -= 1
        self.finalize(X=X, y=y, sample_weight=sample_weight)
        return self.coef


class NewtonCholeskySolver(NewtonSolver):
    """Cholesky based Newton solver.

    Inner solver for finding the Newton step H w_newton = -g uses Cholesky based linear
    solver.
    """

    def __init__(
        self,
        *,
        coef,
        linear_loss=LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True),
        l2_reg_strength=0.0,
        tol=1e-4,
        max_iter=100,
        n_threads=1,
        verbose=0,
    ):
        super().__init__(
            coef=coef,
            linear_loss=linear_loss,
            l1_reg_strength=0,
            l2_reg_strength=l2_reg_strength,
            tol=tol,
            max_iter=max_iter,
            n_threads=n_threads,
            verbose=verbose,
        )

    def setup(self, X, y, sample_weight):
        super().setup(X=X, y=y, sample_weight=sample_weight)
        if self.linear_loss.base_loss.is_multiclass:
            # Easier with ravelled arrays, e.g., for scipy.linalg.solve.
            # As with LinearModelLoss, we always are contiguous in n_classes.
            self.coef = self.coef.ravel(order="F")
            n_classes = self.linear_loss.base_loss.n_classes
        # Note that the computation of gradient in LinearModelLoss follows the shape of
        # coef.
        self.gradient = np.empty_like(self.coef)
        # But the hessian is always 2d.
        n = self.coef.size
        self.hessian = np.empty_like(self.coef, shape=(n, n))
        # To help case distinctions.
        self.is_multinomial_with_intercept = (
            self.linear_loss.base_loss.is_multiclass and self.linear_loss.fit_intercept
        )
        self.is_multinomial_no_penalty = (
            self.linear_loss.base_loss.is_multiclass and self.l2_reg_strength == 0
        )
        if self.is_multinomial_no_penalty:
            # See inner_solve. The provided coef might not adhere to the convention
            # that the last class is set to zero.
            # This is done by the usual freedom of a (overparametrized) multinomial to
            # add a constant to all classes which doesn't change predictions.
            coef = self.coef.reshape(n_classes, -1, order="F")  # easier as 2d
            coef -= coef[-1, :]  # coef -= coef of last class
        elif self.is_multinomial_with_intercept:
            # See inner_solve. Same as above, but only for the intercept.
            # intercept -= intercept of last class
            self.coef[-n_classes:] -= self.coef[-1]

        # Memory buffers for pointwise gradient and hessian.
        n_samples = X.shape[0]
        if self.linear_loss.base_loss.is_multiclass:
            self.grad_pointwise = np.empty_like(
                y, shape=(n_samples, n_classes), order="F"
            )
            self.hess_pointwise = np.empty_like(self.grad_pointwise)
        else:
            self.grad_pointwise = np.empty_like(y, order="F")
            self.hess_pointwise = np.empty_like(self.grad_pointwise)

    def update_gradient_hessian(self, X, y, sample_weight):
        _, _, self.hessian_warning = self.linear_loss.gradient_hessian(
            coef=self.coef,
            X=X,
            y=y,
            sample_weight=sample_weight,
            l2_reg_strength=self.l2_reg_strength,
            n_threads=self.n_threads,
            gradient_out=self.gradient,
            hessian_out=self.hessian,
            raw_prediction=self.raw_prediction,  # this was updated in line_search
            grad_pointwise_out=self.grad_pointwise,
            hess_pointwise_out=self.hess_pointwise,
        )

    def prepare_gradient_hessian(self):
        """Prepare gradient and hessian, in particular for multiclass case.

        This can't go into update_gradient_hessian because we need to keep the original
        self.gradient and self.hessian for line search and convergence checks.
        """
        # Note: The following case distinction could also be shifted to the
        # implementation of HalfMultinomialLoss instead of here within the solver.
        if self.is_multinomial_no_penalty:
            # The multinomial loss is overparametrized for each unpenalized feature, so
            # at least the intercepts. This can be seen by noting that predicted
            # probabilities are invariant under shifting all coefficients of a single
            # feature j for all classes by the same amount c:
            #   coef[k, :] -> coef[k, :] + c    =>    proba stays the same
            # where we have assumed coef.shape = (n_classes, n_features).
            # Therefore, also the loss (-log-likelihood), gradient and hessian stay the
            # same, see
            # Noah Simon and Jerome Friedman and Trevor Hastie. (2013) "A Blockwise
            # Descent Algorithm for Group-penalized Multiresponse and Multinomial
            # Regression". https://doi.org/10.48550/arXiv.1311.6529
            #
            # We choose the standard approach and set all the coefficients of the last
            # class to zero, for all features including the intercept.
            # Note that coef was already dealt with in setup.
            n_classes = self.linear_loss.base_loss.n_classes
            n_dof = self.coef.size // n_classes  # degree of freedom per class
            n = self.coef.size - n_dof  # effective size
            self.gradient[n_classes - 1 :: n_classes] = 0
            self.hessian[n_classes - 1 :: n_classes, :] = 0
            self.hessian[:, n_classes - 1 :: n_classes] = 0
            # We also need the reduced variants of gradient and hessian where the
            # entries set to zero are removed. For 2 features and 3 classes with
            # arbitrary values, "x" means removed:
            #   gradient = [0, 1, x, 3, 4, x]
            #
            #   hessian = [0,  1, x,  3,  4, x]
            #             [1,  7, x,  9, 10, x]
            #             [x,  x, x,  x,  x, x]
            #             [3,  9, x, 21, 22, x]
            #             [4, 10, x, 22, 28, x]
            #             [x,  x, x,  x, x,  x]
            # The following slicing triggers copies of gradient and hessian.
            gradient = self.gradient.reshape(-1, n_classes)[:, :-1].flatten()
            hessian = self.hessian.reshape(n_dof, n_classes, n_dof, n_classes)[
                :, :-1, :, :-1
            ].reshape(n, n)
        elif self.is_multinomial_with_intercept:
            # Here, only intercepts are unpenalized. We again choose the last class and
            # set its intercept to zero.
            # Note that coef was already dealt with in setup.
            self.gradient[-1] = 0
            self.hessian[-1, :] = 0
            self.hessian[:, -1] = 0
            gradient, hessian = self.gradient[:-1], self.hessian[:-1, :-1]
        else:
            gradient, hessian = self.gradient, self.hessian
        return gradient, hessian

    def inner_solve(self, X, y, sample_weight):
        if self.hessian_warning:
            warnings.warn(
                (
                    f"The inner solver of {self.__class__.__name__} detected a "
                    "pointwise hessian with many negative values at iteration "
                    f"#{self.iteration}. It will now resort to lbfgs instead."
                ),
                ConvergenceWarning,
            )
            if self.verbose:
                print(
                    "  The inner solver detected a pointwise Hessian with many "
                    "negative values and resorts to lbfgs instead."
                )
            self.use_fallback_lbfgs_solve = True
            return

        gradient, hessian = self.prepare_gradient_hessian()

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", scipy.linalg.LinAlgWarning)
                self.coef_newton = scipy.linalg.solve(
                    hessian, -gradient, check_finite=False, assume_a="sym"
                )
                if self.is_multinomial_no_penalty:
                    n_classes = self.linear_loss.base_loss.n_classes
                    n_dof = self.coef.size // n_classes  # degree of freedom per class
                    self.coef_newton = np.c_[
                        self.coef_newton.reshape(n_dof, n_classes - 1), np.zeros(n_dof)
                    ].reshape(-1)
                    assert self.coef_newton.flags.f_contiguous
                elif self.is_multinomial_with_intercept:
                    self.coef_newton = np.r_[self.coef_newton, 0]
                self.gradient_times_newton = self.gradient @ self.coef_newton
                if self.gradient_times_newton > 0:
                    if self.verbose:
                        print(
                            "  The inner solver found a Newton step that is not a "
                            "descent direction and resorts to LBFGS steps instead."
                        )
                    self.use_fallback_lbfgs_solve = True
                    return
        except (np.linalg.LinAlgError, scipy.linalg.LinAlgWarning) as e:
            warnings.warn(
                f"The inner solver of {self.__class__.__name__} stumbled upon a "
                "singular or very ill-conditioned Hessian matrix at iteration "
                f"{self.iteration}. It will now resort to lbfgs instead.\n"
                "Further options are to use another solver or to avoid such situation "
                "in the first place. Possible remedies are removing collinear features"
                " of X or increasing the penalization strengths.\n"
                "The original Linear Algebra message was:\n" + str(e),
                scipy.linalg.LinAlgWarning,
            )
            # Possible causes:
            # 1. hess_pointwise is negative. But this is already taken care in
            #    LinearModelLoss.gradient_hessian.
            # 2. X is singular or ill-conditioned
            #    This might be the most probable cause.
            #
            # There are many possible ways to deal with this situation. Most of them
            # add, explicitly or implicitly, a matrix to the hessian to make it
            # positive definite, confer to Chapter 3.4 of Nocedal & Wright 2nd ed.
            # Instead, we resort to lbfgs.
            if self.verbose:
                print(
                    "  The inner solver stumbled upon a singular or ill-conditioned "
                    "Hessian matrix and resorts to LBFGS instead."
                )
            self.use_fallback_lbfgs_solve = True
            return

    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""
        return self.coef_newton @ self.hessian @ self.coef_newton

    def finalize(self, X, y, sample_weight):
        if self.is_multinomial_no_penalty:
            # Our convention is usually the symmetric parametrization where
            # sum(coef[classes, features], axis=0) = 0.
            # We convert now to this convention. Note that it does not change
            # the predicted probabilities.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef = self.coef.reshape(n_classes, -1, order="F")
            self.coef -= np.mean(self.coef, axis=0)
        elif self.is_multinomial_with_intercept:
            # Only the intercept needs an update to the symmetric parametrization.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef[-n_classes:] -= np.mean(self.coef[-n_classes:])


class NewtonCDGramSolver(NewtonCholeskySolver):
    """Coordinate Descent Gram inner solver for a Newton solver.

    This solver can deal with L1 and L2 penalties.

    The inner solver for finding the Newton step H w_newton = -g uses coordinate
    descent:

        H @ coef_newton = -G

    With an L1 penalty, it is better to write down the minimization problem and use
    the 2nd order Taylor approximation only on the smooth parts (loss and L2), see
    Eq. 13 of Yuan, Ho, Lin (2011)

        min 1/2 d' H d + G' d + l1_reg ||coef + d||_1 - l1_reg ||coef||_1

    with d = coef_newton and coef = current coefficients.
    To circumvent the norm ||coef + d||_1, we instead minimized for
    c = coef_new = coef + d. This gives, up to constant terms

        min 1/2 c' H c + (G' - coef' H) c + l1_reg ||c||_1

    with
        c = coef_new = coef + coef_newton
        coef = current coefficients
        G = X.T @ g + l2_reg_strength * P @ coef
        H = X.T @ diag(h) @ X + l2_reg_strength * P
        g = loss.gradient = pointwise gradient
        h = loss.hessian = pointwise hessian
        P = penalty matrix in 1/2 w @ P @ w,
            for a pure L2 penalty without intercept it equals the identity matrix.

    This minimization problem is then solved by enet_coordinate_descent_gram.

    Note that this solver can naturally deal with sparse X.
    """

    def __init__(
        self,
        *,
        coef,
        linear_loss=LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True),
        l1_reg_strength=0.0,
        l2_reg_strength=0.0,
        tol=1e-4,
        max_iter=100,
        n_threads=1,
        verbose=0,
    ):
        super(NewtonCholeskySolver, self).__init__(
            coef=coef,
            linear_loss=linear_loss,
            l1_reg_strength=l1_reg_strength,
            l2_reg_strength=l2_reg_strength,
            tol=tol,
            max_iter=max_iter,
            n_threads=n_threads,
            verbose=verbose,
        )
        self.inner_tol = self.tol

    def setup(self, X, y, sample_weight):
        super().setup(X=X, y=y, sample_weight=sample_weight)
        self.is_multinomial_no_penalty = (
            self.linear_loss.base_loss.is_multiclass
            and self.l1_reg_strength == 0
            and self.l2_reg_strength == 0
        )

        # Memory buffers for Q_centered because it can be a large array.
        if self.linear_loss.fit_intercept:
            n_features = X.shape[1]
            n_classes = self.linear_loss.base_loss.n_classes
            if self.is_multinomial_no_penalty:
                n = (n_classes - 1) * n_features
            elif self.linear_loss.base_loss.is_multiclass:
                n = n_classes * n_features
            else:
                n = n_features
            self.Q_centered = np.empty(shape=(n, n), dtype=X.dtype, order="C")

    def update_gradient_hessian(self, X, y, sample_weight):
        # Same as NewtonCholesky but without the L2 penalty which is directly passed to
        # the CD solver.
        l2_reg = self.l2_reg_strength
        self.l2_reg_strength = 0
        super().update_gradient_hessian(X=X, y=y, sample_weight=sample_weight)
        self.l2_reg_strength = l2_reg

        if not self.linear_loss.base_loss.is_multiclass:
            # For non-canonical link functions and far away from the optimum, the
            # pointwise hessian can be negative.
            # We need non-negative hessians as we take the square root.
            # TODO: This should be done LinearModelLoss.gradient_hessian, as now we
            # have a slight inconsistency between hess_pointwise and the full hessian.
            np.maximum(0, self.hess_pointwise, out=self.hess_pointwise)

    def fallback_lbfgs_solve(self, X, y, sample_weight):
        if self.l1_reg_strength == 0:
            super().fallback_lbfgs_solve(X, y, sample_weight)
        else:
            # TODO(newton-cd): For the time being, we just honestly fail.
            cname = self.__class__.__name__
            msg = (
                f"This solver, {cname}, does not have a fallback for non-zero L1 "
                "penalties."
            )
            raise ConvergenceWarning(msg)

    def inner_solve(self, X, y, sample_weight):
        if self.hessian_warning:
            # TODO(newton-cd): Think about raising an error as there is no fallback.
            raise warnings.warn(
                (
                    f"The inner solver of {self.__class__.__name__} detected a "
                    "pointwise hessian with many negative values at iteration "
                    f"#{self.iteration}."
                ),
                ConvergenceWarning,
            )

        n_samples, n_features = X.shape
        gradient, hessian = self.prepare_gradient_hessian()

        if self.linear_loss.base_loss.is_multiclass:
            # Often needed variables for the multinomial.
            n_classes = self.linear_loss.base_loss.n_classes
            n_dof = self.coef.size // n_classes  # degree of freedom per class

        # Set w (coefficient passed to enet-cd) and Hcoef = H @ coef.
        if self.is_multinomial_no_penalty:
            # This is similar to penalized multinomial, but on top the last class per
            # feature was set to zero together with the corresponding elements in
            # gradient and hessian.
            coef = self.coef.reshape(-1, n_classes)[:, :-1].flatten()
            Hcoef = hessian @ coef
            if self.linear_loss.fit_intercept:
                # n_pc = number of effective coefficients, we keep the name from below.
                n_pc = (n_classes - 1) * (n_dof - 1)
            else:
                n_pc = (n_classes - 1) * n_dof
            w = coef[:n_pc].copy()
        elif self.is_multinomial_with_intercept:
            Hcoef = hessian @ self.coef[:-1]
            n_pc = n_classes * (n_dof - 1)  # number of penalized coefficients
            w = self.coef[:n_pc].copy()
        elif self.linear_loss.fit_intercept:
            Hcoef = hessian @ self.coef
            w = self.coef[:-1].copy()
        else:
            Hcoef = hessian @ self.coef
            w = self.coef.copy()

        # We minimize the 2. order Taylor approximation of the loss:
        #   1/2 c' H c + (G' - coef' H) c + l1_reg ||c||_1
        # c = coef_new = coef + coef_newton
        # Expressed with c = coef_new instead of coef_newton, such that L1 penalty fits
        # the formulation of our CD solver, which minimizes
        #   1/2 w' Q w - q' w + alpha ||w||_1
        # Neglecting intercepts and multiclass, this amounts to setting
        #   q = coef' H - G'= hessian @ coef - gradient
        #   Q = H
        if not self.linear_loss.fit_intercept:
            Q_centered = hessian
            q_centered = Hcoef - gradient
        elif not self.linear_loss.base_loss.is_multiclass:
            # We need to separate the intercept c0. CD solver solves
            #   min 1/2 w' Q w - q' w + alpha ||w||_1
            # Write
            #   Q_full = (Q  Q0')
            #            (Q0 Q00)
            # Optimizing only the intercept term w0 gives (unpenalized)
            # w0_optimal  = argmin_{w0} 1/2 w0' Q00 w0 - (q0 - Q0' w) w0
            #             = (q0 - Q0' w) / Q00
            # Inserting w0_optimal back into the objective (min instead of argmin)
            #   min_{w0} ... = -1/2 (q0 - Q0' w)' (q0 - Q0' w) / Q00
            # Added to the objective without intercept, we get
            #   obj = 1/2 w (Q - Q0 Q0' / Q00) w - (q - q0 Q0' / Q00) w + const
            Hcoef = hessian @ self.coef
            Q = hessian[:-1, :-1]  # shape (n_features, n_features)
            Q0 = hessian[-1, :-1]  # shape (n_features,)
            Q00 = hessian[-1, -1]  # float
            q = Hcoef[:-1] - gradient[:-1]
            q0 = Hcoef[-1] - gradient[-1]
            # Q_centered = Q - np.outer(Q0, Q0 / Q00)
            Q_centered = self.Q_centered  # use allocated memory
            np.copyto(dst=Q_centered, src=Q)
            Q_centered -= np.outer(Q0, Q0 / Q00)
            q_centered = q - q0 / Q00 * Q0
        else:
            # For the general treatment of the multinomial case, see
            # NewtonCholeskySolver.prepare_gradient_hessian.
            # In principle, this is the same case as fit_intercept above. Note:
            # - The intecept of the last class is set to zero and the corresponding
            #   element/row/column has been removed from the arrays gradient and
            #   hessian, see prepare_gradient_hessian.
            # - The step of minimizing the intercept terms alone is more complicated
            #   compared to above, because Q00 is now a 2-d array of shape
            #   (n_classes - 1, n_classes - 1). As it is usually a small array, we dare
            #   to compute its inverse (as it is used several times).
            Q = hessian[:n_pc, :n_pc]  # shape (n_pc, n_pc)
            Q0 = hessian[n_pc:, :n_pc]  # shape (n_classes - 1, n_pc)
            Q00 = hessian[n_pc:, n_pc:]  # shape (n_classes - 1, n_classes - 1)
            Q00_inv = np.linalg.pinv(Q00)
            q = Hcoef[:n_pc] - gradient[:n_pc]  # shape (n_pc,)
            q0 = Hcoef[n_pc:] - gradient[n_pc:]  # shape (n_classes - 1,)
            Q0t_Q00_inv = Q0.T @ Q00_inv
            # Q_centered = Q - Q0t_Q00_inv @ Q0
            Q_centered = self.Q_centered  # use allocated memory
            np.copyto(dst=Q_centered, src=Q)
            Q_centered -= Q0t_Q00_inv @ Q0
            q_centered = q - Q0t_Q00_inv @ q0

        # Which "y" to pass to enet_coordinate_descent_gram? It only effects the
        # stopping tolerance and the dual gap computation but not the coefficient
        # updates. A simple np.ones(n_samples) works. But we want to take advantage of
        # gap safe screening rules. Therefore, a correct dual gap is important.
        # For 1-d targets we need to write the objective as weighted least squares
        #   1/2 c' H c - q' c = 1/2 ||diag(sqrt(h))(y -  X c)||_2^2 + const
        # This gives
        #   y' diag(h) X c = q' c
        #   y = -g / h + X coef
        if self.linear_loss.base_loss.is_multiclass:
            # TODO(newton-cd): This is still the simple solution.
            y_cd = np.ones(shape=n_samples, dtype=X.dtype)
        else:
            y_cd = self.raw_prediction.copy()
            h_zero = self.hess_pointwise != 0
            y_cd[h_zero] -= self.grad_pointwise[h_zero] / self.hess_pointwise[h_zero]
            if self.linear_loss.fit_intercept:
                # As if applying _pre_fit(X, y_cd, ..).
                y_cd -= np.average(y_cd, axis=0, weights=self.hess_pointwise)
            y_cd *= np.sqrt(self.hess_pointwise)

        with warnings.catch_warnings():
            # Ignore warnings that add little information for users.
            warnings.simplefilter("ignore", ConvergenceWarning)
            w, gap, inner_tol, n_inner_iter = enet_coordinate_descent_gram(
                w=w,
                alpha=self.l1_reg_strength,
                beta=self.l2_reg_strength,
                Q=Q_centered,
                q=q_centered,
                y=y_cd,  # used in dual gap
                max_iter=1000,  # TODO(newton-cd): improve
                tol=self.inner_tol,
                rng=np.random.RandomState(0),
                random=False,
                positive=False,
                early_stopping=False,
            )

        # Set self.coef_newton and compute intercept terms.
        if n_inner_iter == 0:
            # Safeguard: if nothing changed nothing should change in this iter.
            # Without it, intercept terms might change.
            self.coef_newton = np.zeros_like(self.coef)
        elif self.is_multinomial_no_penalty:
            if self.linear_loss.fit_intercept:
                # Add intercept tems but for the last class
                intercepts = Q00_inv @ (q0 - Q0 @ w)
                w = np.r_[w, intercepts]
            # Add all zeros for the last class, intercept and coefficients.
            w = w.reshape(-1, n_classes - 1)
            w = np.c_[w, np.zeros(w.shape[0])].flatten()
            self.coef_newton = w - self.coef
        elif self.is_multinomial_with_intercept:
            intercepts = Q00_inv @ (q0 - Q0 @ w)
            self.coef_newton = np.r_[w, intercepts, 0] - self.coef
        elif self.linear_loss.fit_intercept:
            self.coef_newton = np.r_[w, (q0 - Q0 @ w) / Q00] - self.coef
        else:
            self.coef_newton = w - self.coef

        # Tighten inner stopping criterion, see Chapter 6.1 of "An Improved GLMNET".
        if n_inner_iter == 0:
            self.inner_tol *= 0.1
        elif n_inner_iter <= 4:  # "An Improved GLMNET" uses n_inner_iter <= 1
            self.inner_tol *= 0.25
        # elif n_inner_iter <= 8:
        #     # TODO(newton-cd): Check if this improves convergence
        #     self.inner_tol *= 0.5

        if self.verbose >= 2:
            print(
                f"  Inner coordinate descent solver stopped with {n_inner_iter} "
                f"iterations and a dualilty gap={float(gap)}."
            )

        if self.l2_reg_strength:
            # We neglected the l2_reg_strength in update_gradient_hessian. We correct it
            # now for the gradient.
            weights, _ = self.linear_loss.weight_intercept(self.coef)
            if self.linear_loss.base_loss.is_multiclass:
                # weights.shape = (n_classes, n_dof)
                self.gradient.reshape((n_classes, -1), order="F")[:, :n_features] += (
                    self.l2_reg_strength * weights
                )
            else:
                self.gradient[:n_features] += self.l2_reg_strength * weights

        self.gradient_times_newton = self.gradient @ self.coef_newton
        if self.gradient_times_newton > 0:
            if self.verbose:
                # TODO(newton-cd): What to do?
                print(
                    "  The inner solver found a Newton step that is not a "
                    "descent direction."
                )
            return

        return

    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""
        d2 = self.coef_newton @ self.hessian @ self.coef_newton
        if self.l2_reg_strength > 0:
            # We neglected the l2_reg_strength in update_gradient_hessian.
            weights, intercept = self.linear_loss.weight_intercept(self.coef_newton)
            d2 += 2 * self.linear_loss.l2_penalty(weights, self.l2_reg_strength)
        return d2
