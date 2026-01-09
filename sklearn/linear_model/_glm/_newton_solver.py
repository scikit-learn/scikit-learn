# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
Newton solver for Generalized Linear Models
"""

import warnings
from abc import ABC, abstractmethod
from time import perf_counter

import numpy as np
import scipy.linalg
import scipy.optimize
from scipy import sparse

from sklearn._loss.loss import HalfSquaredError
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model._base import _pre_fit
from sklearn.linear_model._cd_fast import (
    enet_coordinate_descent,
    enet_coordinate_descent_gram,
    sparse_enet_coordinate_descent,
)
from sklearn.linear_model._linear_loss import (
    LinearModelLoss,
    Multinomial_LDL_Decomposition,
)
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
        # TODO: Pass correct y=b (see NewtonCDSolver), or at least y=||b||_2^2.

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


class NewtonCDSolver(NewtonSolver):
    """Coordinate Descent inner solver for a Newton solver.

    This solver can deal with L1 and L2 penalties.

    It avoids the explicit computation of hessian H = X' @ diag(h) @ X which makes it
    a good choice for use cases with n_features > n_samples (saves computation and
    saves large memory allocation of H).

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

    Up to constant terms, this is then cast as a Lasso/Enet least squares problem with
    sample weights equal to the pointwise hessian

        min 1/2 ||diag(sqrt(h)) (X c - z)||_2^2 + l1_reg ||c||_1 + 1/2 l2_reg ||c||_2^2

    with

        z = -g/h + X coef

    Note: Defining
        A = diag(sqrt(h)) X
        b = sqrt(h) z = -g / sqrt(h) + diag(sqrt(h)) X coef
    The first term is equivalent to 1/2 ||A c - b||_2^2.
    - A is a square root of H but without L2 penalty: A'A = X' diag(h) X
    - A and b form G - H coef (note: L2 cancels out): -A'b = G - H coef
    - The normal equation of this least squares problem, A'A coef_newton = A'b, is
      again H @ coef_newton = -G (without L2 penalty).

    This minimization problem is then solved by enet_coordinate_descent or
    sparse_enet_coordinate_descent.

    Note that this solver can naturally deal with sparse X.

    Attributes
    ----------
    This solver does not have the attribute
    - self.hessian

    Instead, it has
    - self.grad_pointwise
    - self.hess_pointwise
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
        super().__init__(
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
        if self.linear_loss.base_loss.is_multiclass and self.coef.ndim == 1:
            # NewtonCDSolver prefers 2-dim coef of shape (n_classes, n_dof). The
            # computation of gradient follows the shape of coef.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef = self.coef.reshape(n_classes, -1, order="F")

        self.gradient = np.empty_like(self.coef, order="F")
        self.grad_pointwise = np.empty_like(self.raw_prediction, order="F")
        self.hess_pointwise = np.empty_like(self.raw_prediction, order="F")
        self.sw_sum = X.shape[0] if sample_weight is None else np.sum(sample_weight)
        if sparse.issparse(X):
            if not sparse.isspmatrix_csc(X):
                ValueError(
                    f"X must be a CSC array/matrix for {self.__class__.__name__}"
                )
            if self.linear_loss.base_loss.is_multiclass:
                ValueError(
                    f"Solver {self.__class__.__name__} does not support multiclass "
                    "settings (n_classes >= 3)."
                )
        elif not X.flags.f_contiguous:
            ValueError(f"X must be F-contiguous for {self.__class__.__name__}")

    def update_gradient_hessian(self, X, y, sample_weight):
        """Update gradient and only pointwise hessian.

        Neglect l2_reg_strength because it is passed directly to the CD solver.
        This is later corrected for at the end of inner_solve.

        self.gradient (G)
        self.grad_pointwise (g)
        self.hess_pointwise (h)
        """
        # This duplicates a bit of code from LinearModelLoss.
        n_features = X.shape[1]
        if not self.linear_loss.base_loss.is_multiclass:
            _, _ = self.linear_loss.base_loss.gradient_hessian(
                y_true=y,
                raw_prediction=self.raw_prediction,  # this was updated in line_search
                sample_weight=sample_weight,
                gradient_out=self.grad_pointwise,
                hessian_out=self.hess_pointwise,
                n_threads=self.n_threads,
            )
            # For non-canonical link functions and far away from the optimum, the
            # pointwise hessian can be negative.
            np.maximum(0, self.hess_pointwise, out=self.hess_pointwise)

            self.grad_pointwise /= self.sw_sum
            self.hess_pointwise /= self.sw_sum
            self.gradient[:n_features] = X.T @ self.grad_pointwise
            if self.linear_loss.fit_intercept:
                self.gradient[-1] = self.grad_pointwise.sum()
        else:
            # We use self.hess_pointwise to store the predicted class probabilities.
            _, _ = self.linear_loss.base_loss.gradient_proba(
                y_true=y,
                raw_prediction=self.raw_prediction,  # this was updated in line_search
                sample_weight=sample_weight,
                gradient_out=self.grad_pointwise,
                proba_out=self.hess_pointwise,
                n_threads=self.n_threads,
            )
            self.grad_pointwise /= self.sw_sum
            self.gradient[:, :n_features] = self.grad_pointwise.T @ X
            if self.linear_loss.fit_intercept:
                self.gradient[:, -1] = self.grad_pointwise.sum(axis=0)

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
        n_samples, n_features = X.shape
        if not self.linear_loss.base_loss.is_multiclass:
            # z = self.raw_prediction - self.grad_pointwise / self.hess_pointwise
            z = self.raw_prediction.copy()
            h_zero = self.hess_pointwise != 0
            z[h_zero] -= self.grad_pointwise[h_zero] / self.hess_pointwise[h_zero]

            X, z, X_offset, z_offset, X_scale, _, _ = _pre_fit(
                X=X,
                y=z,
                Xy=None,
                precompute=False,
                fit_intercept=self.linear_loss.fit_intercept,
                copy=True,
                sample_weight=self.hess_pointwise,
            )

            if self.linear_loss.fit_intercept:
                w = self.coef[:-1].copy()
            else:
                w = self.coef.copy()

            with warnings.catch_warnings():
                # Ignore warnings that add little information for users.
                warnings.simplefilter("ignore", ConvergenceWarning)
                if sparse.issparse(X):
                    w, gap, inner_tol, n_inner_iter = sparse_enet_coordinate_descent(
                        w=w,
                        alpha=self.l1_reg_strength,
                        beta=self.l2_reg_strength,
                        X_data=X.data,
                        X_indices=X.indices,
                        X_indptr=X.indptr,
                        y=z,
                        sample_weight=self.hess_pointwise,
                        X_mean=np.asarray(X_offset / X_scale, dtype=X.dtype)
                        if X_offset is not None
                        else np.zeros(n_features, dtype=X.dtype),
                        max_iter=1000,  # TODO(newton-cd): improve
                        tol=self.inner_tol,
                        rng=np.random.RandomState(42),
                        random=False,
                        positive=False,
                        early_stopping=False,
                    )
                else:
                    w, gap, inner_tol, n_inner_iter = enet_coordinate_descent(
                        w=w,
                        alpha=self.l1_reg_strength,
                        beta=self.l2_reg_strength,
                        X=X,
                        y=z,
                        max_iter=1000,  # TODO(newton-cd): improve
                        tol=self.inner_tol,
                        rng=np.random.RandomState(42),
                        random=False,
                        positive=False,
                        early_stopping=False,
                    )
            # Set self.coef_newton and compute intercept terms.
            if n_inner_iter == 0:
                # Safeguard: if nothing changed nothing should change in this iter.
                # Without it, intercept terms might change.
                self.coef_newton = np.zeros_like(self.coef)
            elif self.linear_loss.fit_intercept:
                # Intercept treatment as in class ElasticNet, i.e. mean centering and
                # scaling.
                if X_scale is not None:
                    w /= X_scale
                w0 = z_offset - X_offset @ w
                self.coef_newton = np.r_[w, w0] - self.coef
            else:
                self.coef_newton = w - self.coef
        else:
            # Multinomial multiclass.
            # Unfortunately, the pointwise hessian h is not diagonal and we can't write
            # this as least squares plus penalties:
            #   ||sqrt(h) X w - y||_2^2 + penalties
            # One possible solution is to majorized the hessian by a diagonal matrix t,
            # see Lemma 3.3 in https://arxiv.org/abs/1311.6529:
            #   - one diagonal t_k for each class k
            #     diag(p) - p p' <= 2 diag(p (1 - p)) = t
            #   - identity times constant
            #     diag(p) - p p' <= max(2 diag(p (1 - p))) * identity = t
            # The full Hessian H would then be majorized by
            #     H <= X' diag(t) X  for each class k=1..K
            # Another strategy is pursued in Friedman, Hastie & Tibshirani (2010)
            # https://doi.org/10.18637/JSS.V033.I01. In section 4, they add another
            # middle loop over classes and optimized only for that class. This is
            # the same as using the diagonal majorization above and additionally
            # updating gradient and hessian in this middle loop.
            # Unfortunately, these strategies fail even for simple datasets such as
            #     X, y = make_classification(n_samples=20, n_features=20,
            #         n_informative=10, n_classes=3)
            #     LogisticRegression(C=1).fit(X, y)
            # Therefore, we take a more sophisticated approach:
            # Tanabe & Sagae (1992) https://doi.org/10.1111/J.2517-6161.1992.TB01875.X
            # derive an analytical LDL' decomposition for the matrix
            # h = diag(p) - p p' = LDL'.
            # Defining A = sqrt(D) L' X and b = we can now write
            #     1/2 c' H c + (G' - coef' H) c + alpha ||c||_1 + beta/2 ||c||_2^2
            #   = 1/2 ||A c - b||_2^2 + alpha ||c||_1 + beta/2 ||c||_2^2 + const
            #
            #   A' A = H
            #   G = X' g
            #   -c' A' b = c' (G - H coef)
            #   b = (L sqrt(D))^-1 (LDL' X coef - g)
            proba = self.hess_pointwise
            w = np.copy(self.coef, order="F")

            with warnings.catch_warnings():
                # Ignore warnings that add little information for users.
                warnings.simplefilter("ignore", ConvergenceWarning)
                _, gap, inner_tol, n_inner_iter = (
                    enet_coordinate_descent_multinomial_py(
                        W=w,
                        alpha=self.l1_reg_strength,
                        beta=self.l2_reg_strength,
                        X=X,
                        sample_weight=sample_weight,
                        raw_prediction=self.raw_prediction,
                        grad_pointwise=self.grad_pointwise,
                        proba=proba,
                        fit_intercept=self.linear_loss.fit_intercept,
                        max_iter=1000,  # TODO(newton-cd): improve
                        tol=self.inner_tol,
                        do_screening=True,
                        early_stopping=False,
                        verbose=max(0, self.verbose - 3),  # True for verbose >= 4
                    )
                )

            # Set self.coef_newton.
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
                self.gradient[:, :n_features] += self.l2_reg_strength * weights
            else:
                self.gradient[:n_features] += self.l2_reg_strength * weights

        self.gradient_times_newton = self.gradient.ravel(
            order="F"
        ) @ self.coef_newton.ravel(order="F")
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
        # return self.coef_newton @ self.hessian @ self.coef_newton
        weights, intercept, raw_prediction = self.linear_loss.weight_intercept_raw(
            self.coef_newton, X
        )
        if not self.linear_loss.base_loss.is_multiclass:
            d2 = np.sum(raw_prediction * self.hess_pointwise * raw_prediction)
        else:
            # hess_pointwise is predicted probability
            proba = self.hess_pointwise
            # coef' H coef = raw_prediction' * h * raw_prediction
            # = sum_{i,k,l} raw_{i,k} p_{i,k} (1_{k=l} - p_{i,l}) raw_{i,l}
            # = sum_{i,k} p_{i,k} raw_{i,k}^2
            # - sum_{i} (sum_{k} raw_{i,k} p_{i,k})^2
            proba_raw = proba * raw_prediction
            if sample_weight is None:
                d2 = np.sum(proba_raw * raw_prediction)
                d2 -= np.sum(np.sum(proba_raw, axis=1) ** 2)
            else:
                d2 = np.sum(sample_weight[:, None] * proba_raw * raw_prediction)
                d2 -= np.sum(sample_weight * np.sum(proba_raw, axis=1) ** 2)
            d2 /= self.sw_sum
        if self.l2_reg_strength > 0:
            d2 += 2 * self.linear_loss.l2_penalty(weights, self.l2_reg_strength)
        return d2

    def finalize(self, X, y, sample_weight):
        if (
            self.linear_loss.base_loss.is_multiclass
            and self.l1_reg_strength == 0
            and self.l2_reg_strength == 0
        ):
            # Our convention is usually the symmetric parametrization where
            # sum(coef[classes, features], axis=0) = 0.
            # We convert now to this convention. Note that it does not change
            # the predicted probabilities.
            n_classes = self.linear_loss.base_loss.n_classes
            # self.coef = self.coef.reshape(n_classes, -1, order="F")
            self.coef -= np.mean(self.coef, axis=0)
        elif (
            self.linear_loss.base_loss.is_multiclass and self.linear_loss.fit_intercept
        ):
            # Only the intercept needs an update to the symmetric parametrization.
            self.coef[:, -1] -= np.mean(self.coef[:, -1])


def gap_enet_multinomial(
    w,
    alpha,  # L1 penalty
    beta,  # L2 penalty
    X,
    y,
    R,  # current residual b - A @ w
    LD_R,  # LD_R = LDL.L_sqrt_D_matmul(R.copy()) * sqrt_sw[:, None]
    fit_intercept,
    H00_pinv_H0,
    XtA,  # XtA = X.T @ R - beta * w is calculated inplace
    gap_smaller_eps,
):
    """Compute dual gap for use in enet_coordinate_descent.

    Note that X must always be replaced by A = sqrt(D) L' X
    """
    l1_norm = 0.0
    l2_norm2 = 0.0
    eps = np.finfo(w.dtype).eps
    n_classes = R.shape[1]

    # l2_norm2 = w @ w
    if beta > 0:
        l2_norm2 = np.linalg.norm(w, ord="fro") ** 2
    if alpha > 0:
        l1_norm = np.sum(np.abs(w))
    R_norm2 = np.linalg.norm(R, ord="fro") ** 2
    # Ry = R @ y
    if not (alpha == 0 and beta == 0):
        Ry = R.ravel(order="F") @ y.ravel(order="F")

    primal = 0.5 * R_norm2 + alpha * l1_norm + 0.5 * beta * l2_norm2

    # XtA = X'R -> A'R = X' L sqrt(D) R
    XtA[:, :] = (X.T @ LD_R).T
    if fit_intercept:
        XtA -= (H00_pinv_H0.T @ LD_R[:, :-1].sum(axis=0)).reshape(
            n_classes, -1, order="F"
        )

    if alpha == 0:
        # ||X'R||_2^2
        dual_norm_XtA = np.linalg.norm(XtA, ord="fro") ** 2
        if beta == 0:
            # This is OLS, no dual gap available. Resort to first order condition
            #     X'R = 0
            #     gap = ||X'R||_2^2
            # Compare with stopping criterion of LSQR.
            gap = dual_norm_XtA
            return gap, dual_norm_XtA
        # This is Ridge regression with primal objective
        #     P(w) = 1/2 ||y - X w||_2^2 + beta/2 ||w||_2^2
        # We use the "alternative" dual formulation with alpha=0, see
        # https://github.com/scikit-learn/scikit-learn/issues/22836
        #     D(v) = -1/2 ||v||_2^2 - v'y - 1/(2 beta) |||X'v||_2^2
        # With v = Xw - y = -R (residual), the dual gap reads
        #     G = P(w) - D(-R)
        #       = ||R||_2^2  + beta/2 ||w||_2^2 - R'y + 1/(2 beta) ||X'R||_2^2
        dual = -0.5 * R_norm2 + Ry - 1 / (2 * beta) * dual_norm_XtA
    else:
        # XtA = X.T @ R - beta * w
        XtA -= beta * w
        dual_norm_XtA = np.max(np.abs(XtA))

        if dual_norm_XtA > alpha:
            const_ = alpha / dual_norm_XtA
        else:
            const_ = 1.0
        dual = -0.5 * (const_**2) * (R_norm2 + beta * l2_norm2)
        dual += const_ * Ry

    gap = primal - dual
    if gap_smaller_eps and abs(gap) <= 2 * eps * primal:
        gap = 0.0
    return gap, dual_norm_XtA


def enet_coordinate_descent_multinomial_py(
    W,
    alpha,
    beta,
    X,
    sample_weight,
    raw_prediction,
    grad_pointwise,
    proba,
    fit_intercept,
    max_iter,
    tol,
    do_screening=True,
    early_stopping=True,
    verbose=False,
):
    """Cython coordinate descent algorithm for Elastic-Net multinomial regression.

    See function enet_coordinate_descent.
    We minimize the primal

        P(w) = 1/2 ||y - X w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    But we replace the variables y and X such that P(w) corresponds to the 2nd order
    approximation of the multinomial loss

        1/2 w' H w + (G' - coef' H) w + alpha ||w||_1 + beta/2 ||w||_2^2

    - w = W.ravel(order="F")
    - G = gradient = X.T @ (proba - Y)  with Y_k = (y==k)
    - H = X' * (diag(p) - pp') * X, the full multinomial hessian
      (* is kind of a Kronecker multiplication)

    We use the analytical LDL decomposition of diag(p) - pp' = L D L' of
    Tanabe & Sagae (1992) to replace

        X -> A = sqrt(D) L' * X
        y -> b = (L sqrt(D))^-1 (LDL' * X coef - g)
               = A coef - (L sqrt(D))^-1 g

    Which gives (up to a constant)

        P(w) = 1/2 ||b - A w||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    As the matrix A has n_samples * n_features * n_classes * n_classes entries, we
    avoid to explicitly build it.

    Note:
        - H = A'A
        - A'b = (coef' H - G') w

    Intercepts
    ----------
    If intercepts w0 of shape (n_classes,) are included, we have

        P(w) = 1/2 ||b - A w - A0 w0||_2^2 + alpha ||w||_1 + beta/2 ||w||_2^2

    with A0 = sqrt(D) L' 1 and 1 = 1_{n_samples, n_classes} for the n_classes intercept
    columns.

    Minimizing for w0 gives

        w0 = (A0' A0)^(-1) (A0' b - A0' A w)

        P(w) = 1/2 ||(I - A0 (A0' A0)^(-1) A0') (b - A w)||_2^2 + ...
             = 1/2 ||(b - A0 H00^(-1) q0) - (A - A0 H00^(-1) H0) w||_2^2 + ...

    We therefore replace

        A -> A - A0 H00^(-1) H0 = sqrt(D) L' (X - 1 H00^(-1) H0)

        b -> b - A0 H00^(-1) q0 = sqrt(D) L' (X - 1 H00^(-1) H0) coef
                                - (L sqrt(D))^-1 (I - LDL' 1 H00^(-1) 1') g

    Note:
        - A0 = sqrt(D) L' 1_{n_samples, n_classes}
        - H00 = A0 A0', i.e. the part of the hessian for the the intercept alone
        - H0 = A0' A, i.e. the part of the hessian mixing intercepts with the features:
          H[-n_classes, :-n_classes]
        - q0 = A0' b = 1' (LDL' X coef - g)

    Tracking the Residual
    ---------------------
    Usually, CD tracks the residual R = y - X w -> b - Aw. In the multinomial case, it
    is advantageous to track the rotated residual instead:

        LD_R = L sqrt(D) R

    This makes the coordinate update simple and fast involving just

        X[:, j] @ LD_R[:, k]

    The update of the rotated residual involves basically multiplying by the full
    LDL = diag(p) - pp', i.e.

        L sqrt(D) R -= (w_new_kj - w_old_kj) * LDL[:, k] * X[:, j]

    with LDL[:, k] begin the k-th column of the LDL matrix, having shape
    (n_sampels, n_classes).

    Returns
    -------
    W : ndarray of shape (n_classes, n_features + fit_intercept)
        ElasticNet coefficients.
    gap : float
        Achieved dual gap.
    tol : float
        Equals input `tol` times `np.dot(y, y)`. The tolerance used for the dual gap.
    n_iter : int
        Number of coordinate descent iterations.
    """
    dtype = X.dtype
    time_pre = 0
    time_w = 0
    time_r = 0
    tic = perf_counter()

    # get the data information into easy vars
    n_samples = X.shape[0]
    n_features = X.shape[1]
    n_classes = W.shape[0]
    n_active = np.full(shape=n_classes, fill_value=n_features)
    X_is_sparse = sparse.issparse(X)

    if sample_weight is None:
        sw_sum = n_samples
        sw = np.full(fill_value=1 / sw_sum, shape=n_samples, dtype=dtype)
    else:
        sw_sum = np.sum(sample_weight)
        sw = sample_weight / sw_sum
    sqrt_sw = np.sqrt(sw)
    # Note: Full hessian H of LinearModelLoss.gradient_hessian(coef=W, X=X, y=y)
    # divides by sw_sum. So each LDL.XXX_matmul computation must be multiplied by
    # sqrt(samples_weight / sw_sum) = sqrt_sw, i.e. as if sqrt(D) would have
    # been multiplied by sqrt_sw.
    LDL = Multinomial_LDL_Decomposition(proba=proba)
    XtA = np.zeros((n_classes, n_features), dtype=dtype, order="C")
    R = np.empty((n_samples, n_classes), dtype=dtype, order="F")
    LD_R = np.empty_like(R, order="F")  # memory buffer for LDL.L_sqrt_D_matmul(R)
    sw_proba = sw[:, None] * proba  # precompute this one
    x_k = np.empty(n_samples, dtype=dtype)  # memory buffer for temporaries
    if fit_intercept:
        # memory buffer for temporaries
        xx = np.empty((n_samples, n_classes), dtype=dtype, order="F")

    W_original = W
    H00_pinv_H0 = None
    if fit_intercept:
        # See LinearModelLoss.gradient_hessian and NewtonCDGramSolver.inner_solve.
        # We set the intercept of the last class to zero, loops and shapes often
        # have n_classes - 1 instead of n_classes. The missing entries are all zeros,
        # but we often do not add those zeros explicitly.
        W0 = W[:, -1]
        W0[-1] = 0  # intercept of last class
        W = W[:, :-1]
        # H00 = H[-n_classes:-1, -n_classes:-1] = Hessian part of intercepts
        H00 = np.zeros((n_classes - 1, n_classes - 1), dtype=dtype)
        # H0 = H[-n_classes:-1, :-n_classes] = Hessian part mixing intercepts with
        # features
        H0 = np.zeros((n_classes - 1, n_classes * n_features), dtype=dtype)
        H0_coef = np.zeros((n_classes - 1,), dtype=dtype)  # 1' LDL raw_prediction
        h = x_k
        for k in range(n_classes - 1):
            # Diagonal terms h_kk.
            h[:] = proba[:, k] * (1 - proba[:, k]) * sw
            H00[k, k] = h.sum()
            H0[k, k::n_classes] = X.T @ h
            H0_coef[k] += (h * raw_prediction[:, k]).sum()
            # Off diagonal terms (in classes) hess_kl.
            for l in range(k + 1, n_classes - 1):
                # Upper triangle (in classes).
                h[:] = -proba[:, k] * proba[:, l] * sw
                H00[k, l] = h.sum()
                H00[l, k] = H00[k, l]
                H0[k, l::n_classes] = X.T @ h
                H0[l, k::n_classes] = H0[k, l::n_classes]
                H0_coef[k] += (h * raw_prediction[:, l]).sum()
                H0_coef[l] += (h * raw_prediction[:, k]).sum()
            # So far omitted term for last class where we still need it:
            l = n_classes - 1
            h[:] = -proba[:, k] * proba[:, l] * sw
            H0[k, l::n_classes] = X.T @ h
            H0_coef[k] += (h * raw_prediction[:, l]).sum()

        H00_pinv = np.linalg.pinv(H00, hermitian=True)
        # H0_coef is the same as:
        # H0_coef = H0 @ W.ravel(order="F") + H00 @ W0
        grad_p_sum = grad_pointwise[:, :-1].sum(axis=0)
        q0 = H0_coef - grad_p_sum
        # H00_pinv_H0 = H00_pinv @ H0  # shape (n_classes - 1, n_classes * n_features)
        H00_pinv_H0 = (H0.T @ H00_pinv.T).T  # make it F-contiguous

        # Centering X as
        #   X -= (H00_pinv_H0)[None, :]
        # is not an option because it would mean n_classes^2 copies of X.
        # Center raw_prediction: += -intercepts - 1 H00^(-1) H0) coef
        raw_prediction = raw_prediction - W0  # creates a copy
        raw_prediction[:, :-1] -= (H00_pinv_H0 @ W.ravel(order="F"))[None, :]
        # Center grad_pointwise: g -= LDL' 1 H00^(-1) 1' g
        t = H00_pinv @ grad_p_sum  # H00^(-1) 1' g
        for k in range(n_classes - 1):
            h = proba[:, k] * (1 - proba[:, k]) * sw
            grad_pointwise[:, k] -= h * t[k]
            for l in range(k + 1, n_classes - 1):
                h = -proba[:, k] * proba[:, l] * sw
                grad_pointwise[:, k] -= h * t[l]
                grad_pointwise[:, l] -= h * t[k]

    gap = tol + 1.0
    d_w_tol = tol
    n_iter = 0
    if alpha <= 0:
        do_screening = False
    if do_screening:
        # active_set maps [k, :n_active[k]] -> j
        active_set = np.empty((n_classes, n_features), dtype=np.uint32)
        excluded_set = np.empty((n_classes, n_features), dtype=np.uint8)

    # A w = sqrt(D) L' X w = sqrt(D) L' raw = sqrt_D_Lt_raw
    sqrt_D_Lt_raw = LDL.sqrt_D_Lt_matmul(raw_prediction.copy()) * sqrt_sw[:, None]

    # residual R = b - A w = -(L sqrt(D))^-1 g
    R[:, :] = LDL.inverse_L_sqrt_D_matmul(-grad_pointwise / sqrt_sw[:, None])

    # b = A w - (L sqrt(D))^-1 g
    b = sqrt_D_Lt_raw + R  # shape (n_samples, n_classes)

    # tol = tol * linalg.norm(Y, ord='fro') ** 2
    tol = tol * np.linalg.norm(b, ord="fro") ** 2

    # Rotated residual: L sqrt(D) R sqrt(sw)
    LD_R[:, :] = LDL.L_sqrt_D_matmul(R * sqrt_sw[:, None])
    # IMPORTANT NOTE: With fit_intercept=True, np.sum(LD_R, axis=0) == 0.

    # Check convergence before entering the main loop.
    gap, dual_norm_XtA = gap_enet_multinomial(
        w=W,
        alpha=alpha,
        beta=beta,
        X=X,
        y=b,
        R=R,
        LD_R=LD_R,
        fit_intercept=fit_intercept,
        H00_pinv_H0=H00_pinv_H0,
        XtA=XtA,
        gap_smaller_eps=False,
    )
    if early_stopping and gap <= tol:
        if verbose:
            print(
                f"  inner coordinate descent solver stops early with gap={float(gap)}"
                f"  <= tol={float(tol)}"
            )
        return np.asarray(W), gap, tol, 0

    # Compute squared norms of the columns of X.
    # Same as norm2_cols_X = np.square(X).sum(axis=0)
    # norm2_cols_X = np.einsum(
    #     "ij,ij->j", X, X, dtype=dtype, order="C"
    # )
    # X -> sqrt(D) L' X = A
    # sum_{i} X_{ij} X_{ij} -> sum_i X_{ij} LDL_i X_{ij}
    # These are just the diagonal elements of the full hessian H.
    norm2_cols_X = np.empty((n_classes, n_features), dtype=dtype, order="F")
    for k in range(n_classes):
        h = proba[:, k] * (1 - proba[:, k]) * sw
        if not X_is_sparse:
            norm2_cols_X[k, :] = np.einsum("ij, i, ij -> j", X, h, X)
        else:
            # norm2_cols_X[k, :] = [X[:, j] @ (h * X[:, j]) for j in range(n_features)]

            def Xj_h_Xj(X, h, j):
                startptr = X.indptr[j]
                endptr = X.indptr[j + 1]
                idx = X.indices[startptr:endptr]
                return X.data[startptr:endptr] @ (h[idx] * X.data[startptr:endptr])

            norm2_cols_X[k, :] = [Xj_h_Xj(X, h, j) for j in range(n_features)]
    if fit_intercept:
        # See Q_centered = Q - Q0.T @ Q00_inv @ Q0 in NewtonCDGramSolve.
        norm2_cols_X -= np.einsum("ji, ji -> i", H0, H00_pinv_H0).reshape(
            n_classes, -1, order="F"
        )
        # Guarantee norm2_cols_X >= 0:
        norm2_cols_X[norm2_cols_X < 0] = 0

    # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
    if do_screening:
        n_active[:] = 0
        for k in range(n_classes):
            for j in range(n_features):
                if norm2_cols_X[k, j] == 0:
                    W[k, j] = 0
                    excluded_set[k, j] = 1
                    continue
                Xj_theta = XtA[k, j] / max(alpha, dual_norm_XtA)  # X[:,j] @ dual_theta
                d_j = (1 - np.abs(Xj_theta)) / np.sqrt(norm2_cols_X[k, j] + beta)
                if d_j <= np.sqrt(2 * gap) / alpha:
                    # include feature j of class k
                    active_set[k, n_active[k]] = j
                    excluded_set[k, j] = 0
                    n_active[k] += 1
                else:
                    # R += w[j] * X[:,j] -> w[k, j] A[:, kj]
                    # LD_R += w[k, j] * LDL[:, k] * X[:, j]
                    if not fit_intercept:
                        if not X_is_sparse:
                            x_k[:] = W[k, j] * X[:, j] * sw_proba[:, k]
                        else:
                            x_k[:] = 0
                            idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                            x_k[idx] = (
                                W[k, j]
                                * X.data[X.indptr[j] : X.indptr[j + 1]]
                                * sw_proba[idx, k]
                            )
                        LD_R -= proba * x_k[:, None]
                        LD_R[:, k] += x_k
                    else:
                        jn = n_classes * j
                        xx[:, :-1] = -H00_pinv_H0[:, jn + k][None, :]
                        xx[:, -1] = 0
                        if not X_is_sparse:
                            xx[:, k] += X[:, j]
                        else:
                            idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                            xx[idx, k] += X.data[X.indptr[j] : X.indptr[j + 1]]
                        xx *= W[k, j] * sw_proba
                        LD_R -= proba * np.sum(xx, axis=1)[:, None]
                        LD_R += xx
                    W[k, j] = 0
                    excluded_set[k, j] = 1
        if np.sum(n_active) == 0:
            # We want to guarantee at least one CD step.
            n_active[:] = n_features
            active_set[:, :] = np.arange(n_features, dtype=np.uint32)[None, :]
            excluded_set[:, :] = 0

    time_pre += perf_counter() - tic

    for n_iter in range(max_iter):
        w_max = 0.0
        d_w_max = 0.0
        for k in range(n_classes):  # Loop over coordinates
            for j in range(n_active[k]):  # Loop over coordinates
                tic = perf_counter()
                if do_screening:
                    j = active_set[k, j]

                if norm2_cols_X[k, j] == 0.0:
                    continue
                w_kj = W[k, j]  # Store previous value

                # tmp = X[:,j].T @ (R + w_j * X[:,j])
                #    -> A[:,jk].T @ (R + w_kj * A[:,kj])
                # With precomputed LD_R: A'R = X' L sqrt(D) R
                tmp = X[:, j] @ LD_R[:, k] + w_kj * norm2_cols_X[k, j]
                # Note: With fit_intercept=True, np.sum(LD_R, axis=0) == 0. So the
                # additional term
                #   tmp -= H00_pinv_H0[:, n_classes * j + k].T @ np.sum(LD_R, axis=0)
                # is zero.

                W[k, j] = (
                    np.sign(tmp)
                    * max(np.abs(tmp) - alpha, 0)
                    / (norm2_cols_X[k, j] + beta)
                )
                time_w += perf_counter() - tic
                tic = perf_counter()

                if W[k, j] != w_kj:
                    # Update residual
                    # R -= (w[j] - w_j) * X[:,j] -> (w[kj] - w_kj) * A[:,kj]
                    if not fit_intercept:
                        # Update rotated residual LD_R instead of R
                        # L sqrt(D) R -= (w[kj] - w_kj) * LDL[:, k] * X[:, j]
                        if not X_is_sparse:
                            x_k[:] = (w_kj - W[k, j]) * X[:, j] * sw_proba[:, k]
                        else:
                            x_k[:] = 0
                            idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                            x_k[idx] = (
                                (w_kj - W[k, j])
                                * X.data[X.indptr[j] : X.indptr[j + 1]]
                                * sw_proba[idx, k]
                            )
                        # LD_R[:, k] += (1 - proba[:, k]) * x_k  # diagonal LDL
                        # for l in range(n_classes):
                        #     if l != k:
                        #         LD_R[:, l] -= proba[:, l] * x_k
                        # faster version:
                        LD_R -= proba * x_k[:, None]
                        LD_R[:, k] += x_k
                    else:
                        # Update rotated residual LD_R instead of R
                        xx[:, :-1] = -H00_pinv_H0[:, n_classes * j + k][None, :]
                        xx[:, -1] = 0
                        if not X_is_sparse:
                            xx[:, k] += X[:, j]
                        else:
                            idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                            xx[idx, k] += X.data[X.indptr[j] : X.indptr[j + 1]]
                        xx *= (w_kj - W[k, j]) * sw_proba
                        # L sqrt(D) R += LDL xx
                        # for l in range(n_classes):
                        #     LD_R[:, l] += (1 - proba[:, l]) * xx[:, l]
                        #     for m in range(l + 1, n_classes):
                        #         LD_R[:, l] -= proba[:, l] * xx[:, m]
                        #         LD_R[:, m] -= proba[:, m] * xx[:, l]
                        # faster version:
                        LD_R -= proba * np.sum(xx, axis=1)[:, None]
                        LD_R += xx

                time_r += perf_counter() - tic
                tic = perf_counter()

                # update the maximum absolute coefficient update
                d_w_j = abs(W[k, j] - w_kj)
                d_w_max = max(d_w_max, d_w_j)
                w_max = max(w_max, np.abs(W[k, j]))

        if w_max == 0.0 or d_w_max / w_max <= d_w_tol or n_iter == max_iter - 1:
            # The biggest coordinate update of this iteration was smaller
            # than the tolerance: check the duality gap as ultimate
            # stopping criterion.
            R = LDL.inverse_L_sqrt_D_matmul(LD_R / sqrt_sw[:, None])
            gap, dual_norm_XtA = gap_enet_multinomial(
                w=W,
                alpha=alpha,
                beta=beta,
                X=X,
                y=b,
                R=R,
                LD_R=LD_R,
                fit_intercept=fit_intercept,
                H00_pinv_H0=H00_pinv_H0,
                XtA=XtA,
                gap_smaller_eps=True,
            )
            if verbose:
                print(
                    f"  inner coordinate descent solver at iter {n_iter} "
                    f"checks for gap={float(gap)} <= tol={float(tol)} ? {gap <= tol}"
                )
            if gap <= tol:
                # return if we reached desired tolerance
                break

            # Gap Safe Screening Rules, see https://arxiv.org/abs/1802.07481, Eq. 11
            if do_screening:
                n_active[:] = 0
                for k in range(n_classes):
                    for j in range(n_features):
                        if excluded_set[k, j]:
                            continue
                        Xj_theta = XtA[k, j] / max(
                            alpha, dual_norm_XtA
                        )  # X[:,j] @ dual_theta
                        d_j = (1 - np.abs(Xj_theta)) / np.sqrt(
                            norm2_cols_X[k, j] + beta
                        )
                        if d_j <= np.sqrt(2 * gap) / alpha:
                            # include feature j of class k
                            active_set[k, n_active[k]] = j
                            excluded_set[k, j] = 0
                            n_active[k] += 1
                        else:
                            # R += w[j] * X[:,j] -> w[k, j] A[:, kj]
                            # LD_R += w[k, j] * LDL[:, k] * X[:, j]
                            if not fit_intercept:
                                if not X_is_sparse:
                                    x_k[:] = W[k, j] * X[:, j] * sw_proba[:, k]
                                else:
                                    x_k[:] = 0
                                    idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                                    x_k[idx] = (
                                        W[k, j]
                                        * X.data[X.indptr[j] : X.indptr[j + 1]]
                                        * sw_proba[idx, k]
                                    )
                                LD_R -= proba * x_k[:, None]
                                LD_R[:, k] += x_k
                            else:
                                jn = n_classes * j
                                xx[:, :-1] = -H00_pinv_H0[:, jn + k][None, :]
                                xx[:, -1] = 0
                                if not X_is_sparse:
                                    xx[:, k] += X[:, j]
                                else:
                                    idx = X.indices[X.indptr[j] : X.indptr[j + 1]]
                                    xx[idx, k] += X.data[X.indptr[j] : X.indptr[j + 1]]
                                xx *= W[k, j] * sw_proba
                                LD_R -= proba * np.sum(xx, axis=1)[:, None]
                                LD_R += xx
                            W[k, j] = 0
                            excluded_set[k, j] = 1
    else:
        # for/else, runs if for doesn't end with a `break`
        message = (
            "Objective did not converge. You might want to increase "
            "the number of iterations, check the scale of the "
            "features or consider increasing regularisation."
            f" Duality gap: {gap:.6e}, tolerance: {tol:.3e}"
        )
        warnings.warn(message, ConvergenceWarning)

    if fit_intercept:
        W0[:-1] = H00_pinv @ (q0 - H0 @ W.ravel(order="F"))
        # W0[-1] = 0  # was already set at the beginning.

    # print(
    #     f"{time_pre=:6.4f} {time_w=:6.4f} {time_r=:6.4f} "
    #     f"total={time_pre + time_w + time_r:6.4f}"
    # )
    return np.asarray(W_original), gap, tol, n_iter + 1
