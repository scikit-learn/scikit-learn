"""
Newton solver for Generalized Linear Models
"""

# Author: Christian Lorentzen <lorentzen.ch@gmail.com>
# License: BSD 3 clause

import warnings
from abc import ABC, abstractmethod

import numpy as np
import scipy.linalg
import scipy.optimize

from ..._loss.loss import HalfSquaredError
from ...exceptions import ConvergenceWarning
from ...utils.optimize import _check_optimize_result
from .._linear_loss import LinearModelLoss, Multinomial_LDL_Decomposition


class NewtonSolver(ABC):
    """Newton solver for GLMs.

    This class implements Newton/2nd-order optimization routines for GLMs. Each Newton
    iteration aims at finding the Newton step which is done by the inner solver. With
    Hessian H, gradient g and coefficients coef, one step solves:

        H @ coef_newton = -G

    For our GLM / LinearModelLoss, we have gradient G and Hessian H:

        G = X.T @ g + l2_reg_strength * P @ coef_old
        H = X.T @ diag(h) @ X + l2_reg_strength * P
        g = loss.gradient = pointwise gradient
        h = loss.hessian = pointwise hessian
        P = penalty matrix in 1/2 w @ P @ w,
            for a pure L2 penalty without intercept it equals the identity matrix.

    stemming from the 2nd order Taylor series of the loss:

        loss(coef_old) + g @ X @ coef_newton
        + 1/2 * coef_newton @ X.T @ diag(h) @ X @ coef_newton
        + 1/2 * l2_reg_strength * coef @ P @ coef

    In the last line, we have coef = coef_old + coef_newton.

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

    Parameters
    ----------
    coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
        Initial coefficients of a linear model.
        If shape (n_classes * n_dof,), the classes of one feature are contiguous,
        i.e. one reconstructs the 2d-array via
        coef.reshape((n_classes, -1), order="F").

    linear_loss : LinearModelLoss
        The loss to be minimized.

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
        l2_reg_strength=0.0,
        tol=1e-4,
        max_iter=100,
        n_threads=1,
        verbose=0,
    ):
        self.coef = coef
        self.linear_loss = linear_loss
        self.l2_reg_strength = l2_reg_strength
        self.tol = tol
        self.max_iter = max_iter
        self.n_threads = n_threads
        self.verbose = verbose

    def setup(self, X, y, sample_weight):
        """Precomputations

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

    def fallback_lbfgs_solve(self, X, y, sample_weight):
        """Fallback solver in case of emergency.

        If a solver detects convergence problems, it may fall back to this methods in
        the hope to exit with success instead of raising an error.

        Sets:
            - self.coef
            - self.converged
        """
        # Note that LinearModelLoss expects coef.ravel(order"F") for multiclass case.
        opt_res = scipy.optimize.minimize(
            self.linear_loss.loss_gradient,
            self.coef.ravel(order="F"),
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": self.max_iter,
                "maxls": 50,  # default is 20
                "iprint": self.verbose - 1,
                "gtol": self.tol,
                "ftol": 64 * np.finfo(np.float64).eps,
            },
            args=(X, y, sample_weight, self.l2_reg_strength, self.n_threads),
        )
        self.n_iter_ = _check_optimize_result("lbfgs", opt_res)
        self.coef = opt_res.x
        if self.linear_loss.base_loss.is_multiclass:
            # No test case was found (yet) where NewtonLSMRSolver ends up in this code
            # path. Note that NewtonCholeskySolver does not support muliclass problems.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef = self.coef.reshape((n_classes, -1), order="F")
        self.converged = opt_res.status == 0

    def line_search(self, X, y, sample_weight):
        """Backtracking line search.

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
        eps = 16 * np.finfo(self.loss_value.dtype).eps
        t = 1  # step size

        # gradient_times_newton = self.gradient @ self.coef_newton
        # was computed in inner_solve.
        armijo_term = sigma * self.gradient_times_newton
        _, _, raw_prediction_newton = self.linear_loss.weight_intercept_raw(
            self.coef_newton, X
        )

        self.coef_old = self.coef
        self.loss_value_old = self.loss_value
        self.gradient_old = self.gradient

        # np.sum(np.abs(self.gradient_old))
        sum_abs_grad_old = -1

        is_verbose = self.verbose >= 2
        if is_verbose:
            print("  Backtracking Line Search")
            print(f"    eps=10 * finfo.eps={eps}")

        for i in range(21):  # until and including t = beta**20 ~ 1e-6
            self.coef = self.coef_old + t * self.coef_newton
            raw = self.raw_prediction + t * raw_prediction_newton
            self.loss_value, self.gradient = self.linear_loss.loss_gradient(
                coef=self.coef,
                X=X,
                y=y,
                sample_weight=sample_weight,
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
                    f"    line search iteration={i+1}, step size={t}\n"
                    f"      check loss improvement <= armijo term: {loss_improvement} "
                    f"<= {t * armijo_term} {check}"
                )
            if check:
                break
            # 2. Deal with relative loss differences around machine precision.
            tiny_loss = np.abs(self.loss_value_old * eps)
            check = np.abs(loss_improvement) <= tiny_loss
            if is_verbose:
                print(
                    "      check loss |improvement| <= eps * |loss_old|:"
                    f" {np.abs(loss_improvement)} <= {tiny_loss} {check}"
                )
            if check:
                if sum_abs_grad_old < 0:
                    sum_abs_grad_old = scipy.linalg.norm(self.gradient_old, ord=1)
                # 2.1 Check sum of absolute gradients as alternative condition.
                sum_abs_grad = scipy.linalg.norm(self.gradient, ord=1)
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
                    f" iteration #{self.iteration} did no converge after 21 line search"
                    " refinement iterations. It will now resort to lbfgs instead."
                ),
                ConvergenceWarning,
            )
            if self.verbose:
                print("  Line search did not converge and resorts to lbfgs instead.")
            self.use_fallback_lbfgs_solve = True
            return

        self.raw_prediction = raw

    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""
        return self.coef_newton @ self.hessian @ self.coef_newton

    def check_convergence(self, X, y, sample_weight):
        """Check for convergence.

        Sets self.converged.
        """
        if self.verbose:
            print("  Check Convergence")
            loss_value = self.linear_loss.loss(
                coef=self.coef,
                X=X,
                y=y,
                sample_weight=sample_weight,
                l2_reg_strength=self.l2_reg_strength,
                n_threads=self.n_threads,
            )
            print(f"    loss = {loss_value}.")

        # Note: Checking maximum relative change of coefficient <= tol is a bad
        # convergence criterion because even a large step could have brought us close
        # to the true minimum.
        # coef_step = self.coef - self.coef_old
        # check = np.max(np.abs(coef_step) / np.maximum(1, np.abs(self.coef_old)))

        # 1. Criterion: maximum |gradient| <= tol
        #    The gradient was already updated in line_search()
        check = np.max(np.abs(self.gradient))
        if self.verbose:
            print(f"    1. max |gradient| {check} <= {self.tol}")
        if check > self.tol:
            return

        # 2. Criterion: For Newton decrement d, check 1/2 * d^2 <= tol
        #       d = sqrt(grad @ hessian^-1 @ grad)
        #         = sqrt(coef_newton @ hessian @ coef_newton)
        #    See Boyd, Vanderberghe (2009) "Convex Optimization" Chapter 9.5.1.
        d2 = self.compute_d2(X, sample_weight=sample_weight)
        if self.verbose:
            print(f"    2. Newton decrement {0.5 * d2} <= {self.tol}")
        if 0.5 * d2 > self.tol:
            return

        if self.verbose:
            loss_value = self.linear_loss.loss(
                coef=self.coef,
                X=X,
                y=y,
                sample_weight=sample_weight,
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
        # setup usually:
        #   - initializes self.coef if needed
        #   - initializes and calculates self.raw_predictions, self.loss_value
        self.setup(X=X, y=y, sample_weight=sample_weight)

        self.iteration = 1
        self.converged = False

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

    def setup(self, X, y, sample_weight):
        super().setup(X=X, y=y, sample_weight=sample_weight)
        n_dof = X.shape[1]
        if self.linear_loss.fit_intercept:
            n_dof += 1
        self.gradient = np.empty_like(self.coef)
        self.hessian = np.empty_like(self.coef, shape=(n_dof, n_dof))

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
        )

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

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", scipy.linalg.LinAlgWarning)
                self.coef_newton = scipy.linalg.solve(
                    self.hessian, -self.gradient, check_finite=False, assume_a="sym"
                )
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
                f"#{self.iteration}. It will now resort to lbfgs instead.\n"
                "Further options are to use another solver or to avoid such situation "
                "in the first place. Possible remedies are removing collinear features"
                " of X or increasing the penalization strengths.\n"
                "The original Linear Algebra message was:\n"
                + str(e),
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
                    "  The inner solver stumbled upon an singular or ill-conditioned "
                    "Hessian matrix and resorts to LBFGS instead."
                )
            self.use_fallback_lbfgs_solve = True
            return


class NewtonLSMRSolver(NewtonSolver):
    """LSMR based inexact Newton solver.

    The inner solver uses LSMR [1] after the Newton update is cast into the iteratively
    reweighted least squares (IRLS) formulation. This means

        H @ coef_newton = -G

    with

        G = X.T @ g + l2_reg_strength * P @ coef
        H = X.T @ diag(h) @ X + l2_reg_strength * P
        g = loss.gradient = pointwise gradient
        h = loss.hessian = pointwise hessian
        P = penalty matrix in 1/2 w @ P @ w,
            for a pure L2 penalty without intercept it equals the identity matrix.

    is cast as a least squares problem

        min ||A @ coef_newton - b||_2^2

    with

        A = [                    sqrt(h) * X]
            [sqrt(l2_reg_strength) * sqrt(P)]
        b = [                                - g / sqrt(h)]
            [- sqrt(l2_reg_strength) * sqrt(P) @ self.coef]

    Notes:
    - A is a square root of H: H = A'A
    - A and b form G: G = -A'b
    - The normal equation of this least squares problem, A'A coef_newton = A'b, is
      again H @ coef_newton = -G.

    After the Newton iteration, and neglecting line search, the new coefficient is
    self.coef += coef_newton.

    Note that this solver can naturally deal with sparse X.

    References
    ----------
    .. [1] :arxiv:`Fong & Saunders "LSMR: An iterative algorithm for sparse
           least-squares problems" <1006.0758>`
           See also
           https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html

    """  # noqa: E501

    def setup(self, X, y, sample_weight):
        """Setup.

        Additionally to super().setup(), also sets:
            - self.g = pointwise gradient
            - self.h = pointwise hessian
            - self.gradient
            - self.sqrt_P = sqrt(l2_reg_strength) * sqrt(P)
            - self.A_norm
            - self.r_norm
            - self.lsmr_iter
        """
        super().setup(X=X, y=y, sample_weight=sample_weight)
        (
            self.g,
            self.h,
        ) = self.linear_loss.base_loss.init_gradient_and_hessian(
            n_samples=X.shape[0], dtype=X.dtype
        )
        n_samples, n_features = X.shape
        # For a general (symmetric) penalty matrix P, we can use any square root of it,
        # e.g. the Cholesky decomposition. For the simple L2-penalty, we can use the
        # identity matrix and handle the intercept separately.
        # We use a 1-d array of shape (n_features,) instead of a 2-d diagonal array to
        # save memory. We also simply omit the zero at the end for the intercept.
        # For the multiclass case, we use the same 1-d array of shape (n_features,) and
        # do not store n_classes copies of it.
        self.sqrt_P = np.full(
            shape=n_features, fill_value=np.sqrt(self.l2_reg_strength), dtype=X.dtype
        )

        # In the inner_solve with LSMR, we set atol ~ 1 / (||A|| * ||r||) with the
        # Frobenius norm ||A||, see below. We track this in every iteration with
        # self.A_norm and self.r_norm. For the first call of inner_solve, we need an
        # initial estimation of ||A||. We assume h = 1 and get
        #   ||A||^2 = ||X||^2 + sqrt(l2_reg_strength) ||sqrt(P)||^2
        # if scipy.sparse.issparse(X):
        #     self.A_norm = scipy.sparse.linalg.norm(X) ** 2
        # else:
        #     self.A_norm = scipy.linalg.norm(X) ** 2
        # self.A_norm += scipy.linalg.norm(self.sqrt_P) ** 2
        # self.A_norm = np.sqrt(self.A_norm)
        #
        # Another, even simpler choice is ||A|| = 1 which likely understimates
        # ||A|| because ||X|| > 1 is very likely for typical (e.g. standardized) X.
        # Underestimating ||A|| means overestimating atol ~ 1 / ||A||, which means the
        # first iteration stops earlier. This is usually a good thing!
        # Note: For n_featues > n_samples and no regularization, it seems that a higher
        # starting A_norm is better which forces more initial LSMR iterations. We could
        # even condider to make it dependent on log(l2_reg_strength).
        self.A_norm = max(1, n_features / n_samples)
        if n_features > n_samples and self.l2_reg_strength < 1e-7:
            # We are in or close to an unpenalized underparametrized setting, where a
            # large first LSMR iteration is beneficial.
            # The multiplicative term is 1 <= term <= 50, log(1e-22) ~ -50
            self.A_norm *= 1 - np.log((self.l2_reg_strength + 1e-22) * 1e7)
        self.r_norm = 1
        self.lsmr_iter = 0  # number of total LSMR iterations

    def update_gradient_hessian(self, X, y, sample_weight):
        """Update gradient and hessian.

        Update pointwise gradient and hessian, self.g and self.h,
        as well as the full gradient, self.gradient.
        """
        n_features = X.shape[1]
        if not self.linear_loss.base_loss.is_multiclass:
            self.linear_loss.base_loss.gradient_hessian(
                y_true=y,
                raw_prediction=self.raw_prediction,  # this was updated in line_search
                sample_weight=sample_weight,
                gradient_out=self.g,
                hessian_out=self.h,
                n_threads=self.n_threads,
            )

            # See LinearModelLoss.gradient_hessian
            # For non-canonical link functions and far away from the optimum, the
            # pointwise hessian can be negative. We take care that 75% ot the hessian
            # entries are positive.
            self.hessian_warning = np.mean(self.h <= 0) > 0.25
            self.h = np.abs(self.h)

            # This duplicates a bit of code from LinearModelLoss.gradient.
            weights, _ = self.linear_loss.weight_intercept(self.coef)
            self.gradient = np.empty_like(self.coef, dtype=self.coef.dtype)
            self.gradient[:n_features] = X.T @ self.g + self.l2_reg_strength * weights
            if self.linear_loss.fit_intercept:
                self.gradient[-1] = self.g.sum()
        else:
            # Here we may safely assume HalfMultinomialLoss.
            # We use self.h to store the predicted probabilities instead of the hessian
            self.linear_loss.base_loss.gradient_proba(
                y_true=y,
                raw_prediction=self.raw_prediction,  # this was updated in line_search
                sample_weight=sample_weight,
                gradient_out=self.g,
                proba_out=self.h,
                n_threads=self.n_threads,
            )
            # Note: Multinomial hessian is always non-negative.

            # This duplicates a bit of code from LinearModelLoss.gradient.
            n_classes = self.linear_loss.base_loss.n_classes
            n_dof = n_features + int(self.linear_loss.fit_intercept)
            weights, _ = self.linear_loss.weight_intercept(self.coef)
            self.gradient = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
            # self.g.shape = (n_samples, n_classes)
            self.gradient[:, :n_features] = (
                self.g.T @ X + self.l2_reg_strength * weights
            )
            if self.linear_loss.fit_intercept:
                self.gradient[:, -1] = self.g.sum(axis=0)

    def compute_A_b(self, X, y, sample_weight):
        """Compute A and b for IRLS formulation.

        Returns
        -------
        A : LinearOperator of shape (n_samples + n_features, n_dof) or \
            (n_samples + n_features) * n_classes, n_dof * n_classes)
            For the multiclass case, `A @ x` expects `x = coef.ravel(order="C")` for
            `coef.shape = (n_classes, n_dof)`, and `A.T @ x` expects
            `x = b.ravel(order="F")` for `b.shape=(n_samples + n_features, n_classes)`.

        b : ndarray of shape (n_samples + n_features) or \
            ((n_samples + n_features) * n_classes)
        """
        n_samples, n_features = X.shape

        if not self.linear_loss.base_loss.is_multiclass:
            n_classes = 1
            sqrt_h = np.sqrt(self.h)

            b = np.r_[-self.g / sqrt_h, -self.sqrt_P * self.coef[:n_features]]

            if self.linear_loss.fit_intercept:
                n_dof = n_features + 1

                def matvec(x):
                    # A @ x with intercept
                    # We assume self.sqrt_P to be 1-d array of shape (n_features,),
                    # representing a diagonal matrix.
                    return np.r_[sqrt_h * (X @ x[:-1] + x[-1]), self.sqrt_P * x[:-1]]

                def rmatvec(x):
                    # A.T @ x with intercept
                    return np.r_[
                        X.T @ (sqrt_h * x[:n_samples]) + self.sqrt_P * x[n_samples:],
                        sqrt_h @ x[:n_samples],
                    ]

            else:
                n_dof = n_features

                def matvec(x):
                    # A @ x without intercept
                    return np.r_[sqrt_h * (X @ x), self.sqrt_P * x]

                def rmatvec(x):
                    # A.T @ x without intercept
                    return X.T @ (sqrt_h * x[:n_samples]) + self.sqrt_P * x[n_samples:]

        else:
            # Here we may safely assume HalfMultinomialLoss.
            n_classes = self.linear_loss.base_loss.n_classes
            p = self.h  # probability array of shape (n_samples, n_classes)
            # We need a square root of X' h X with pointwise hessian h. Note that h is
            # a priori a 4-dimensional matrix (n_samples ** 2, n_classes ** 2) and is
            # diagonal in n_samples ** 2, i.e. effectively a 3-dimensional matrix.
            # To accomplish this, we need the Cholesky or LDL' decomposition of h
            #   h = diag(p) - p' p
            # or with indices i and j for classes
            #   h_ij = p_i * delta_ij - p_i * p_j
            # This holds for every single point (n_sample). To tackle this, we use the
            # class Multinomial_LDL_Decomposition, which provides further details.
            #
            # We have h = L D L' with lower triangular L having unit diagonal and
            # diagonal D. This gives
            #   C = sqrt(D) L' X
            #   X' h X = C' C
            # and C is our searched for square root in
            #   A = [                              C]
            #       [sqrt(l2_reg_strength) * sqrt(P)]
            #   b = [                          -(L sqrt(D))^(-1) g]
            #       [- sqrt(l2_reg_strength) * sqrt(P) @ self.coef]

            # Note again
            # self.g.shape = p.shape = (n_samples, n_classes)
            # self.coef.shape = (n_classes, n_dof)
            # sqrt(h) becomes sqrt(D) L'
            LDL = Multinomial_LDL_Decomposition(proba=p)
            self.LDL = LDL  # store it for compute_d2

            # We need "g/sqrt(h)" such that g = L sqrt(D) "g/sqrt(h)", i.e.
            # "g/sqrt(h)" = 1/sqrt(D) L^-1 @ g.
            # Note that L^-1 is again lower triangular.
            g_over_h = self.g.copy()  # F-contiguous, shape (n_samples, n_classes)
            LDL.inverse_L_sqrt_D_matmul(g_over_h)
            if sample_weight is not None:
                # Sample weights sw enter the hessian multiplicatively, i.e. sw * h.
                # They can be dealt with by adding them to D, which means sqrt(sw)
                # to sqrt(D) in every place.
                sqrt_sw = np.sqrt(sample_weight)
                g_over_h /= sqrt_sw[:, None]

            # For ravelled results we use the convention that all values for the same
            # class are in sequence. For n_classes = 2, the ravelled b looks like
            #   [            -g_over_h[:, 0]]    n_samples elements of class 0
            #   [-self.sqrt_P * coef.T[:, 0]]    n_features elements of class 0
            #   [            -g_over_h[:, 1]]    n_samples elements of class 1
            #   [-self.sqrt_P * coef.T[:, 1]]    n_features elements of class 1
            # Note that this convention is different from
            # LinearModelLoss.gradient_hessian_product which expects raveled coef to
            # to have all classes in sequence.
            b = np.r_[
                -g_over_h,
                -self.sqrt_P[:, None] * self.coef[:, :n_features].T,
            ].ravel(order="F")

            # Note on performance:
            # X @ coef.T returns a C-contiguous array, but due to slicing along the
            # first axis (along samples like x[:, i]), sqrt_D_Lt_matmul is more
            # efficient for F-contiguous arrays. The output of X @ c can be made
            # F-contigous by (c.T @ X.T).T.

            if self.linear_loss.fit_intercept:
                n_dof = n_features + 1

                def matvec(x):
                    # A @ x with or without intercept
                    # We assume self.sqrt_P to be 1-d array of shape (n_features,),
                    # representing a diagonal matrix.
                    coef = x.reshape((n_classes, -1), order="C")
                    # C_coef = sqrt(D) L' X coef
                    # C_coef = X @ coef[:, :-1].T + coef[:, -1]  # C-contiguous
                    # F-contiguous version, see note above:
                    C_coef = (coef[:, :-1] @ X.T + coef[:, -1].T[:, None]).T
                    LDL.sqrt_D_Lt_matmul(C_coef)
                    # C_coef.shape = (n_samples, n_classes)
                    if sample_weight is not None:
                        C_coef *= sqrt_sw[:, None]
                    # For n_classes = 2, the ravelled result looks like
                    #   [         C_Coef[:, 0]]    n_samples elements of class 0
                    #   [sqrt_P * coef.T[:, 0]]    n_features elements of class 0
                    #   [         C_Coef[:, 1]]    n_samples elements of class 1
                    #   [sqrt_P * coef.T[:, 1]]    n_features elements of class 1
                    return np.r_[
                        C_coef,
                        self.sqrt_P[:, None] * coef[:, :n_features].T,
                    ].ravel(order="F")

                def rmatvec(x):
                    # A.T @ x with intercept
                    # Ct_y = X' L sqrt(D) y
                    x = x.reshape((-1, n_classes), order="F")
                    y = x[:n_samples, :]  # y is neither C nor F-contiguous
                    y = y.copy(order="F")
                    if sample_weight is not None:
                        y *= sqrt_sw[:, None]
                    L_sqrtD_y = LDL.L_sqrt_D_matmul(y)
                    Ct_y = X.T @ L_sqrtD_y  # shape = (n_features, n_classes)
                    return np.r_[
                        Ct_y + self.sqrt_P[:, None] * x[n_samples:, :],
                        np.sum(L_sqrtD_y, axis=0)[None, :],
                    ].ravel(order="F")

            else:
                n_dof = n_features

                def matvec(x):
                    # A @ x without intercept
                    coef = x.reshape((n_classes, -1), order="C")
                    # C_coef = sqrt(D) L' X coef
                    # C_coef = X @ coef.T  # C-contiguous
                    # F-contiguous version, see note above:
                    C_coef = (coef @ X.T).T
                    LDL.sqrt_D_Lt_matmul(C_coef)
                    if sample_weight is not None:
                        C_coef *= sqrt_sw[:, None]
                    return np.r_[
                        C_coef,
                        self.sqrt_P[:, None] * coef.T,
                    ].ravel(order="F")

                def rmatvec(x):
                    # A.T @ x without intercept
                    # Ct_y = X' L sqrt(D) y
                    x = x.reshape((-1, n_classes), order="F")
                    y = x[:n_samples, :]  # y is neither C nor F-contiguous
                    y = y.copy(order="F")
                    if sample_weight is not None:
                        y *= sqrt_sw[:, None]
                    L_sqrtD_y = LDL.L_sqrt_D_matmul(y)
                    Ct_y = X.T @ L_sqrtD_y  # shape = (n_features, n_classes)
                    return (Ct_y + self.sqrt_P[:, None] * x[n_samples:, :]).ravel(
                        order="F"
                    )

        # Note that initializing LinearOperator seems to have some surprisingly sizable
        # overhead.
        A = scipy.sparse.linalg.LinearOperator(
            shape=((n_samples + n_features) * n_classes, n_dof * n_classes),
            matvec=matvec,
            rmatvec=rmatvec,
        )

        return A, b

    def inner_solve(self, X, y, sample_weight):
        """Compute Newton step.

        Sets:
            - self.coef_newton via LSMR
              As LSMR is an iterative method, NewtonLSMRSolver is an inexact Newton
              method.
            - self.gradient_times_newton
            - self.A_norm
            - self.r_norm

        A_norm and r_norm are used for the inner tolerance used in LSMR.
        """
        n_samples, n_features = X.shape
        if not self.linear_loss.base_loss.is_multiclass:
            n_classes = 1
        else:
            n_classes = self.linear_loss.base_loss.n_classes
        A, b = self.compute_A_b(X, y, sample_weight)
        # The choice of atol in LSMR is essential for stability and for computation
        # time. For n_samples > n_features, we most certainly have a least squares
        # problem (no solution to the linear equation A x = b), such that the following
        # stopping criterion with residual r = b - A x and x = coef_newton applies:
        #   ||A' r|| = ||H x + G|| <= atol * ||A|| * ||r||.
        # For inexact Newton solvers, res = H x + G is called the residual and one
        # usually chooses, see Eq. 7.3 Nocedal & Wright 2nd ed, a stopping criterion
        #   ||res|| = ||A' r|| <= eta * ||G||
        # with a forcing sequence 0 < eta < 1 (eta_k for iteration k).
        # As for our Newton conjugate gradient solver "_newton_cg", we set
        #   eta = min(0.5, np.sqrt(||G||))
        # which establishes a superlinear rate of convergence.
        # Fortunately, we get approximations of the Frobenius norm of A and the norm of
        # r, ||A|| and ||r|| respectively, for free by LSMR such that we can set
        #   atol = eta * ||G|| / (||A|| * ||r||)
        # This results in our desired stopping criterion
        #   ||res|| = ||A' r|| <= eta * ||G||
        # at least approximately, as we have to use ||A|| and ||r|| from the last and
        # not the current iteration.
        #
        # In the case of perfect interpolation, n_features >= n_samples, we might have
        # an exact solution to A x = b. LSMR then uses the stopping criterion
        #   ||r|| <= atol * ||A|| * ||x|| + btol * ||b||.
        # Note that
        #   ||b||^2 = ||g / sqrt(h)||^2 + l2_reg_strength * x0 @ P @ x0
        # and
        #   1/2 * ||r||^2 = 1/2 * x X' diag(h) X x + g X x
        #                 + 1/2 * l2_reg_strength * (x + x0) P (x + x0)
        #                 + 1/2 *||g/sqrt(h)||^2
        # Here, x is the solution of Ax=b, i.e. x=coef_newton, and, neglecting line
        # search, coef = x + x0. Thus, 1/2 ||r||^2 is just the Taylor series of the
        # objective with the zero order term loss(x0) replaced by
        # 1/2 * ||g/sqrt(h)||^2.
        #
        # We set btol to the strict self.tol to help cases of collinearity in X when
        # there is (almost) no L2 penalty.
        norm_G = scipy.linalg.norm(self.gradient)
        eta = min(0.5, np.sqrt(norm_G))
        if self.verbose >= 3:
            print(f"    norm(gradient) = {norm_G}")
        result = scipy.sparse.linalg.lsmr(
            A,
            b,
            damp=0,
            atol=eta * norm_G / (self.A_norm * self.r_norm),
            btol=self.tol,
            maxiter=max(n_samples, n_features) * n_classes,  # default is min(A.shape)
            # default conlim = 1e8, for compatible systems 1e12 is still reasonable,
            # see LSMR documentation
            conlim=1e12,
            show=self.verbose >= 3,
        )
        # We store the estimated Frobenius norm of A and norm of residual r in
        # self.A_norm and self.r_norm for tolerance of next iteration.
        (
            self.coef_newton,
            istop,
            itn,
            self.r_norm,
            normar,
            self.A_norm,
            conda,
            normx,
        ) = result
        self.lsmr_iter += itn
        if self.verbose >= 2:
            print(f"  Inner iterations in LSMR = {itn}, total = {self.lsmr_iter}")
        if self.coef.dtype == np.float32:
            self.coef_newton = self.coef_newton.astype(np.float32)
        if not self.linear_loss.base_loss.is_multiclass:
            self.gradient_times_newton = self.gradient @ self.coef_newton
        else:
            # Center result such that sum over classes equals zero.
            self.coef_newton = self.coef_newton.reshape((n_classes, -1), order="C")
            # Note: The L2 penalty guarantees that coef_newton step is centered, i.e.
            # np.mean(self.coef_newton, axis=0) = 0.
            # But for the intercept term self.coef_newton[:, -1] this is less clear.
            # Nevertheless, all tests pass without the following line.
            # self.coef_newton -= np.mean(self.coef_newton, axis=0)

            self.gradient_times_newton = self.gradient.ravel(
                order="F"
            ) @ self.coef_newton.ravel(order="F")

        # In the first Newton iteraton, we tolerate istop == 7 and other things. We
        # just got started and are therefore forgiving.
        if self.iteration == 1:
            return
        # Note: We could detect too large steps by comparing norm(coef_newton) = normx
        # with norm(gradient) or with the already available condition number of A, e.g.
        # conda.
        if istop == 7:
            self.use_fallback_lbfgs_solve = True
            msg = (
                f"The inner solver of {self.__class__.__name__} reached "
                f"maxiter={itn} before the other stopping conditions were "
                f"satisfied at iteration #{self.iteration}. "
            )
        elif istop in (3, 6):
            self.use_fallback_lbfgs_solve = True
            msg = (
                f"The inner solver of {self.__class__.__name__} complained that the "
                f"condition number of A (A'A = Hessian), conda={conda}, seems to be "
                "greater than the given limit conlim=1e8 at iteration "
                f"#{self.iteration}."
            )

        if self.use_fallback_lbfgs_solve:
            warnings.warn(
                msg
                + "It will now resort to lbfgs instead.\n"
                "This may be caused by singular or very ill-conditioned Hessian "
                "matrix. "
                "Further options are to use another solver or to avoid such situation "
                "in the first place. Possible remedies are removing collinear features"
                "of X or increasing the penalization strengths.",
                ConvergenceWarning,
            )
            if self.verbose:
                print(
                    "  The inner solver had problems to converge and resorts to lbfgs."
                )
            self.use_fallback_lbfgs_solve = True
            return

    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""
        weights, intercept, raw_prediction = self.linear_loss.weight_intercept_raw(
            self.coef_newton, X
        )
        if not self.linear_loss.base_loss.is_multiclass:
            d2 = np.sum(raw_prediction * self.h * raw_prediction)
        else:
            d = self.LDL.sqrt_D_Lt_matmul(raw_prediction)
            if sample_weight is not None:
                d *= sample_weight[:, None]
            d = d.ravel(order="F")
            d2 = d @ d
        d2 += 2 * self.linear_loss.l2_penalty(weights, self.l2_reg_strength)
        return d2

    def finalize(self, X, y, sample_weight):
        if hasattr(self, "LDL"):
            delattr(self, "LDL")
