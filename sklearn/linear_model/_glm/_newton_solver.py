# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
Newton solver for Generalized Linear Models
"""

import math
import warnings
from abc import ABC, abstractmethod

import numpy as np
import scipy.linalg
import scipy.optimize

from sklearn._loss.loss import HalfSquaredError
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.utils.fixes import _get_additional_lbfgs_options_dict
from sklearn.utils.optimize import _check_optimize_result


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
            args=(X, y, sample_weight, self.l2_reg_strength, self.n_threads),
        )
        self.iteration += _check_optimize_result("lbfgs", opt_res, max_iter=max_iter)
        self.coef = opt_res.x
        self.converged = opt_res.status == 0
        if len(coef_shape) > 1:
            self.coef = self.coef.reshape(coef_shape, order="F")

    def line_search(self, X, y, sample_weight):
        """Backtracking line search.

        Define

            phi(alpha) = loss(coef_old + alpha * coef_newton)

        Searches for a step length alpha satisfying the sufficient decrease (Armijo)
        condition:

            phi(alpha) <= phi(0) + sigma * alpha * phi'(0)

        with sigma = 1/2 ** 11 ~ 5e-4.

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
        sigma = 0.00048828125  # 1/2**11, sometimes called c1
        min_step_reduction = 1e-2  # minimum factor of decrease of alpha per step
        min_step_length = 1e-12  # absolute minimum value of alpha
        # Remember: dtype follows X, also the dtype of self.loss_value. For Array API
        # support, self.loss_value might be float instead of np.floatXX.
        eps = 16 * np.finfo(X.dtype).eps
        alpha = 1  # initial step size, Newton methods should always try 1 first.
        alpha_old = 1

        # gradient_times_newton = self.gradient @ self.coef_newton
        # was computed in inner_solve.
        armijo_term = sigma * self.gradient_times_newton
        _, _, raw_prediction_newton = self.linear_loss.weight_intercept_raw(
            self.coef_newton, X
        )

        self.coef_old = self.coef
        self.loss_value_old = self.loss_value
        self.gradient_old = self.gradient
        phi_0 = phi_old = self.loss_value_old  # phi(0)
        phi_prime_0 = self.gradient_times_newton  # phi'(0) = gradient @ coef_newton

        # np.sum(np.abs(self.gradient_old))
        sum_abs_grad_old = -1

        is_verbose = self.verbose >= 2
        if is_verbose:
            print("  Backtracking Line Search")
            print(f"    eps=16 * finfo.eps={eps}")

        for i in range(21):
            self.coef = self.coef_old + alpha * self.coef_newton
            raw = self.raw_prediction + alpha * raw_prediction_newton
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
            check = loss_improvement <= alpha * armijo_term
            if is_verbose:
                print(
                    f"    line search iteration={i + 1}, step size={alpha}\n"
                    f"      check loss improvement <= armijo term: {loss_improvement} "
                    f"<= {alpha * armijo_term} {check}"
                )
            if check:
                break
            # 2. Tiny gradient / Armijo term.
            # If we are already close to the minimum, gradient and Armijo term are
            # tiny. It is best to use a Newton step length alpha = 1.
            if i == 0 and np.abs(armijo_term) <= self.tol:
                # Note that final convergence is checked with the infinity norm.
                g_max_abs = np.max(np.abs(self.gradient))
                check = g_max_abs <= self.tol
                if is_verbose:
                    print(
                        "      check max |gradient| <= tol: "
                        f"{g_max_abs} <= {self.tol} {check}"
                    )
                if check:
                    break
            # 3. Deal with differences around machine precision.
            # 3.1 Check relative loss difference ~ machine precision
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
                # 3.2 Check sum of absolute gradients as alternative condition.
                sum_abs_grad = scipy.linalg.norm(self.gradient, ord=1)
                check = sum_abs_grad < sum_abs_grad_old
                if is_verbose:
                    print(
                        "      check sum(|gradient|) < sum(|gradient_old|): "
                        f"{sum_abs_grad} < {sum_abs_grad_old} {check}"
                    )
                if check:
                    break

            # Set a smart new value of alpha, smaller than previous one, larger than 0.
            # We know that
            #   - phi(0) = phi_0
            #   - phi(alpha) = phi_1
            #   - phi'(0) = phi_prime_0
            # We fit a quadratic polynomial through those 3 points and take the minimum
            # alpha as next trial step length.
            # See Nocedal & Wright 2nd ed. Chapter 3.5, page 58, Eq 3.58.
            phi_1 = self.loss_value
            alpha_trial = (
                -phi_prime_0 * alpha**2 / (2 * (phi_1 - phi_0 - phi_prime_0 * alpha))
            )
            if i > 0:
                # We additionally know phi(alpha_old) = phi_old and use cubic
                # interpolation.
                # See Nocedal & Wright 2nd ed. Chapter 3.5, page 58, below Eq 3.58.
                denom = (alpha_old * alpha) ** 2 * (alpha - alpha_old)
                vec0 = phi_1 - phi_0 - phi_prime_0 * alpha
                vec1 = phi_old - phi_0 - phi_prime_0 * alpha_old
                a = (alpha_old**2 * vec0 - alpha**2 * vec1) / denom
                b = (alpha**3 * vec1 - alpha_old**3 * vec0) / denom
                if a != 0 and b**2 - 3 * a * phi_prime_0 >= 0:
                    alpha_trial = (-b + math.sqrt(b**2 - 3 * a * phi_prime_0)) / (3 * a)
                # else we keep the quadratic alpha_trial from above

            alpha_old = alpha
            phi_old = phi_1
            # Safeguards
            if alpha_trial >= alpha or alpha_trial <= min_step_length:
                alpha_trial = 0.5 * alpha
            # Avoid too large a reduction of alpha.
            alpha = max(alpha_trial, min_step_reduction * alpha, min_step_length)
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
        #    The gradient was already updated in line_search()
        g_max_abs = np.max(np.abs(self.gradient))
        check = g_max_abs <= self.tol
        if self.verbose:
            print(f"    1. max |gradient| {g_max_abs} <= {self.tol} {check}")
        if not check:
            return

        # 2. Criterion: For Newton decrement d, check 1/2 * d^2 <= tol
        #       d = sqrt(grad @ hessian^-1 @ grad)
        #         = sqrt(coef_newton @ hessian @ coef_newton)
        #    See Boyd, Vandenberghe (2009) "Convex Optimization" Chapter 9.5.1. and
        #    Eq. 20 of Yuan, Ho, Lin (2011).
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

    def setup(self, X, y, sample_weight):
        super().setup(X=X, y=y, sample_weight=sample_weight)
        if self.linear_loss.base_loss.is_multiclass:
            # Easier with ravelled arrays, e.g., for scipy.linalg.solve.
            # As with LinearModelLoss, we always are contiguous in n_classes.
            self.coef = self.coef.ravel(order="F")
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
            n_classes = self.linear_loss.base_loss.n_classes
            coef = self.coef.reshape(n_classes, -1, order="F")  # easier as 2d view
            coef -= coef[-1, :]  # coef -= coef of last class
        elif self.is_multinomial_with_intercept:
            # See inner_solve. Same as above, but only for the intercept.
            n_classes = self.linear_loss.base_loss.n_classes
            # intercept -= intercept of last class
            self.coef[-n_classes:] -= self.coef[-1]
        # Track linear algebra warnings
        self.has_already_warned = False

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

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", scipy.linalg.LinAlgWarning)
                self.coef_newton = scipy.linalg.solve(
                    hessian, -gradient, check_finite=False, assume_a="sym"
                )
                if self.is_multinomial_no_penalty:
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
            # Possible causes:
            # 1. hess_pointwise is negative. But this is already taken care of in
            #    LinearModelLoss.gradient_hessian.
            # 2. X is singular or ill-conditioned
            #    This might be the most probable cause.
            if not self.has_already_warned:
                # We only warn once.
                warnings.warn(
                    f"The inner solver of {self.__class__.__name__} stumbled upon a "
                    "singular or very ill-conditioned Hessian matrix at iteration "
                    f"{self.iteration}. It will now resort to using an "
                    "eigendecomposition of the Hessian, which is expensive to compute,"
                    " and only take large enough eigenvalues into consideration. For "
                    "tiny or negative eigenvalues, it will resort to lbfgs instead.\n"
                    "Further options are to use another solver or to avoid such a "
                    "situation in the first place. Possible remedies are removing "
                    "collinear features of X or increasing the penalization strengths."
                    "\n"
                    "The original Linear Algebra message was:\n" + str(e),
                    scipy.linalg.LinAlgWarning,
                )
            self.has_already_warned = True
            # There are many possible ways to deal with this situation. Most of them
            # add, explicitly or implicitly, a matrix to the hessian to make it
            # positive definite, confer to Chapter 3.4 of Nocedal & Wright 2nd ed.
            # We use an eigenvalue decomposition of the hessian and only use the
            # eigenvectors corresponding to large enough eigenvalues. The eigenvalues
            # correspond to curvature, negative curvature is bad (non-convex problem),
            # tiny eigenvalues mean flat space.
            eval, evec = scipy.linalg.eigh(self.hessian)
            max_eval = eval[-1]
            eps = np.sqrt(np.finfo(X.dtype).eps)
            if max_eval <= eps:
                # For flat space and negative curvature, we resort to lbfgs instead.
                if self.verbose:
                    print(
                        "  The inner solver stumbled upon a singular or ill-"
                        "conditioned Hessian matrix and resorts to LBFGS instead."
                    )
                self.use_fallback_lbfgs_solve = True
                return

            if self.verbose:
                print(
                    "  The inner solver stumbled upon a singular or ill-conditioned "
                    "Hessian matrix and resorts to an eigendecomposition of the "
                    "Hessian."
                )
            idx = np.flatnonzero(eval >= eps * np.abs(max_eval))
            # Eigendecomposition Hessian = Q D Q', Newton step = -Q D^(-1) Q' grad
            # As Hessian is symmetric, Q is an orthonormal matrix.
            Q = evec[:, idx]
            self.coef_newton = -Q @ ((1 / eval[idx]) * (Q.T @ self.gradient))
            self.gradient_times_newton = self.gradient @ self.coef_newton
            return

    def compute_d2(self, X, sample_weight):
        """Compute square of Newton decrement."""
        return self.coef_newton @ self.hessian @ self.coef_newton

    def finalize(self, X, y, sample_weight):
        if self.has_already_warned and not self.linear_loss.base_loss.is_multiclass:
            # X might be singular, try to find minimum norm solution.
            if self.linear_loss.fit_intercept:
                Xe = np.concatenate(
                    (X, np.ones((X.shape[0], 1), dtype=X.dtype)), axis=1
                )
            else:
                Xe = X
            eps = np.finfo(X.dtype).eps
            # As long as Xe @ coef remains constant, we are allowed to change coef.
            # This corresponds to the nullspace Z of Xe.
            Z = scipy.linalg.null_space(Xe, rcond=np.sqrt(eps))
            if (k := Z.shape[1]) > 0:
                # For any x of length Z.shape[1], we are free to add Z @ x to the
                # coefficients: coef_new = coef + Z @ x. We want to minimize
                # ||coef_new||_2^2 = ||coef + Z @ x||_2^2 which is a least squares.
                res = scipy.linalg.lstsq(Z, -self.coef)
                self.coef += Z @ res[0]

        if self.is_multinomial_no_penalty:
            # Our convention is usually the symmetric parametrization where
            # sum(coef[classes, features], axis=0) = 0.
            # We convert now to this convention. Note that it does not change
            # the predicted probabilities.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef = self.coef.reshape(n_classes, -1, order="F")  # easier as 2d
            self.coef -= np.mean(self.coef, axis=0)
        elif self.is_multinomial_with_intercept:
            # Only the intercept needs an update to the symmetric parametrization.
            n_classes = self.linear_loss.base_loss.n_classes
            self.coef[-n_classes:] -= np.mean(self.coef[-n_classes:])
