"""
Our own implementation of the Newton algorithm

Unlike the scipy.optimize version, this version of the Newton conjugate
gradient solver uses only one function call to retrieve the
func value, the gradient value and a callable for the Hessian matvec
product. If the function call is very expensive (e.g. for logistic
regression with large design matrix), this approach gives very
significant speedups.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# This is a modified file from scipy.optimize
# Original authors: Travis Oliphant, Eric Jones

import warnings

import scipy
from scipy.optimize._linesearch import (
    line_search_wolfe2,
    scalar_search_wolfe1,
)

from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._array_api import get_namespace_and_device, size


class _LineSearchError(RuntimeError):
    pass


# Copied from scipy
# https://github.com/scipy/scipy/blob/7a7fbca0b9baa1b709e4a5e0afaf9f94bd34941c/scipy/optimize/_linesearch.py#L37
# Modified for array API compliance: np.dot(a, b) -> a @ b
# TODO: use the `line_search_wolfe1` from `scipy` when it is array API compliant.
# Reference: https://github.com/scipy/scipy/pull/25022
def _line_search_wolfe1(
    f,
    fprime,
    xk,
    pk,
    gfk=None,
    old_fval=None,
    old_old_fval=None,
    args=(),
    c1=1e-4,
    c2=0.9,
    amax=50,
    amin=1e-8,
    xtol=1e-14,
):
    """
    Same as `scalar_search_wolfe1` but do a line search to direction `pk`
    """
    if gfk is None:
        gfk = fprime(xk, *args)

    gval = [gfk]
    gc = [0]
    fc = [0]

    def phi(s):
        fc[0] += 1
        return f(xk + s * pk, *args)

    def derphi(s):
        gval[0] = fprime(xk + s * pk, *args)
        gc[0] += 1
        return gval[0] @ pk

    derphi0 = gfk @ pk

    stp, fval, old_fval = scalar_search_wolfe1(
        phi,
        derphi,
        old_fval,
        old_old_fval,
        derphi0,
        c1=c1,
        c2=c2,
        amax=amax,
        amin=amin,
        xtol=xtol,
    )

    return stp, fc[0], gc[0], fval, old_fval, gval[0]


def _line_search_wolfe12(
    f, fprime, xk, pk, gfk, old_fval, old_old_fval, xp, device, verbose=0, **kwargs
):
    """
    Same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.

    Raises
    ------
    _LineSearchError
        If no suitable step size is found.

    """
    is_verbose = verbose >= 2
    eps = 16 * xp.finfo(xk.dtype).eps
    if is_verbose:
        print("  Line Search")
        print(f"    eps=16 * finfo.eps={eps}")
        print("    try line search wolfe1")

    ret = _line_search_wolfe1(f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs)

    if is_verbose:
        _not_ = "not " if ret[0] is None else ""
        print("    wolfe1 line search was " + _not_ + "successful")

    if ret[0] is None:
        # Have a look at the line_search method of our NewtonSolver class. We borrow
        # the logic from there
        # Deal with relative loss differences around machine precision.
        args = kwargs.get("args", tuple())
        fval = f(xk + pk, *args)
        tiny_loss = xp.abs(old_fval * eps)
        loss_improvement = fval - old_fval
        check = xp.abs(loss_improvement) <= tiny_loss
        if is_verbose:
            print(
                "    check loss |improvement| <= eps * |loss_old|:"
                f" {xp.abs(loss_improvement)} <= {tiny_loss} {check}"
            )
        if check:
            # 2.1 Check sum of absolute gradients as alternative condition.
            sum_abs_grad_old = scipy.linalg.norm(gfk, ord=1)
            grad = fprime(xk + pk, *args)
            sum_abs_grad = scipy.linalg.norm(grad, ord=1)
            check = sum_abs_grad < sum_abs_grad_old
            if is_verbose:
                print(
                    "    check sum(|gradient|) < sum(|gradient_old|): "
                    f"{sum_abs_grad} < {sum_abs_grad_old} {check}"
                )
            if check:
                ret = (
                    1.0,  # step size
                    ret[1] + 1,  # number of function evaluations
                    ret[2] + 1,  # number of gradient evaluations
                    fval,
                    old_fval,
                    grad,
                )

    if ret[0] is None:
        # line search failed: try different one.
        # TODO: It seems that the new check for the sum of absolute gradients above
        # catches all cases that, earlier, ended up here. In fact, our tests never
        # trigger this "if branch" here and we can consider to remove it.
        if is_verbose:
            print("    last resort: try line search wolfe2")
        ret = line_search_wolfe2(
            f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs
        )
        if is_verbose:
            _not_ = "not " if ret[0] is None else ""
            print("    wolfe2 line search was " + _not_ + "successful")

    if ret[0] is None:
        raise _LineSearchError()

    return ret


def _cg(fhess_p, fgrad, maxiter, tol, xp, device, verbose=0):
    """
    Solve iteratively the linear system 'fhess_p . xsupi = fgrad'
    with a conjugate gradient descent.

    Parameters
    ----------
    fhess_p : callable
        Function that takes the gradient as a parameter and returns the
        matrix product of the Hessian and gradient.

    fgrad : ndarray of shape (n_features,) or (n_features + 1,)
        Gradient vector.

    maxiter : int
        Number of CG iterations.

    tol : float
        Stopping criterion.

    Returns
    -------
    xsupi : ndarray of shape (n_features,) or (n_features + 1,)
        Estimated solution.
    """
    eps = 16 * xp.finfo(fgrad.dtype).eps
    xsupi = xp.zeros(size(fgrad), dtype=fgrad.dtype, device=device)
    ri = xp.asarray(fgrad, copy=True)  # residual = fgrad - fhess_p @ xsupi
    psupi = -ri
    i = 0
    dri0 = ri @ ri
    # We also keep track of |p_i|^2.
    psupi_norm2 = dri0
    is_verbose = verbose >= 2

    while i <= maxiter:
        if (norm1_re := xp.sum(xp.abs(ri))) <= tol:
            if is_verbose:
                print(
                    f"  Inner CG solver iteration {i} stopped with\n"
                    f"    sum(|residuals|) <= tol: {norm1_re} <= {tol}"
                )
            break

        Ap = fhess_p(psupi)
        # check curvature
        curv = psupi @ Ap
        if 0 <= curv <= eps * psupi_norm2:
            # See https://arxiv.org/abs/1803.02924, Algo 1 Capped Conjugate Gradient.
            if is_verbose:
                print(
                    f"  Inner CG solver iteration {i} stopped with\n"
                    f"    tiny_|p| = eps * ||p||^2, eps = {eps}, "
                    f"squared L2 norm ||p||^2 = {psupi_norm2}\n"
                    f"    curvature <= tiny_|p|: {curv} <= {eps * psupi_norm2}"
                )
            break
        elif curv < 0:
            if i > 0:
                if is_verbose:
                    print(
                        f"  Inner CG solver iteration {i} stopped with negative "
                        f"curvature, curvature = {curv}"
                    )
                break
            else:
                # fall back to steepest descent direction
                xsupi += dri0 / curv * psupi
                if is_verbose:
                    print("  Inner CG solver iteration 0 fell back to steepest descent")
                break
        alphai = dri0 / curv
        xsupi += alphai * psupi
        ri += alphai * Ap
        dri1 = ri @ ri
        betai = dri1 / dri0
        psupi = -ri + betai * psupi
        # We use  |p_i|^2 = |r_i|^2 + beta_i^2 |p_{i-1}|^2
        psupi_norm2 = dri1 + betai**2 * psupi_norm2
        i = i + 1
        dri0 = dri1  # update ri @ri for next time.
    if is_verbose and i > maxiter:
        print(
            f"  Inner CG solver stopped reaching maxiter={i - 1} with "
            f"sum(|residuals|) = {xp.sum(xp.abs(ri))}"
        )
    return xsupi


def _newton_cg(
    grad_hess,
    func,
    grad,
    x0,
    args=(),
    tol=1e-4,
    maxiter=100,
    maxinner=200,
    line_search=True,
    warn=True,
    verbose=0,
):
    """
    Minimization of scalar function of one or more variables using the
    Newton-CG algorithm.

    Parameters
    ----------
    grad_hess : callable
        Should return the gradient and a callable returning the matvec product
        of the Hessian.

    func : callable
        Should return the value of the function.

    grad : callable
        Should return the function value and the gradient. This is used
        by the linesearch functions.

    x0 : array-like of float
        Initial guess.

    args : tuple, default=()
        Arguments passed to func_grad_hess, func and grad.

    tol : float, default=1e-4
        Stopping criterion. The iteration will stop when
        ``max{|g_i | i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient.

    maxiter : int, default=100
        Number of Newton iterations.

    maxinner : int, default=200
        Number of CG iterations.

    line_search : bool, default=True
        Whether to use a line search or not.

    warn : bool, default=True
        Whether to warn when didn't converge.

    Returns
    -------
    xk : array-like of float
        Estimated minimum.
    """
    xp, _, device = get_namespace_and_device(x0)
    x0 = xp.asarray(x0, device=device)
    if x0.ndim != 1:
        msg = f"x0 must be 1-dimensional; got {x0.ndim=}"
        raise ValueError(msg)
    xk = xp.asarray(x0, copy=True)  # np.copy(x0)
    k = 0

    if line_search:
        old_fval = func(x0, *args)
        old_old_fval = None
    else:
        old_fval = 0

    is_verbose = verbose > 0

    # Outer loop: our Newton iteration
    while k < maxiter:
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - fgrad f(xk) starting from 0.
        fgrad, fhess_p = grad_hess(xk, *args)

        absgrad = xp.abs(fgrad)
        max_absgrad = xp.max(absgrad)
        check = max_absgrad <= tol
        if is_verbose:
            print(f"Newton-CG iter = {k}")
            print("  Check Convergence")
            print(f"    max |gradient| <= tol: {max_absgrad} <= {tol} {check}")
        if check:
            break

        maggrad = xp.sum(absgrad)
        eta = min([0.5, xp.sqrt(maggrad)])
        termcond = eta * maggrad

        # Inner loop: solve the Newton update by conjugate gradient, to
        # avoid inverting the Hessian
        xsupi = _cg(
            fhess_p,
            fgrad,
            maxiter=maxinner,
            tol=termcond,
            xp=xp,
            device=device,
            verbose=verbose,
        )

        alphak = 1.0

        if line_search:
            try:
                alphak, fc, gc, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                    func,
                    grad,
                    xk,
                    xsupi,
                    fgrad,
                    old_fval,
                    old_old_fval,
                    xp=xp,
                    device=device,
                    verbose=verbose,
                    args=args,
                )
            except _LineSearchError:
                warnings.warn("Line Search failed")
                break

        xk += alphak * xsupi  # upcast if necessary
        k += 1

    if warn and k >= maxiter:
        warnings.warn(
            (
                f"newton-cg failed to converge at loss = {old_fval}. Increase the"
                " number of iterations."
            ),
            ConvergenceWarning,
        )
    elif is_verbose:
        print(f"  Solver did converge at loss = {old_fval}.")
    return xk, k


def _check_optimize_result(solver, result, max_iter=None, extra_warning_msg=None):
    """Check the OptimizeResult for successful convergence

    Parameters
    ----------
    solver : str
       Solver name. Currently only `lbfgs` is supported.

    result : OptimizeResult
       Result of the scipy.optimize.minimize function.

    max_iter : int, default=None
       Expected maximum number of iterations.

    extra_warning_msg : str, default=None
        Extra warning message.

    Returns
    -------
    n_iter : int
       Number of iterations.
    """
    # handle both scipy and scikit-learn solver names
    if solver == "lbfgs":
        if max_iter is not None:
            # In scipy <= 1.0.0, nit may exceed maxiter for lbfgs.
            # See https://github.com/scipy/scipy/issues/7854
            n_iter_i = min(result.nit, max_iter)
        else:
            n_iter_i = result.nit

        if result.status != 0:
            warning_msg = (
                f"{solver} failed to converge after {n_iter_i} iteration(s) "
                f"(status={result.status}):\n"
                f"{result.message}\n"
            )
            # Append a recommendation to increase iterations only when the
            # number of iterations reaches the maximum allowed (max_iter),
            # as this suggests the optimization may have been prematurely
            # terminated due to the iteration limit.
            if max_iter is not None and n_iter_i == max_iter:
                warning_msg += (
                    f"\nIncrease the number of iterations to improve the "
                    f"convergence (max_iter={max_iter})."
                )
            warning_msg += (
                "\nYou might also want to scale the data as shown in:\n"
                "    https://scikit-learn.org/stable/modules/"
                "preprocessing.html"
            )
            if extra_warning_msg is not None:
                warning_msg += "\n" + extra_warning_msg
            warnings.warn(warning_msg, ConvergenceWarning, stacklevel=2)

    else:
        raise NotImplementedError

    return n_iter_i
