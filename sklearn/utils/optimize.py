"""
Our own implementation of the Newton algorithm

Unlike the scipy.optimize version, this version of the Newton conjugate
gradient solver uses only one function call to retrieve the
func value, the gradient value and a callable for the Hessian matvec
product. If the function call is very expensive (e.g. for logistic
regression with large design matrix), this approach gives very
significant speedups.
"""
# This is a modified file from scipy.optimize
# Original authors: Travis Oliphant, Eric Jones
# Modifications by Gael Varoquaux, Mathieu Blondel and Tom Dupre la Tour
# License: BSD

import numpy as np
import warnings
from scipy.optimize.linesearch import line_search_wolfe2, line_search_wolfe1

from ..exceptions import ConvergenceWarning


class _LineSearchError(RuntimeError):
    pass


def _line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval,
                         **kwargs):
    """
    Same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.

    Raises
    ------
    _LineSearchError
        If no suitable step size is found

    """
    ret = line_search_wolfe1(f, fprime, xk, pk, gfk,
                             old_fval, old_old_fval,
                             **kwargs)

    if ret[0] is None:
        # line search failed: try different one.
        ret = line_search_wolfe2(f, fprime, xk, pk, gfk,
                                 old_fval, old_old_fval, **kwargs)

    if ret[0] is None:
        raise _LineSearchError()

    return ret


def _cg(fhess_p, fgrad, maxiter, tol):
    """
    Solve iteratively the linear system 'fhess_p . xsupi = fgrad'
    with a conjugate gradient descent.

    Parameters
    ----------
    fhess_p : callable
        Function that takes the gradient as a parameter and returns the
        matrix product of the Hessian and gradient

    fgrad : ndarray, shape (n_features,) or (n_features + 1,)
        Gradient vector

    maxiter : int
        Number of CG iterations.

    tol : float
        Stopping criterion.

    Returns
    -------
    xsupi : ndarray, shape (n_features,) or (n_features + 1,)
        Estimated solution
    """
    xsupi = np.zeros(len(fgrad), dtype=fgrad.dtype)
    ri = fgrad
    psupi = -ri
    i = 0
    dri0 = np.dot(ri, ri)

    while i <= maxiter:
        if np.sum(np.abs(ri)) <= tol:
            break

        Ap = fhess_p(psupi)
        # check curvature
        curv = np.dot(psupi, Ap)
        if 0 <= curv <= 3 * np.finfo(np.float64).eps:
            break
        elif curv < 0:
            if i > 0:
                break
            else:
                # fall back to steepest descent direction
                xsupi += dri0 / curv * psupi
                break
        alphai = dri0 / curv
        xsupi += alphai * psupi
        ri = ri + alphai * Ap
        dri1 = np.dot(ri, ri)
        betai = dri1 / dri0
        psupi = -ri + betai * psupi
        i = i + 1
        dri0 = dri1          # update np.dot(ri,ri) for next time.

    return xsupi


def newton_cg(grad_hess, func, grad, x0, args=(), tol=1e-4,
              maxiter=100, maxinner=200, line_search=True, warn=True):
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

    x0 : array of float
        Initial guess.

    args: tuple, optional
        Arguments passed to func_grad_hess, func and grad.

    tol : float
        Stopping criterion. The iteration will stop when
        ``max{|g_i | i = 1, ..., n} <= tol``
        where ``g_i`` is the i-th component of the gradient.

    maxiter : int
        Number of Newton iterations.

    maxinner : int
        Number of CG iterations.

    line_search: boolean
        Whether to use a line search or not.

    warn: boolean
        Whether to warn when didn't converge.

    Returns
    -------
    xk : ndarray of float
        Estimated minimum.
    """
    x0 = np.asarray(x0).flatten()
    xk = x0
    k = 0

    if line_search:
        old_fval = func(x0, *args)
        old_old_fval = None

    # Outer loop: our Newton iteration
    while k < maxiter:
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - fgrad f(xk) starting from 0.
        fgrad, fhess_p = grad_hess(xk, *args)

        absgrad = np.abs(fgrad)
        if np.max(absgrad) < tol:
            break

        maggrad = np.sum(absgrad)
        eta = min([0.5, np.sqrt(maggrad)])
        termcond = eta * maggrad

        # Inner loop: solve the Newton update by conjugate gradient, to
        # avoid inverting the Hessian
        xsupi = _cg(fhess_p, fgrad, maxiter=maxinner, tol=termcond)

        alphak = 1.0

        if line_search:
            try:
                alphak, fc, gc, old_fval, old_old_fval, gfkp1 = \
                    _line_search_wolfe12(func, grad, xk, xsupi, fgrad,
                                         old_fval, old_old_fval, args=args)
            except _LineSearchError:
                warnings.warn('Line Search failed')
                break

        xk = xk + alphak * xsupi        # upcast if necessary
        k += 1

    if warn and k >= maxiter:
        warnings.warn("newton-cg failed to converge. Increase the "
                      "number of iterations.", ConvergenceWarning)
    return xk, k
