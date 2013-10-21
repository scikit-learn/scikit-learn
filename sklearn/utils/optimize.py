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
# Modifications by Gael Varoquaux
# License: BSD

import numpy as np
import warnings
from scipy.optimize.linesearch import line_search_wolfe2, line_search_wolfe1

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

def newton_cg(func_grad_hess, func, grad, x0, args=(), xtol=1e-5, eps=1e-4,
              maxiter=100, disp=False):
    """
    Minimization of scalar function of one or more variables using the
    Newton-CG algorithm.

    func: callable
        Should return the value of the function, the gradient, and a
        callable returning the matvec product of the Hessian
    """
    avextol = xtol

    x0 = np.asarray(x0).flatten()
    xtol = len(x0) * avextol
    update = [2 * xtol]
    xk = x0
    k = 0
    old_fval = func(x0, *args)
    old_old_fval = None

    # Outer loop: our Newton iteration
    while (np.sum(np.abs(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - fgrad f(xk) starting from 0.
        fval, fgrad, fhess_p = func_grad_hess(xk, *args)
        maggrad = np.sum(np.abs(fgrad))
        eta = min([0.5, np.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = np.zeros(len(x0), dtype=x0.dtype)
        ri = fgrad
        psupi = -ri
        i = 0
        dri0 = np.dot(ri, ri)

        # Inner loop: solve the Newton update by conjugate gradient, to
        # avoid inverting the Hessian
        while np.sum(np.abs(ri)) > termcond:
            Ap = fhess_p(psupi)
            # check curvature
            curv = np.dot(psupi, Ap)
            if 0 <= curv <= 3*np.finfo(np.float64).eps:
                break
            elif curv < 0:
                if (i > 0):
                    break
                else:
                    # fall back to steepest descent direction
                    xsupi = xsupi + dri0 / curv * psupi
                    break
            alphai = dri0 / curv
            xsupi = xsupi + alphai * psupi
            ri = ri + alphai * Ap
            dri1 = np.dot(ri, ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i = i + 1
            dri0 = dri1          # update np.dot(ri,ri) for next time.

        try:
            alphak, fc, gc, old_fval, old_old_fval, gfkp1 = \
                _line_search_wolfe12(func, grad, xk, xsupi, fgrad,
                                     old_fval, old_old_fval, args=args)
        except _LineSearchError:
            warnings.warn('Line Search failed')
            break

        update = alphak * xsupi
        xk = xk + update        # upcast if necessary
        k += 1

    return xk


###############################################################################
# Tests

if __name__ == "__main__":
    A = np.random.normal(size=(10, 10))

    def func(x):
        print 'Call to f: x %r' % x
        Ax = A.dot(x)
        return .5*(Ax).dot(Ax)

    def func_grad_hess(x):
        print 'Call to f_g_h: x %r' % x
        return func(x), A.T.dot(A.dot(x)), lambda x: A.T.dot(A.dot(x))

    x0 = np.ones(10)
    out = newton_cg(func_grad_hess, func, x0)