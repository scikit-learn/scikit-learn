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
from scipy.optimize.linesearch import line_search_BFGS, line_search_wolfe2

def newton_cg(func_grad_hess, func, x0, args=(), xtol=1e-5, eps=1e-4,
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
    xtol = len(x0)*avextol
    update = [2*xtol]
    xk = x0
    k = 0
    old_fval = None

    # Outer loop: our Newton iteration
    while (np.sum(np.abs(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        fval, grad, fhess_p = func_grad_hess(xk, *args)
        maggrad = np.sum(np.abs(grad))
        eta = min([0.5, np.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = np.zeros(len(x0), dtype=x0.dtype)
        ri = grad
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

        if old_fval is None:
            old_fval = fval

        alphak, fc, gc, old_fval = line_search_BFGS(func, xk, xsupi, grad,
                                                    old_fval, args=args)
        if alphak is None:
            # line search failed
            out = line_search_wolfe2(
                lambda x: func(x, *args),
                lambda x: func_grad_hess(x, *args)[1], xk, xsupi)
            alphak, fc, gc = out[0], out[1], out[2]
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