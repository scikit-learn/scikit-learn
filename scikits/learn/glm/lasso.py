# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr> 
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg
from lasso_cd import lasso_coordinate_descent as lasso_coordinate_descent_slow

"""Trying to improve speed with cython... no success for now :(
"""
try:
    from lasso_cd_fast import lasso_coordinate_descent as lasso_coordinate_descent_fast
    lasso_coordinate_descent = lasso_coordinate_descent_fast
except ImportError, e:
    lasso_coordinate_descent = lasso_coordinate_descent_slow
    print "Using Python version of coordinate descent"

def enet_dual_gap(X, y, w, alpha, beta=0):
    """Compute dual gap for Elastic-Net model
        to check KKT optimality conditions
        returns gap, primal_objective, dual_objective
        gap should be positive
        gap = primal_objective - dual_objective
        gap < 1e-6 means convergence in practice
    """
    Xw = np.dot(X,w)
    A = (y - Xw)
    if beta > 0:
        B = - np.sqrt(beta*n)*w
    XtA = np.dot(X.T,A)
    if beta > 0:
        XtA += np.sqrt(beta*n) * B
    dual_norm_XtA = np.max(XtA)
    if (dual_norm_XtA > alpha):
        A *= alpha / dual_norm_XtA
        if beta > 0:
            B *= alpha / dual_norm_XtA
    pobj = 0.5 * linalg.norm(y - Xw)**2 + alpha * np.abs(w).sum() + 0.5 * beta * linalg.norm(w)**2
    dobj = - 0.5 * linalg.norm(A)**2 + np.dot(A.T, y)
    if beta > 0:
        dobj += - 0.5 * linalg.norm(B)**2
    gap = pobj - dobj
    return gap, pobj, dobj

class LassoCD(object):
    """docstring for LassoCD"""

    learner = lasso_coordinate_descent

    def __init__(self, alpha=None, w0=None):
        self.alpha = alpha
        self.w = w0
        self.E = None

    def fit(self, X, y, maxit=10):
        """fit Lasso model with coordinate descent"""
        nsamples, nfeatures = X.shape

        if self.w is None:
            self.w = np.zeros(nfeatures)

        self.w, self.E = lasso_coordinate_descent(X, y, self.alpha, self.w, maxit=10)

        # Check convergence
        self.gap, _, _ = enet_dual_gap(X, y, self.w, self.alpha, beta=0)
        return self

    def predict(self, X):
        """Predict with Linear Model
        """
        y = np.dot(X,self.w)
        return y

if __name__ == '__main__':
    # N, P, maxit = 5, 10, 30
    N, P, maxit = 100, 10000, 30
    np.random.seed(0)
    y = np.random.randn(N)
    X = np.random.randn(N,P)

    import time
    t0 = time.time()
    model_slow = LassoCD(alpha=1)
    model_slow.learner = lasso_coordinate_descent_slow
    model_slow.fit(X, y, maxit=maxit)
    print time.time() - t0

    t0 = time.time()
    model_fast = LassoCD(alpha=1)
    model_fast.learner = lasso_coordinate_descent_fast
    model_fast.fit(X, y, maxit=maxit)
    print time.time() - t0

    import pylab as pl
    pl.close('all')
    pl.plot(model_fast.E,"rx-")
    pl.plot(model_slow.E,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Slow', 'Fast'])
    pl.show()
