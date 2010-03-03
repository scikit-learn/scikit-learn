# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$
"""Implementation of regularized linear regression with Coordinate Descent

We focus the implementation on regularizer that lead to sparse parameters (many
zeros) such as the laplacian (L1) and Elastic Net (L1 + L2) priors.
"""

import numpy as np
import scipy.linalg as linalg
from lasso_cd import lasso_coordinate_descent as lasso_coordinate_descent_slow
from enet_cd import enet_coordinate_descent as enet_coordinate_descent_slow

# Attempt to improve speed with cython
try:
    from lasso_cd_fast import lasso_coordinate_descent as lasso_coordinate_descent_fast
    from enet_cd_fast import enet_coordinate_descent as enet_coordinate_descent_fast
    lasso_coordinate_descent = lasso_coordinate_descent_fast
    enet_coordinate_descent = enet_coordinate_descent_fast
except ImportError, e:
    lasso_coordinate_descent = lasso_coordinate_descent_slow
    enet_coordinate_descent = enet_coordinate_descent_slow
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
        B = - np.sqrt(beta) * w
    XtA = np.dot(X.T,A)
    if beta > 0:
        XtA += np.sqrt(beta) * B
    dual_norm_XtA = np.max(XtA)
    if (dual_norm_XtA > alpha):
        A *= alpha / dual_norm_XtA
        if beta > 0:
            B *= alpha / dual_norm_XtA
    pobj = 0.5 * linalg.norm(y - Xw)**2 + alpha * np.abs(w).sum() \
           + 0.5 * beta * linalg.norm(w)**2
    dobj = - 0.5 * linalg.norm(A)**2 + np.dot(A.T, y)
    if beta > 0:
        dobj += - 0.5 * linalg.norm(B)**2
    gap = pobj - dobj
    return gap, pobj, dobj

class LinearModel(object):
    """Generic class for Linear Model optimized
    coordinate descent
    """

    def __init__(self, w0=None):
        self.w = w0
        self.E = None

    def predict(self, X):
        """Predict with Linear Model
        """
        y = np.dot(X,self.w)
        return y

class Lasso(LinearModel):
    """Class Lasso"""

    def __init__(self, alpha=None, w0=None):
        super(Lasso, self).__init__(w0)
        self.alpha = alpha
        self.learner = lasso_coordinate_descent

    def fit(self, X, y, maxit=10):
        """Fit Lasso model with coordinate descent"""
        nsamples, nfeatures = X.shape

        if self.w is None:
            self.w = np.zeros(nfeatures)

        self.w, self.E = self.learner(X, y, self.alpha, self.w, maxit=maxit)

        # Check convergence
        self.gap, _, _ = enet_dual_gap(X, y, self.w, self.alpha, beta=0)
        return self

class ElasticNet(LinearModel):
    """Class ElasticNet"""

    def __init__(self, alpha=None, beta=None, w0=None):
        super(ElasticNet, self).__init__(w0)
        self.alpha = alpha
        self.beta = beta
        self.learner = enet_coordinate_descent

    def fit(self, X, y, maxit=10):
        """Fit Elastic Net model with coordinate descent"""
        nsamples, nfeatures = X.shape

        if self.w is None:
            self.w = np.zeros(nfeatures)

        self.w, self.E = self.learner(X, y, self.alpha, self.beta, \
                                                     self.w, maxit=maxit)

        # Check convergence
        self.gap, _, _ = enet_dual_gap(X, y, self.w, self.alpha, beta=self.beta)
        return self

if __name__ == '__main__':
    N, P, maxit = 5, 10, 30
    np.random.seed(0)
    y = np.random.randn(N)
    X = np.random.randn(N,P)

    """Tests Lasso implementations (python and cython)
    """

    alpha = 1.0

    import time
    t0 = time.time()
    lasso_slow = Lasso(alpha=alpha)
    lasso_slow.learner = lasso_coordinate_descent_slow
    lasso_slow.fit(X, y, maxit=maxit)
    print time.time() - t0

    t0 = time.time()
    lasso_fast = Lasso(alpha=alpha)
    lasso_fast.learner = lasso_coordinate_descent_fast
    lasso_fast.fit(X, y, maxit=maxit)
    print time.time() - t0

    print "Duality gap Lasso (should be small): %f"%lasso_fast.gap

    import pylab as pl
    pl.close('all')
    pl.plot(lasso_fast.E,"rx-")
    pl.plot(lasso_slow.E,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Slow', 'Fast'])
    pl.title('Lasso')
    pl.show()

    """Tests Elastic-Net implementations (python and cython)
    """

    alpha = 1.0
    beta = 1.0

    import time
    t0 = time.time()
    enet_slow = ElasticNet(alpha=alpha, beta=beta)
    enet_slow.learner = enet_coordinate_descent_slow
    enet_slow.fit(X, y, maxit=maxit)
    print time.time() - t0

    t0 = time.time()
    enet_fast = ElasticNet(alpha=alpha, beta=beta)
    enet_fast.learner = enet_coordinate_descent_fast
    enet_fast.fit(X, y, maxit=maxit)
    print time.time() - t0

    print "Duality gap (should be small): %f"%enet_fast.gap

    pl.figure()
    pl.plot(enet_fast.E,"rx-")
    pl.plot(enet_slow.E,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Slow', 'Fast'])
    pl.title('Elastic-Net')
    pl.show()
