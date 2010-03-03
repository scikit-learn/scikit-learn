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
    """Compute dual gap for Elastic-Net model to check KKT optimality conditions

    Returns
    -------
    gap : the difference  primal_objective - dual_objective (should be positive)
        A value less that 1e-5 means convergence in practice
    primal_objective : the value of the objective function of the primal problem
    dual_objective : the value of the objective function of the dual problem

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
    """Generic class for Linear Model optimized with coordinate descent"""

    def __init__(self, w0=None):
        self.w = w0
        self.E = None

    def predict(self, X):
        """Predict with Linear Model
        """
        y = np.dot(X,self.w)
        return y

    def density(self):
        """Ratio of non-zero weights in the model"""
        return 0 if self.w is None else float((self.w != 0).sum()) / self.w.size


class Lasso(LinearModel):
    """Linear Model trained with L1 prior as regularizer (a.k.a. the Lasso)"""

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
    """Linear Model trained with L1 and L2 prior as regularizer"""

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

def lasso_path(X, y, factor=0.95, n_alphas = 10):
    """Compute Lasso path with coordinate descent"""
    alpha_max = np.dot(X.T, y).max()
    alpha = alpha_max
    model = Lasso(alpha=alpha)
    weights = []
    alphas = []
    for _ in range(n_alphas):
        model.alpha *= factor
        model.fit(X, y)

        alphas.append(model.alpha)
        weights.append(model.w.copy())

    alphas = np.asarray(alphas)
    weights = np.asarray(weights)
    return alphas, weights

def enet_path(X, y, factor=0.95, n_alphas = 10, beta=1.0):
    """Compute Elastic-Net path with coordinate descent"""
    alpha_max = np.dot(X.T, y).max()
    alpha = alpha_max
    model = ElasticNet(alpha=alpha, beta=beta)
    weights = []
    alphas = []
    for _ in range(n_alphas):
        model.alpha *= factor
        model.fit(X, y)

        alphas.append(model.alpha)
        weights.append(model.w.copy())

    alphas = np.asarray(alphas)
    weights = np.asarray(weights)
    return alphas, weights

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

    print "Duality gap Lasso (should be small): %f" % lasso_fast.gap

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

    print "Duality gap (should be small): %f" % enet_fast.gap

    pl.figure()
    pl.plot(enet_fast.E,"rx-")
    pl.plot(enet_slow.E,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Slow', 'Fast'])
    pl.title('Elastic-Net')
    pl.show()

    """Test path functions
    """

    alphas_lasso, weights_lasso = lasso_path(X, y, factor=0.97, n_alphas = 100)
    alphas_enet, weights_enet = enet_path(X, y, factor=0.97, n_alphas = 100, beta=0.1)

    from itertools import cycle
    color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

    import pylab as pl
    pl.close('all')
    for color, weight_lasso, weight_enet in zip(color_iter, weights_lasso.T, weights_enet.T):
        pl.plot(-np.log(alphas_lasso), weight_lasso, color)
        pl.plot(-np.log(alphas_enet), weight_enet, color+'x')
    pl.xlabel('-log(lambda)')
    pl.ylabel('weights')
    pl.title('Lasso and Elastic-Net Paths')
    pl.legend(['Lasso','Elastic-Net'])
    pl.show()