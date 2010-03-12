# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$
"""Implementation of regularized linear regression with Coordinate Descent

This implementation is focused on regularizers that lead to sparse parameters
(many zeros) such as the laplacian (L1) and Elastic Net (L1 + L2) priors:

  http://en.wikipedia.org/wiki/Generalized_linear_model

The objective function to minimize is for the Lasso::

        0.5 * ||R||_2 ^ 2 + alpha * ||w||_1

and for the Elastic Network::

        0.5 * ||R||_2 ^ 2 + alpha * ||w||_1 + beta * ||w||_2 ^ 2

Where R are the residuals between the output of the model and the expected
value and w is the vector of weights to fit.
"""

import numpy as np
import scipy.linalg as linalg
from lasso_cd import lasso_coordinate_descent as lasso_coordinate_descent_slow
from enet_cd import enet_coordinate_descent as enet_coordinate_descent_slow

# Attempt to improve speed with cython
try:
    from lasso_cd_fast import lasso_coordinate_descent \
            as lasso_coordinate_descent_fast
    from enet_cd_fast import enet_coordinate_descent \
            as enet_coordinate_descent_fast
    lasso_coordinate_descent = lasso_coordinate_descent_fast
    enet_coordinate_descent = enet_coordinate_descent_fast
except ImportError:
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
    Xw = np.dot(X, w)
    A = (y - Xw)
    if beta > 0:
        B = - np.sqrt(beta) * w
    XtA = np.dot(X.T, A)
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


class BaseIterationCallback(object):
    """Base callback to be called at the end of each iteration of CD

    - record the value of the current objective cost
    - record the density of the model
    - record and check the duality gap for early stop of the optim
      (before maxiter)

    To be subclassed if more monitoring is required.
    """

    def __init__(self, linear_model, gap_tolerance=1e-4):
        self.linear_model = linear_model
        self.gap_tolerance = gap_tolerance

    def __call__(self, X, y, R, alpha, w, iter):
        # TODO: check the last time stamp to avoid computing the stats too often
        lm = self.linear_model
        lm.compute_objective(X, y, R, record=True)
        lm.compute_density(record=True)
        gap = lm.compute_gap(X, y, record=True)

        # should go on?
        if len(lm.gap) > 1 and lm.gap[-1] > lm.gap[-2]:
            # if the gap increases it means that it means we reached convergence
            # this is a consequence of the way we compute dual_objective
            return False
        return gap > self.gap_tolerance


class LinearModel(object):
    """Base class for Linear Model optimized with coordinate descent"""

    def __init__(self, w0=None):
        # weights of the model (can be lazily initialized by the ``fit`` method)
        self.w = w0

        # recorded historic data at each training iteration, suitable for
        # plotting and monitoring of the convergence
        self.objective = []
        self.gap = []
        self.density = []

        # callback that handles recording of the historic data
        self.callback = BaseIterationCallback(self)

    def predict(self, X):
        """Linear model prediction: compute the dot product with the weights"""
        X = np.asanyarray(X)
        y = np.dot(X, self.w)
        return y

    def compute_density(self, record=False):
        """Ratio of non-zero weights in the model"""
        d = 0 if self.w is None else float((self.w != 0).sum()) / self.w.size
        if record:
            self.density.append(d)
        return d


class Lasso(LinearModel):
    """Linear Model trained with L1 prior as regularizer (a.k.a. the Lasso)"""

    def __init__(self, alpha=1.0, w0=None):
        super(Lasso, self).__init__(w0)
        self.alpha = alpha
        self.learner = lasso_coordinate_descent

    def __repr__(self):
	return "Lasso cd"

    def fit(self, X, y, maxit=10):
        """Fit Lasso model with coordinate descent"""
        X, y = np.asanyarray(X), np.asanyarray(y)
        nsamples, nfeatures = X.shape

        if self.w is None:
            self.w = np.zeros(nfeatures)

        self.w = self.learner(X, y, self.alpha, self.w, maxit=maxit,
                              callback=self.callback)

        # return self for chaining fit and predict calls
        return self

    def compute_gap(self, X, y, record=False):
        """Evaluate the duality gap of the current state of the model"""
        gap, _, _ = enet_dual_gap(X, y, self.w, self.alpha, beta=0)
        if record:
            self.gap.append(gap)
        return gap

    def compute_objective(self, X, y, R=None, record=False):
        """Evaluate the cost function to minimize"""
        if R is None:
            R = y - np.dot(X, self.w)
        cost = 0.5 * linalg.norm(R) ** 2 + self.alpha * np.abs(self.w).sum()
        if record:
            self.objective.append(cost)
        return cost



class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer"""

    def __init__(self, alpha=1.0, beta=1.0, w0=None):
        super(ElasticNet, self).__init__(w0)
        self.alpha = alpha
        self.beta = beta
        self.learner = enet_coordinate_descent

    def __repr__(self):
	return "ElasticNet cd"

    def fit(self, X, y, maxit=10):
        """Fit Elastic Net model with coordinate descent"""
        X, y = np.asanyarray(X), np.asanyarray(y)
        nsamples, nfeatures = X.shape

        if self.w is None:
            self.w = np.zeros(nfeatures)

        self.w = self.learner(X, y, self.alpha, self.beta, self.w, maxit=maxit,
                              callback=self.callback)

        # return self for chaining fit and predict calls
        return self

    def compute_gap(self, X, y, record=False):
        gap, _, _ = enet_dual_gap(X, y, self.w, self.alpha, self.beta)
        if record:
            self.gap.append(gap)
        return gap

    def compute_objective(self, X, y, R=None, record=False):
        """Evaluate the cost function to minimize"""
        if R is None:
            nsamples, nfeatures = X.shape
            R = np.empty(nfeatures + nsamples)
            R[:nsamples] = y - np.dot(X, self.w)
            R[nsamples:] = - sqrt(self.beta) * self.w
        cost = 0.5 * linalg.norm(R) ** 2 + self.alpha * np.abs(self.w).sum() + \
                0.5 * self.beta * linalg.norm(self.w) ** 2
        if record:
            self.objective.append(cost)
        return cost


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
    X = np.random.randn(N, P)

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

    print "Duality gap Lasso (should be small): %f" % lasso_fast.gap[-1]

    import pylab as pl
    pl.close('all')
    pl.plot(lasso_fast.objective,"rx-")
    pl.plot(lasso_slow.objective,"bo--")
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

    print "Duality gap (should be small): %f" % enet_fast.gap[-1]

    pl.figure()
    pl.plot(enet_fast.objective,"rx-")
    pl.plot(enet_slow.objective,"bo--")
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

    pl.figure()
    for color, weight_lasso, weight_enet in zip(color_iter, weights_lasso.T, weights_enet.T):
        pl.plot(-np.log(alphas_lasso), weight_lasso, color)
        pl.plot(-np.log(alphas_enet), weight_enet, color+'x')
    pl.xlabel('-log(lambda)')
    pl.ylabel('weights')
    pl.title('Lasso and Elastic-Net Paths')
    pl.legend(['Lasso','Elastic-Net'])
    pl.show()

