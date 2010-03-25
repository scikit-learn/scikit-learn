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

        0.5 * ||R||_2 ^ 2 + alpha * ||w||_1 + beta * 0.5 * ||w||_2 ^ 2

Where R are the residuals between the output of the model and the expected
value and w is the vector of weights to fit.
"""

import numpy as np
import scipy.linalg as linalg
from lasso_cd import lasso_coordinate_descent as lasso_coordinate_descent_slow
from enet_cd import enet_coordinate_descent as enet_coordinate_descent_slow
from iteration_callbacks import IterationCallbackMaxit, IterationCallbackFunc
from utils import enet_dual_gap, lasso_dual_gap, lasso_objective, \
                  enet_objective, density

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

class LinearModel(object):
    """Base class for Linear Model optimized with coordinate descent"""

    def __init__(self, w0=None, callbacks=None):
        # weights of the model (can be lazily initialized by the ``fit`` method)
        self.w = w0

        # callbacks that handles recording of the historic data
        # and can stop iterations
        self.callbacks = []
        if callbacks is not None:
            for callback in callbacks:
                self.callbacks.append(callback)

        self.learner = None
        self.dual_gap_func = None

    def fit(self, X, y, maxit=100, tol=1e-4):
        """Fit Lasso model with coordinate descent"""
        X, y = np.asanyarray(X), np.asanyarray(y)
        n_samples, n_features = X.shape

        if tol is not None:
            cb_dual_gap = IterationCallbackFunc(self._dual_gap_func, tol=tol)
            self.callbacks.append(cb_dual_gap)

        if self.w is None:
            self.w = np.zeros(n_features)

        self.w = self.learner(self, X, y, maxit)

        # return self for chaining fit and predict calls
        return self

    def predict(self, X):
        """Linear model prediction: compute the dot product with the weights"""
        X = np.asanyarray(X)
        y = np.dot(X, self.w)
        return y

    def compute_density(self):
        """Ratio of non-zero weights in the model"""
        return density(self.w)

class Lasso(LinearModel):
    """Linear Model trained with L1 prior as regularizer (a.k.a. the Lasso)"""

    def __init__(self, alpha=1.0, w0=None, callbacks=None):
        super(Lasso, self).__init__(w0, callbacks)
        self.alpha = alpha
        self.learner = lasso_coordinate_descent

    def _dual_gap_func(self, X, y, w, **kw):
        return lasso_dual_gap(X, y, w, kw['alpha'])[0]

    def __repr__(self):
        return "Lasso cd"


class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer"""

    def __init__(self, alpha=1.0, beta=1.0, w0=None, callbacks=None):
        super(ElasticNet, self).__init__(w0, callbacks)
        self.alpha = alpha
        self.beta = beta
        self.learner = enet_coordinate_descent

    def _dual_gap_func(self, X, y, w, **kw):
        return enet_dual_gap(X, y, w, kw['alpha'], kw['beta'])[0]

    def __repr__(self):
        return "ElasticNet cd"


def lasso_path(X, y, factor=0.95, n_alphas = 10, **kwargs):
    """Compute Lasso path with coordinate descent"""
    alpha_max = np.abs(np.dot(X.T, y)).max()
    alpha = alpha_max
    model = Lasso(alpha=alpha)
    weights = []
    alphas = []
    for _ in range(n_alphas):
        model.alpha *= factor
        model.fit(X, y, **kwargs)

        alphas.append(model.alpha)
        weights.append(model.w.copy())

    alphas = np.asarray(alphas)
    weights = np.asarray(weights)
    return alphas, weights

def enet_path(X, y, factor=0.95, n_alphas=10, beta=1.0, **kwargs):
    """Compute Elastic-Net path with coordinate descent"""
    alpha_max = np.abs(np.dot(X.T, y)).max()
    alpha = alpha_max
    model = ElasticNet(alpha=alpha, beta=beta)
    weights = []
    alphas = []
    for _ in range(n_alphas):
        model.alpha *= factor
        model.fit(X, y, **kwargs)

        alphas.append(model.alpha)
        weights.append(model.w.copy())

    alphas = np.asarray(alphas)
    weights = np.asarray(weights)
    return alphas, weights

if __name__ == '__main__':
    import time
    import pylab as pl

    n_samples, n_features, maxit = 5, 10, 30
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    """Tests Lasso implementations (python and cython)
    """

    alpha = 1.0

    tol = 1e-5

    # Callbacks to store objective values and densities
    callback_objective = IterationCallbackFunc(lasso_objective)
    callback_density = IterationCallbackFunc(density)

    t0 = time.time()
    lasso_slow = Lasso(alpha=alpha, callbacks=[callback_objective,
                                               callback_density])
    lasso_slow.learner = lasso_coordinate_descent_slow
    lasso_slow.fit(X, y, maxit=maxit, tol=tol)
    print time.time() - t0

    objective_convergence_slow = callback_objective.values
    density_slow = callback_density.values

    print "Duality gap Lasso (should be small): %f" % \
            lasso_dual_gap(X, y, lasso_slow.w, alpha)[0]

    t0 = time.time()
    lasso_fast = Lasso(alpha=alpha, callbacks=[callback_objective,
                                               callback_density])
    lasso_fast.learner = lasso_coordinate_descent_fast
    lasso_fast.fit(X, y, maxit=maxit, tol=tol)
    print time.time() - t0

    print "Duality gap Lasso (should be small): %f" % \
            lasso_dual_gap(X, y, lasso_slow.w, alpha)[0]

    objective_convergence_fast = callback_objective.values
    density_fast = callback_density.values

    pl.close('all')
    pl.plot(objective_convergence_fast,"rx-")
    pl.plot(objective_convergence_slow,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Fast', 'Slow'])
    pl.title('Lasso')

    pl.figure()
    pl.plot(density_fast,"rx-")
    pl.plot(density_slow,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Density')
    pl.legend(['Fast', 'Slow'])
    pl.title('Lasso')

    """Tests Elastic-Net implementations (python and cython)
    """

    alpha = 1.0
    beta = 1.0

    callback_objective = IterationCallbackFunc(enet_objective)

    import time
    t0 = time.time()
    enet_slow = ElasticNet(alpha=alpha, beta=beta, callbacks=[callback_objective,
                                                              callback_density])
    enet_slow.learner = enet_coordinate_descent_slow
    enet_slow.fit(X, y, maxit=maxit)
    print time.time() - t0

    print "Duality gap (should be small): %f" % \
            enet_dual_gap(X, y, enet_slow.w, alpha)[0]

    objective_convergence_slow = callback_objective.values
    density_slow = callback_density.values

    t0 = time.time()
    enet_fast = ElasticNet(alpha=alpha, beta=beta, callbacks=[callback_objective,
                                                              callback_density])

    enet_fast.learner = enet_coordinate_descent_fast
    enet_fast.fit(X, y, maxit=maxit)
    print time.time() - t0

    print "Duality gap (should be small): %f" % \
            enet_dual_gap(X, y, enet_fast.w, alpha)[0]

    objective_convergence_fast = callback_objective.values
    density_fast = callback_density.values

    pl.figure()
    pl.plot(objective_convergence_fast,"rx-")
    pl.plot(objective_convergence_slow,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.legend(['Fast', 'Slow'])
    pl.title('Elastic-Net')

    pl.figure()
    pl.plot(density_fast,"rx-")
    pl.plot(density_slow,"bo--")
    pl.xlabel('Iteration')
    pl.ylabel('Density')
    pl.legend(['Fast', 'Slow'])
    pl.title('Elastic-Net')

    """Test path functions
    """

    alphas_lasso, weights_lasso = lasso_path(X, y, factor=0.97, n_alphas = 100, tol=1-2)
    alphas_enet, weights_enet = enet_path(X, y, factor=0.97, n_alphas = 100, beta=0.1, tol=1-2)

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

