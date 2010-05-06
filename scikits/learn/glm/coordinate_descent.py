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

        0.5 * ||R||_2 ^ 2 + alpha * rho * ||w||_1 + alpha * (1-rho) * 0.5 * ||w||_2 ^ 2

Where R are the residuals between the output of the model and the expected
value and w is the vector of weights to fit.
"""

import warnings
import numpy as np

from cd_fast import lasso_coordinate_descent, enet_coordinate_descent
from utils import lasso_objective, enet_objective, density

class LinearModel(object):
    """Base class for Linear Model optimized with coordinate descent"""

    def __init__(self, w0=None):
        # weights of the model (can be lazily initialized by the ``fit`` method)
        self.coef_ = w0

    def predict(self, X):
        """Linear model prediction: compute the dot product with the weights"""
        X = np.asanyarray(X)
        return np.dot(X, self.coef_)

    def compute_density(self):
        """Ratio of non-zero weights in the model"""
        return density(self.coef_)


class Lasso(LinearModel):
    """Linear Model trained with L1 prior as regularizer (a.k.a. the Lasso)"""

    def __init__(self, alpha=1.0, w0=None):
        super(Lasso, self).__init__(w0)
        self.alpha = float(alpha)

    def fit(self, X, Y, maxit=100, tol=1e-4):
        """Fit Lasso model with coordinate descent"""
        X = np.asanyarray(X, dtype=np.float64)
        Y = np.asanyarray(Y, dtype=np.float64)
        nsamples =  X.shape[0]

        if self.coef_ is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)

        self.coef_, self.dual_gap_, self.eps_ = \
                    lasso_coordinate_descent(self.coef_, self.alpha, X, Y, maxit, 10, tol)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want to increase the number of interations')

        # return self for chaining fit and predict calls
        return self


    def __repr__(self):
        return "Lasso cd"


class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer"""

    def __init__(self, alpha=1.0, rho=0.5, w0=None):
        super(ElasticNet, self).__init__(w0)
        self.alpha = alpha
        self.rho = rho

    def fit(self, X, Y, maxit=100, tol=1e-4):
        """Fit Elastic Net model with coordinate descent"""
        X = np.asanyarray(X, dtype=np.float64)
        Y = np.asanyarray(Y, dtype=np.float64)

        if self.coef_ is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)

        nsamples = X.shape[0]
        alpha = self.alpha * self.rho * nsamples
        beta = self.alpha * (1.0 - self.rho) * nsamples
        self.coef_, self.dual_gap_, self.eps_ = \
                    enet_coordinate_descent(self.coef_, alpha, beta, X, Y, maxit, 10, tol)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want to increase the number of interations')

        # return self for chaining fit and predict calls
        return self

    def __repr__(self):
        return "ElasticNet cd"

def lasso_path(X, y, eps=1e-3, n_alphas=100, **kwargs):
    """Compute Lasso path with coordinate descent"""
    nsamples = X.shape[0]
    alpha_max = np.abs(np.dot(X.T, y)).max()
    model = Lasso(alpha=alpha_max)
    weights = []
    alphas = np.linspace(np.log(alpha_max), np.log(eps * alpha_max), n_alphas)
    alphas = np.exp(alphas)
    for alpha in alphas:
        model.alpha = alpha
        model.fit(X, y, **kwargs)
        weights.append(model.coef_.copy())

    weights = np.asarray(weights)
    return alphas, weights

def enet_path(X, y, eps=1e-3, n_alphas=100, rho=0.5, **kwargs):
    """Compute Elastic-Net path with coordinate descent"""
    nsamples = X.shape[0]
    alpha_max = np.abs(np.dot(X.T, y)).max()
    model = ElasticNet(alpha=alpha_max, rho=rho)
    weights = []
    alphas = np.linspace(np.log(alpha_max), np.log(eps * alpha_max), n_alphas)
    alphas = np.exp(alphas)
    for alpha in alphas:
        model.alpha = alpha
        model.fit(X, y, **kwargs)
        weights.append(model.coef_.copy())

    weights = np.asarray(weights)
    return alphas, weights
