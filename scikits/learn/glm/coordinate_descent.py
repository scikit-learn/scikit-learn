# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

"""Implementation of regularized linear regression with Coordinate Descent

This implementation is focused on regularizers that lead to sparse parameters
(many zeros) such as the laplacian (L1) and Elastic Net (L1 + L2) priors:

  http://en.wikipedia.org/wiki/Generalized_linear_model

The objective function to minimize is for the Lasso::

        0.5 * ||y - X w||_2 ^ 2 + alpha * ||w||_1

and for the Elastic Network::

        0.5 * ||y - X w||_2 ^ 2 + alpha * rho * ||w||_1
        + alpha * (1-rho) * 0.5 * ||w||_2 ^ 2

Where R are the residuals between the output of the model and the expected
value and w is the vector of weights to fit.
"""

import warnings
import numpy as np

from .cd_fast import lasso_coordinate_descent, enet_coordinate_descent
from .utils import density
from ..cross_val import KFold

class LinearModel(object):
    """Base class for Linear Model optimized with coordinate descent"""

    def __init__(self, coef=None):
        # weights of the model (can be lazily initialized by the ``fit`` method)
        self.coef_ = coef

    def predict(self, X):
        """
        Predict using the linear model

        Parameters
        ----------
        X : numpy array of shape [nsamples,nfeatures]

        Returns
        -------
        C : array, shape = [nsample]
            Returns predicted values.
        """
        X = np.asanyarray(X)
        return np.dot(X, self.coef_) + self.intercept_

    def compute_density(self):
        """Ratio of non-zero weights in the model"""
        return density(self.coef_)

    def compute_rsquared(self, X, Y):
        """Compute explained variance a.k.a. r^2"""
        self.rsquared_ = 1 - np.linalg.norm(Y - np.dot(X, self.coef_))**2 \
                         / np.linalg.norm(Y)**2

    def __str__(self):
        if self.coef_ is not None:
            n_non_zeros = (np.abs(self.coef_) != 0).sum()
            return ("%s with %d non-zero coefficients (%.2f%%)\n" + \
                    " * Regularisation parameter = %.7f\n" +\
                    " * Training r^2: %.4f") % \
                    (self.__class__.__name__, n_non_zeros,
                     n_non_zeros / float(len(self.coef_)) * 100,
                     self.alpha, self.rsquared_)
        else:
            return ("%s\n" + \
                    " * Regularisation parameter = %.7f\n" +\
                    " * No fit") % \
                    (self.__class__.__name__, self.alpha)

class Lasso(LinearModel):
    """
    Linear Model trained with L1 prior as regularizer (a.k.a. the
    lasso).

    The lasso estimate solves the minization of the least-squares
    penalty with alpha * ||beta||_1 added, where alpha is a constant and
    ||beta||_1 is the L1-norm of the parameter vector.

    This formulation is useful in some context due to its tendency to
    prefer solutions with fewer parameter values, effectively reducing
    the number of variables upon which the given solution is
    dependent. For this reason, the LASSO and its variants are
    fundamental to the field of compressed sensing.

    Parameters
    ----------
    alpha : double
        Constant that multiplies the L1 term.

    Attributes
    ----------
    `coef_` : array, shape = [nfeatures]
        parameter vector (w in the fomulation formula)

    Examples
    --------
    >>> from scikits.learn import glm
    >>> clf = glm.Lasso()
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso Coordinate Descent
    >>> print clf.coef_
    [ 0.4  0. ]

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.x
    """

    def __init__(self, alpha=1.0, coef=None, tol=1e-4):
        super(Lasso, self).__init__(coef)
        self.alpha = float(alpha)
        self.tol = tol

    def fit(self, X, Y, intercept=True, maxit=1000):
        """
        Fit Lasso model.

        Parameters
        ----------
        X : numpy array of shape [nsamples,nfeatures]
            Training data
        Y : numpy array of shape [nsamples]
            Target values
        intercept : boolean
            whether to calculate the intercept for this model. If set
            to false, no intercept will be used in calculations
            (e.g. data is expected to be already centered).

        Returns
        -------
        self : returns an instance of self.
        """
        X = np.asanyarray(X, dtype=np.float64)
        Y = np.asanyarray(Y, dtype=np.float64)

        self._intercept = intercept
        if self._intercept:
            self._xmean = X.mean(axis=0)
            self._ymean = Y.mean(axis=0)
            X = X - self._xmean
            Y = Y - self._ymean
        else:
            self._xmean = 0.
            self._ymean = 0.

        nsamples = X.shape[0]
        alpha = self.alpha * nsamples

        if self.coef_ is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)

        self.coef_, self.dual_gap_, self.eps_ = \
                    lasso_coordinate_descent(self.coef_, alpha, X, Y, maxit, \
                    10, self.tol)

        self.intercept_ = self._ymean - np.dot(self._xmean, self.coef_)

        # TODO: why not define a method rsquared that computes this ?
        self.compute_rsquared(X, Y)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want to increase the number of interations')

        # return self for chaining fit and predict calls
        return self

    def __repr__(self):
        return "Lasso Coordinate Descent"


class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer

    rho=1 is the lasso penalty. Currently, rho <= 0.01 is not
    reliable, unless you supply your own sequence of alpha.

    Parameters
    ----------
    alpha : double
        TODO
    rho : double
        The ElasticNet mixing parameter, with 0 < rho <= 1.
    """

    def __init__(self, alpha=1.0, rho=0.5, coef=None, tol=1e-4):
        super(ElasticNet, self).__init__(coef)
        self.alpha = alpha
        self.rho = rho
        self.tol = tol

    def fit(self, X, Y, intercept=True, maxit=1000):
        """Fit Elastic Net model with coordinate descent"""
        X = np.asanyarray(X, dtype=np.float64)
        Y = np.asanyarray(Y, dtype=np.float64)

        self._intercept = intercept
        if self._intercept:
            self._xmean = X.mean(axis=0)
            self._ymean = Y.mean(axis=0)
            X = X - self._xmean
            Y = Y - self._ymean
        else:
            self._xmean = 0
            self._ymean = 0

        if self.coef_ is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)

        nsamples = X.shape[0]
        alpha = self.alpha * self.rho * nsamples
        beta = self.alpha * (1.0 - self.rho) * nsamples
        self.coef_, self.dual_gap_, self.eps_ = \
                enet_coordinate_descent(self.coef_, alpha, beta, X, Y,
                                        maxit, 10, self.tol)

        self.intercept_ = self._ymean - np.dot(self._xmean, self.coef_)

        self.compute_rsquared(X, Y)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want to increase the number of interations')

        # return self for chaining fit and predict calls
        return self

    def __repr__(self):
        return "ElasticNet cd"


#########################################################################
#                                                                       #
# The following classes store linear models along a regularization path #
#                                                                       #
#########################################################################

def lasso_path(X, y, eps=1e-3, n_alphas=100, alphas=None, **fit_kwargs):
    """
    Compute Lasso path with coordinate descent

    See examples/plot_lasso_coordinate_descent_path.py for an example.
    """
    nsamples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(np.dot(X.T, y)).max() / nsamples
        alphas = np.linspace(np.log(alpha_max), np.log(eps * alpha_max), n_alphas)
        alphas = np.exp(alphas)
    else:
        alphas = np.sort(alphas)[::-1] # make sure alphas are properly ordered
    coef = None # init coef_
    models = []
    for alpha in alphas:
        model = Lasso(coef=coef, alpha=alpha)
        model.fit(X, y, **fit_kwargs)
        coef = model.coef_.copy()
        models.append(model)
    return models

def enet_path(X, y, rho=0.5, eps=1e-3, n_alphas=100, alphas=None, **fit_kwargs):
    """Compute Elastic-Net path with coordinate descent"""
    nsamples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(np.dot(X.T, y)).max() / (nsamples*rho)
        alphas = np.linspace(np.log(alpha_max), np.log(eps * alpha_max), n_alphas)
        alphas = np.exp(alphas)
    else:
        alphas = np.sort(alphas)[::-1] # make sure alphas are properly ordered
    coef = None # init coef_
    models = []
    for alpha in alphas:
        model = ElasticNet(coef=coef, alpha=alpha, rho=rho)
        model.fit(X, y, **fit_kwargs)
        coef = model.coef_.copy()
        models.append(model)
    return models

def optimized_lasso(X, y, cv=None, n_alphas=100, alphas=None,
                                eps=1e-3, intercept=True, **fit_kwargs):
    """Returns an optimized Lasso instance
    """
    # Start to compute path on full data
    models = lasso_path(X, y, eps=eps, n_alphas=n_alphas, alphas=alphas,
                                intercept=intercept, **fit_kwargs)

    n_samples = y.size
    # init cross-validation generator
    cv = cv if cv else KFold(n_samples, 5)

    alphas = [model.alpha for model in models]
    n_alphas = len(alphas)
    # Compute path for all folds and compute MSE to get the best alpha
    mse_alphas = np.zeros(n_alphas)
    for train, test in cv:
        models_train = lasso_path(X[train], y[train], eps, n_alphas,
                                    alphas=alphas,
                                    intercept=intercept, **fit_kwargs)
        for i_alpha, model in enumerate(models_train):
            y_ = model.predict(X[test])
            mse_alphas[i_alpha] += ((y_ - y[test]) ** 2).mean()

    i_best_alpha = np.argmin(mse_alphas)
    return models[i_best_alpha]

def optimized_enet(X, y, rho=0.5, cv=None, n_alphas=100, alphas=None,
                                 eps=1e-3, intercept=True, **fit_kwargs):
    """Returns an optimized ElasticNet instance
    """
    # Start to compute path on full data
    models = enet_path(X, y, rho=rho, eps=eps, n_alphas=n_alphas,
                                alphas=alphas, intercept=intercept, **fit_kwargs)

    n_samples = y.size
    # init cross-validation generator
    cv = cv if cv else KFold(n_samples, 5)

    alphas = [model.alpha for model in models]
    n_alphas = len(alphas)
    # Compute path for all folds and compute MSE to get the best alpha
    mse_alphas = np.zeros(n_alphas)
    for train, test in cv:
        models_train = enet_path(X[train], y[train], rho=rho,
                                    alphas=alphas, eps=eps, n_alphas=n_alphas,
                                    intercept=intercept, **fit_kwargs)
        for i_alpha, model in enumerate(models_train):
            y_ = model.predict(X[test])
            mse_alphas[i_alpha] += ((y_ - y[test]) ** 2).mean()

    i_best_alpha = np.argmin(mse_alphas)
    return models[i_best_alpha]

class LinearModelPath(LinearModel):
    """Base class for iterative model fitting along a regularization path"""

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, intercept=True):
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        # self.path = None
        self.intercept = intercept

    def fit(self, X, y, cv=None, **kwargs):
        """Fit linear model with coordinate descent along decreasing alphas

        The same model is reused with warm restarts. Early stopping can happen
        before reaching n_alphas if the cross validation detects overfitting
        when decreasing the strength of the regularization.
        """
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        self.path_ = []
        n_samples = X.shape[0]

        model = self.path(X, y, cv=cv, eps=self.eps, n_alphas=self.n_alphas,
                                    intercept=self.intercept, **kwargs)

        self.__dict__.update(model.__dict__)
        return self

class LassoPath(LinearModelPath):
    """Lasso linear model with iterative fitting along a regularization path"""

    @property
    def path(self):
        return optimized_lasso

class ElasticNetPath(LinearModelPath):
    """Elastic Net model with iterative fitting along a regularization path"""

    @property
    def path(self):
        return optimized_enet

    def __init__(self, rho=0.5, **kwargs):
        super(ElasticNetPath, self).__init__(**kwargs)
        self.rho = rho