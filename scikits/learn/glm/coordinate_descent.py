# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

import warnings
import numpy as np

from .base import LinearModel
from ..cross_val import KFold
from . import cd_fast


###############################################################################
# ElasticNet model

class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer

    rho=1 is the lasso penalty. Currently, rho <= 0.01 is not
    reliable, unless you supply your own sequence of alpha.

    Parameters
    ----------
    alpha : float
        Constant that multiplies the L1 term. Defaults to 1.0
    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1.
    coef: ndarray of shape n_features
        The initial coeffients to warm-start the optimization
    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.
    """

    def __init__(self, alpha=1.0, rho=0.5, coef_=None, 
                fit_intercept=True):
        self.alpha = alpha
        self.rho = rho
        self.coef_ = coef_
        self.fit_intercept = fit_intercept

    def fit(self, X, Y, maxit=1000, tol=1e-4, **params):
        """Fit Elastic Net model with coordinate descent"""
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64)
        Y = np.asanyarray(Y, dtype=np.float64)

        X, Y, Xmean, Ymean = self._center_data (X, Y)

        if self.coef_ is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)

        n_samples = X.shape[0]
        alpha = self.alpha * self.rho * n_samples
        beta = self.alpha * (1.0 - self.rho) * n_samples

        X = np.asfortranarray(X) # make data contiguous in memory

        self.coef_, self.dual_gap_, self.eps_ = \
                cd_fast.enet_coordinate_descent(self.coef_, alpha, beta, X, Y,
                                        maxit, tol)

        self._set_intercept(Xmean, Ymean)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want'
                                'to increase the number of interations')

        # Store explained variance for __str__
        self.explained_variance_ = self._explained_variance(X, Y)

        # return self for chaining fit and predict calls
        return self


###############################################################################
# Lasso model

class Lasso(ElasticNet):
    """Linear Model trained with L1 prior as regularizer (aka the Lasso)

    Technically the Lasso model is optimizing the same objective function as
    the Elastic Net with rho=1.0 (no L2 penalty).

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1 term. Defaults to 1.0

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    Examples
    --------
    >>> from scikits.learn import glm
    >>> clf = glm.Lasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso(alpha=0.1, coef_=array([ 0.85,  0.  ]), fit_intercept=True)
    >>> print clf.coef_
    [ 0.85  0.  ]
    >>> print clf.intercept_
    0.15

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.
    """

    def __init__(self, alpha=1.0, fit_intercept=True, coef_=None):
        super(Lasso, self).__init__(alpha=alpha, rho=1.0,
                                    fit_intercept=fit_intercept, coef_=coef_)

###############################################################################
# Classes to store linear models along a regularization path 

def lasso_path(X, y, eps=1e-3, n_alphas=100, alphas=None,
               verbose=False, fit_params=dict()):
    """
    Compute Lasso path with coordinate descent

    Parameters
    ----------
    X : numpy array of shape [n_samples,n_features]
        Training data

    Y : numpy array of shape [n_samples]
        Target values

    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    fit_params : dict, optional
        keyword arguments passed to the Lasso fit method

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.
    """
    n_samples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(np.dot(X.T, y)).max() / n_samples
        alphas = np.logspace(np.log10(alpha_max*eps), np.log10(alpha_max),
                             num=n_alphas)[::-1]
    else:
        # XXX: Maybe should reorder the models when outputing them, so
        # that they are ordered in the order of the initial alphas
        alphas = np.sort(alphas)[::-1] # make sure alphas are properly ordered
    coef_ = None # init coef_
    models = []
    for alpha in alphas:
        model = Lasso(coef_=coef_, alpha=alpha)
        model.fit(X, y, **fit_params)
        if verbose: 
            print model
        coef_ = model.coef_.copy()
        models.append(model)
    return models


def enet_path(X, y, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
              verbose=False, fit_params=dict()):
    """Compute Elastic-Net path with coordinate descent

    Parameters
    ----------
    X : numpy array of shape [n_samples,n_features]
        Training data

    Y : numpy array of shape [n_samples]
        Target values

    eps : float
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    fit_params : dict, optional
        keyword arguments passed to the ElasticNet fit method

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.
    """
    n_samples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(np.dot(X.T, y)).max() / (n_samples*rho)
        alphas = np.logspace(np.log10(alpha_max*eps), np.log10(alpha_max),
                             num=n_alphas)[::-1]
    else:
        alphas = np.sort(alphas)[::-1] # make sure alphas are properly ordered
    coef_ = None # init coef_
    models = []
    for alpha in alphas:
        model = ElasticNet(coef_=coef_, alpha=alpha, rho=rho)
        model.fit(X, y, **fit_params)
        if verbose: print model
        coef_ = model.coef_.copy()
        models.append(model)
    return models


class LinearModelCV(LinearModel):
    """Base class for iterative model fitting along a regularization path"""

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None):
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas

    def fit(self, X, y, cv=None, **fit_params):
        """Fit linear model with coordinate descent along decreasing alphas
        using cross-validation

        Parameters
        ----------

        X : numpy array of shape [n_samples,n_features]
            Training data

        Y : numpy array of shape [n_samples]
            Target values

        cv : cross-validation generator, optional
             If None, KFold will be used.

        fit_params : kwargs
            keyword arguments passed to the Lasso fit method

        """

        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        n_samples = X.shape[0]

        # Start to compute path on full data
        models = self.path(X, y, fit_params=fit_params, **self._get_params())

        alphas = [model.alpha for model in models]
        n_alphas = len(alphas)

        # init cross-validation generator
        cv = cv if cv else KFold(n_samples, 5)

        params = self._get_params()
        params['alphas'] = alphas
        params['n_alphas'] = n_alphas

        # Compute path for all folds and compute MSE to get the best alpha
        mse_alphas = np.zeros(n_alphas)
        for train, test in cv:
            models_train = self.path(X[train], y[train], fit_params=fit_params,
                                        **params)
            for i_alpha, model in enumerate(models_train):
                y_ = model.predict(X[test])
                mse_alphas[i_alpha] += ((y_ - y[test]) ** 2).mean()

        i_best_alpha = np.argmin(mse_alphas)
        model = models[i_best_alpha]

        self.coef_ = model.coef_
        self.intercept_ = model.intercept_
        self.explained_variance_ = model.explained_variance_
        self.alpha = model.alpha
        self.alphas = np.asarray(alphas)
        return self


class LassoCV(LinearModelCV):
    """Lasso linear model with iterative fitting along a regularization path

    The best model is selected by cross-validation.

    Parameters
    ----------
    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    Notes
    -----
    See examples/glm/lasso_path_with_crossvalidation.py for an example.
    """

    path = staticmethod(lasso_path)


class ElasticNetCV(LinearModelCV):
    """Elastic Net model with iterative fitting along a regularization path

    The best model is selected by cross-validation.

    Parameters
    ----------
    rho : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties)

    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    Notes
    -----
    See examples/glm/lasso_path_with_crossvalidation.py for an example.
    """

    path = staticmethod(enet_path)

    def __init__(self, rho=0.5, eps=1e-3, n_alphas=100, alphas=None):
        self.rho = rho
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas

