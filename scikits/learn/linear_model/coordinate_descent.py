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

    coef_: ndarray of shape n_features
        The initial coeffients to warm-start the optimization

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    Notes
    -----
    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """

    def __init__(self, alpha=1.0, rho=0.5, fit_intercept=True):
        self.alpha = alpha
        self.rho = rho
        self.coef_ = None
        self.fit_intercept = fit_intercept

    # @profile
    def fit(self, X, y, precompute='auto', Xy=None, max_iter=1000, tol=1e-4,
            coef_init=None, **params):
        """Fit Elastic Net model with coordinate descent

        Parameters
        -----------
        X: ndarray, (n_samples, n_features)
            Data
        y: ndarray, (n_samples)
            Target
        precompute : True | False | 'auto' | array-like
            Whether to use a precomputed Gram matrix to speed up
            calculations. If set to 'auto' let us decide. The Gram
            matrix can also be passed as argument.
        Xy : array-like, optional
            Xy = np.dot(X.T, y) that can be precomputed. It is useful
            only when the Gram matrix is precomuted.
        max_iter: int, optional
            The maximum number of iterations
        tol: float, optional
            The tolerance for the optimization: if the updates are
            smaller than 'tol', the optimization code checks the
            dual gap for optimality and continues until it is smaller
            than tol.

        Notes
        -----

        Coordinate descent is an algorithm that considers each column of
        data at a time hence it will automatically convert the X input
        as a fortran contiguous numpy array if necessary.

        To avoid memory re-allocation it is advised to allocate the
        initial data in memory directly using that format.
        """
        self._set_params(**params)
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

        if coef_init is None:
            self.coef_ = np.zeros(X.shape[1], dtype=np.float64)
        else:
            self.coef_ = coef_init

        n_samples = X.shape[0]
        alpha = self.alpha * self.rho * n_samples
        beta = self.alpha * (1.0 - self.rho) * n_samples

        X = np.asfortranarray(X) # make data contiguous in memory

        # precompute if n_samples > n_features
        if hasattr(precompute, '__array__'):
            Gram = precompute
        elif precompute == True or \
               (precompute == 'auto' and X.shape[0] > X.shape[1]):
            Gram = np.dot(X.T, X)
        else:
            Gram = None

        if Gram is None:
            self.coef_, self.dual_gap_, self.eps_ = \
                    cd_fast.enet_coordinate_descent(self.coef_, alpha, beta,
                                                    X, y, max_iter, tol)
        else:
            if Xy is None:
                Xy = np.dot(X.T, y)
            self.coef_, self.dual_gap_, self.eps_ = \
                    cd_fast.enet_coordinate_descent_gram(self.coef_, alpha,
                                beta, Gram, Xy, y, max_iter, tol)

        self._set_intercept(Xmean, ymean)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want'
                          ' to increase the number of interations')

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
    >>> from scikits.learn import linear_model
    >>> clf = linear_model.Lasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso(alpha=0.1, fit_intercept=True)
    >>> print clf.coef_
    [ 0.85  0.  ]
    >>> print clf.intercept_
    0.15

    See also
    --------
    LassoLARS

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """

    def __init__(self, alpha=1.0, fit_intercept=True):
        super(Lasso, self).__init__(alpha=alpha, rho=1.0,
                                    fit_intercept=fit_intercept)

###############################################################################
# Classes to store linear models along a regularization path

def lasso_path(X, y, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
               verbose=False, **fit_params):
    """Compute Lasso path with coordinate descent

    Parameters
    ----------
    X : numpy array of shape [n_samples,n_features]
        Training data. Pass directly as fortran contiguous data to avoid
        unnecessary memory duplication

    y : numpy array of shape [n_samples]
        Target values

    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    fit_params : kwargs
        keyword arguments passed to the Lasso fit method

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """
    return enet_path(X, y, rho=1., eps=eps, n_alphas=n_alphas, alphas=alphas,
                  fit_intercept=fit_intercept, verbose=verbose, **fit_params)


def enet_path(X, y, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
              fit_intercept=True, verbose=False, **fit_params):
    """Compute Elastic-Net path with coordinate descent

    Parameters
    ----------
    X : numpy array of shape [n_samples, n_features]
        Training data. Pass directly as fortran contiguous data to avoid
        unnecessary memory duplication

    y : numpy array of shape [n_samples]
        Target values

    rho : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). rho=1 corresponds to the Lasso

    eps : float
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    fit_params : kwargs
        keyword arguments passed to the Lasso fit method

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.
    """
    X, y, Xmean, ymean = LinearModel._center_data(X, y, fit_intercept)
    X = np.asfortranarray(X) # make data contiguous in memory

    n_samples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(np.dot(X.T, y)).max() / (n_samples * rho)
        alphas = np.logspace(np.log10(alpha_max*eps), np.log10(alpha_max),
                             num=n_alphas)[::-1]
    else:
        alphas = np.sort(alphas)[::-1] # make sure alphas are properly ordered
    coef_ = None # init coef_
    models = []

    if not 'precompute' in fit_params \
        or fit_params['precompute'] is True \
        or (fit_intercept and hasattr(fit_params['precompute'], '__array__')):
        fit_params['precompute'] = np.dot(X.T, X)
        if not 'Xy' in fit_params or fit_params['Xy'] is None:
            fit_params['Xy'] = np.dot(X.T, y)

    for alpha in alphas:
        model = ElasticNet(alpha=alpha, rho=rho, fit_intercept=False)
        model.fit(X, y, coef_init=coef_, **fit_params)
        if fit_intercept:
            model.fit_intercept = True
            model._set_intercept(Xmean, ymean)
        if verbose:
            print model
        coef_ = model.coef_.copy()
        models.append(model)
    return models


class LinearModelCV(LinearModel):
    """Base class for iterative model fitting along a regularization path"""

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True):
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept

    def fit(self, X, y, cv=None, **fit_params):
        """Fit linear model with coordinate descent along decreasing alphas
        using cross-validation

        Parameters
        ----------

        X : numpy array of shape [n_samples,n_features]
            Training data. Pass directly as fortran contiguous data to avoid
            unnecessary memory duplication

        y : numpy array of shape [n_samples]
            Target values

        cv : cross-validation generator, optional
            If None, KFold will be used.

        fit_params : kwargs
            keyword arguments passed to the Lasso fit method

        """

        X = np.asfortranarray(X, dtype=np.float64)
        y = np.asanyarray(y, dtype=np.float64)

        n_samples = X.shape[0]

        # Start to compute path on full data
        path_params = fit_params.copy()
        path_params.update(self._get_params())
        models = self.path(X, y, **path_params)

        alphas = [model.alpha for model in models]
        n_alphas = len(alphas)

        # init cross-validation generator
        cv = cv if cv else KFold(n_samples, 5)

        params = self._get_params()
        params['alphas'] = alphas
        params['n_alphas'] = n_alphas

        # Compute path for all folds and compute MSE to get the best alpha
        folds = list(cv)
        mse_alphas = np.zeros((len(folds), n_alphas))
        fit_params.update(params)
        for i, (train, test) in enumerate(folds):
            models_train = self.path(X[train], y[train], **fit_params)
            for i_alpha, model in enumerate(models_train):
                y_ = model.predict(X[test])
                mse_alphas[i, i_alpha] += ((y_ - y[test]) ** 2).mean()

        i_best_alpha = np.argmin(np.mean(mse_alphas, axis=0))
        model = models[i_best_alpha]

        self.coef_ = model.coef_
        self.intercept_ = model.intercept_
        self.alpha = model.alpha
        self.alphas = np.asarray(alphas)
        self.coef_path_ = np.asarray([model.coef_ for model in models])
        self.mse_path_ = mse_alphas.T
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
    See examples/linear_model/lasso_path_with_crossvalidation.py
    for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
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
    See examples/linear_model/lasso_path_with_crossvalidation.py
    for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """

    path = staticmethod(enet_path)

    def __init__(self, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True):
        self.rho = rho
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
