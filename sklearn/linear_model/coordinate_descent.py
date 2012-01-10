# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

import sys
import warnings
import numpy as np

from .base import LinearModel
from ..utils import as_float_array
from ..cross_validation import check_cv
from . import cd_fast


###############################################################################
# ElasticNet model

class ElasticNet(LinearModel):
    """Linear Model trained with L1 and L2 prior as regularizer

    Minimizes the objective function::

            1 / (2 * n_samples) * ||y - Xw||^2_2 +
            + alpha * rho * ||w||_1 + 0.5 * alpha * (1 - rho) * ||w||^2_2

    If you are interested in controlling the L1 and L2 penalty
    separately, keep in mind that this is equivalent to::

            a * L1 + b * L2

    where::

            alpha = a + b and rho = a / (a + b)

    The parameter rho corresponds to alpha in the glmnet R package while
    alpha corresponds to the lambda parameter in glmnet. Specifically, rho =
    1 is the lasso penalty. Currently, rho <= 0.01 is not reliable, unless
    you supply your own sequence of alpha.

    Parameters
    ----------
    alpha : float
        Constant that multiplies the penalty terms. Defaults to 1.0
        See the notes for the exact mathematical meaning of this
        parameter

    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1. For rho = 0
        the penalty is an L1 penalty. For rho = 1 it is an L2 penalty.
        For 0 < rho < 1, the penalty is a combination of L1 and L2

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: int, optional
        The maximum number of iterations

    copy_X : boolean, optional, default False
        If True, X will be copied; else, it may be overwritten.

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    Notes
    -----
    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """
    def __init__(self, alpha=1.0, rho=0.5, fit_intercept=True,
                 normalize=False, precompute='auto', max_iter=1000,
                 copy_X=True, tol=1e-4):
        self.alpha = alpha
        self.rho = rho
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol

    def fit(self, X, y, Xy=None, coef_init=None):
        """Fit Elastic Net model with coordinate descent

        Parameters
        -----------
        X: ndarray, (n_samples, n_features)
            Data
        y: ndarray, (n_samples)
            Target
        Xy : array-like, optional
            Xy = np.dot(X.T, y) that can be precomputed. It is useful
            only when the Gram matrix is precomputed.
        coef_init: ndarray of shape n_features
            The initial coeffients to warm-start the optimization

        Notes
        -----

        Coordinate descent is an algorithm that considers each column of
        data at a time hence it will automatically convert the X input
        as a fortran contiguous numpy array if necessary.

        To avoid memory re-allocation it is advised to allocate the
        initial data in memory directly using that format.
        """
        # X and y must be of type float64
        X = np.asanyarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)

        n_samples, n_features = X.shape

        X_init = X
        X, y, X_mean, y_mean, X_std = self._center_data(X, y,
                self.fit_intercept, self.normalize, copy=self.copy_X)
        precompute = self.precompute
        if X_init is not X and hasattr(precompute, '__array__'):
            # recompute Gram
            # FIXME: it could be updated from precompute and X_mean
            # instead of recomputed
            precompute = 'auto'
        if X_init is not X and Xy is not None:
            Xy = None  # recompute Xy

        if coef_init is None:
            self.coef_ = np.zeros(n_features, dtype=np.float64)
        else:
            self.coef_ = coef_init

        alpha = self.alpha * self.rho * n_samples
        beta = self.alpha * (1.0 - self.rho) * n_samples

        X = np.asfortranarray(X)  # make data contiguous in memory

        # precompute if n_samples > n_features
        if hasattr(precompute, '__array__'):
            Gram = precompute
        elif precompute == True or \
               (precompute == 'auto' and n_samples > n_features):
            Gram = np.dot(X.T, X)
        else:
            Gram = None

        if Gram is None:
            self.coef_, self.dual_gap_, self.eps_ = \
                    cd_fast.enet_coordinate_descent(self.coef_, alpha, beta,
                                                    X, y, self.max_iter,
                                                    self.tol)
        else:
            if Xy is None:
                Xy = np.dot(X.T, y)
            self.coef_, self.dual_gap_, self.eps_ = \
                    cd_fast.enet_coordinate_descent_gram(self.coef_, alpha,
                                beta, Gram, Xy, y, self.max_iter, self.tol)

        self._set_intercept(X_mean, y_mean, X_std)

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want'
                          ' to increase the number of iterations')

        # return self for chaining fit and predict calls
        return self


###############################################################################
# Lasso model

class Lasso(ElasticNet):
    """Linear Model trained with L1 prior as regularizer (aka the Lasso)

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

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

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: int, optional
        The maximum number of iterations

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.


    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.Lasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
       normalize=False, precompute='auto', tol=0.0001)
    >>> print clf.coef_
    [ 0.85  0.  ]
    >>> print clf.intercept_
    0.15

    See also
    --------
    lars_path
    lasso_path
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """

    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 precompute='auto', copy_X=True, max_iter=1000,
                 tol=1e-4):
        super(Lasso, self).__init__(alpha=alpha, rho=1.0,
                            fit_intercept=fit_intercept, normalize=normalize,
                            precompute=precompute, copy_X=copy_X,
                            max_iter=max_iter, tol=tol)


###############################################################################
# Classes to store linear models along a regularization path

def lasso_path(X, y, eps=1e-3, n_alphas=100, alphas=None,
               precompute='auto', Xy=None, fit_intercept=True,
               normalize=False, copy_X=True, verbose=False,
               **params):
    """Compute Lasso path with coordinate descent

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

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

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    Xy : array-like, optional
        Xy = np.dot(X.T, y) that can be precomputed. It is useful
        only when the Gram matrix is precomputed.

    fit_intercept : bool
        Fit or not an intercept

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    verbose : bool or integer
        Amount of verbosity

    params : kwargs
        keyword arguments passed to the Lasso objects

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.

    See also
    --------
    lars_path
    Lasso
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode
    """
    return enet_path(X, y, rho=1., eps=eps, n_alphas=n_alphas, alphas=alphas,
                     precompute=precompute, Xy=Xy,
                     fit_intercept=fit_intercept, normalize=normalize,
                     copy_X=copy_X, verbose=verbose, **params)


def enet_path(X, y, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
              precompute='auto', Xy=None, fit_intercept=True,
              normalize=False, copy_X=True, verbose=False,
              **params):
    """Compute Elastic-Net path with coordinate descent

    The Elastic Net optimization function is::

        1 / (2 * n_samples) * ||y - Xw||^2_2 +
        + alpha * rho * ||w||_1 + 0.5 * alpha * (1 - rho) * ||w||^2_2

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

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    Xy : array-like, optional
        Xy = np.dot(X.T, y) that can be precomputed. It is useful
        only when the Gram matrix is precomputed.

    fit_intercept : bool
        Fit or not an intercept

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    verbose : bool or integer
        Amount of verbosity

    params : kwargs
        keyword arguments passed to the Lasso objects

    Returns
    -------
    models : a list of models along the regularization path

    Notes
    -----
    See examples/plot_lasso_coordinate_descent_path.py for an example.

    See also
    --------
    ElasticNet
    ElasticNetCV
    """
    X = as_float_array(X, copy_X)

    X_init = X
    X, y, X_mean, y_mean, X_std = LinearModel._center_data(X, y,
                                                           fit_intercept,
                                                           normalize,
                                                           copy=False)
    X = np.asfortranarray(X)  # make data contiguous in memory
    n_samples, n_features = X.shape

    if X_init is not X and hasattr(precompute, '__array__'):
        precompute = 'auto'
    if X_init is not X and Xy is not None:
        Xy = None

    if 'precompute' is True or \
                ((precompute == 'auto') and (n_samples > n_features)):
        precompute = np.dot(X.T, X)

    if Xy is None:
        Xy = np.dot(X.T, y)

    n_samples = X.shape[0]
    if alphas is None:
        alpha_max = np.abs(Xy).max() / (n_samples * rho)
        alphas = np.logspace(np.log10(alpha_max * eps), np.log10(alpha_max),
                             num=n_alphas)[::-1]
    else:
        alphas = np.sort(alphas)[::-1]  # make sure alphas are properly ordered
    coef_ = None  # init coef_
    models = []

    n_alphas = len(alphas)
    for i, alpha in enumerate(alphas):
        model = ElasticNet(alpha=alpha, rho=rho, fit_intercept=False,
                           precompute=precompute)
        model.set_params(**params)
        model.fit(X, y, coef_init=coef_, Xy=Xy)
        if fit_intercept:
            model.fit_intercept = True
            model._set_intercept(X_mean, y_mean, X_std)
        if verbose:
            if verbose > 2:
                print model
            elif verbose > 1:
                print 'Path: %03i out of %03i' % (i, n_alphas)
            else:
                sys.stderr.write('.')
        coef_ = model.coef_.copy()
        models.append(model)
    return models


class LinearModelCV(LinearModel):
    """Base class for iterative model fitting along a regularization path"""

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
            normalize=False, precompute='auto', max_iter=1000, tol=1e-4,
            copy_X=True, cv=None, verbose=False):
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.tol = tol
        self.copy_X = copy_X
        self.cv = cv
        self.verbose = verbose

    def fit(self, X, y):
        """Fit linear model with coordinate descent along decreasing alphas
        using cross-validation

        Parameters
        ----------

        X : numpy array of shape [n_samples,n_features]
            Training data. Pass directly as fortran contiguous data to avoid
            unnecessary memory duplication

        y : numpy array of shape [n_samples]
            Target values

        fit_params : kwargs
            keyword arguments passed to the Lasso fit method

        """
        X = np.asfortranarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)

        # All LinearModelCV parameters except 'cv' are acceptable
        path_params = self._get_params()
        del path_params['cv']

        # Start to compute path on full data
        models = self.path(X, y, **path_params)

        # Update the alphas list
        alphas = [model.alpha for model in models]
        n_alphas = len(alphas)
        path_params.update({'alphas': alphas, 'n_alphas': n_alphas})

        # init cross-validation generator
        cv = check_cv(self.cv, X)

        # Compute path for all folds and compute MSE to get the best alpha
        folds = list(cv)
        n_folds = len(folds)
        mse_alphas = np.zeros((n_folds, n_alphas))
        for i, (train, test) in enumerate(folds):
            if self.verbose:
                print '%s: fold % 2i out of % 2i' % (
                        self.__class__.__name__, i, n_folds),
                sys.stdout.flush()
            models_train = self.path(X[train], y[train], **path_params)
            for i_alpha, model in enumerate(models_train):
                y_ = model.predict(X[test])
                mse_alphas[i, i_alpha] += ((y_ - y[test]) ** 2).mean()
            if self.verbose == 1:
                print ''

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

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

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

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: int, optional
        The maximum number of iterations

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    verbose : bool or integer
        amount of verbosity

    Notes
    -----
    See examples/linear_model/lasso_path_with_crossvalidation.py
    for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.

    See also
    --------
    lars_path
    lasso_path
    LassoLars
    Lasso
    LassoLarsCV
    """
    path = staticmethod(lasso_path)
    estimator = Lasso


class ElasticNetCV(LinearModelCV):
    """Elastic Net model with iterative fitting along a regularization path

    The best model is selected by cross-validation.

    Parameters
    ----------
    rho : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). For rho = 0
        the penalty is an L1 penalty. For rho = 1 it is an L2 penalty.
        For 0 < rho < 1, the penalty is a combination of L1 and L2

    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3.

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : numpy array, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: int, optional
        The maximum number of iterations

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    verbose : bool or integer
        amount of verbosity

    Notes
    -----
    See examples/linear_model/lasso_path_with_crossvalidation.py
    for an example.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.

    The parameter rho corresponds to alpha in the glmnet R package
    while alpha corresponds to the lambda parameter in glmnet.
    More specifically, the optimization objective is::

        1 / (2 * n_samples) * ||y - Xw||^2_2 +
        + alpha * rho * ||w||_1 + 0.5 * alpha * (1 - rho) * ||w||^2_2

    If you are interested in controlling the L1 and L2 penalty
    separately, keep in mind that this is equivalent to::

        a * L1 + b * L2

    for::

        alpha = a + b and rho = a / (a + b)

    See also
    --------
    enet_path
    ElasticNet

    """
    path = staticmethod(enet_path)
    estimator = ElasticNet

    def __init__(self, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True, normalize=False, precompute='auto',
                 max_iter=1000, tol=1e-4, cv=None, copy_X=True,
                 verbose=0):
        self.rho = rho
        self.eps = eps
        self.n_alphas = n_alphas
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.tol = tol
        self.cv = cv
        self.copy_X = copy_X
        self.verbose = verbose
