# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Gael Varoquaux <gael.varoquaux@inria.fr>
#
# License: BSD Style.

import sys
import warnings
import itertools
import operator
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import sparse

from .base import LinearModel
from ..base import RegressorMixin
from .base import sparse_center_data, center_data
from ..utils import array2d, atleast2d_or_csc
from ..cross_validation import check_cv
from ..externals.joblib import Parallel, delayed
from ..utils.extmath import safe_sparse_dot

from . import cd_fast


###############################################################################
# ElasticNet model


class ElasticNet(LinearModel, RegressorMixin):
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
        the penalty is an L2 penalty. For rho = 1 it is an L1 penalty.
        For 0 < rho < 1, the penalty is a combination of L1 and L2

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered.

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument. For sparse input
        this option is always True to preserve sparsity.

    max_iter: int, optional
        The maximum number of iterations

    copy_X : boolean, optional, default False
        If True, X will be copied; else, it may be overwritten.

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    positive: bool, optional
        When set to True, forces the coefficients to be positive.

    Attributes
    ----------
    `coef_` : array, shape = (n_features,)
        parameter vector (w in the cost function formula)

    `sparse_coef_` : scipy.sparse matrix, shape = (n_features, 1)
        `sparse_coef_` is a readonly property derived from `coef_`

    `intercept_` : float | array, shape = (n_targets,)
        independent term in decision function.

    Notes
    -----
    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """
    def __init__(self, alpha=1.0, rho=0.5, fit_intercept=True,
                 normalize=False, precompute='auto', max_iter=1000,
                 copy_X=True, tol=1e-4, warm_start=False, positive=False):
        self.alpha = alpha
        self.rho = rho
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.positive = positive
        self.intercept_ = 0.0

    def fit(self, X, y, Xy=None, coef_init=None):
        """Fit Elastic Net model with coordinate descent

        Parameters
        -----------
        X: ndarray or scipy.sparse matrix, (n_samples, n_features)
            Data
        y: ndarray, shape = (n_samples,) or (n_samples, n_targets)
            Target
        Xy : array-like, optional
            Xy = np.dot(X.T, y) that can be precomputed. It is useful
            only when the Gram matrix is precomputed.
        coef_init: ndarray of shape n_features or (n_targets, n_features)
            The initial coeffients to warm-start the optimization

        Notes
        -----

        Coordinate descent is an algorithm that considers each column of
        data at a time hence it will automatically convert the X input
        as a fortran contiguous numpy array if necessary.

        To avoid memory re-allocation it is advised to allocate the
        initial data in memory directly using that format.
        """
        X = atleast2d_or_csc(X, dtype=np.float64, order='F',
                             copy=self.copy_X and self.fit_intercept)
        # From now on X can be touched inplace
        y = np.asarray(y, dtype=np.float64)
        # now all computation with X can be done inplace
        fit = self._sparse_fit if sparse.isspmatrix(X) else self._dense_fit
        fit(X, y, Xy, coef_init)
        return self

    def _dense_fit(self, X, y, Xy=None, coef_init=None):

        X, y, X_mean, y_mean, X_std = center_data(X, y,
                self.fit_intercept, self.normalize,
                copy=False)  # copy was done in fit if necessary

        if y.ndim == 1:
            y = y[:, np.newaxis]
        if Xy is not None and Xy.ndim == 1:
            Xy = Xy[:, np.newaxis]

        n_samples, n_features = X.shape
        n_targets = y.shape[1]

        precompute = self.precompute
        if hasattr(precompute, '__array__') \
                and not np.allclose(X_mean, np.zeros(n_features)) \
                and not np.allclose(X_std, np.ones(n_features)):
            # recompute Gram
            precompute = 'auto'
            Xy = None

        coef_ = self._init_coef(coef_init, n_features, n_targets)
        dual_gap_ = np.empty(n_targets)
        eps_ = np.empty(n_targets)

        l1_reg = self.alpha * self.rho * n_samples
        l2_reg = self.alpha * (1.0 - self.rho) * n_samples

        # precompute if n_samples > n_features
        if hasattr(precompute, '__array__'):
            Gram = precompute
        elif precompute == True or \
               (precompute == 'auto' and n_samples > n_features):
            Gram = np.dot(X.T, X)
        else:
            Gram = None

        for k in xrange(n_targets):
            if Gram is None:
                coef_[k, :], dual_gap_[k], eps_[k] = \
                        cd_fast.enet_coordinate_descent(coef_[k, :],
                        l1_reg, l2_reg, X, y[:, k], self.max_iter, self.tol,
                        self.positive)
            else:
                Gram = Gram.copy()
                if Xy is None:
                    this_Xy = np.dot(X.T, y[:, k])
                else:
                    this_Xy = Xy[:, k]
                coef_[k, :], dual_gap_[k], eps_[k] = \
                        cd_fast.enet_coordinate_descent_gram(coef_[k, :],
                        l1_reg, l2_reg, Gram, this_Xy, y[:, k], self.max_iter,
                        self.tol, self.positive)

            if dual_gap_[k] > eps_[k]:
                warnings.warn('Objective did not converge for ' +
                              'target %d, you might want' % k +
                              ' to increase the number of iterations')

        self.coef_, self.dual_gap_, self.eps_ = (np.squeeze(a) for a in (
                                                    coef_, dual_gap_, eps_))
        self._set_intercept(X_mean, y_mean, X_std)

        # return self for chaining fit and predict calls
        return self

    def _sparse_fit(self, X, y, Xy=None, coef_init=None):

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "Note: Sparse matrices cannot be indexed w/" +
                             "boolean masks (use `indices=True` in CV).")

        # NOTE: we are explicitly not centering the data the naive way to
        # avoid breaking the sparsity of X
        X_data, y, X_mean, y_mean, X_std = sparse_center_data(X, y,
                                                       self.fit_intercept,
                                                       self.normalize)

        if y.ndim == 1:
            y = y[:, np.newaxis]

        n_samples, n_features = X.shape[0], X.shape[1]
        n_targets = y.shape[1]

        coef_ = self._init_coef(coef_init, n_features, n_targets)
        dual_gap_ = np.empty(n_targets)
        eps_ = np.empty(n_targets)

        l1_reg = self.alpha * self.rho * n_samples
        l2_reg = self.alpha * (1.0 - self.rho) * n_samples

        for k in xrange(n_targets):
            coef_[k, :], dual_gap_[k], eps_[k] = \
                    cd_fast.sparse_enet_coordinate_descent(
                        coef_[k, :], l1_reg, l2_reg, X_data, X.indices,
                        X.indptr, y[:, k], X_mean / X_std,
                        self.max_iter, self.tol, self.positive)

            if dual_gap_[k] > eps_[k]:
                warnings.warn('Objective did not converge for ' +
                              'target %d, you might want' % k +
                              ' to increase the number of iterations')

        self.coef_, self.dual_gap_, self.eps_ = (np.squeeze(a) for a in (
                                                    coef_, dual_gap_, eps_))
        self._set_intercept(X_mean, y_mean, X_std)

        # return self for chaining fit and predict calls
        return self

    def _init_coef(self, coef_init, n_features, n_targets):
        if coef_init is None:
            if not self.warm_start or self.coef_ is None:
                coef_ = np.zeros((n_targets, n_features), dtype=np.float64)
            else:
                coef_ = self.coef_
        else:
            coef_ = coef_init

        if coef_.ndim == 1:
            coef_ = coef_[np.newaxis, :]
        if coef_.shape != (n_targets, n_features):
            raise ValueError("X and coef_init have incompatible " +
                 "shapes (%s != %s)." % (coef_.shape, (n_targets, n_features)))

        return coef_

    @property
    def sparse_coef_(self):
        """ sparse representation of the fitted coef """
        return sparse.csr_matrix(self.coef_)

    def decision_function(self, X):
        """Decision function of the linear model

        Parameters
        ----------
        X : numpy array or scipy.sparse matrix of shape (n_samples, n_features)

        Returns
        -------
        T : array, shape = (n_samples,)
            The predicted decision function
        """
        if sparse.isspmatrix(X):
            return np.ravel(safe_sparse_dot(self.coef_, X.T, \
                                        dense_output=True) + self.intercept_)
        else:
            return super(ElasticNet, self).decision_function(X)


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
        matrix can also be passed as argument. For sparse input
        this option is always True to preserve sparsity.

    max_iter: int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    positive : bool, optional
        When set to True, forces the coefficients to be positive.


    Attributes
    ----------
    `coef_` : array, shape = (n_features,)
        parameter vector (w in the cost function formula)

    `sparse_coef_` : scipy.sparse matrix, shape = (n_features, 1)
        `sparse_coef_` is a readonly property derived from `coef_`

    `intercept_` : float
        independent term in decision function.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.Lasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
    Lasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
       normalize=False, positive=False, precompute='auto', tol=0.0001,
       warm_start=False)
    >>> print(clf.coef_)
    [ 0.85  0.  ]
    >>> print(clf.intercept_)
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
                 tol=1e-4, warm_start=False, positive=False):
        super(Lasso, self).__init__(alpha=alpha, rho=1.0,
                            fit_intercept=fit_intercept, normalize=normalize,
                            precompute=precompute, copy_X=copy_X,
                            max_iter=max_iter, tol=tol, warm_start=warm_start,
                            positive=positive)


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
    X : ndarray, shape = (n_samples, n_features)
        Training data. Pass directly as fortran contiguous data to avoid
        unnecessary memory duplication

    y : ndarray, shape = (n_samples,)
        Target values

    eps : float, optional
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : ndarray, optional
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
    See examples/linear_model/plot_lasso_coordinate_descent_path.py
    for an example.

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
    X : ndarray, shape = (n_samples, n_features)
        Training data. Pass directly as fortran contiguous data to avoid
        unnecessary memory duplication

    y : ndarray, shape = (n_samples,)
        Target values

    rho : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). rho=1 corresponds to the Lasso

    eps : float
        Length of the path. eps=1e-3 means that
        alpha_min / alpha_max = 1e-3

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : ndarray, optional
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
    X = atleast2d_or_csc(X, dtype=np.float64, order='F',
                        copy=copy_X and fit_intercept)
    # From now on X can be touched inplace
    if not sparse.isspmatrix(X):
        X, y, X_mean, y_mean, X_std = center_data(X, y, fit_intercept,
                                                  normalize, copy=False)
        # XXX : in the sparse case the data will be centered
        # at each fit...

    n_samples, n_features = X.shape

    if hasattr(precompute, '__array__') \
            and not np.allclose(X_mean, np.zeros(n_features)) \
            and not np.allclose(X_std, np.ones(n_features)):
        # recompute Gram
        precompute = 'auto'
        Xy = None

    if precompute is True or \
                ((precompute == 'auto') and (n_samples > n_features)):
        if sparse.isspmatrix(X):
            warnings.warn("precompute is ignored for sparse data")
            precompute = False
        else:
            precompute = np.dot(X.T, X)

    if Xy is None:
        Xy = safe_sparse_dot(X.T, y, dense_output=True)

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
        model = ElasticNet(alpha=alpha, rho=rho,
            fit_intercept=fit_intercept if sparse.isspmatrix(X) else False,
            precompute=precompute)
        model.set_params(**params)
        model.fit(X, y, coef_init=coef_, Xy=Xy)
        if fit_intercept and not sparse.isspmatrix(X):
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


def _path_residuals(X, y, train, test, path, path_params, rho=1):
    this_mses = list()
    if 'rho' in path_params:
        path_params['rho'] = rho
    models_train = path(X[train], y[train], **path_params)
    this_mses = np.empty(len(models_train))
    for i_model, model in enumerate(models_train):
        y_ = model.predict(X[test])
        this_mses[i_model] = ((y_ - y[test]) ** 2).mean()
    return this_mses, rho


class LinearModelCV(LinearModel):
    """Base class for iterative model fitting along a regularization path"""
    __metaclass__ = ABCMeta

    @abstractmethod
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
        """Fit linear model with coordinate descent

        Fit is on grid of alphas and best alpha estimated by cross-validation.

        Parameters
        ----------

        X : array-like, shape (n_samples, n_features)
            Training data. Pass directly as fortran contiguous data to avoid
            unnecessary memory duplication

        y : narray, shape (n_samples,) or (n_samples, n_targets)
            Target values

        """
        X = atleast2d_or_csc(X, dtype=np.float64, order='F',
                             copy=self.copy_X and self.fit_intercept)
        # From now on X can be touched inplace
        y = np.asarray(y, dtype=np.float64)
        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have inconsistent dimensions (%d != %d)"
                              % (X.shape[0], y.shape[0]))

        # All LinearModelCV parameters except 'cv' are acceptable
        path_params = self.get_params()
        if 'rho' in path_params:
            rhos = np.atleast_1d(path_params['rho'])
            # For the first path, we need to set rho
            path_params['rho'] = rhos[0]
        else:
            rhos = [1, ]
        path_params.pop('cv', None)
        path_params.pop('n_jobs', None)

        # Start to compute path on full data
        # XXX: is this really useful: we are fitting models that we won't
        # use later
        models = self.path(X, y, **path_params)

        # Update the alphas list
        alphas = [model.alpha for model in models]
        n_alphas = len(alphas)
        path_params.update({'alphas': alphas, 'n_alphas': n_alphas})

        # init cross-validation generator
        cv = check_cv(self.cv, X)

        # Compute path for all folds and compute MSE to get the best alpha
        folds = list(cv)
        best_mse = np.inf
        all_mse_paths = list()

        # We do a double for loop folded in one, in order to be able to
        # iterate in parallel on rho and folds
        for rho, mse_alphas in itertools.groupby(
                    Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                        delayed(_path_residuals)(X, y, train, test,
                                    self.path, path_params, rho=rho)
                            for rho in rhos for train, test in folds
                    ), operator.itemgetter(1)):

            mse_alphas = [m[0] for m in mse_alphas]
            mse_alphas = np.array(mse_alphas)
            mse = np.mean(mse_alphas, axis=0)
            i_best_alpha = np.argmin(mse)
            this_best_mse = mse[i_best_alpha]
            all_mse_paths.append(mse_alphas.T)
            if this_best_mse < best_mse:
                model = models[i_best_alpha]
                best_rho = rho

        if hasattr(model, 'rho'):
            if model.rho != best_rho:
                # Need to refit the model
                model.rho = best_rho
                model.fit(X, y)
            self.rho_ = model.rho
        self.coef_ = model.coef_
        self.intercept_ = model.intercept_
        self.alpha_ = model.alpha
        self.alphas_ = np.asarray(alphas)
        self.coef_path_ = np.asarray([model.coef_ for model in models])
        self.mse_path_ = np.squeeze(all_mse_paths)
        return self

    @property
    def alpha(self):
        warnings.warn("Use alpha_. Using alpha is deprecated"
                "since version 0.12, and backward compatibility "
                "won't be maintained from version 0.14 onward. ",
                DeprecationWarning, stacklevel=1)
        return self.alpha_


class LassoCV(LinearModelCV, RegressorMixin):
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

    Attributes
    ----------
    `alpha_`: float
        The amount of penalization choosen by cross validation

    `coef_` : array, shape = (n_features,)
        parameter vector (w in the cost function formula)

    `intercept_` : float
        independent term in decision function.

    `mse_path_`: array, shape = (n_alphas, n_folds)
        mean square error for the test set on each fold, varying alpha

    `alphas_`: numpy array
        The grid of alphas used for fitting

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
    n_jobs = 1

    def __init__(self, eps=1e-3, n_alphas=100, alphas=None, fit_intercept=True,
            normalize=False, precompute='auto', max_iter=1000, tol=1e-4,
            copy_X=True, cv=None, verbose=False):
        super(LassoCV, self).__init__(eps=eps, n_alphas=n_alphas,
                alphas=alphas, fit_intercept=fit_intercept,
                normalize=normalize, precompute=precompute, max_iter=max_iter,
                tol=tol, copy_X=copy_X, cv=cv, verbose=verbose)


class ElasticNetCV(LinearModelCV, RegressorMixin):
    """Elastic Net model with iterative fitting along a regularization path

    The best model is selected by cross-validation.

    Parameters
    ----------
    rho : float, optional
        float between 0 and 1 passed to ElasticNet (scaling between
        l1 and l2 penalties). For rho = 0
        the penalty is an L2 penalty. For rho = 1 it is an L1 penalty.
        For 0 < rho < 1, the penalty is a combination of L1 and L2
        This parameter can be a list, in which case the different
        values are tested by cross-validation and the one giving the best
        prediction score is used. Note that a good choice of list of
        values for rho is often to put more values close to 1
        (i.e. Lasso) and less close to 0 (i.e. Ridge), as in [.1, .5, .7,
        .9, .95, .99, 1]

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

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
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

    n_jobs : integer, optional
        Number of CPUs to use during the cross validation. If '-1', use
        all the CPUs. Note that this is used only if multiple values for
        rho are given.

    Attributes
    ----------
    `alpha_` : float
        The amount of penalization choosen by cross validation

    `rho_` : float
        The compromise between l1 and l2 penalization choosen by
        cross validation

    `coef_` : array, shape = (n_features,)
        parameter vector (w in the cost function formula)

    `intercept_` : float
        independent term in decision function.

    `mse_path_` : array, shape = (n_rho, n_alpha, n_folds)
        mean square error for the test set on each fold, varying rho and
        alpha

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

    def __init__(self, rho=0.5, eps=1e-3, n_alphas=100, alphas=None,
                 fit_intercept=True, normalize=False, precompute='auto',
                 max_iter=1000, tol=1e-4, cv=None, copy_X=True,
                 verbose=0, n_jobs=1):
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
        self.n_jobs = n_jobs


###############################################################################
# Multi Task ElasticNet and Lasso models (with joint feature selection)

class MultiTaskElasticNet(Lasso):
    """Multi-task ElasticNet model trained with L1/L2 mixed-norm as regularizer

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||Y - XW||^Fro_2
        + alpha * rho * ||W||_21 + 0.5 * alpha * (1 - rho) * ||W||_Fro^2

    Where::

        ||W||_21 = \sum_i \sqrt{\sum_j w_{ij}^2}

    i.e. the sum of norm of earch row.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1/L2 term. Defaults to 1.0

    rho : float
        The ElasticNet mixing parameter, with 0 < rho <= 1. For rho = 0
        the penalty is an L1/L2 penalty. For rho = 1 it is an L1 penalty.
        For 0 < rho < 1, the penalty is a combination of L1/L2 and L2

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    `intercept_` : array, shape = (n_tasks,)
        Independent term in decision function.

    `coef_` : array, shape = (n_tasks, n_features)
        Parameter vector (W in the cost function formula). If a 1D y is \
        passed in at fit (non multi-task usage), `coef_` is then a 1D array

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.MultiTaskElasticNet(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [[0, 0], [1, 1], [2, 2]])
    ... #doctest: +NORMALIZE_WHITESPACE
    MultiTaskElasticNet(alpha=0.1, copy_X=True, fit_intercept=True,
            max_iter=1000, normalize=False, rho=0.5, tol=0.0001,
            warm_start=False)
    >>> print clf.coef_
    [[ 0.45663524  0.45612256]
     [ 0.45663524  0.45612256]]
    >>> print clf.intercept_
    [ 0.0872422  0.0872422]

    See also
    --------
    ElasticNet, MultiTaskLasso

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """
    def __init__(self, alpha=1.0, rho=0.5, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=1000, tol=1e-4, warm_start=False):
        self.alpha = alpha
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.rho = rho

    def fit(self, X, y, Xy=None, coef_init=None):
        """Fit MultiTaskLasso model with coordinate descent

        Parameters
        -----------
        X: ndarray, shape = (n_samples, n_features)
            Data
        y: ndarray, shape = (n_samples, n_tasks)
            Target
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
        X = array2d(X, dtype=np.float64, order='F',
                    copy=self.copy_X and self.fit_intercept)
        y = np.asarray(y, dtype=np.float64)

        squeeze_me = False
        if y.ndim == 1:
            squeeze_me = True
            y = y[:, np.newaxis]

        n_samples, n_features = X.shape
        _, n_tasks = y.shape

        X, y, X_mean, y_mean, X_std = center_data(X, y,
                self.fit_intercept, self.normalize, copy=False)

        if coef_init is None:
            if not self.warm_start or self.coef_ is None:
                self.coef_ = np.zeros((n_tasks, n_features), dtype=np.float64,
                                       order='F')
        else:
            self.coef_ = coef_init

        l1_reg = self.alpha * self.rho * n_samples
        l2_reg = self.alpha * (1.0 - self.rho) * n_samples

        self.coef_ = np.asfortranarray(self.coef_)  # coef contiguous in memory

        self.coef_, self.dual_gap_, self.eps_ = \
                cd_fast.enet_coordinate_descent_multi_task(self.coef_, l1_reg,
                        l2_reg, X, y, self.max_iter, self.tol)

        self._set_intercept(X_mean, y_mean, X_std)

        # Make sure that the coef_ have the same shape as the given 'y',
        # to predict with the same shape
        if squeeze_me:
            self.coef_ = self.coef_.squeeze()

        if self.dual_gap_ > self.eps_:
            warnings.warn('Objective did not converge, you might want'
                          ' to increase the number of iterations')

        # return self for chaining fit and predict calls
        return self


class MultiTaskLasso(MultiTaskElasticNet):
    """Multi-task Lasso model trained with L1/L2 mixed-norm as regularizer

    The optimization objective for Lasso is::

        (1 / (2 * n_samples)) * ||Y - XW||^2_Fro + alpha * ||W||_21

    Where::

        ||W||_21 = \sum_i \sqrt{\sum_j w_{ij}^2}

    i.e. the sum of norm of earch row.

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1/L2 term. Defaults to 1.0

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        The maximum number of iterations

    tol : float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

    warm_start : bool, optional
        When set to True, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.

    Attributes
    ----------
    `coef_` : array, shape = (n_tasks, n_features)
        parameter vector (W in the cost function formula)

    `intercept_` : array, shape = (n_tasks,)
        independent term in decision function.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.MultiTaskLasso(alpha=0.1)
    >>> clf.fit([[0,0], [1, 1], [2, 2]], [[0, 0], [1, 1], [2, 2]])
    MultiTaskLasso(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
            normalize=False, tol=0.0001, warm_start=False)
    >>> print clf.coef_
    [[ 0.89393398  0.        ]
     [ 0.89393398  0.        ]]
    >>> print clf.intercept_
    [ 0.10606602  0.10606602]

    See also
    --------
    Lasso, MultiTaskElasticNet

    Notes
    -----
    The algorithm used to fit the model is coordinate descent.

    To avoid unnecessary memory duplication the X argument of the fit method
    should be directly passed as a fortran contiguous numpy array.
    """
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=1000, tol=1e-4, warm_start=False):
        self.alpha = alpha
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.max_iter = max_iter
        self.copy_X = copy_X
        self.tol = tol
        self.warm_start = warm_start
        self.rho = 1.0
