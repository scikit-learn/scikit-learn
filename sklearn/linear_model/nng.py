"""
Non-Negative Garrote
"""

# Author: Jaques Grobler <jaques.grobler@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#
# License: BSD Style.


import numpy as np

from ..utils import array2d

from .base import LinearModel
from ..linear_model import LinearRegression, Lasso, lars_path, lasso_path


def non_negative_garrote(X, y, alpha, tol=0.001, fit_intercept=False,
                         normalize=True, max_iter=1000, precompute='auto'):
    """Function that implements the Non-negative garrote method

    Parameters
    -----------
    X : array, shape: (n_samples, n_features)
        Input data.

    y : array, shape: (n_samples)
        Input targets.


    tol: float, optional
        The tolerance for the optimization

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional, default False
        If ``True``, the regressors X will be normalized before regression.

    max_iter : integer, optional
        Maximum number of iterations to perform, set to infinity for no limit.

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to ``'auto'`` let us decide. The Gram
        matrix can also be passed as argument.

    Returns
    --------
    coefs: array, shape (n_features, max_features + 1)
        Coefficients along the path

    shrink_coef: array, shape (n_features, max_features + 1)
        Shrinkage coefficients obtained using Lasso

    See also
    --------
    lasso_path
    LassoLars
    Lars
    LassoLarsCV
    LarsCV
    sklearn.decomposition.sparse_encode

    Notes
    ------
    * http://libra.msra.cn/Publication/2083384/better-subset-regression-using-the-nonnegative-garrote

    """
    # Obtain the ordinary least squares coefficients from our data
    coef_ols = LinearRegression(fit_intercept=fit_intercept).fit(X, y).coef_

    X = X * coef_ols[np.newaxis, :]
    # find the shrinkage factor by minimising the sum of square residuals
    # under the restriction that it is positive (positive=True)
    shrink_coef = Lasso(alpha=alpha, fit_intercept=fit_intercept,
                        max_iter=max_iter,
                        positive=True, precompute=precompute,
                        normalize=normalize,
                        tol=tol).fit(X, y).coef_
    # Shrunken betas
    coef = coef_ols * shrink_coef
    return coef, shrink_coef


def non_negative_garrote_path(X, y, n_alphas=100, alphas=None,
                              precompute='auto', fit_intercept=False,
                              **params):
    """
    Compute the Non-negative Garrote path

    Parameters
    ----------
    X : ndarray, shape = (n_samples, n_features)
        Training data. Pass directly as fortran contiguous data to avoid
        unnecessary memory duplication

    y : ndarray, shape = (n_samples,)
        Target values

    n_alphas : int, optional
        Number of alphas along the regularization path

    alphas : ndarray, optional
        List of alphas where to compute the models.
        If None alphas are set automatically

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    fit_intercept : bool
        Fit or not an intercept

    params : kwargs
        keyword arguments passed to the Lasso objects

    Returns
    -------
    coef_path: array, shape (n_features, max_features + 1)
        Coefficients along the path

    shrink_coef_path: array, shape (n_features, max_features + 1)
        Shrinkage coefficients along the path using lars_path

    See also
    --------
    lars_path
    Lasso
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode
    """

    # Obtain the ordinary least squares coefficients from our data
    coef_ols = LinearRegression(fit_intercept=fit_intercept).fit(X, y).coef_

    X = X * coef_ols[np.newaxis, :]
    _, _, shrink_coef_path = lars_path(X, y, method='lasso')

    # Lasso_path version
    #shrink_coef_path = lasso_path(X, y, positive=True, n_alphas=n_alphas)[1].coef_
    #coef_path = coef_ols * shrink_coef_path

    # Shrunken betas
    coef_path = shrink_coef_path * coef_ols[:, np.newaxis]
    return coef_path, shrink_coef_path


class NonNegativeGarrote(LinearModel):
    """
    A method for subset regression, it minimizes the same equation as
    the Lasso, under the constraint that the coefficients always be positive.

    Ref:
    Breiman, L. (1995), "Better Subset Regression Using the Nonnegative
    Garrote," Technometrics, 37, 373-384. [349,351]

    Parameters
    ----------
    alpha : float, optional
        Constant that multiplies the L1 term. Defaults to 1.0
        alpha = 0 is equivalent to an ordinary least square, solved
        by the LinearRegression object in the scikit.

    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    tol: float, optional
        The tolerance for the optimization: if the updates are
        smaller than 'tol', the optimization code checks the
        dual gap for optimality and continues until it is smaller
        than tol.

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

    Attributes
    ----------
    ``coef_`` : array, shape = (n_features,)
        The parameter vector

    ``shrink_coef_`` : array, shape = (n_features,)
        The shrinkage coefficients that must always be positive.

    Examples
    --------
    >>> from sklearn import linear_model
    >>> clf = linear_model.NonNegativeGarrote(alpha=0.1)
    >>> clf.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    NonNegativeGarrote(alpha=0.1, copy_X=True, fit_intercept=True, max_iter=1000,
       normalize=False, precompute='auto', tol=0.0001)
    >>> clf.coef_ # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([ 0.0026...,  0.        ])
    >>> clf.intercept_ # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    0.36...

    See also
    --------
    lars_path
    Lasso
    LassoLars
    LassoCV
    LassoLarsCV
    sklearn.decomposition.sparse_encode
    """
    def __init__(self, alpha=0.35, fit_intercept=True, max_iter=1000,
                 tol=1e-4, normalize=False,
                 copy_X=True, precompute='auto'):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.tol = tol
        self.normalize = normalize
        self.copy_X = copy_X
        self.precompute = precompute

    def fit(self, X, y):
        """Fit the model using X, y as training data.

        parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            training data.

        y : array-like, shape = [n_samples]
            target values.

        returns
        -------
        self : object
            returns an instance of self.
        """
        # X and y must be of type float64
        X = array2d(X, dtype=np.float64, order='F',
                    copy=self.copy_X and self.fit_intercept)
        y = np.asarray(y, dtype=np.float64)

        X, y, X_mean, y_mean, X_std =\
            LinearModel._center_data(X, y, self.fit_intercept,
                                     self.normalize, self.copy_X)

        self.coef_, self.shrink_coef_ = \
            non_negative_garrote(X, y, self.alpha,
                                 self.tol, self.fit_intercept, self.normalize,
                                 self.max_iter, self.precompute)
        self._set_intercept(X_mean, y_mean, X_std)

        return self
