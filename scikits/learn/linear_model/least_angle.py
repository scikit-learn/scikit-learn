"""
Least Angle Regression algorithm. See the documentation on the
Generalized Linear Model for a complete discussion.
"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#
# License: BSD Style.

import numpy as np
from scipy import linalg, interpolate
from scipy.linalg.lapack import get_lapack_funcs

from .base import LinearModel
from ..utils import arrayfuncs
from ..utils import deprecated
from ..cross_val import KFold
from ..externals.joblib import Parallel, delayed


def lars_path(X, y, Xy=None, Gram=None, max_iter=500,
              alpha_min=0, method='lar', overwrite_X=False,
              overwrite_Gram=False, verbose=False):
    """Compute Least Angle Regression and LASSO path

    Parameters
    -----------
    X: array, shape: (n_samples, n_features)
        Input data

    y: array, shape: (n_samples)
        Input targets

    max_iter: integer, optional
        Maximum number of iterations to perform, set to infinity for no limit.

    Gram: None, 'auto', array, shape: (n_features, n_features), optional
        Precomputed Gram matrix (X' * X), if 'auto', the Gram
        matrix is precomputed from the given X, if there are more samples
        than features

    alpha_min: float, optional
        Minimum correlation along the path. It corresponds to the
        regularization parameter alpha parameter in the Lasso.

    method: {'lar', 'lasso'}
        Specifies the returned model. Select 'lar' for Least Angle
        Regression, 'lasso' for the Lasso.

    Returns
    --------
    alphas: array, shape: (max_iter,)
        Maximum of covariances (in absolute value) at each
        iteration.

    active: array, shape (max_features,)
        Indices of active variables at the end of the path.

    coefs: array, shape (n_features, max_iter)
        Coefficients along the path

    See also
    --------
    :ref:`LassoLars`, :ref:`Lars`

    Notes
    ------
    * http://en.wikipedia.org/wiki/Least-angle_regression

    * http://en.wikipedia.org/wiki/Lasso_(statistics)#LASSO_method
    """

    n_features = X.shape[1]
    n_samples = y.size
    max_features = min(max_iter, n_features)

    coefs = np.zeros((max_features + 1, n_features))
    alphas = np.zeros(max_features + 1)
    n_iter, n_active = 0, 0
    active, indices = list(), np.arange(n_features)
    # holds the sign of covariance
    sign_active = np.empty(max_features, dtype=np.int8)
    drop = False
    eps = np.finfo(X.dtype).eps

    # will hold the cholesky factorization. Only lower part is
    # referenced.
    L = np.empty((max_features, max_features), dtype=X.dtype)
    swap, nrm2 = linalg.get_blas_funcs(('swap', 'nrm2'), (X,))
    potrs, = get_lapack_funcs(('potrs',), (X,))

    if Gram is None:
        if not overwrite_X:
            # force copy. setting the array to be fortran-ordered
            # speeds up the calculation of the (partial) Gram matrix
            # and allows to easily swap columns
            X = X.copy('F')
    elif Gram == 'auto':
        Gram = None
        if X.shape[0] > X.shape[1]:
            Gram = np.dot(X.T, X)
    else:
        if not overwrite_Gram:
            Gram = Gram.copy()

    if Xy is None:
        Cov = np.dot(X.T, y)
    else:
        Cov = Xy.copy()

    if verbose:
        print "Step\t\tAdded\t\tDropped\t\tActive set size\t\tC"

    while 1:

        if Cov.size:
            C_idx = np.argmax(np.abs(Cov))
            C_ = Cov[C_idx]
            C = np.fabs(C_)
        else:
            C = 0.

        alphas[n_iter] = C / n_samples
        if alphas[n_iter] < alpha_min:  # early stopping
            # interpolation factor 0 <= ss < 1
            if n_iter > 0:
                # In the first iteration, all alphas are zero, the formula
                # below would make ss a NaN
                ss = (alphas[n_iter - 1] - alpha_min) / (alphas[n_iter - 1] -
                                                    alphas[n_iter])
                coefs[n_iter] = coefs[n_iter - 1] + ss * (coefs[n_iter] -
                                coefs[n_iter - 1])
            alphas[n_iter] = alpha_min
            break

        if n_iter >= max_iter or n_active >= n_features:
            break

        if not drop:

            # Update the Cholesky factorization of (Xa * Xa') #
            #                                                 #
            #            ( L   0 )                            #
            #     L  ->  (       )  , where L * w = b         #
            #            ( w   z )    z = 1 - ||w||           #
            #                                                 #
            #   where u is the last added to the active set   #

            sign_active[n_active] = np.sign(C_)
            m, n = n_active, C_idx + n_active

            Cov[C_idx], Cov[0] = swap(Cov[C_idx], Cov[0])
            indices[n], indices[m] = indices[m], indices[n]
            Cov = Cov[1:]  # remove Cov[0]

            if Gram is None:
                X.T[n], X.T[m] = swap(X.T[n], X.T[m])
                c = nrm2(X.T[n_active]) ** 2
                L[n_active, :n_active] = \
                    np.dot(X.T[n_active], X.T[:n_active].T)
            else:
                # swap does only work inplace if matrix is fortran
                # contiguous ...
                Gram[m], Gram[n] = swap(Gram[m], Gram[n])
                Gram[:, m], Gram[:, n] = swap(Gram[:, m], Gram[:, n])
                c = Gram[n_active, n_active]
                L[n_active, :n_active] = Gram[n_active, :n_active]

            # Update the cholesky decomposition for the Gram matrix
            arrayfuncs.solve_triangular(L[:n_active, :n_active],
                                        L[n_active, :n_active])
            v = np.dot(L[n_active, :n_active], L[n_active, :n_active])
            diag = max(np.sqrt(np.abs(c - v)), eps)
            L[n_active,  n_active] = diag

            active.append(indices[n_active])
            n_active += 1

            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, active[-1], '',
                                                            n_active, C)

        # least squares solution
        least_squares, info = potrs(L[:n_active, :n_active],
                               sign_active[:n_active], lower=True)

        # is this really needed ?
        AA = 1. / np.sqrt(np.sum(least_squares * sign_active[:n_active]))
        least_squares *= AA

        if Gram is None:
            # equiangular direction of variables in the active set
            eq_dir = np.dot(X.T[:n_active].T, least_squares)
            # correlation between each unactive variables and
            # eqiangular vector
            corr_eq_dir = np.dot(X.T[n_active:], eq_dir)
        else:
            # if huge number of features, this takes 50% of time, I
            # think could be avoided if we just update it using an
            # orthogonal (QR) decomposition of X
            corr_eq_dir = np.dot(Gram[:n_active, n_active:].T,
                                 least_squares)

        g1 = arrayfuncs.min_pos((C - Cov) / (AA - corr_eq_dir))
        g2 = arrayfuncs.min_pos((C + Cov) / (AA + corr_eq_dir))
        gamma_ = min(g1, g2, C / AA)

        # TODO: better names for these variables: z
        drop = False
        z = - coefs[n_iter, active] / least_squares
        z_pos = arrayfuncs.min_pos(z)
        if z_pos < gamma_:

            # some coefficients have changed sign
            idx = np.where(z == z_pos)[0]

            # update the sign, important for LAR
            sign_active[idx] = -sign_active[idx]

            if method == 'lasso':
                gamma_ = z_pos
            drop = True

        n_iter += 1

        if n_iter >= coefs.shape[0]:
            # resize the coefs and alphas array
            add_features = 2 * max(1, (max_features - n_active))
            coefs.resize((n_iter + add_features, n_features))
            alphas.resize(n_iter + add_features)

        coefs[n_iter, active] = coefs[n_iter - 1, active] + \
                                gamma_ * least_squares

        # update correlations
        Cov -= gamma_ * corr_eq_dir

        # See if any coefficient has changed sign
        if drop and method == 'lasso':

            arrayfuncs.cholesky_delete(L[:n_active, :n_active], idx)

            n_active -= 1
            m, n = idx, n_active
            drop_idx = active.pop(idx)

            if Gram is None:
                # propagate dropped variable
                for i in range(idx, n_active):
                    X.T[i], X.T[i + 1] = swap(X.T[i], X.T[i + 1])
                    indices[i], indices[i + 1] =  \
                            indices[i + 1], indices[i]  # yeah this is stupid

                # TODO: this could be updated
                residual = y - np.dot(X[:, :n_active],
                                      coefs[n_iter, active])
                temp = np.dot(X.T[n_active], residual)

                Cov = np.r_[temp, Cov]
            else:
                for i in range(idx, n_active):
                    indices[i], indices[i + 1] =  \
                                indices[i + 1], indices[i]
                    Gram[i], Gram[i + 1] = swap(Gram[i], Gram[i + 1])
                    Gram[:, i], Gram[:, i + 1] = swap(Gram[:, i],
                                                      Gram[:, i + 1])

                # Cov_n = Cov_j + x_j * X + increment(betas) TODO:
                # will this still work with multiple drops ?

                # recompute covariance. Probably could be done better
                # wrong as Xy is not swapped with the rest of variables

                # TODO: this could be updated
                residual = y - np.dot(X, coefs[n_iter])
                temp = np.dot(X.T[drop_idx], residual)
                Cov = np.r_[temp, Cov]

            sign_active = np.delete(sign_active, idx)
            sign_active = np.append(sign_active, 0.)  # just to maintain size
            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, '', drop_idx,
                                                      n_active, abs(temp))

    # resize coefs in case of early stop
    alphas = alphas[:n_iter + 1]
    coefs = coefs[:n_iter + 1]

    return alphas, active, coefs.T


###############################################################################
# Estimator classes

class Lars(LinearModel):
    """Least Angle Regression model a.k.a. LAR

    Parameters
    ----------
    n_nonzero_coefs : int, optional
        Target number of non-zero coefficients. Use np.inf for no limit.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.


    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    Examples
    --------
    >>> from scikits.learn import linear_model
    >>> clf = linear_model.Lars(n_nonzero_coefs=1)
    >>> clf.fit([[-1,1], [0, 0], [1, 1]], [-1, 0, -1])
    Lars(normalize=True, precompute='auto', n_nonzero_coefs=1, verbose=False,
       fit_intercept=True)
    >>> print clf.coef_
    [ 0. -1.]

    References
    ----------
    http://en.wikipedia.org/wiki/Least_angle_regression

    See also
    --------
    lars_path, LassoLARS, LarsCV, LassoLarsCV
    """
    def __init__(self, fit_intercept=True, verbose=False, normalize=True,
                 precompute='auto', n_nonzero_coefs=500):
        self.fit_intercept = fit_intercept
        self.verbose = verbose
        self.normalize = normalize
        self.method = 'lar'
        self.precompute = precompute
        self.n_nonzero_coefs = n_nonzero_coefs

    def fit(self, X, y, overwrite_X=False, **params):
        """Fit the model using X, y as training data.

        parameters
        ----------
        x : array-like, shape = [n_samples, n_features]
            training data.

        y : array-like, shape = [n_samples]
            target values.

        returns
        -------
        self : object
            returns an instance of self.
        """
        self._set_params(**params)

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)
        alpha = getattr(self, 'alpha', 0.)
        if hasattr(self, 'n_nonzero_coefs'):
            alpha = 0. # n_nonzero_coefs parametrization takes priority
            max_iter = self.n_nonzero_coefs
        else:
            max_iter = self.max_iter

        if self.normalize:
            norms = np.sqrt(np.sum(X ** 2, axis=0))
            nonzeros = np.flatnonzero(norms)
            if not overwrite_X:
                X = X.copy()
                overwrite_X = True
            X[:, nonzeros] /= norms[nonzeros]

        # precompute if n_samples > n_features
        precompute = self.precompute
        if hasattr(precompute, '__array__'):
            # copy as it's going to be modified
            Gram = precompute.copy()
        elif precompute == 'auto':
            Gram = 'auto'
        else:
            Gram = None

        self.alphas_, self.active_, self.coef_path_ = lars_path(X, y,
                  Gram=Gram, overwrite_X=overwrite_X,
                  overwrite_Gram=True, alpha_min=alpha,
                  method=self.method, verbose=self.verbose,
                  max_iter=max_iter)

        if self.normalize:
            self.coef_path_ /= norms[:, np.newaxis]
        self.coef_ = self.coef_path_[:, -1]

        self._set_intercept(Xmean, ymean)

        return self


class LassoLars(Lars):
    """Lasso model fit with Least Angle Regression a.k.a. Lars

    It is a Linear Model trained with an L1 prior as regularizer.
    lasso).

    Parameters
    ----------
    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: integer, optional
        Maximum number of iterations to perform.


    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    Examples
    --------
    >>> from scikits.learn import linear_model
    >>> clf = linear_model.LassoLars(alpha=0.01)
    >>> clf.fit([[-1,1], [0, 0], [1, 1]], [-1, 0, -1])
    LassoLars(normalize=True, verbose=False, fit_intercept=True, max_iter=500,
         precompute='auto', alpha=0.01)
    >>> print clf.coef_
    [ 0.         -0.96325765]

    References
    ----------
    http://en.wikipedia.org/wiki/Least_angle_regression

    See also
    --------
    lars_path, Lasso
    """

    def __init__(self, alpha=1.0, fit_intercept=True, verbose=False,
                 normalize=True, precompute='auto', max_iter=500):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.verbose = verbose
        self.normalize = normalize
        self.method = 'lasso'
        self.precompute = precompute


# Deprecated classes
class LARS(Lars):
    pass
LARS = deprecated("Use Lars instead")(LARS)


class LassoLARS(LassoLars):
    pass
LassoLARS = deprecated("Use LassoLars instead")(LassoLARS)

################################################################################
# Cross-validated estimator classes

def _lars_path_residues(X_train, y_train, X_test, y_test, Gram=None,
                     overwrite_data=False, method='lars', verbose=False, 
                     fit_intercept=True, max_iter=500):
    """Compute the residues on left-out data for a full LARS path

    Returns
    --------
    alphas: array, shape: (max_features + 1,)
        Maximum of covariances (in absolute value) at each
        iteration.

    active: array, shape (max_features,)
        Indices of active variables at the end of the path.

    coefs: array, shape (n_features, max_features+1)
        Coefficients along the path

    residues: array, shape (n_features, max_features+1)
        Residues of the prediction on the test data
    """
    if not overwrite_data:
        X_train = X_train.copy()
        y_train = y_train.copy()
        X_test = X_test.copy()
        y_test = y_test.copy()
    if fit_intercept:
        X_mean = X_train.mean(axis=0)
        X_train -= X_mean
        X_test -= X_mean
        y_mean = y_train.mean(axis=0)
        y_train -= y_mean
        y_test -= y_mean
    alphas, active, coefs = lars_path(X_train, y_train, Gram=Gram, 
                            overwrite_X=True, overwrite_Gram=True,
                            method=method, verbose=verbose,
                            max_iter=max_iter)
    residues = np.array([(np.dot(X_test, coef) - y_test)
                         for coef in coefs.T])
    return alphas, active, coefs, residues


class LarsCV(LARS):
    """Cross-validated Least Angle Regression model

    Parameters
    ----------
    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: integer, optional
        Maximum number of iterations to perform.

    cv : crossvalidation generator, optional
        see scikits.learn.cross_val module. If None is passed, default to
        a 5-fold strategy

    n_jobs : integer, optional
        Number of CPUs to use during the cross validation. If '-1', use
        all the CPUs

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    `coef_path`: array, shape = [n_features, n_alpha]
        the varying values of the coefficients along the path

    See also
    --------
    lars_path, LassoLARS, LassoLarsCV
    """

    method = 'lar'

    def __init__(self, fit_intercept=True, verbose=False, max_iter=500, 
                 normalize=True, precompute='auto', cv=None, n_jobs=1):
        self.fit_intercept = fit_intercept
        self.max_iter = max_iter
        self.verbose = verbose
        self.normalize = normalize
        self.precompute = precompute
        self.cv = cv
        self.n_jobs = n_jobs

    def fit(self, X, y, **params):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training data.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
            returns an instance of self.
        """
        self._set_params(**params)

        n_samples = len(X)
        # init cross-validation generator
        cv = self.cv if self.cv else KFold(n_samples, 5)

        Gram = 'auto' if self.precompute else None

        cv_paths = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                    delayed(_lars_path_residues)(X[train], y[train], 
                            X[test], y[test], Gram=Gram, 
                            overwrite_data=True, method=self.method, 
                            verbose=max(0, self.verbose-1),
                            fit_intercept=self.fit_intercept,
                            max_iter=self.max_iter)
                    for train, test in cv)
        all_alphas = np.concatenate(zip(*cv_paths)[0])
        all_alphas.sort()

        mse_path = list()
        #coef_path = list()
        for alphas, active, coefs, residues in cv_paths:
            #this_coefs = interpolate.interp1d(alphas, coefs)(all_alphas)
            #coef_path.append(this_coefs)
            this_residues = interpolate.interp1d(alphas[::-1], 
                                                 residues[::-1], 
                                                 bounds_error=False,
                                                 fill_value=residues.max(),
                                                 axis=0)(all_alphas)
            this_residues **= 2
            this_residues = this_residues.sum(axis=-1)
            mse_path.append(this_residues)

        mse_path = np.array(mse_path).T
        mask = np.all(np.isfinite(mse_path), axis=-1)
        all_alphas = all_alphas[mask]
        mse_path = mse_path[mask]
        # Select the alpha that minimizes left-out error
        i_best_alpha = np.argmin(mse_path.mean(axis=-1))
        best_alpha = all_alphas[i_best_alpha]

        # Store our parameters
        self.alpha = best_alpha/n_samples
        self.cv_alphas = all_alphas
        self.cv_mse_path_ = mse_path

        # Now compute the full model
        LARS.fit(self, X, y)
        return self


class LassoLarsCV(LarsCV):
    """Cross-validated Lasso, using the LARS algorithm 

    Parameters
    ----------
    fit_intercept : boolean
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    verbose : boolean or integer, optional
        Sets the verbosity amount

    normalize : boolean, optional
        If True, the regressors X are normalized

    precompute : True | False | 'auto' | array-like
        Whether to use a precomputed Gram matrix to speed up
        calculations. If set to 'auto' let us decide. The Gram
        matrix can also be passed as argument.

    max_iter: integer, optional
        Maximum number of iterations to perform.

    cv : crossvalidation generator, optional
        see scikits.learn.cross_val module. If None is passed, default to
        a 5-fold strategy

    n_jobs : integer, optional
        Number of CPUs to use during the cross validation. If '-1', use
        all the CPUs

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    `coef_path`: array, shape = [n_features, n_alpha]
        the varying values of the coefficients along the path

    `alphas_`: array, shape = [n_alpha]
        the different values of alpha along the path

    `cv_alphas`: array, shape = [n_cv_alphas]
        all the values of alpha along the path for the different folds

    `cv_mse_path_`: array, shape = [n_folds, n_cv_alphas]
        the mean square error on left-out for each fold along the path
        (alpha values given by cv_alphas)

    See also
    --------
    lars_path, LassoLARS, LarsCV
    """


    method = 'lasso'

