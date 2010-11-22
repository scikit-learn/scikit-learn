"""
Least Angle Regression algorithm. See the documentation on the
Generalized Linear Model for a complete discussion.
"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#
# License: BSD Style.

import warnings
import numpy as np
from scipy import linalg
from scipy.linalg.lapack import get_lapack_funcs

from .base import LinearModel
from ..utils import arrayfuncs

def lars_path(X, y, Xy=None, Gram=None, max_features=None,
              alpha_min=0, method="lar", overwrite_X=False,
              overwrite_Gram=False, overwrite_y=False, verbose=False):

    """ Compute Least Angle Regression and LASSO path

        Parameters
        -----------
        X: array, shape: (n_samples, n_features)
            Input data

        y: array, shape: (n_samples)
            Input targets

        max_features: integer, optional
            Maximum number of selected features.

        Gram: array, shape: (n_features, n_features), optional
            Precomputed Gram matrix (X' * X)

        alpha_min: float, optional
            Minimum correlation along the path. It corresponds to the
            regularization parameter alpha parameter in the Lasso.

        method: 'lar' or 'lasso'
            Specifies the problem solved: the LAR or its variant the
            LASSO-LARS that gives the solution of the LASSO problem
            for any regularization parameter.

        Returns
        --------
        alphas: array, shape: (k)
            The alphas along the path

        active: array, shape (?)
            Indices of active variables at the end of the path.

        coefs: array, shape (p,k)
            Coefficients along the path

        Notes
        ------
        * http://en.wikipedia.org/wiki/Least-angle_regression

        * http://en.wikipedia.org/wiki/Lasso_(statistics)#LASSO_method
    """
    # : make sure it works with non-normalized columns of X

    n_samples, n_features = X.shape

    if max_features is None:
        max_features = min(n_samples, n_features)

    coefs = np.zeros((max_features + 1, n_features))
    alphas = np.zeros(max_features + 1)
    n_iter, n_active = 0, 0
    active, indices = list(), np.arange(n_features)
    # holds the sign of covariance
    sign_active = np.empty(max_features, dtype=np.int8)
    drop = False

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
    else:
        if not overwrite_Gram:
            Gram = Gram.copy()

    if Xy is None:
        Cov = np.dot(X.T, y)
    else:
        Cov = Xy.copy()

    if (Gram is None) and (not overwrite_y):
        y = y.copy()

    if verbose:
        print "Step\t\tAdded\t\tDropped\t\tActive set size\t\tC"

    while 1:

        if Cov.size:
            imax = np.argmax(np.abs(Cov)) # TODO: rename
            C_ = Cov[imax]
            C = np.fabs(C_)
            # to match a for computing gamma_
        else:
            if Gram is None:
                C -= gamma_ * np.abs(np.dot(X.T[0], eq_dir))
            else:
                C -= gamma_ * np.abs(np.dot(Gram[0], least_squares))

        alphas[n_iter] = C

        if (C < alpha_min) or (n_active == max_features): break

        if not drop:

            # Update the Cholesky factorization of (Xa * Xa') #
            #                                                 #
            #            ( L   0 )                            #
            #     L  ->  (       )  , where L * w = b         #
            #            ( w   z )    z = 1 - ||w||           #
            #                                                 #
            #   where u is the last added to the active set   #

            sign_active[n_active] = np.sign(C_)
            m, n = n_active, imax+n_active

            Cov[imax], Cov[0] = swap(Cov[imax], Cov[0])
            indices[n], indices[m] = indices[m], indices[n]
            Cov = Cov[1:] # remove Cov[0]

            if Gram is None:
                X.T[n], X.T[m] = swap(X.T[n], X.T[m])
                c = nrm2(X.T[n_active])**2
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
            arrayfuncs.solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            v = np.dot(L[n_active, :n_active], L[n_active, :n_active])
            L[n_active,  n_active] = np.sqrt(c - v)

            active.append(indices[n_active])
            n_active += 1

            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, active[-1], '',
                                                            n_active, C)

        # least squares solution

        least_squares, info = potrs(L[:n_active, :n_active],
                               sign_active[:n_active], lower=True)

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

        # equation 2.13
        g1 = arrayfuncs.min_pos((C - Cov) / (AA - corr_eq_dir))
        g2 = arrayfuncs.min_pos((C + Cov) / (AA + corr_eq_dir))
        gamma_ = min(g1, g2, C/AA)

        if method == 'lasso':
            drop = False
            z = - coefs[n_iter, active] / least_squares
            z_pos = arrayfuncs.min_pos(z)
            if z_pos < gamma_:
                idx = np.where(z == z_pos)[0]
                gamma_ = z_pos
                drop = True

        n_iter += 1

        if n_iter >= coefs.shape[0]: # resize
            add_features = 2 * (max_features - n_active) # heuristic
            coefs.resize(n_iter + add_features, n_features)
            alphas.resize(n_iter + add_features)

        coefs[n_iter, active] = coefs[n_iter-1, active] + \
                                gamma_ * least_squares

        # update correlations
        Cov -= gamma_ * corr_eq_dir

        if n_active > n_features:
            break

        if drop:
            # got to add corresponding covariance 
            arrayfuncs.cholesky_delete(L[:n_active, :n_active], idx)

            n_active -= 1
            m, n = idx, n_active
            drop_idx = active.pop(idx)

            if Gram is None:
                # propagate dropped variable
                for i in range(idx, n_active):
                    X.T[i], X.T[i+1] = swap(X.T[i], X.T[i+1])
                    indices[i], indices[i+1] =  \
                                indices[i+1], indices[i] # yeah this is stupid

                # TODO: this could be updated
                residual = y - np.dot(X[:, :n_active],
                                      coefs[n_iter, active])
                temp = np.dot(X.T[n_active], residual)
                Cov = np.r_[temp, Cov]
            else:
                for i in range(idx, n_active):
                    indices[i], indices[i+1] =  \
                                indices[i+1], indices[i]
                    Gram[i], Gram[i+1] = swap(Gram[i], Gram[i+1])
                    Gram[:, i], Gram[:, i+1] = swap(Gram[:, i], Gram[:, i+1])

                # Cov_n = Cov_j + x_j * X + increment(betas) TODO:
                # will this still work with multiple drops ?

                # recompute covariance. Probably could be done better
                # wrong as Xy is not swapped with the rest of variables

                # TODO: this could be updated
                residual = y - np.dot(X, coefs[n_iter])
                temp = np.dot(X.T[drop_idx], residual)
                Cov = np.r_[temp, Cov]

            # do an append to maintain size
            sign_active = np.delete(sign_active, idx)
            sign_active = np.append(sign_active, 0.)
            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, '', drop_idx,

                                                            n_active, C)
    if C < alpha_min: # interpolate
        # interpolation factor 0 <= ss < 1
        ss = (alphas[n_iter-1] - alpha_min) / (alphas[n_iter-1] - alphas[n_iter])
        coefs[n_iter] = coefs[n_iter-1] + ss*(coefs[n_iter] - coefs[n_iter-1])
        alphas[n_iter] = alpha_min


    # resize coefs in case of early stop
    alphas = alphas[:n_iter+1]
    coefs = coefs[:n_iter+1]

    return alphas, active, coefs.T


class LARS(LinearModel):
    """Least Angle Regression model a.k.a. LAR

    Parameters
    ----------
    n_features : int, optional
        Number of selected active features

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
    >>> clf = glm.LARS(n_features=1)
    >>> clf.fit([[-1,1], [0, 0], [1, 1]], [-1, 0, -1])
    LARS(normalize=True, n_features=1, verbose=False, fit_intercept=True)
    >>> print clf.coef_
    [ 0.         -0.81649658]

    Notes
    -----
    See also scikits.learn.glm.LassoLARS that fits a LASSO model
    using a variant of Least Angle Regression

    http://en.wikipedia.org/wiki/Least_angle_regression

    See examples/glm/plot_lar.py for an example.
    """
    def __init__(self, n_features, normalize=True, fit_intercept=True, verbose=False):
        self.n_features = n_features
        self.normalize = normalize
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.verbose = verbose
        if not normalize:
            warnings.warn('LARS can only be fit with normalized regressors. Please set normalize to True')

    def fit (self, X, y, Gram=None, overwrite_X=False, **params):
        self._set_params(**params)

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

        if self.normalize:
            norms = np.sqrt(np.sum(X**2, axis=0))
            nonzeros = np.flatnonzero(norms)
            X[:, nonzeros] /= norms[nonzeros]

        alphas_, active, coef_path_ = lars_path(X, y, method='lar',
                                Gram=Gram, overwrite_X=overwrite_X,
                                max_features=self.n_features)

        self.coef_ = coef_path_[:,-1]

        self._set_intercept(Xmean, ymean)

        return self


class LassoLARS (LinearModel):
    """ Lasso model fit with Least Angle Regression a.k.a. LARS

    It is a Linear Model trained with an L1 prior as regularizer.
    lasso).

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
    >>> clf = glm.LassoLARS(alpha=0.1)
    >>> clf.fit([[-1,1], [0, 0], [1, 1]], [-1, 0, -1])
    LassoLARS(max_features=None, alpha=0.1, normalize=True, verbose=False,
         fit_intercept=True)
    >>> print clf.coef_
    [ 0.         -0.51649658]

    Notes
    -----
    See examples/glm/plot_lasso_lars.py for an example.

    See also scikits.learn.glm.Lasso that fits the same model using
    an alternative optimization strategy called 'coordinate descent.'
    """

    def __init__(self, alpha=1.0, max_features=None, normalize=True,
                        fit_intercept=True, verbose=False):
        """ XXX : add doc
                # will only normalize non-zero columns
        """
        self.alpha = alpha
        self.normalize = normalize
        self.coef_ = None
        self.max_features = max_features
        self.fit_intercept = fit_intercept
        self.verbose = verbose

    def fit (self, X, y, Gram=None, overwrite_X=False,
             overwrite_Gram=False, **params):
        """ XXX : add doc
        """
        self._set_params(**params)

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

        n_samples = X.shape[0]
        alpha = self.alpha * n_samples # scale alpha with number of samples

        # XXX : should handle also unnormalized datasets
        if self.normalize:
            norms = np.sqrt(np.sum(X**2, axis=0))
            nonzeros = np.flatnonzero(norms)
            X[:, nonzeros] /= norms[nonzeros]

        alphas_, active, coef_path_ = lars_path(X, y, Gram=Gram,
                  overwrite_X=overwrite_X,
                  overwrite_Gram=overwrite_Gram, alpha_min=alpha,
                  method='lasso', verbose=self.verbose,
                  max_features=self.max_features)

        self.coef_ = coef_path_[:,-1]

        self._set_intercept(Xmean, ymean)

        return self

