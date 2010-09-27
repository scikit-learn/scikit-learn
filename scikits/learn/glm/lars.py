"""
Least Angle Regression algorithm. See the documentation on the
Generalized Linear Model for a complete discussion.
"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#
# License: BSD Style.

from math import fabs, sqrt
import numpy as np
from scipy import linalg

from .base import LinearModel
from ..utils import arrayfuncs

def lars_path(X, y, Gram=None, max_features=None, alpha_min=0,
              method="lar", verbose=False):
    """ Compute Least Angle Regression and LASSO path

        Parameters
        -----------
        X: array, shape: (n, p)
            Input data

        y: array, shape: (n)
            Input targets

        max_features: integer, optional
            The number of selected features

        Gram: array, shape: (p, p), optional
            Precomputed Gram matrix (X' * X)

        alpha_min: float, optional
            The minimum correlation along the path. It corresponds to
            the regularization parameter alpha parameter in the Lasso.

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
        http://en.wikipedia.org/wiki/Least-angle_regression
        http://en.wikipedia.org/wiki/Lasso_(statistics)#LASSO_method
        XXX : add reference papers
        
        XXX : make sure it works with non-normalized columns of X

    """
    # TODO: detect stationary points.

    X = np.atleast_2d(X)
    y = np.atleast_1d(y)

    n_samples, n_features = X.shape

    if max_features is None:
        max_features = min(n_samples, n_features)

    max_iter = max_features # OK for now but can be expanded dynomically
                            # for the lasso case

    # because of some restrictions in Cython, boolean values are
    # simulated using np.int8
    coefs = np.zeros((max_iter + 1, n_features))
    alphas = np.zeros(max_iter + 1)
    n_iter, n_active = 0, 0
    active = list()
    n_inactive = n_features
    active_mask = np.zeros(n_features, dtype=np.uint8)
    # holds the sign of covariance
    sign_active = np.empty(max_features, dtype=np.int8)
    Cov = np.empty(n_features)
    a = np.empty(n_features)
    drop = False

    # will hold the cholesky factorization
    # only lower part is referenced. We do not create it as
    # empty array because chol_solve calls chkfinite on the
    # whole array, which can cause problems.
    L = np.zeros((max_features, max_features), dtype=np.float64)

    if Gram is None:
        # setting the array to be fortran-ordered speeds up the
        # calculation of the (partial) Gram matrix, but not very
        # useful if the Gram matrix is already precomputed.
        X = np.asfortranarray(X)
    Xt  = X.T

    if Gram is not None:
        Xty = np.dot(Xt, y)
    else:
        res = y.copy() # Residual to be kept up to date

    if verbose:
        print "Step\t\tAdded\t\tDropped\t\tActive set size\t\tC"

    while 1:

        n_inactive = n_features - n_active # number of inactive elements
        inactive_mask = np.logical_not(active_mask)
        inactive = np.where(inactive_mask)[0]

        # Calculate covariance matrix and get maximum
        if n_inactive:
            if Gram is None:
                # Compute X[:,inactive].T * res where res = y - X beta
                # To get the most correlated variable not already in the active set
                arrayfuncs.dot_over(Xt, res, active_mask, np.False_, Cov)
            else:
                arrayfuncs.dot_over(Gram, coefs[n_iter], active_mask, np.False_, a)
                Cov = Xty[inactive_mask] - a[:n_inactive]

            imax = np.argmax(np.abs(Cov[:n_inactive])) # rename
            C_ = Cov[imax]
        else: # special case when all elements are in the active set
            if Gram is None:
                C_ = np.dot(Xt[0], res)
            else:
                C_ = np.dot(Gram[0], coefs[n_iter]) - Xty[0]

        C = fabs(C_)
        alphas[n_iter] = C

        if n_active >= max_features:
            break

        if (C < alpha_min): break

        if not drop:
            imax = inactive[imax] # needs to be sorted for this to work

            # Update the Cholesky factorization of (Xa * Xa') #
            #                                                 #
            #            ( L   0 )                            #
            #     L  ->  (       )  , where L * w = b         #
            #            ( w   z )    z = 1 - ||w||           #
            #                                                 #
            #   where u is the last added to the active set   #

            sign_active[n_active] = np.sign(C_)

            if Gram is None:
                X_max = Xt[imax]
                c = linalg.norm(X_max)**2
                b = np.dot(X_max, X[:, active])
            else:
                c = Gram[imax, imax]
                b = Gram[imax, active]

            # Do cholesky update of the Gram matrix of the active set
            L[n_active, n_active] = c
            active.append(imax)
            if n_active > 0:
                arrayfuncs.solve_triangular(L[:n_active, :n_active], b)
                L[n_active, :n_active] = b[:]
                v = np.dot(L[n_active, :n_active], L[n_active, :n_active])
                L[n_active,  n_active] = np.sqrt(c - v)

            n_active += 1

            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, imax+1, '',
                                                            n_active, C)

        # Now we go into the normal equations dance.
        # (Golub & Van Loan, 1996)

        # compute eqiangular vector
        b = linalg.cho_solve((L[:n_active, :n_active], True),
                             sign_active[:n_active])
        AA = 1. / sqrt(np.sum(b * sign_active[:n_active]))
        b *= AA

        if Gram is None:
            eqdir = np.dot(X[:,active], b) # equiangular direction (unit vector)
            # correlation between active variables and eqiangular vector
            arrayfuncs.dot_over(Xt, eqdir, active_mask, np.False_, a)
        else:
            # in the case of huge number of features, this takes 50% of time
            # to avoid the copies, we could reorder the rows of Gram ...
            arrayfuncs.dot_over (Gram[active].T, b, active_mask, np.False_, a)

        if n_active >= n_features:
            gamma_ = C / AA
        else:
            # equation 2.13
            g1 = (C - Cov[:n_inactive]) / (AA - a[:n_inactive])
            g2 = (C + Cov[:n_inactive]) / (AA + a[:n_inactive])
            gamma_ = np.r_[g1[g1 > 0], g2[g2 > 0], C / AA].min()

        if not drop:
            # Quickfix
            active_mask[imax] = np.True_

        if method == 'lasso':
            drop = False
            z = - coefs[n_iter, active] / b
            z_pos = z[z > 0.]
            if z_pos.size > 0:
                gamma_tilde_ = np.r_[z_pos, gamma_].min()
                if gamma_tilde_ < gamma_:
                    idx = np.where(z == gamma_tilde_)[0]
                    gamma_ = gamma_tilde_
                    drop = True

        n_iter += 1

        if n_iter > max_iter: # resize
            coefs = np.r_[coefs, np.zeros((max_iter + 1, n_features))]
            alphas = np.r_[alphas, np.zeros(max_iter + 1)]
            max_iter += max_iter

        coefs[n_iter, active] = coefs[n_iter - 1, active] + gamma_ * b

        if Gram is None:
            res -= gamma_ * eqdir # update residual

        if n_active > n_features:
            break

        if drop:
            arrayfuncs.cholesky_delete(L[:n_active, :n_active], idx)
            n_active -= 1
            drop_idx = active.pop(idx)
            active_mask[drop_idx] = False
            # do an append to maintain size
            sign_active = np.delete(sign_active, idx)
            sign_active = np.append(sign_active, 0.)
            if verbose:
                print "%s\t\t%s\t\t%s\t\t%s\t\t%s" % (n_iter, '', drop_idx+1,
                                                            n_active, C)

    if C < alpha_min: # interpolate
        # interpolation factor 0 <= ss < 1
        ss = (alphas[n_iter-1] - alpha_min) / (alphas[n_iter-1] - alphas[n_iter])
        coefs[n_iter] = coefs[n_iter-1] + ss*(coefs[n_iter] - coefs[n_iter-1]);
        alphas[n_iter] = alpha_min

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
    LARS(normalize=True, n_features=1)
    >>> print clf.coef_
    [ 0.         -0.81649658]

    Notes
    -----
    See also scikits.learn.glm.LassoLARS that fits a LASSO model
    using a variant of Least Angle Regression

    http://en.wikipedia.org/wiki/Least_angle_regression

    See examples. XXX : add examples names
    """
    def __init__(self, n_features, normalize=True):
        self.n_features = n_features
        self.normalize = normalize
        self.coef_ = None
        self.fit_intercept = True

    def fit (self, X, y, Gram=None, **params):
        self._set_params(**params)
                # will only normalize non-zero columns

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, Ymean = self._center_data(X, y)

        if self.normalize:
            norms = np.sqrt(np.sum(X**2, axis=0))
            nonzeros = np.flatnonzero(norms)
            X[:, nonzeros] /= norms[nonzeros]

        method = 'lar'
        alphas_, active, coef_path_ = lars_path(X, y, Gram=Gram,
                                max_features=self.n_features, method=method)
        self.coef_ = coef_path_[:,-1]
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
    LassoLARS(max_features=None, alpha=0.1, normalize=True, fit_intercept=True)
    >>> print clf.coef_
    [ 0.         -0.51649658]

    Notes
    -----
    See also scikits.learn.glm.Lasso that fits the same model using
    an alternative optimization strategy called 'coordinate descent.'
    """

    def __init__(self, alpha=1.0, max_features=None, normalize=True,
                        fit_intercept=True):
        """ XXX : add doc
                # will only normalize non-zero columns
        """
        self.alpha = alpha
        self.normalize = normalize
        self.coef_ = None
        self.max_features = max_features
        self.fit_intercept = fit_intercept

    def fit (self, X, y, Gram=None, **params):
        """ XXX : add doc
        """
        self._set_params(**params)

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, Ymean = self._center_data(X, y)

        n_samples = X.shape[0]
        alpha = self.alpha * n_samples # scale alpha with number of samples

        # XXX : should handle also unnormalized datasets
        if self.normalize:
            norms = np.sqrt(np.sum(X**2, axis=0))
            nonzeros = np.flatnonzero(norms)
            X[:, nonzeros] /= norms[nonzeros]

        method = 'lasso'
        alphas_, active, coef_path_ = lars_path(X, y, Gram=Gram,
                                            alpha_min=alpha, method=method,
                                            max_features=self.max_features)

        self.coef_ = coef_path_[:,-1]

        self._set_intercept(Xmean, Ymean)

        return self

