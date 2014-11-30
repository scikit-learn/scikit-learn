"""
The :mod:`sklearn.pls` module implements Partial Least Squares (PLS).
"""

# Author: Edouard Duchesnay <edouard.duchesnay@cea.fr>
# License: BSD 3 clause

from ..base import BaseEstimator, RegressorMixin, TransformerMixin
from ..utils import check_array, check_consistent_length
from ..externals import six

import warnings
from abc import ABCMeta, abstractmethod
import numpy as np
from scipy import linalg
from ..utils import arpack

__all__ = ['PLSCanonical', 'PLSRegression', 'PLSSVD']


def _nipals_twoblocks_inner_loop(X, Y, mode="A", max_iter=500, tol=1e-06,
                                 norm_y_weights=False):
    """Inner loop of the iterative NIPALS algorithm.

    Provides an alternative to the svd(X'Y); returns the first left and right
    singular vectors of X'Y.  See PLS for the meaning of the parameters.  It is
    similar to the Power method for determining the eigenvectors and
    eigenvalues of a X'Y.
    """
    y_score = Y[:, [0]]
    x_weights_old = 0
    ite = 1
    X_pinv = Y_pinv = None
    # Inner loop of the Wold algo.
    while True:
        # 1.1 Update u: the X weights
        if mode == "B":
            if X_pinv is None:
                X_pinv = linalg.pinv(X)   # compute once pinv(X)
            x_weights = np.dot(X_pinv, y_score)
        else:  # mode A
        # Mode A regress each X column on y_score
            x_weights = np.dot(X.T, y_score) / np.dot(y_score.T, y_score)
        # 1.2 Normalize u
        x_weights /= np.sqrt(np.dot(x_weights.T, x_weights))
        # 1.3 Update x_score: the X latent scores
        x_score = np.dot(X, x_weights)
        # 2.1 Update y_weights
        if mode == "B":
            if Y_pinv is None:
                Y_pinv = linalg.pinv(Y)    # compute once pinv(Y)
            y_weights = np.dot(Y_pinv, x_score)
        else:
            # Mode A regress each Y column on x_score
            y_weights = np.dot(Y.T, x_score) / np.dot(x_score.T, x_score)
        ## 2.2 Normalize y_weights
        if norm_y_weights:
            y_weights /= np.sqrt(np.dot(y_weights.T, y_weights))
        # 2.3 Update y_score: the Y latent scores
        y_score = np.dot(Y, y_weights) / np.dot(y_weights.T, y_weights)
        ## y_score = np.dot(Y, y_weights) / np.dot(y_score.T, y_score) ## BUG
        x_weights_diff = x_weights - x_weights_old
        if np.dot(x_weights_diff.T, x_weights_diff) < tol or Y.shape[1] == 1:
            break
        if ite == max_iter:
            warnings.warn('Maximum number of iterations reached')
            break
        x_weights_old = x_weights
        ite += 1
    return x_weights, y_weights, ite


def _svd_cross_product(X, Y):
    C = np.dot(X.T, Y)
    U, s, Vh = linalg.svd(C, full_matrices=False)
    u = U[:, [0]]
    v = Vh.T[:, [0]]
    return u, v


def _center_scale_xy(X, Y, scale=True):
    """ Center X, Y and scale if the scale parameter==True

    Returns
    -------
        X, Y, x_mean, y_mean, x_std, y_std
    """
    # center
    x_mean = X.mean(axis=0)
    X -= x_mean
    y_mean = Y.mean(axis=0)
    Y -= y_mean
    # scale
    if scale:
        x_std = X.std(axis=0, ddof=1)
        x_std[x_std == 0.0] = 1.0
        X /= x_std
        y_std = Y.std(axis=0, ddof=1)
        y_std[y_std == 0.0] = 1.0
        Y /= y_std
    else:
        x_std = np.ones(X.shape[1])
        y_std = np.ones(Y.shape[1])
    return X, Y, x_mean, y_mean, x_std, y_std


class _PLS(six.with_metaclass(ABCMeta), BaseEstimator, TransformerMixin,
           RegressorMixin):
    """Partial Least Squares (PLS)

    This class implements the generic PLS algorithm, constructors' parameters
    allow to obtain a specific implementation such as:

    - PLS2 regression, i.e., PLS 2 blocks, mode A, with asymmetric deflation
      and unnormalized y weights such as defined by [Tenenhaus 1998] p. 132.
      With univariate response it implements PLS1.

    - PLS canonical, i.e., PLS 2 blocks, mode A, with symmetric deflation and
      normalized y weights such as defined by [Tenenhaus 1998] (p. 132) and
      [Wegelin et al. 2000]. This parametrization implements the original Wold
      algorithm.

    We use the terminology defined by [Wegelin et al. 2000].
    This implementation uses the PLS Wold 2 blocks algorithm based on two
    nested loops:
        (i) The outer loop iterate over components.
        (ii) The inner loop estimates the weights vectors. This can be done
        with two algo. (a) the inner loop of the original NIPALS algo. or (b) a
        SVD on residuals cross-covariance matrices.

    n_components : int, number of components to keep. (default 2).

    scale : boolean, scale data? (default True)

    deflation_mode : str, "canonical" or "regression". See notes.

    mode : "A" classical PLS and "B" CCA. See notes.

    norm_y_weights: boolean, normalize Y weights to one? (default False)

    algorithm : string, "nipals" or "svd"
        The algorithm used to estimate the weights. It will be called
        n_components times, i.e. once for each iteration of the outer loop.

    max_iter : an integer, the maximum number of iterations (default 500)
        of the NIPALS inner loop (used only if algorithm="nipals")

    tol : non-negative real, default 1e-06
        The tolerance used in the iterative algorithm.

    copy : boolean
        Whether the deflation should be done on a copy. Let the default
        value to True unless you don't care about side effects.

    Attributes
    ----------
    x_weights_ : array, [p, n_components]
        X block weights vectors.

    y_weights_ : array, [q, n_components]
        Y block weights vectors.

    x_loadings_ : array, [p, n_components]
        X block loadings vectors.

    y_loadings_ : array, [q, n_components]
        Y block loadings vectors.

    x_scores_ : array, [n_samples, n_components]
        X scores.

    y_scores_ : array, [n_samples, n_components]
        Y scores.

    x_rotations_ : array, [p, n_components]
        X block to latents rotations.

    y_rotations_ : array, [q, n_components]
        Y block to latents rotations.

    coef_: array, [p, q]
        The coefficients of the linear model: Y = X coef_ + Err

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component. Not useful if the algorithm given is "svd".

    References
    ----------

    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In French but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    PLSCanonical
    PLSRegression
    CCA
    PLS_SVD
    """

    @abstractmethod
    def __init__(self, n_components=2, scale=True, deflation_mode="regression",
                 mode="A", algorithm="nipals", norm_y_weights=False,
                 max_iter=500, tol=1e-06, copy=True):
        self.n_components = n_components
        self.deflation_mode = deflation_mode
        self.mode = mode
        self.norm_y_weights = norm_y_weights
        self.scale = scale
        self.algorithm = algorithm
        self.max_iter = max_iter
        self.tol = tol
        self.copy = copy

    def fit(self, X, Y):
        """Fit model to data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples in the number of samples and
            n_features is the number of predictors.

        Y : array-like of response, shape = [n_samples, n_targets]
            Target vectors, where n_samples in the number of samples and
            n_targets is the number of response variables.
        """

        # copy since this will contains the residuals (deflated) matrices
        check_consistent_length(X, Y)
        X = check_array(X, dtype=np.float, copy=self.copy)
        Y = check_array(Y, dtype=np.float, copy=self.copy)

        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))
        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')
        if self.algorithm not in ("svd", "nipals"):
            raise ValueError("Got algorithm %s when only 'svd' "
                             "and 'nipals' are known" % self.algorithm)
        if self.algorithm == "svd" and self.mode == "B":
            raise ValueError('Incompatible configuration: mode B is not '
                             'implemented with svd algorithm')
        if not self.deflation_mode in ["canonical", "regression"]:
            raise ValueError('The deflation mode is unknown')
        # Scale (in place)
        X, Y, self.x_mean_, self.y_mean_, self.x_std_, self.y_std_\
            = _center_scale_xy(X, Y, self.scale)
        # Residuals (deflated) matrices
        Xk = X
        Yk = Y
        # Results matrices
        self.x_scores_ = np.zeros((n, self.n_components))
        self.y_scores_ = np.zeros((n, self.n_components))
        self.x_weights_ = np.zeros((p, self.n_components))
        self.y_weights_ = np.zeros((q, self.n_components))
        self.x_loadings_ = np.zeros((p, self.n_components))
        self.y_loadings_ = np.zeros((q, self.n_components))
        self.n_iter_ = []

        # NIPALS algo: outer loop, over components
        for k in range(self.n_components):
            #1) weights estimation (inner loop)
            # -----------------------------------
            if self.algorithm == "nipals":
                x_weights, y_weights, n_iter_ = \
                    _nipals_twoblocks_inner_loop(
                        X=Xk, Y=Yk, mode=self.mode, max_iter=self.max_iter,
                        tol=self.tol, norm_y_weights=self.norm_y_weights)
                self.n_iter_.append(n_iter_)
            elif self.algorithm == "svd":
                x_weights, y_weights = _svd_cross_product(X=Xk, Y=Yk)
            # compute scores
            x_scores = np.dot(Xk, x_weights)
            if self.norm_y_weights:
                y_ss = 1
            else:
                y_ss = np.dot(y_weights.T, y_weights)
            y_scores = np.dot(Yk, y_weights) / y_ss
            # test for null variance
            if np.dot(x_scores.T, x_scores) < np.finfo(np.double).eps:
                warnings.warn('X scores are null at iteration %s' % k)
            #2) Deflation (in place)
            # ----------------------
            # Possible memory footprint reduction may done here: in order to
            # avoid the allocation of a data chunk for the rank-one
            # approximations matrix which is then subtracted to Xk, we suggest
            # to perform a column-wise deflation.
            #
            # - regress Xk's on x_score
            x_loadings = np.dot(Xk.T, x_scores) / np.dot(x_scores.T, x_scores)
            # - subtract rank-one approximations to obtain remainder matrix
            Xk -= np.dot(x_scores, x_loadings.T)
            if self.deflation_mode == "canonical":
                # - regress Yk's on y_score, then subtract rank-one approx.
                y_loadings = (np.dot(Yk.T, y_scores)
                              / np.dot(y_scores.T, y_scores))
                Yk -= np.dot(y_scores, y_loadings.T)
            if self.deflation_mode == "regression":
                # - regress Yk's on x_score, then subtract rank-one approx.
                y_loadings = (np.dot(Yk.T, x_scores)
                              / np.dot(x_scores.T, x_scores))
                Yk -= np.dot(x_scores, y_loadings.T)
            # 3) Store weights, scores and loadings # Notation:
            self.x_scores_[:, k] = x_scores.ravel()  # T
            self.y_scores_[:, k] = y_scores.ravel()  # U
            self.x_weights_[:, k] = x_weights.ravel()  # W
            self.y_weights_[:, k] = y_weights.ravel()  # C
            self.x_loadings_[:, k] = x_loadings.ravel()  # P
            self.y_loadings_[:, k] = y_loadings.ravel()  # Q
        # Such that: X = TP' + Err and Y = UQ' + Err

        # 4) rotations from input space to transformed space (scores)
        # T = X W(P'W)^-1 = XW* (W* : p x k matrix)
        # U = Y C(Q'C)^-1 = YC* (W* : q x k matrix)
        self.x_rotations_ = np.dot(
            self.x_weights_,
            linalg.pinv(np.dot(self.x_loadings_.T, self.x_weights_)))
        if Y.shape[1] > 1:
            self.y_rotations_ = np.dot(
                self.y_weights_,
                linalg.pinv(np.dot(self.y_loadings_.T, self.y_weights_)))
        else:
            self.y_rotations_ = np.ones(1)

        if True or self.deflation_mode == "regression":
            # Estimate regression coefficient
            # Regress Y on T
            # Y = TQ' + Err,
            # Then express in function of X
            # Y = X W(P'W)^-1Q' + Err = XB + Err
            # => B = W*Q' (p x q)
            self.coef_ = np.dot(self.x_rotations_, self.y_loadings_.T)
            self.coef_ = (1. / self.x_std_.reshape((p, 1)) * self.coef_ *
                          self.y_std_)
        return self

    def transform(self, X, Y=None, copy=True):
        """Apply the dimension reduction learned on the train data.

        Parameters
        ----------
        X : array-like of predictors, shape = [n_samples, p]
            Training vectors, where n_samples in the number of samples and
            p is the number of predictors.

        Y : array-like of response, shape = [n_samples, q], optional
            Training vectors, where n_samples in the number of samples and
            q is the number of response variables.

        copy : boolean
            Whether to copy X and Y, or perform in-place normalization.

        Returns
        -------
        x_scores if Y is not given, (x_scores, y_scores) otherwise.
        """
        # Normalize
        if copy:
            Xc = (np.asarray(X) - self.x_mean_) / self.x_std_
            if Y is not None:
                Yc = (np.asarray(Y) - self.y_mean_) / self.y_std_
        else:
            X = np.asarray(X)
            Xc -= self.x_mean_
            Xc /= self.x_std_
            if Y is not None:
                Y = np.asarray(Y)
                Yc -= self.y_mean_
                Yc /= self.y_std_
        # Apply rotation
        x_scores = np.dot(Xc, self.x_rotations_)
        if Y is not None:
            y_scores = np.dot(Yc, self.y_rotations_)
            return x_scores, y_scores

        return x_scores

    def predict(self, X, copy=True):
        """Apply the dimension reduction learned on the train data.

        Parameters
        ----------
        X : array-like of predictors, shape = [n_samples, p]
            Training vectors, where n_samples in the number of samples and
            p is the number of predictors.

        copy : boolean
            Whether to copy X and Y, or perform in-place normalization.

        Notes
        -----
        This call requires the estimation of a p x q matrix, which may
        be an issue in high dimensional space.
        """
        # Normalize
        if copy:
            Xc = (np.asarray(X) - self.x_mean_)
        else:
            X = np.asarray(X)
            Xc -= self.x_mean_
            Xc /= self.x_std_
        Ypred = np.dot(Xc, self.coef_)
        return Ypred + self.y_mean_

    def fit_transform(self, X, y=None, **fit_params):
        """Learn and apply the dimension reduction on the train data.

        Parameters
        ----------
        X : array-like of predictors, shape = [n_samples, p]
            Training vectors, where n_samples in the number of samples and
            p is the number of predictors.

        Y : array-like of response, shape = [n_samples, q], optional
            Training vectors, where n_samples in the number of samples and
            q is the number of response variables.

        copy : boolean
            Whether to copy X and Y, or perform in-place normalization.

        Returns
        -------
        x_scores if Y is not given, (x_scores, y_scores) otherwise.
        """
        return self.fit(X, y, **fit_params).transform(X, y)


class PLSRegression(_PLS):
    """PLS regression

    PLSRegression implements the PLS 2 blocks regression known as PLS2 or PLS1
    in case of one dimensional response.
    This class inherits from _PLS with mode="A", deflation_mode="regression",
    norm_y_weights=False and algorithm="nipals".

    Parameters
    ----------
    X : array-like of predictors, shape = [n_samples, p]
        Training vectors, where n_samples in the number of samples and
        p is the number of predictors.

    Y : array-like of response, shape = [n_samples, q]
        Training vectors, where n_samples in the number of samples and
        q is the number of response variables.

    n_components : int, (default 2)
        Number of components to keep.

    scale : boolean, (default True)
        whether to scale the data

    max_iter : an integer, (default 500)
        the maximum number of iterations of the NIPALS inner loop (used
        only if algorithm="nipals")

    tol : non-negative real
        Tolerance used in the iterative algorithm default 1e-06.

    copy : boolean, default True
        Whether the deflation should be done on a copy. Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_ : array, [p, n_components]
        X block weights vectors.

    y_weights_ : array, [q, n_components]
        Y block weights vectors.

    x_loadings_ : array, [p, n_components]
        X block loadings vectors.

    y_loadings_ : array, [q, n_components]
        Y block loadings vectors.

    x_scores_ : array, [n_samples, n_components]
        X scores.

    y_scores_ : array, [n_samples, n_components]
        Y scores.

    x_rotations_ : array, [p, n_components]
        X block to latents rotations.

    y_rotations_ : array, [q, n_components]
        Y block to latents rotations.

    coef_: array, [p, q]
        The coefficients of the linear model: Y = X coef_ + Err

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component.

    Notes
    -----
    For each component k, find weights u, v that optimizes:
    ``max corr(Xk u, Yk v) * var(Xk u) var(Yk u)``, such that ``|u| = 1``

    Note that it maximizes both the correlations between the scores and the
    intra-block variances.

    The residual matrix of X (Xk+1) block is obtained by the deflation on
    the current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current X score. This performs the PLS regression known as PLS2. This
    mode is prediction oriented.

    This implementation provides the same results that 3 PLS packages
    provided in the R language (R-project):

        - "mixOmics" with function pls(X, Y, mode = "regression")
        - "plspm " with function plsreg2(X, Y)
        - "pls" with function oscorespls.fit(X, Y)

    Examples
    --------
    >>> from sklearn.cross_decomposition import PLSRegression
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> pls2 = PLSRegression(n_components=2)
    >>> pls2.fit(X, Y)
    ... # doctest: +NORMALIZE_WHITESPACE
    PLSRegression(copy=True, max_iter=500, n_components=2, scale=True,
            tol=1e-06)
    >>> Y_pred = pls2.predict(X)

    References
    ----------

    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In french but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.
    """

    def __init__(self, n_components=2, scale=True,
                 max_iter=500, tol=1e-06, copy=True):
        _PLS.__init__(self, n_components=n_components, scale=scale,
                      deflation_mode="regression", mode="A",
                      norm_y_weights=False, max_iter=max_iter, tol=tol,
                      copy=copy)


class PLSCanonical(_PLS):
    """ PLSCanonical implements the 2 blocks canonical PLS of the original Wold
    algorithm [Tenenhaus 1998] p.204, referred as PLS-C2A in [Wegelin 2000].

    This class inherits from PLS with mode="A" and deflation_mode="canonical",
    norm_y_weights=True and algorithm="nipals", but svd should provide similar
    results up to numerical errors.

    Parameters
    ----------
    X : array-like of predictors, shape = [n_samples, p]
        Training vectors, where n_samples is the number of samples and
        p is the number of predictors.

    Y : array-like of response, shape = [n_samples, q]
        Training vectors, where n_samples is the number of samples and
        q is the number of response variables.

    n_components : int, number of components to keep. (default 2).

    scale : boolean, scale data? (default True)

    algorithm : string, "nipals" or "svd"
        The algorithm used to estimate the weights. It will be called
        n_components times, i.e. once for each iteration of the outer loop.

    max_iter : an integer, (default 500)
        the maximum number of iterations of the NIPALS inner loop (used
        only if algorithm="nipals")

    tol : non-negative real, default 1e-06
        the tolerance used in the iterative algorithm

    copy : boolean, default True
        Whether the deflation should be done on a copy. Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_ : array, shape = [p, n_components]
        X block weights vectors.

    y_weights_ : array, shape = [q, n_components]
        Y block weights vectors.

    x_loadings_ : array, shape = [p, n_components]
        X block loadings vectors.

    y_loadings_ : array, shape = [q, n_components]
        Y block loadings vectors.

    x_scores_ : array, shape = [n_samples, n_components]
        X scores.

    y_scores_ : array, shape = [n_samples, n_components]
        Y scores.

    x_rotations_ : array, shape = [p, n_components]
        X block to latents rotations.

    y_rotations_ : array, shape = [q, n_components]
        Y block to latents rotations.

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component. Not useful if the algorithm provided is "svd".

    Notes
    -----
    For each component k, find weights u, v that optimize::
    max corr(Xk u, Yk v) * var(Xk u) var(Yk u), such that ``|u| = |v| = 1``

    Note that it maximizes both the correlations between the scores and the
    intra-block variances.

    The residual matrix of X (Xk+1) block is obtained by the deflation on the
    current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current Y score. This performs a canonical symmetric version of the PLS
    regression. But slightly different than the CCA. This is mostly used
    for modeling.

    This implementation provides the same results that the "plspm" package
    provided in the R language (R-project), using the function plsca(X, Y).
    Results are equal or collinear with the function
    ``pls(..., mode = "canonical")`` of the "mixOmics" package. The difference
    relies in the fact that mixOmics implementation does not exactly implement
    the Wold algorithm since it does not normalize y_weights to one.

    Examples
    --------
    >>> from sklearn.cross_decomposition import PLSCanonical
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> plsca = PLSCanonical(n_components=2)
    >>> plsca.fit(X, Y)
    ... # doctest: +NORMALIZE_WHITESPACE
    PLSCanonical(algorithm='nipals', copy=True, max_iter=500, n_components=2,
                 scale=True, tol=1e-06)
    >>> X_c, Y_c = plsca.transform(X, Y)

    References
    ----------

    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    CCA
    PLSSVD
    """

    def __init__(self, n_components=2, scale=True, algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        _PLS.__init__(self, n_components=n_components, scale=scale,
                      deflation_mode="canonical", mode="A",
                      norm_y_weights=True, algorithm=algorithm,
                      max_iter=max_iter, tol=tol, copy=copy)

    @property
    def coefs(self):
        DeprecationWarning("'coefs' attribute has been deprecated. Use 'coef_' "
                           "instead")
        return self._coef


class PLSSVD(BaseEstimator, TransformerMixin):
    """Partial Least Square SVD

    Simply perform a svd on the crosscovariance matrix: X'Y
    There are no iterative deflation here.

    Parameters
    ----------
    X : array-like of predictors, shape = [n_samples, p]
        Training vector, where n_samples is the number of samples and
        p is the number of predictors. X will be centered before any analysis.

    Y : array-like of response, shape = [n_samples, q]
        Training vector, where n_samples is the number of samples and
        q is the number of response variables. X will be centered before any
        analysis.

    n_components : int, (default 2).
        number of components to keep.

    scale : boolean, (default True)
        whether to scale X and Y.

    Attributes
    ----------
    x_weights_ : array, [p, n_components]
        X block weights vectors.

    y_weights_ : array, [q, n_components]
        Y block weights vectors.

    x_scores_ : array, [n_samples, n_components]
        X scores.

    y_scores_ : array, [n_samples, n_components]
        Y scores.

    See also
    --------
    PLSCanonical
    CCA
    """

    def __init__(self, n_components=2, scale=True, copy=True):
        self.n_components = n_components
        self.scale = scale
        self.copy = copy

    def fit(self, X, Y):
        # copy since this will contains the centered data
        check_consistent_length(X, Y)
        X = check_array(X, dtype=np.float, copy=self.copy)
        Y = check_array(Y, dtype=np.float, copy=self.copy)

        p = X.shape[1]

        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')

        # Scale (in place)
        X, Y, self.x_mean_, self.y_mean_, self.x_std_, self.y_std_ =\
            _center_scale_xy(X, Y, self.scale)
        # svd(X'Y)
        C = np.dot(X.T, Y)

        # The arpack svds solver only works if the number of extracted
        # components is smaller than rank(X) - 1. Hence, if we want to extract
        # all the components (C.shape[1]), we have to use another one. Else,
        # let's use arpacks to compute only the interesting components.
        if self.n_components == C.shape[1]:
            U, s, V = linalg.svd(C, full_matrices=False)
        else:
            U, s, V = arpack.svds(C, k=self.n_components)
        V = V.T
        self.x_scores_ = np.dot(X, U)
        self.y_scores_ = np.dot(Y, V)
        self.x_weights_ = U
        self.y_weights_ = V
        return self

    def transform(self, X, Y=None):
        """Apply the dimension reduction learned on the train data."""
        Xr = (X - self.x_mean_) / self.x_std_
        x_scores = np.dot(Xr, self.x_weights_)
        if Y is not None:
            Yr = (Y - self.y_mean_) / self.y_std_
            y_scores = np.dot(Yr, self.y_weights_)
            return x_scores, y_scores
        return x_scores

    def fit_transform(self, X, y=None, **fit_params):
        """Learn and apply the dimension reduction on the train data.

        Parameters
        ----------
        X : array-like of predictors, shape = [n_samples, p]
            Training vectors, where n_samples in the number of samples and
            p is the number of predictors.

        Y : array-like of response, shape = [n_samples, q], optional
            Training vectors, where n_samples in the number of samples and
            q is the number of response variables.

        Returns
        -------
        x_scores if Y is not given, (x_scores, y_scores) otherwise.
        """
        return self.fit(X, y, **fit_params).transform(X, y)
