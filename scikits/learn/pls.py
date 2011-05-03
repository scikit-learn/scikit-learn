""" Partial Least Square
"""

# Author: Edouard Duchesnay <edouard.duchesnay@cea.fr>
# License: BSD Style.

from .base import BaseEstimator
#from scikits.learn.base import BaseEstimator
import warnings
import numpy as np
from scipy import linalg


def _nipals_twoblocks_inner_loop(X, Y, mode="A", max_iter=500, tol=1e-06):
    """Inner loop of the iterative NIPALS algorithm. provide an alternative
    of the svd(X'Y) ie. return the first left and rigth singular vectors of X'Y
    See PLS for the meaning of the parameters.
    It is similar to the Power method for determining the eigenvectors and
    eigenvalues of a X'Y
    """
    y_score = Y[:, [0]]
    u_old = 0
    ite = 1
    X_pinv = Y_pinv = None
    # Inner loop of the Wold algo.
    while True:
        # 1.1 Update u: the X weights
        if mode is "B":
            if X_pinv is None:
                X_pinv = linalg.pinv(X)   # compute once pinv(X)
            u = np.dot(X_pinv, y_score)
        else:  # mode A
        # Mode A regress each X column on y_score
            u = np.dot(X.T, y_score) / np.dot(y_score.T, y_score)
        # 1.2 Normalize u
        u /= np.sqrt(np.dot(u.T, u))
        # 1.3 Update x_score: the X latent scores
        x_score = np.dot(X, u)

        # 2.1 Update v: the Y weights
        if mode is "B":
            if Y_pinv is None:
                Y_pinv = linalg.pinv(Y)    # compute once pinv(Y)
            v = np.dot(Y_pinv, x_score)
        else:
            # Mode A regress each X column on y_score
            v = np.dot(Y.T, x_score) / np.dot(x_score.T, x_score)
        # 2.2 Normalize v
        v /= np.sqrt(np.dot(v.T, v))
        # 2.3 Update y_score: the Y latent scores
        y_score = np.dot(Y, v)

        u_diff = u - u_old
        if np.dot(u_diff.T, u_diff) < tol or Y.shape[1] == 1:
            break
        if ite == max_iter:
            warnings.warn('Maximum number of iterations reached')
            break
        u_old = u
        ite += 1
    return u, v


def _svd_cross_product(X, Y):
    C = np.dot(X.T, Y)
    U, s, Vh = linalg.svd(C, full_matrices=False)
    u = U[:, [0]]
    v = Vh.T[:, [0]]
    return u, v


def _center_scale_xy(X, Y, scale=True):
    """ Center X, Y and scale if the scale parameter==True
    Return
    ------
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
        X /= x_std
        y_std = Y.std(axis=0, ddof=1)
        Y /= y_std
    else:
        x_std = np.ones(X.shape[1])
        y_std = np.ones(Y.shape[1])
    return X, Y, x_mean, y_mean, x_std, y_std


class _PLS(BaseEstimator):
    """Partial Least Square (PLS)

    We use the therminology defined by [Wegelin et al. 2000].
    This implementation uses the PLS Wold 2 blocks algorithm or NIPALS which is
    based on two nested loops:
    (i) The outer loop iterate over compoments.
        (ii) The inner loop estimates the loading vectors. This can be done
        with two algo. (a) the inner loop of the original NIPALS algo or (b) a
        SVD on residuals cross-covariance matrices.

    This implementation provides:
    - PLS regression, ie., PLS 2 blocks, mode A, with asymetric deflation.
      A.k.a. PLS2, with multivariate response or PLS1 with univariate response.
    - PLS canonical, ie., PLS 2 blocks, mode A, with symetric deflation.
    - CCA, ie.,  PLS 2 blocks, mode B, with symetric deflation.

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vectors, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vectors, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).

    deflation_mode: str, "canonical" or "regression". See notes.

    mode: "A" classical PLS and "B" CCA. See notes.

    scale: boolean, scale data? (default True)

    algorithm: str "nipals" or "svd" the algorithm used to estimate the
        weights, it will be called "n_components" time ie.: for each iteration
        of the outer loop.

    max_iter: an integer, the maximum number of iterations (default 500) of the
        NIPALS inner loop (used only if algorithm="nipals")

    tol: a not negative real, the tolerance used in the iterative algorithm
         default 1e-06.

    copy: boolean, should the deflation been made on a copy? Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_: array, [p, n_components]
        X block weights vectors.

    y_weights_: array, [q, n_components]
        Y block weights vectors.

    x_loadings_: array, [p, n_components]
        X block loadings vectors.

    y_loadings_: array, [q, n_components]
        Y block loadings vectors.

    x_scores_: array, [n_samples, n_components]
        X scores.

    y_scores_: array, [n_samples, n_components]
        Y scores.

    x_rotations_: array, [p, n_components]
        X block to latents rotations.

    y_rotations_: array, [q, n_components]
        Y block to latents rotations.

    coefs: array, [p, q]
        The coeficients of the linear model: Y = X coefs + Err

    References
    ----------
    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In french but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    PLSCanonical
    PLSRegression
    CCA
    PLS_SVD
    """

    def __init__(self, n_components=2, deflation_mode="canonical", mode="A",
                 scale=True,
                 algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        self.n_components = n_components
        self.deflation_mode = deflation_mode
        self.mode = mode
        self.scale = scale
        self.algorithm = algorithm
        self.max_iter = max_iter
        self.tol = tol
        self.copy = copy

    def fit(self, X, Y, **params):
        self._set_params(**params)
        # copy since this will contains the residuals (deflated) matrices
        if self.copy:
            X = np.asanyarray(X, dtype=np.float).copy()
            Y = np.asanyarray(Y, dtype=np.float).copy()
        else:
            X = np.asanyarray(X, dtype=np.float)
            Y = np.asanyarray(Y, dtype=np.float)

        if X.ndim != 2:
            raise ValueError('X must be a 2D array')
        if Y.ndim == 1:
            Y = Y.reshape((Y.size, 1))
        if Y.ndim != 2:
            raise ValueError('Y must be a 1D or a 2D array')

        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))
        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')
        if self.algorithm is "svd" and self.mode is "B":
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

        # NIPALS algo: outer loop, over components
        for k in xrange(self.n_components):
            #1) weights estimation (inner loop)
            # -----------------------------------
            if self.algorithm is "nipals":
                u, v = _nipals_twoblocks_inner_loop(
                        X=Xk, Y=Yk, mode=self.mode,
                        max_iter=self.max_iter, tol=self.tol)
            if self.algorithm is "svd":
                u, v = _svd_cross_product(X=Xk, Y=Yk)
            # compute scores
            x_score = np.dot(Xk, u)
            y_score = np.dot(Yk, v)
            # test for null variance
            if np.dot(x_score.T, x_score) < np.finfo(np.double).eps:
                warnings.warn('X scores are null at iteration %s' % k)
            #2) Deflation (in place)
            # ----------------------
            # Possible memory footprint reduction may done here: in order to
            # avoid the allocation of a data chunk for the rank-one
            # approximations matrix which is then substracted to Xk, we suggest
            # to perform a column-wise deflation.
            #
            # - regress Xk's on x_score
            x_loadings = np.dot(Xk.T, x_score) / np.dot(x_score.T, x_score)
            # - substract rank-one approximations to obtain remainder matrix
            Xk -= np.dot(x_score, x_loadings.T)
            if self.deflation_mode is "canonical":
                # - regress Yk's on y_score, then substract rank-one approx.
                y_loadings = np.dot(Yk.T, y_score) / np.dot(y_score.T, y_score)
                Yk -= np.dot(y_score, y_loadings.T)
            if self.deflation_mode is "regression":
                # - regress Yk's on x_score, then substract rank-one approx.
                y_loadings = np.dot(Yk.T, x_score) / np.dot(x_score.T, x_score)
                Yk -= np.dot(x_score, y_loadings.T)
            # 3) Store weights, scores and loadings # Notation:
            self.x_scores_[:, k] = x_score.ravel()  # T
            self.y_scores_[:, k] = y_score.ravel()  # U
            self.x_weights_[:, k] = u.ravel()  # W
            self.y_weights_[:, k] = v.ravel()  # C
            self.x_loadings_[:, k] = x_loadings.ravel()  # P
            self.y_loadings_[:, k] = y_loadings.ravel()  # Q
        # Such that: X = TP' + Err and Y = UQ' + Err

        # 4) rotations from input space to transformed space (scores)
        # T = X W(P'W)^-1 = XW* (W* : p x k matrix)
        # U = Y C(Q'C)^-1 = YC* (W* : q x k matrix)
        self.x_rotations_ = np.dot(self.x_weights_,
            linalg.inv(np.dot(self.x_loadings_.T, self.x_weights_)))
        if Y.shape[1] > 1:
            self.y_rotations_ = np.dot(self.y_weights_,
                linalg.inv(np.dot(self.y_loadings_.T, self.y_weights_)))
        else:
            self.y_rotations_ = np.ones(1)

        if True or self.deflation_mode is "regression":
            # Estimate regression coeficient
            # Regress Y on T
            # Y = TQ' + Err,
            # Then express in function of X
            # Y = X W(P'W)^-1Q' + Err = XB + Err
            # => B = W*Q' (p x q)
            self.coefs = np.dot(self.x_rotations_, self.y_loadings_.T)
            self.coefs = 1. / self.x_std_.reshape((p, 1)) * \
                    self.coefs * self.y_std_
        return self

    def transform(self, X, Y=None, copy=True):
        """Apply the dimension reduction learned on the train data.
            Parameters
            ----------
            X: array-like of predictors, shape (n_samples, p)
                Training vectors, where n_samples in the number of samples and
                p is the number of predictors.

            Y: array-like of response, shape (n_samples, q), optional
                Training vectors, where n_samples in the number of samples and
                q is the number of response variables.

            copy: X and Y have to be normalize, do it on a copy or in place
                with side effect!

            Returns
            -------
            x_scores if Y is not given, (x_scores, y_scores) otherwise.

        """
        # Normalize
        if copy:
            Xc = (np.asanyarray(X) - self.x_mean_) / self.x_std_
            if Y is not None:
                Yc = (np.asanyarray(Y) - self.y_mean_) / self.y_std_
        else:
            X = np.asanyarray(X)
            Xc -= self.x_mean_
            Xc /= self.x_std_
            if Y is not None:
                Y = np.asanyarray(Y)
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
            X: array-like of predictors, shape (n_samples, p)
                Training vectors, where n_samples in the number of samples and
                p is the number of predictors.

            copy: X has to be normalize, do it on a copy or in place
                with side effect!

            Notes
            -----
            This call require the estimation of a p x q matrix, which may
            be an issue in high dimensional space.
        """
        # Normalize
        if copy:
            Xc = (np.asanyarray(X) - self.x_mean_)
        else:
            X = np.asanyarray(X)
            Xc -= self.x_mean_
            Xc /= self.x_std_
        Ypred = np.dot(Xc, self.coefs)
        return Ypred + self.y_mean_


class PLSRegression(_PLS):
    """PLS regression (Also known PLS2 or PLS in case of one dimensional
    response). PLSregression inherits from PLS with mode="A" and
    deflation_mode="regression".

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vectors, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vectors, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).

    scale: boolean, scale data? (default True)

    algorithm: str "nipals" or "svd" the algorithm used to estimate the
        weights, it will be called "n_components" time ie.: for each iteration
        of the outer loop.

    max_iter: an integer, the maximum number of iterations (default 500) of the
        NIPALS inner loop (used only if algorithm="nipals")

    tol: a not negative real, the tolerance used in the iterative algorithm
         default 1e-06.

    copy: boolean, should the deflation been made on a copy? Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_: array, [p, n_components]
        X block weights vectors.

    y_weights_: array, [q, n_components]
        Y block weights vectors.

    x_loadings_: array, [p, n_components]
        X block loadings vectors.

    y_loadings_: array, [q, n_components]
        Y block loadings vectors.

    x_scores_: array, [n_samples, n_components]
        X scores.

    y_scores_: array, [n_samples, n_components]
        Y scores.

    x_rotations_: array, [p, n_components]
        X block to latents rotations.

    y_rotations_: array, [q, n_components]
        Y block to latents rotations.

    coefs: array, [p, q]
        The coeficients of the linear model: Y = X coefs + Err

    Notes
    -----
    For each component k, find weights u, v that optimizes:
    max corr(Xk u, Yk v) * var(Xk u) var(Yk u), such that |u| = |v| = 1

    Note that it maximizes both the correlations between the scores and the
    intra-block variances.

    The residual matrix of X (Xk+1) block is obtained by the deflation on the
    current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current X score. This performs the PLS regression known as PLS2. This
    mode is prediction oriented.

    Examples
    --------
    >>> from scikits.learn.pls import PLSCanonical, PLSRegression, CCA
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> pls2 = PLSRegression()
    >>> pls2.fit(X, Y, n_components=2)
    PLSRegression(scale=True, algorithm='nipals', max_iter=500, n_components=2,
           tol=1e-06, copy=True)
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

    def __init__(self, n_components=2, scale=True, algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        _PLS.__init__(self, n_components=n_components,
                        deflation_mode="regression", mode="A",
                        scale=scale, algorithm=algorithm,
                        max_iter=max_iter, tol=tol, copy=copy)


class PLSCanonical(_PLS):
    """PLS canonical. PLSCanonical inherits from PLS with mode="A" and
    deflation_mode="canonical".

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vectors, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vectors, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).

    scale: boolean, scale data? (default True)

    algorithm: str "nipals" or "svd" the algorithm used to estimate the
        weights, it will be called "n_components" time ie.: for each iteration
        of the outer loop.

    max_iter: an integer, the maximum number of iterations (default 500) of the
        NIPALS inner loop (used only if algorithm="nipals")

    tol: a not negative real, the tolerance used in the iterative algorithm
         default 1e-06.

    copy: boolean, should the deflation been made on a copy? Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_: array, [p, n_components]
        X block weights vectors.

    y_weights_: array, [q, n_components]
        Y block weights vectors.

    x_loadings_: array, [p, n_components]
        X block loadings vectors.

    y_loadings_: array, [q, n_components]
        Y block loadings vectors.

    x_scores_: array, [n_samples, n_components]
        X scores.

    y_scores_: array, [n_samples, n_components]
        Y scores.

    x_rotations_: array, [p, n_components]
        X block to latents rotations.

    y_rotations_: array, [q, n_components]
        Y block to latents rotations.

    Notes
    -----
    For each component k, find weights u, v that optimizes:
    max corr(Xk u, Yk v) * var(Xk u) var(Yk u), such that |u| = |v| = 1

    Note that it maximizes both the correlations between the scores and the
    intra-block variances.

    The residual matrix of X (Xk+1) block is obtained by the deflation on the
    current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current Y score. This performs a canonical symetric version of the PLS
    regression. But slightly different than the CCA. This is mode mostly used
    for modeling

    Examples
    --------
    >>> from scikits.learn.pls import PLSCanonical, PLSRegression, CCA
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> plsca = PLSCanonical()
    >>> plsca.fit(X, Y, n_components=2)
    PLSCanonical(scale=True, algorithm='nipals', max_iter=500, n_components=2,
           tol=1e-06, copy=True)
    >>> X_c, Y_c = plsca.transform(X, Y)

    References
    ----------
    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In french but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    CCA
    PLSSVD
    """

    def __init__(self, n_components=2, scale=True, algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        _PLS.__init__(self, n_components=n_components,
                        deflation_mode="canonical", mode="A",
                        scale=scale, algorithm=algorithm,
                        max_iter=max_iter, tol=tol, copy=copy)


class CCA(_PLS):
    """CCA Canonical Correlation Analysis. CCA inherits from PLS with
    mode="B" and deflation_mode="canonical".

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vectors, where n_samples in the number of samples and
        p is the number of predictors.

    Y: array-like of response, shape (n_samples, q)
        Training vectors, where n_samples in the number of samples and
        q is the number of response variables.

    n_components: int, number of components to keep. (default 2).

    scale: boolean, scale data? (default True)

    algorithm: str "nipals" or "svd" the algorithm used to estimate the
        weights, it will be called "n_components" time ie.: for each iteration
        of the outer loop.

    max_iter: an integer, the maximum number of iterations (default 500) of the
        NIPALS inner loop (used only if algorithm="nipals")

    tol: a not negative real, the tolerance used in the iterative algorithm
         default 1e-06.

    copy: boolean, should the deflation been made on a copy? Let the default
        value to True unless you don't care about side effect

    Attributes
    ----------
    x_weights_: array, [p, n_components]
        X block weights vectors.

    y_weights_: array, [q, n_components]
        Y block weights vectors.

    x_loadings_: array, [p, n_components]
        X block loadings vectors.

    y_loadings_: array, [q, n_components]
        Y block loadings vectors.

    x_scores_: array, [n_samples, n_components]
        X scores.

    y_scores_: array, [n_samples, n_components]
        Y scores.

    x_rotations_: array, [p, n_components]
        X block to latents rotations.

    y_rotations_: array, [q, n_components]
        Y block to latents rotations.

    Notes
    -----
    For each component k, find the weights u, v that maximizes
    max corr(Xk u, Yk v), such that |u| = |v| = 1

    Note that it maximizes only the correlations between the scores.

    The residual matrix of X (Xk+1) block is obtained by the deflation on the
    current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current Y score.

    Examples
    --------
    >>> from scikits.learn.pls import PLSCanonical, PLSRegression, CCA
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> cca = CCA()
    >>> cca.fit(X, Y, n_components=2)
    CCA(scale=True, algorithm='nipals', max_iter=500, n_components=2, tol=1e-06,
      copy=True)
    >>> X_c, Y_c = cca.transform(X, Y)

    References
    ----------
    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In french but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    PLSCanonical
    PLSSVD
    """

    def __init__(self, n_components=2, scale=True, algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        _PLS.__init__(self, n_components=n_components,
                        deflation_mode="canonical", mode="B",
                        scale=scale, algorithm=algorithm,
                        max_iter=max_iter, tol=tol, copy=copy)


class PLSSVD(BaseEstimator):
    """Partial Least Square SVD

    Simply perform a svd on the crosscovariance matrix: X'Y
    The are no iterative deflation here.

    Parameters
    ----------
    X: array-like of predictors, shape (n_samples, p)
        Training vector, where n_samples in the number of samples and
        p is the number of predictors. X will be centered before any analysis.

    Y: array-like of response, shape (n_samples, q)
        Training vector, where n_samples in the number of samples and
        q is the number of response variables. X will be centered before any
        analysis.

    n_components: int, number of components to keep. (default 2).

    scale: boolean, scale X and Y (default True)

    Attributes
    ----------
    x_weights_: array, [p, n_components]
        X block weights vectors.

    y_weights_: array, [q, n_components]
        Y block weights vectors.

    x_scores_: array, [n_samples, n_components]
        X scores.

    y_scores_: array, [n_samples, n_components]
        Y scores.

    See also
    --------
    PLSCanonical
    CCA
    """

    def __init__(self, n_components=2, scale=True, copy=True):
        self.n_components = n_components
        self.scale = scale

    def fit(self, X, Y, **params):
        self._set_params(**params)
        # copy since this will contains the centered data
        if self.copy:
            X = np.asanyarray(X).copy()
            Y = np.asanyarray(Y).copy()
        else:
            X = np.asanyarray(X)
            Y = np.asanyarray(Y)

        n = X.shape[0]
        p = X.shape[1]

        if X.ndim != 2:
            raise ValueError('X must be a 2D array')

        if n != Y.shape[0]:
            raise ValueError(
                'Incompatible shapes: X has %s samples, while Y '
                'has %s' % (X.shape[0], Y.shape[0]))

        if self.n_components < 1 or self.n_components > p:
            raise ValueError('invalid number of components')

        # Scale (in place)
        X, Y, self.x_mean_, self.y_mean_, self.x_std_, self.y_std_ =\
            _center_scale_xy(X, Y, self.scale)
        # svd(X'Y)
        C = np.dot(X.T, Y)
        U, s, V = linalg.svd(C, full_matrices=False)
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
