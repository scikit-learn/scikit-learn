"""
The :mod:`sklearn.pls` module implements Partial Least Squares (PLS).
"""

# Author: Edouard Duchesnay <edouard.duchesnay@cea.fr>
# License: BSD 3 clause

import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.linalg import pinv2, svd
from scipy.sparse.linalg import svds

from ..base import BaseEstimator, RegressorMixin, TransformerMixin
from ..base import MultiOutputMixin
from ..utils import check_array, check_consistent_length
from ..utils.extmath import svd_flip
from ..utils.validation import check_is_fitted, FLOAT_DTYPES
from ..utils.validation import _deprecate_positional_args
from ..exceptions import ConvergenceWarning
from ..utils.deprecation import deprecated

__all__ = ['PLSCanonical', 'PLSRegression', 'PLSSVD']


def _get_first_singular_vectors_power_method(X, Y, mode="A", max_iter=500,
                                             tol=1e-06, norm_y_weights=False):
    """Return the first left and right singular vectors of X'Y.

    Provides an alternative to the svd(X'Y) and uses the power method instead.
    With norm_y_weights to True and in mode A, this corresponds to the
    algorithm section 11.3 of the Wegelin's review, except this starts at the
    "update saliences" part.
    """

    eps = np.finfo(X.dtype).eps
    y_score = next(col for col in Y.T if np.any(np.abs(col) > eps))
    x_weights_old = 100  # init to big value for first convergence check

    if mode == 'B':
        # Precompute pseudo inverse matrices
        # Basically: X_pinv = (X.T X)^-1 X.T
        # Which requires inverting a (n_features, n_features) matrix.
        # As a result, and as detailed in the Wegelin's review, CCA (i.e. mode
        # B) will be unstable if n_features > n_samples or n_targets >
        # n_samples
        X_pinv = pinv2(X, check_finite=False)
        Y_pinv = pinv2(Y, check_finite=False)

    for i in range(max_iter):
        if mode == "B":
            x_weights = np.dot(X_pinv, y_score)
        else:
            x_weights = np.dot(X.T, y_score) / np.dot(y_score, y_score)

        x_weights /= np.sqrt(np.dot(x_weights, x_weights)) + eps
        x_score = np.dot(X, x_weights)

        if mode == "B":
            y_weights = np.dot(Y_pinv, x_score)
        else:
            y_weights = np.dot(Y.T, x_score) / np.dot(x_score.T, x_score)

        if norm_y_weights:
            y_weights /= np.sqrt(np.dot(y_weights, y_weights)) + eps

        # FIXME: what's with the division here??
        y_score = np.dot(Y, y_weights) / (np.dot(y_weights, y_weights) + eps)

        x_weights_diff = x_weights - x_weights_old
        if np.dot(x_weights_diff, x_weights_diff) < tol or Y.shape[1] == 1:
            break
        x_weights_old = x_weights

    n_iter = i + 1
    if n_iter == max_iter:
        warnings.warn('Maximum number of iterations reached',
                      ConvergenceWarning)

    return x_weights, y_weights, n_iter


def _get_first_singular_vectors_svd(X, Y):
    """Return the first left and right singular vectors of X'Y.

    Here the whole SVD is computed.
    """
    C = np.dot(X.T, Y)
    U, _, Vt = svd(C, full_matrices=False)
    return U[:, 0], Vt[0, :]


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


def _svd_flip_1d(u, v):
    """Same as svd_flip but works on 1d arrays, and is inplace"""
    # svd_flip would force us to convert to 2d array and would also return 2d
    # arrays. We don't want that.
    biggest_abs_val_idx = np.argmax(np.abs(u))
    sign = np.sign(u[biggest_abs_val_idx])
    u *= sign
    v *= sign


class _PLS(TransformerMixin, RegressorMixin, MultiOutputMixin, BaseEstimator,
           metaclass=ABCMeta):
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

    norm_y_weights : boolean, normalize Y weights to one? (default False)

    algorithm : string, "nipals" or "svd"
        The algorithm used to estimate the weights. It will be called
        n_components times, i.e. once for each iteration of the outer loop.

    max_iter : int (default 500)
        The maximum number of iterations
        of the NIPALS inner loop (used only if algorithm="nipals")

    tol : non-negative real, default 1e-06
        The tolerance used in the iterative algorithm.

    copy : boolean, default True
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

    x_mean_ : array, [p]
        X mean for each predictor.

    y_mean_ : array, [q]
        Y mean for each response variable.

    x_std_ : array, [p]
        X standard deviation for each predictor.

    y_std_ : array, [q]
        Y standard deviation for each response variable.

    coef_ : array, [p, q]
        The coefficients of the linear model: ``Y = X coef_ + Err``

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component. Not useful if the algorithm given is "svd".

    References
    ----------

    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    See also
    --------
    PLSCanonical
    PLSRegression
    CCA
    PLS_SVD
    """

    @abstractmethod
    def __init__(self, n_components=2, *, scale=True,
                 deflation_mode="regression",
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
        X : array-like of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of predictors.

        Y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target vectors, where n_samples is the number of samples and
            n_targets is the number of response variables.
        """

        check_consistent_length(X, Y)
        X = self._validate_data(X, dtype=np.float64, copy=self.copy,
                                ensure_min_samples=2)
        Y = check_array(Y, dtype=np.float64, copy=self.copy, ensure_2d=False)
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)

        n = X.shape[0]
        p = X.shape[1]
        q = Y.shape[1]

        max_n_components = min(n, p, q)  # see Wegelin page 12
        # if not 1 <= self.n_components <= max_n_components:
        if not 1 <= self.n_components <= p:
            raise ValueError(
                f"n_components({self.n_components}) should be no lower than "
                f"1 and no greater than n_features={p}."
            )
        if self.algorithm not in ("svd", "nipals"):
            raise ValueError("Got algorithm %s when only 'svd' "
                             "and 'nipals' are known" % self.algorithm)
        if self.algorithm == "svd" and self.mode == "B":
            raise ValueError('Incompatible configuration: mode B is not '
                             'implemented with svd algorithm')
        if self.deflation_mode not in ["canonical", "regression"]:
            raise ValueError('The deflation mode is unknown')

        assert self.norm_y_weights == (self.deflation_mode == 'canonical')

        # Scale (in place)
        Xk, Yk, self.x_mean_, self.y_mean_, self.x_std_, self.y_std_ = (
            _center_scale_xy(X, Y, self.scale))

        # Results matrices
        self.x_scores_ = np.zeros((n, self.n_components))
        self.y_scores_ = np.zeros((n, self.n_components))
        self.x_weights_ = np.zeros((p, self.n_components))
        self.y_weights_ = np.zeros((q, self.n_components))
        self.x_loadings_ = np.zeros((p, self.n_components))
        self.y_loadings_ = np.zeros((q, self.n_components))
        self.n_iter_ = []

        # This whole thing corresponds to the algorithm in section 4.1 of the
        # review from Wegelin. See below for a notation mapping from code to
        # paper.
        Y_eps = np.finfo(Yk.dtype).eps
        for k in range(self.n_components):
            if np.all(np.dot(Yk.T, Yk) < np.finfo(np.double).eps):
                # Yk constant
                # That thing is raised with W2A and PLSSVD whenever target is
                # single output??????
                warnings.warn('Y residual constant at iteration %s' % k)
                break
            # Find first left and right singular vectors of the X.T.dot(Y)
            # cross-covariance matrix.
            if self.algorithm == "nipals":
                # Replace columns that are all close to zero with zeros
                Yk_mask = np.all(np.abs(Yk) < 10 * Y_eps, axis=0)
                Yk[:, Yk_mask] = 0.0

                x_weights, y_weights, n_iter_ = \
                    _get_first_singular_vectors_power_method(
                        Xk, Yk, mode=self.mode, max_iter=self.max_iter,
                        tol=self.tol, norm_y_weights=self.norm_y_weights)
                self.n_iter_.append(n_iter_)

            elif self.algorithm == "svd":
                x_weights, y_weights = _get_first_singular_vectors_svd(Xk, Yk)

            # inplace sign flip for consistency across solvers and archs
            _svd_flip_1d(x_weights, y_weights)

            # compute scores, i.e. the projections of X and Y
            x_scores = np.dot(Xk, x_weights)
            if self.norm_y_weights:
                y_ss = 1
            else:
                y_ss = np.dot(y_weights, y_weights)
            y_scores = np.dot(Yk, y_weights) / y_ss

            # Deflation: subtract rank-one approx to obtain Xk+1 and Yk+1
            x_loadings = np.dot(x_scores, Xk) / np.dot(x_scores, x_scores)
            Xk -= np.outer(x_scores, x_loadings)

            if self.deflation_mode == "canonical":
                # regress Yk on y_score
                y_loadings = np.dot(y_scores, Yk) / np.dot(y_scores, y_scores)
                Yk -= np.outer(y_scores, y_loadings)
            if self.deflation_mode == "regression":
                # regress Yk on x_score
                y_loadings = np.dot(x_scores, Yk) / np.dot(x_scores, x_scores)
                Yk -= np.outer(x_scores, y_loadings)

            self.x_weights_[:, k] = x_weights  # U
            self.y_weights_[:, k] = y_weights  # V
            self.x_scores_[:, k] = x_scores  # Xi
            self.y_scores_[:, k] = y_scores  # Omega
            self.x_loadings_[:, k] = x_loadings  # Gamma
            self.y_loadings_[:, k] = y_loadings  # Delta

        # X was approximated as Xi . Gamma.T + X_(R+1)
        # Xi . Gamma.T is a sum of n_components rank-1 matrices. X_(R+1) is
        # whatever is left to fully reconstruct X, and can be 0 if X is of rank
        # n_components.
        # Similiarly, Y was approximated as Omega . Delta.T + Y_(R+1)

        # Transformation matrices (rotations_)
        # X_rot = U (Gamma.T . U)^-1
        # Y_rot = V (Omega.T . V)^-1
        # Note that, with the training data and ommitting X_(R+1), we have
        # transform(X) =def= X . X_rot == Xi . Gamma.T . X_rot = Xi
        # Same for tranform(Y_train)
        self.x_rotations_ = np.dot(
            self.x_weights_,
            pinv2(np.dot(self.x_loadings_.T, self.x_weights_),
                  check_finite=False))
        self.y_rotations_ = np.dot(
            self.y_weights_, pinv2(np.dot(self.y_loadings_.T, self.y_weights_),
                                   check_finite=False))

        if True or self.deflation_mode == "regression":
            # FIXME what's with the if?
            # We regress Y on the transformed training data, i.e. Xi
            # Which means the model is Y_pred = Xi . Delta.T + err,
            # Xi is by def X . X_rot
            # which leads to Y_pred = X . X_rot . Delta.T + err
            # and so the coefficients matrix is X_rot . Delta.T
            self.coef_ = np.dot(self.x_rotations_, self.y_loadings_.T)
            self.coef_ = self.coef_ * self.y_std_
        return self

    def transform(self, X, Y=None, copy=True):
        """Apply the dimension reduction learned on the train data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of predictors.

        Y : array-like of shape (n_samples, n_targets)
            Target vectors, where n_samples is the number of samples and
            n_targets is the number of response variables.

        copy : boolean, default True
            Whether to copy X and Y, or perform in-place normalization.

        Returns
        -------
        x_scores if Y is not given, (x_scores, y_scores) otherwise.
        """
        check_is_fitted(self)
        X = check_array(X, copy=copy, dtype=FLOAT_DTYPES)
        # Normalize
        X -= self.x_mean_
        X /= self.x_std_
        # Apply rotation
        x_scores = np.dot(X, self.x_rotations_)
        if Y is not None:
            Y = check_array(Y, ensure_2d=False, copy=copy, dtype=FLOAT_DTYPES)
            if Y.ndim == 1:
                Y = Y.reshape(-1, 1)
            Y -= self.y_mean_
            Y /= self.y_std_
            y_scores = np.dot(Y, self.y_rotations_)
            return x_scores, y_scores

        return x_scores

    def inverse_transform(self, X):
        """Transform data back to its original space.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_components)
            New data, where n_samples is the number of samples
            and n_components is the number of pls components.

        Returns
        -------
        x_reconstructed : array-like of shape (n_samples, n_features)

        Notes
        -----
        This transformation will only be exact if n_components=n_features
        """
        check_is_fitted(self)
        X = check_array(X, dtype=FLOAT_DTYPES)
        # From pls space to original space
        X_reconstructed = np.matmul(X, self.x_loadings_.T)

        # Denormalize
        X_reconstructed *= self.x_std_
        X_reconstructed += self.x_mean_
        return X_reconstructed

    def predict(self, X, copy=True):
        """Apply the dimension reduction learned on the train data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of predictors.

        copy : boolean, default True
            Whether to copy X and Y, or perform in-place normalization.

        Notes
        -----
        This call requires the estimation of a p x q matrix, which may
        be an issue in high dimensional space.
        """
        check_is_fitted(self)
        X = check_array(X, copy=copy, dtype=FLOAT_DTYPES)
        # Normalize
        X -= self.x_mean_
        X /= self.x_std_
        Ypred = np.dot(X, self.coef_)
        return Ypred + self.y_mean_

    def fit_transform(self, X, y=None):
        """Learn and apply the dimension reduction on the train data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples and
            n_features is the number of predictors.

        y : array-like of shape (n_samples, n_targets)
            Target vectors, where n_samples is the number of samples and
            n_targets is the number of response variables.

        Returns
        -------
        x_scores if Y is not given, (x_scores, y_scores) otherwise.
        """
        return self.fit(X, y).transform(X, y)

    def _more_tags(self):
        return {'poor_score': True,
                'requires_y': False}


class PLSRegression(_PLS):
    """PLS regression

    PLSRegression implements the PLS 2 blocks regression known as PLS2 or PLS1
    in case of one dimensional response.
    This class inherits from _PLS with mode="A", deflation_mode="regression",
    norm_y_weights=False and algorithm="nipals".

    Read more in the :ref:`User Guide <cross_decomposition>`.

    .. versionadded:: 0.8

    Parameters
    ----------
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

    coef_ : array, [p, q]
        The coefficients of the linear model: ``Y = X coef_ + Err``

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component.

    Notes
    -----
    Matrices::

        T: x_scores_
        U: y_scores_
        W: x_weights_
        C: y_weights_
        P: x_loadings_
        Q: y_loadings_

    Are computed such that::

        X = T P.T + Err and Y = U Q.T + Err
        T[:, k] = Xk W[:, k] for k in range(n_components)
        U[:, k] = Yk C[:, k] for k in range(n_components)
        x_rotations_ = W (P.T W)^(-1)
        y_rotations_ = C (Q.T C)^(-1)

    where Xk and Yk are residual matrices at iteration k.

    `Slides explaining
    PLS <http://www.eigenvector.com/Docs/Wise_pls_properties.pdf>`_


    For each component k, find weights u, v that optimizes:
    ``max corr(Xk u, Yk v) * std(Xk u) std(Yk u)``, such that ``|u| = 1``

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
    PLSRegression()
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
    @_deprecate_positional_args
    def __init__(self, n_components=2, *, scale=True,
                 max_iter=500, tol=1e-06, copy=True):
        super().__init__(
            n_components=n_components, scale=scale,
            deflation_mode="regression", mode="A",
            norm_y_weights=False, algorithm='nipals', max_iter=max_iter,
            tol=tol, copy=copy)


class PLSCanonical(_PLS):
    """ PLSCanonical implements the 2 blocks canonical PLS of the original Wold
    algorithm [Tenenhaus 1998] p.204, referred as PLS-W2A in [Wegelin 2000].

    This class inherits from PLS with mode="A" and deflation_mode="canonical",
    norm_y_weights=True and algorithm="nipals", but svd should provide similar
    results up to numerical errors.

    Read more in the :ref:`User Guide <cross_decomposition>`.

    .. versionadded:: 0.8

    Parameters
    ----------
    n_components : int, (default 2).
        Number of components to keep

    scale : boolean, (default True)
        Option to scale data

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

    coef_ : array of shape (p, q)
        The coefficients of the linear model: ``Y = X coef_ + Err``

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component. Not useful if the algorithm provided is "svd".

    Notes
    -----
    Matrices::

        T: x_scores_
        U: y_scores_
        W: x_weights_
        C: y_weights_
        P: x_loadings_
        Q: y_loadings__

    Are computed such that::

        X = T P.T + Err and Y = U Q.T + Err
        T[:, k] = Xk W[:, k] for k in range(n_components)
        U[:, k] = Yk C[:, k] for k in range(n_components)
        x_rotations_ = W (P.T W)^(-1)
        y_rotations_ = C (Q.T C)^(-1)

    where Xk and Yk are residual matrices at iteration k.

    `Slides explaining PLS
    <http://www.eigenvector.com/Docs/Wise_pls_properties.pdf>`_

    For each component k, find weights u, v that optimize::

        max corr(Xk u, Yk v) * std(Xk u) std(Yk u), such that ``|u| = |v| = 1``

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
    PLSCanonical()
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
    @_deprecate_positional_args
    def __init__(self, n_components=2, *, scale=True, algorithm="nipals",
                 max_iter=500, tol=1e-06, copy=True):
        super().__init__(
            n_components=n_components, scale=scale,
            deflation_mode="canonical", mode="A",
            norm_y_weights=True, algorithm=algorithm,
            max_iter=max_iter, tol=tol, copy=copy)


class PLSSVD(TransformerMixin, BaseEstimator):
    """Partial Least Square SVD.

    This transformer simply performs a SVD on the crosscovariance matrix X'Y.
    It is able to project both the training data `X` and the targets `Y`. The
    training data X is projected on the left singular vectors, while the
    targets are projected on the right singular vectors.

    Read more in the :ref:`User Guide <cross_decomposition>`.

    .. versionadded:: 0.8

    Parameters
    ----------
    n_components : int, default 2
        The number of components to keep. Should be in `[1,
        min(n_features, n_targets)]`.

    scale : boolean, default True
        Whether to scale `X` and `Y`.

    copy : boolean, default True
        Whether to copy `X` and `Y` in fit before applying centering, and
        potentially scaling. If False, these operations will be done inplace,
        modifying both arrays.

    Attributes
    ----------
    x_weights_ : ndarray of shape (n_features, n_components)
        The left singular vectors of the SVD of the cross-covariance matrix.
        Used to project `X` in `transform`.

    y_weights_ : ndarray of (n_targets, n_components)
        The right singular vectors of the SVD of the cross-covariance matrix.
        Used to project `X` in `transform`.

    x_scores_ : ndarray of shape (n_samples, n_components)
        The transformed training samples.

        .. deprecated:: 0.24
           `x_scores_` is deprecated in 0.24 and will be removed in 0.26. You
           can just call `transform` on the training data instead.

    y_scores_ : ndarray of shape (n_samples, n_components)
        The transformed training targets.

        .. deprecated:: 0.24
           `y_scores_` is deprecated in 0.24 and will be removed in 0.26. You
           can just call `transform` on the training data instead.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.cross_decomposition import PLSSVD
    >>> X = np.array([[0., 0., 1.],
    ...               [1., 0., 0.],
    ...               [2., 2., 2.],
    ...               [2., 5., 4.]])
    >>> Y = np.array([[0.1, -0.2],
    ...               [0.9, 1.1],
    ...               [6.2, 5.9],
    ...               [11.9, 12.3]])
    >>> pls = PLSSVD(n_components=2).fit(X, Y)
    >>> X_c, Y_c = pls.transform(X, Y)
    >>> X_c.shape, Y_c.shape
    ((4, 2), (4, 2))

    See also
    --------
    PLSCanonical
    CCA
    """
    @_deprecate_positional_args
    def __init__(self, n_components=2, *, scale=True, copy=True):
        self.n_components = n_components
        self.scale = scale
        self.copy = copy

    def fit(self, X, Y):
        """Fit model to data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training samples

        Y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Targets
        """
        check_consistent_length(X, Y)
        X = self._validate_data(X, dtype=np.float64, copy=self.copy,
                                ensure_min_samples=2)
        Y = check_array(Y, dtype=np.float64, copy=self.copy, ensure_2d=False)
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)

        # we'll compute the SVD of the cross-covariance matrix = X.T.dot(Y)
        # This matrix rank is at most min(n_features, n_targets) so
        # n_components cannot be bigger than that.
        n_components = self.n_components
        rank_upper_bound = min(X.shape[1], Y.shape[1])
        if not 1 <= n_components <= rank_upper_bound:
            # TODO: raise an error in 0.26
            warnings.warn(
                f"As of version 0.24, n_components({n_components}) should be "
                f"in [1, min(n_features, n_targets)] = "
                f"[1, {rank_upper_bound}]. "
                f"n_components={rank_upper_bound} will be used instead. "
                f"In version 0.26, an error will be raised.",
                FutureWarning
            )
            n_components = rank_upper_bound

        X, Y, self.x_mean_, self.y_mean_, self.x_std_, self.y_std_ = (
            _center_scale_xy(X, Y, self.scale))

        # Compute SVD of cross-covariance matrix
        C = np.dot(X.T, Y)
        U, s, Vt = svd(C, full_matrices=False)
        U = U[:, :n_components]
        Vt = Vt[:n_components]
        U, Vt = svd_flip(U, Vt)
        V = Vt.T

        self._x_scores = np.dot(X, U)  # TODO: remove in 0.26
        self._y_scores = np.dot(Y, V)  # TODO: remove in 0.26
        self.x_weights_ = U
        self.y_weights_ = V
        return self

    @deprecated("Attribute x_scores_ was deprecated in version 0.24 and "
                "will be removed in 0.26. Use est.transform(X) on the "
                "training data instead.")
    @property
    def x_scores_(self):
        return self._x_scores

    @deprecated("Attribute y_scores_ was deprecated in version 0.24 and "
                "will be removed in 0.26. Use est.transform(X, Y) on the "
                "training data instead.")
    @property
    def y_scores_(self):
        return self._y_scores

    def transform(self, X, Y=None):
        """
        Apply the dimensionality reduction learned on the training data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training samples

        Y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Targets

        Returns
        -------
        out : array-like or tuple of array-like:
            The transformed data X_tranformed if Y is not None,
            (X_transformed, Y_transformed) otherwise
        """
        check_is_fitted(self)
        X = check_array(X, dtype=np.float64)
        Xr = (X - self.x_mean_) / self.x_std_
        x_scores = np.dot(Xr, self.x_weights_)
        if Y is not None:
            Y = check_array(Y, ensure_2d=False, dtype=np.float64)
            if Y.ndim == 1:
                Y = Y.reshape(-1, 1)
            Yr = (Y - self.y_mean_) / self.y_std_
            y_scores = np.dot(Yr, self.y_weights_)
            return x_scores, y_scores
        return x_scores

    def fit_transform(self, X, y=None):
        """Learn and apply the dimensionality reduction on the training data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training samples

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Targets

        Returns
        -------
        out : array-like or tuple of array-like:
            The transformed data X_tranformed if Y is not None,
            (X_transformed, Y_transformed) otherwise
        """
        return self.fit(X, y).transform(X, y)
