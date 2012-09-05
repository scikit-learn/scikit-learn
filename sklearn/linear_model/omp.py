"""Orthogonal matching pursuit algorithms
"""

# Author: Vlad Niculae
#
# License: BSD Style.

import warnings

import numpy as np
from scipy import linalg
from scipy.linalg.lapack import get_lapack_funcs

from .base import LinearModel
from ..base import RegressorMixin
from ..utils import array2d
from ..utils.arrayfuncs import solve_triangular

premature = """ Orthogonal matching pursuit ended prematurely due to linear
dependence in the dictionary. The requested precision might not have been met.
"""


def _cholesky_omp(X, y, n_nonzero_coefs, tol=None, copy_X=True):
    """Orthogonal Matching Pursuit step using the Cholesky decomposition.

    Parameters:
    -----------
    X: array, shape = (n_samples, n_features)
        Input dictionary. Columns are assumed to have unit norm.

    y: array, shape = (n_samples,)
        Input targets

    n_nonzero_coefs: int
        Targeted number of non-zero elements

    tol: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    copy_X: bool, optional
        Whether the design matrix X must be copied by the algorithm. A false
        value is only helpful if X is already Fortran-ordered, otherwise a
        copy is made anyway.

    Returns:
    --------
    gamma: array, shape = (n_nonzero_coefs,)
        Non-zero elements of the solution

    idx: array, shape = (n_nonzero_coefs,)
        Indices of the positions of the elements in gamma within the solution
        vector

    """
    if copy_X:
        X = X.copy('F')
    else:  # even if we are allowed to overwrite, still copy it if bad order
        X = np.asfortranarray(X)

    min_float = np.finfo(X.dtype).eps
    nrm2, swap = linalg.get_blas_funcs(('nrm2', 'swap'), (X,))
    potrs, = get_lapack_funcs(('potrs',), (X,))

    alpha = np.dot(X.T, y)
    residual = y
    gamma = np.empty(0)
    n_active = 0
    indices = range(X.shape[1])  # keeping track of swapping

    max_features = X.shape[1] if tol is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=X.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.argmax(np.abs(np.dot(X.T, residual)))
        if lam < n_active or alpha[lam] ** 2 < min_float:
            # atom already selected or inner product too small
            warnings.warn(premature, RuntimeWarning, stacklevel=2)
            break
        if n_active > 0:
            # Updates the Cholesky decomposition of X' X
            L[n_active, :n_active] = np.dot(X[:, :n_active].T, X[:, lam])
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            v = nrm2(L[n_active, :n_active]) ** 2
            if 1 - v <= min_float:  # selected atoms are dependent
                warnings.warn(premature, RuntimeWarning, stacklevel=2)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        X.T[n_active], X.T[lam] = swap(X.T[n_active], X.T[lam])
        alpha[n_active], alpha[lam] = alpha[lam], alpha[n_active]
        indices[n_active], indices[lam] = indices[lam], indices[n_active]
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], alpha[:n_active], lower=True,
                         overwrite_b=False)

        residual = y - np.dot(X[:, :n_active], gamma)
        if tol is not None and nrm2(residual) ** 2 <= tol:
            break
        elif n_active == max_features:
            break

    return gamma, indices[:n_active]


def _gram_omp(Gram, Xy, n_nonzero_coefs, tol_0=None, tol=None,
              copy_Gram=True, copy_Xy=True):
    """Orthogonal Matching Pursuit step on a precomputed Gram matrix.

    This function uses the the Cholesky decomposition method.

    Parameters:
    -----------
    Gram: array, shape = (n_features, n_features)
        Gram matrix of the input data matrix

    Xy: array, shape = (n_features,)
        Input targets

    n_nonzero_coefs: int
        Targeted number of non-zero elements

    tol_0: float
        Squared norm of y, required if tol is not None.

    tol: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    copy_Gram: bool, optional
        Whether the gram matrix must be copied by the algorithm. A false
        value is only helpful if it is already Fortran-ordered, otherwise a
        copy is made anyway.

    copy_Xy: bool, optional
        Whether the covariance vector Xy must be copied by the algorithm.
        If False, it may be overwritten.

    Returns:
    --------
    gamma: array, shape = (n_nonzero_coefs,)
        Non-zero elements of the solution

    idx: array, shape = (n_nonzero_coefs,)
        Indices of the positions of the elements in gamma within the solution
        vector

    """
    Gram = Gram.copy('F') if copy_Gram else np.asfortranarray(Gram)

    if copy_Xy:
        Xy = Xy.copy()

    min_float = np.finfo(Gram.dtype).eps
    nrm2, swap = linalg.get_blas_funcs(('nrm2', 'swap'), (Gram,))
    potrs, = get_lapack_funcs(('potrs',), (Gram,))

    indices = range(len(Gram))  # keeping track of swapping
    alpha = Xy
    tol_curr = tol_0
    delta = 0
    gamma = np.empty(0)
    n_active = 0

    max_features = len(Gram) if tol is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=Gram.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.argmax(np.abs(alpha))
        if lam < n_active or alpha[lam] ** 2 < min_float:
            # selected same atom twice, or inner product too small
            warnings.warn(premature, RuntimeWarning, stacklevel=2)
            break
        if n_active > 0:
            L[n_active, :n_active] = Gram[lam, :n_active]
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            v = nrm2(L[n_active, :n_active]) ** 2
            if 1 - v <= min_float:  # selected atoms are dependent
                warnings.warn(premature, RuntimeWarning, stacklevel=2)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        Gram[n_active], Gram[lam] = swap(Gram[n_active], Gram[lam])
        Gram.T[n_active], Gram.T[lam] = swap(Gram.T[n_active], Gram.T[lam])
        indices[n_active], indices[lam] = indices[lam], indices[n_active]
        Xy[n_active], Xy[lam] = Xy[lam], Xy[n_active]
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], Xy[:n_active], lower=True,
                         overwrite_b=False)

        beta = np.dot(Gram[:, :n_active], gamma)
        alpha = Xy - beta
        if tol is not None:
            tol_curr += delta
            delta = np.inner(gamma, beta[:n_active])
            tol_curr -= delta
            if tol_curr <= tol:
                break
        elif n_active == max_features:
            break

    return gamma, indices[:n_active]


def orthogonal_mp(X, y, n_nonzero_coefs=None, tol=None, precompute_gram=False,
                  copy_X=True):
    """Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems.
    An instance of the problem has the form:

    When parametrized by the number of non-zero coefficients using
    `n_nonzero_coefs`:
    argmin ||y - X\gamma||^2 subject to ||\gamma||_0 <= n_{nonzero coefs}

    When parametrized by error using the parameter `tol`:
    argmin ||\gamma||_0 subject to ||y - X\gamma||^2 <= tol

    Parameters
    ----------
    X: array, shape = (n_samples, n_features)
        Input data. Columns are assumed to have unit norm.

    y: array, shape = (n_samples,) or (n_samples, n_targets)
        Input targets

    n_nonzero_coefs: int
        Desired number of non-zero entries in the solution. If None (by
        default) this value is set to 10% of n_features.

    tol: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    precompute_gram: {True, False, 'auto'},
        Whether to perform precomputations. Improves performance when n_targets
        or n_samples is very large.

    copy_X: bool, optional
        Whether the design matrix X must be copied by the algorithm. A false
        value is only helpful if X is already Fortran-ordered, otherwise a
        copy is made anyway.

    Returns
    -------
    coef: array, shape = (n_features,) or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    OrthogonalMatchingPursuit
    orthogonal_mp_gram
    lars_path
    decomposition.sparse_encode

    Notes
    -----
    Orthogonal matching pursuit was introduced in G. Mallat, Z. Zhang,
    Matching pursuits with time-frequency dictionaries, IEEE Transactions on
    Signal Processing, Vol. 41, No. 12. (December 1993), pp. 3397-3415.
    (http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf)

    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

    """
    X = array2d(X, order='F', copy=copy_X)
    copy_X = False
    y = np.asarray(y)
    if y.ndim == 1:
        y = y[:, np.newaxis]
    if y.shape[1] > 1:  # subsequent targets will be affected
        copy_X = True
    if n_nonzero_coefs == None and tol == None:
        n_nonzero_coefs = int(0.1 * X.shape[1])
    if tol is not None and tol < 0:
        raise ValueError("Epsilon cannot be negative")
    if tol is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if tol is None and n_nonzero_coefs > X.shape[1]:
        raise ValueError("The number of atoms cannot be more than the number "
                         "of features")
    if precompute_gram == 'auto':
        precompute_gram = X.shape[0] > X.shape[1]
    if precompute_gram:
        G = np.dot(X.T, X)
        G = np.asfortranarray(G)
        Xy = np.dot(X.T, y)
        if tol is not None:
            norms_squared = np.sum((y ** 2), axis=0)
        else:
            norms_squared = None
        return orthogonal_mp_gram(G, Xy, n_nonzero_coefs, tol, norms_squared,
                                  copy_Gram=copy_X, copy_Xy=False)

    coef = np.zeros((X.shape[1], y.shape[1]))
    for k in xrange(y.shape[1]):
        x, idx = _cholesky_omp(X, y[:, k], n_nonzero_coefs, tol,
                               copy_X=copy_X)
        coef[idx, k] = x
    return np.squeeze(coef)


def orthogonal_mp_gram(Gram, Xy, n_nonzero_coefs=None, tol=None,
                       norms_squared=None, copy_Gram=True,
                       copy_Xy=True):
    """Gram Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems using only
    the Gram matrix X.T * X and the product X.T * y.

    Parameters
    ----------
    Gram: array, shape = (n_features, n_features)
        Gram matrix of the input data: X.T * X

    Xy: array, shape = (n_features,) or (n_features, n_targets)
        Input targets multiplied by X: X.T * y

    n_nonzero_coefs: int
        Desired number of non-zero entries in the solution. If None (by
        default) this value is set to 10% of n_features.

    tol: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    norms_squared: array-like, shape = (n_targets,)
        Squared L2 norms of the lines of y. Required if tol is not None.

    copy_Gram: bool, optional
        Whether the gram matrix must be copied by the algorithm. A false
        value is only helpful if it is already Fortran-ordered, otherwise a
        copy is made anyway.

    copy_Xy: bool, optional
        Whether the covariance vector Xy must be copied by the algorithm.
        If False, it may be overwritten.

    Returns
    -------
    coef: array, shape = (n_features,) or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    OrthogonalMatchingPursuit
    orthogonal_mp
    lars_path
    decomposition.sparse_encode

    Notes
    -----
    Orthogonal matching pursuit was introduced in G. Mallat, Z. Zhang,
    Matching pursuits with time-frequency dictionaries, IEEE Transactions on
    Signal Processing, Vol. 41, No. 12. (December 1993), pp. 3397-3415.
    (http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf)

    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

    """
    Gram = array2d(Gram, order='F', copy=copy_Gram)
    Xy = np.asarray(Xy)
    if Xy.ndim > 1 and Xy.shape[1] > 1:  # or subsequent target will be affected
        copy_Gram = True
    if Xy.ndim == 1:
        Xy = Xy[:, np.newaxis]
        if tol is not None:
            norms_squared = [norms_squared]

    if n_nonzero_coefs == None and tol is None:
        n_nonzero_coefs = int(0.1 * len(Gram))
    if tol is not None and norms_squared == None:
        raise ValueError('Gram OMP needs the precomputed norms in order '
                         'to evaluate the error sum of squares.')
    if tol is not None and tol < 0:
        raise ValueError("Epsilon cannot be negative")
    if tol is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if tol is None and n_nonzero_coefs > len(Gram):
        raise ValueError("The number of atoms cannot be more than the number "
                          "of features")
    coef = np.zeros((len(Gram), Xy.shape[1]))
    for k in range(Xy.shape[1]):
        x, idx = _gram_omp(Gram, Xy[:, k], n_nonzero_coefs,
                           norms_squared[k] if tol is not None else None, tol,
                           copy_Gram=copy_Gram, copy_Xy=copy_Xy)
        coef[idx, k] = x
    return np.squeeze(coef)


class OrthogonalMatchingPursuit(LinearModel, RegressorMixin):
    """Orthogonal Mathching Pursuit model (OMP)

    Parameters
    ----------
    n_nonzero_coefs : int, optional
        Desired number of non-zero entries in the solution. If None (by
        default) this value is set to 10% of n_features.

    tol : float, optional
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    fit_intercept : boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If False, the regressors X are assumed to be already normalized.

    precompute_gram : {True, False, 'auto'},
        Whether to use a precomputed Gram and Xy matrix to speed up
        calculations. Improves performance when `n_targets` or `n_samples` is
        very large. Note that if you already have such matrices, you can pass
        them directly to the fit method.

    copy_X : bool, optional
        Whether the design matrix X must be copied by the algorithm. A false
        value is only helpful if X is already Fortran-ordered, otherwise a
        copy is made anyway.

    copy_Gram : bool, optional
        Whether the gram matrix must be copied by the algorithm. A false
        value is only helpful if X is already Fortran-ordered, otherwise a
        copy is made anyway.

    copy_Xy : bool, optional
        Whether the covariance vector Xy must be copied by the algorithm.
        If False, it may be overwritten.


    Attributes
    ----------
    `coef_` : array, shape = (n_features,) or (n_features, n_targets)
        parameter vector (w in the fomulation formula)

    `intercept_` : float or array, shape =(n_targets,)
        independent term in decision function.

    Notes
    -----
    Orthogonal matching pursuit was introduced in G. Mallat, Z. Zhang,
    Matching pursuits with time-frequency dictionaries, IEEE Transactions on
    Signal Processing, Vol. 41, No. 12. (December 1993), pp. 3397-3415.
    (http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf)

    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

    See also
    --------
    orthogonal_mp
    orthogonal_mp_gram
    lars_path
    Lars
    LassoLars
    decomposition.sparse_encode

    """
    def __init__(self, copy_X=True, copy_Gram=True,
            copy_Xy=True, n_nonzero_coefs=None, tol=None,
            fit_intercept=True, normalize=True, precompute_gram=False):
        self.n_nonzero_coefs = n_nonzero_coefs
        self.tol = tol
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute_gram = precompute_gram
        self.copy_Gram = copy_Gram
        self.copy_Xy = copy_Xy
        self.copy_X = copy_X

    def fit(self, X, y, Gram=None, Xy=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X: array-like, shape = (n_samples, n_features)
            Training data.

        y: array-like, shape = (n_samples,) or (n_samples, n_targets)
            Target values.

        Gram: array-like, shape = (n_features, n_features) (optional)
            Gram matrix of the input data: X.T * X

        Xy: array-like, shape = (n_features,) or (n_features, n_targets)
            (optional)
            Input targets multiplied by X: X.T * y


        Returns
        -------
        self: object
            returns an instance of self.
        """
        X = array2d(X)
        y = np.asarray(y)
        n_features = X.shape[1]

        X, y, X_mean, y_mean, X_std = self._center_data(X, y,
                                                        self.fit_intercept,
                                                        self.normalize,
                                                        self.copy_X)

        if y.ndim == 1:
            y = y[:, np.newaxis]

        if self.n_nonzero_coefs == None and self.tol is None:
            self.n_nonzero_coefs = int(0.1 * n_features)
        if (Gram is not None or Xy is not None) and (self.fit_intercept is True
                                                 or self.normalize is True):
            warnings.warn('Mean subtraction (fit_intercept) and '
                 'normalization cannot be applied on precomputed Gram '
                 'and Xy matrices. Your precomputed values are ignored '
                 'and recomputed. To avoid this, do the scaling yourself '
                 'and call with fit_intercept and normalize set to False.',
                 RuntimeWarning, stacklevel=2)
            Gram, Xy = None, None

        if Gram is not None:
            if Xy is None:
                Xy = np.dot(X.T, y)
            else:
                if self.copy_Xy:
                    Xy = Xy.copy()
                if self.normalize:
                    if len(Xy.shape) == 1:
                        Xy /= X_std
                    else:
                        Xy /= X_std[:, np.newaxis]

            if self.normalize:
                Gram /= X_std
                Gram /= X_std[:, np.newaxis]

            norms_sq = np.sum(y ** 2, axis=0) if self.tol is not None else None
            self.coef_ = orthogonal_mp_gram(Gram, Xy, self.n_nonzero_coefs,
                                            self.tol, norms_sq,
                                            self.copy_Gram, True).T
        else:
            precompute_gram = self.precompute_gram
            if precompute_gram == 'auto':
                precompute_gram = X.shape[0] > X.shape[1]
            self.coef_ = orthogonal_mp(X, y, self.n_nonzero_coefs, self.tol,
                                       precompute_gram=self.precompute_gram,
                                       copy_X=self.copy_X).T

        self._set_intercept(X_mean, y_mean, X_std)
        return self
