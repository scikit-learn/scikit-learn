"""Orthogonal matching pursuit algorithms
"""

# Author: Vlad Niculae
#
# License: BSD Style.

from warnings import warn

import numpy as np
from scipy import linalg
from scipy.linalg.lapack import get_lapack_funcs

from .base import LinearModel
from ..utils.arrayfuncs import solve_triangular
from ..utils import as_float_array

premature = """ Orthogonal matching pursuit ended prematurely due to linear
dependence in the dictionary. The requested precision might not have been met.
"""


def _cholesky_omp(X, y, n_nonzero_coefs, eps=None, overwrite_X=False):
    """Orthogonal Matching Pursuit step using the Cholesky decomposition.

    Parameters:
    -----------
    X: array, shape = (n_samples, n_features)
        Input dictionary. Columns are assumed to have unit norm.

    y: array, shape = (n_samples,)
        Input targets

    n_nonzero_coefs: int
        Targeted number of non-zero elements

    eps: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    overwrite_X: bool,
        Whether the design matrix X can be overwritten by the algorithm. This
        is only helpful if X is already Fortran-ordered, otherwise a copy
        is made anyway.

    Returns:
    --------
    gamma: array, shape = (n_nonzero_coefs,)
        Non-zero elements of the solution

    idx: array, shape = (n_nonzero_coefs,)
        Indices of the positions of the elements in gamma within the solution
        vector

    """
    if not overwrite_X:
        X = X.copy('F')
    else:  # even if we are allowed to overwrite, still copy it if bad order
        X = np.asfortranarray(X)

    min_float = np.finfo(X.dtype).eps
    nrm2, swap = linalg.get_blas_funcs(('nrm2', 'swap'), (X,))
    potrs, = get_lapack_funcs(('potrs',), (X,))

    alpha = np.dot(X.T, y)
    residual = y
    n_active = 0
    idx = []

    max_features = X.shape[1] if eps is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=X.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.argmax(np.abs(np.dot(X.T, residual)))
        if lam < n_active or alpha[lam] ** 2 < min_float:
            # atom already selected or inner product too small
            warn(premature)
            break
        if n_active > 0:
            # Updates the Cholesky decomposition of X' X
            L[n_active, :n_active] = np.dot(X[:, :n_active].T, X[:, lam])
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            v = nrm2(L[n_active, :n_active]) ** 2
            if 1 - v <= min_float:  # selected atoms are dependent
                warn(premature)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        idx.append(lam)
        X.T[n_active], X.T[lam] = swap(X.T[n_active], X.T[lam])
        alpha[n_active], alpha[lam] = alpha[lam], alpha[n_active]
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], alpha[:n_active], lower=True,
                         overwrite_b=False)

        residual = y - np.dot(X[:, :n_active], gamma)
        if eps is not None and nrm2(residual) ** 2 <= eps:
            break
        elif n_active == max_features:
            break

    return gamma, idx


def _gram_omp(Gram, Xy, n_nonzero_coefs, eps_0=None, eps=None,
              overwrite_gram=False, overwrite_Xy=False):
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

    eps_0: float
        Squared norm of y, required if eps is not None.

    eps: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    overwrite_gram: bool,
        Whether the gram matrix can be overwritten by the algorithm. This
        is only helpful if it is already Fortran-ordered, otherwise a copy
        is made anyway.

    overwrite_Xy: bool,
        Whether the covariance vector Xy can be overwritten by the algorithm.

    Returns:
    --------
    gamma: array, shape = (n_nonzero_coefs,)
        Non-zero elements of the solution

    idx: array, shape = (n_nonzero_coefs,)
        Indices of the positions of the elements in gamma within the solution
        vector

    """
    if not overwrite_gram:
        Gram = Gram.copy('F')
    else:
        Gram = np.asfortranarray(Gram)

    if not overwrite_Xy:
        Xy = Xy.copy()

    min_float = np.finfo(Gram.dtype).eps
    nrm2, swap = linalg.get_blas_funcs(('nrm2', 'swap'), (Gram,))
    potrs, = get_lapack_funcs(('potrs',), (Gram,))

    idx = []
    alpha = Xy
    eps_curr = eps_0
    delta = 0
    n_active = 0

    max_features = len(Gram) if eps is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=Gram.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.argmax(np.abs(alpha))
        if lam < n_active or alpha[lam] ** 2 < min_float:
            # selected same atom twice, or inner product too small
            warn(premature)
            break
        if n_active > 0:
            L[n_active, :n_active] = Gram[lam, :n_active]
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            v = nrm2(L[n_active, :n_active]) ** 2
            if 1 - v <= min_float:  # selected atoms are dependent
                warn(premature)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        idx.append(lam)
        Gram[n_active], Gram[lam] = swap(Gram[n_active], Gram[lam])
        Gram.T[n_active], Gram.T[lam] = swap(Gram.T[n_active], Gram.T[lam])
        Xy[n_active], Xy[lam] = Xy[lam], Xy[n_active]
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], Xy[:n_active], lower=True,
                         overwrite_b=False)

        beta = np.dot(Gram[:, :n_active], gamma)
        alpha = Xy - beta
        if eps is not None:
            eps_curr += delta
            delta = np.inner(gamma, beta[:n_active])
            eps_curr -= delta
            if eps_curr <= eps:
                break
        elif n_active == max_features:
            break

    return gamma, idx


def orthogonal_mp(X, y, n_nonzero_coefs=None, eps=None, precompute_gram=False,
                  overwrite_X=False):
    """Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems.
    An instance of the problem has the form:

    When parametrized by the number of non-zero coefficients using
    `n_nonzero_coefs`:
    argmin ||y - X\gamma||^2 subject to ||\gamma||_0 <= n_{nonzero coefs}

    When parametrized by error using the parameter `eps`:
    argmin ||\gamma||_0 subject to ||y - X\gamma||^2 <= \epsilon

    Parameters
    ----------
    X: array, shape = (n_samples, n_features)
        Input data. Columns are assumed to have unit norm.

    y: array, shape = (n_samples,) or (n_samples, n_targets)
        Input targets

    n_nonzero_coefs: int
        Desired number of non-zero entries in the solution. If None (by
        default) this value is set to 10% of n_features.

    eps: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    precompute_gram: {True, False, 'auto'},
        Whether to perform precomputations. Improves performance when n_targets
        or n_samples is very large.

    overwrite_X: bool,
        Whether the design matrix X can be overwritten by the algorithm. This
        is only helpful if X is already Fortran-ordered, otherwise a copy
        is made anyway.

    Returns
    -------
    coef: array, shape = (n_features,) or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    OrthogonalMatchingPursuit
    orthogonal_mp_gram
    lars_path

    Notes
    -----
    Orthogonal matching pursuit was introduced in G. Mallat, Z. Zhang,
    Matching pursuits with time-frequency dictionaries, IEEE Transactions on
    Signal Processing, Vol. 41, No. 12. (December 1993), pp. 3397-3415.
    (http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf)

    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVX-OMP-v2.pdf

    """
    X, y = map(np.asanyarray, (X, y))
    if y.ndim == 1:
        y = y[:, np.newaxis]
    if not overwrite_X:
        X = X.copy('F')
        overwrite_X = True
    else:
        X = np.asfortranarray(X)
    if y.shape[1] > 1:  # subsequent targets will be affected
        overwrite_X = False
    if n_nonzero_coefs == None and eps == None:
        n_nonzero_coefs = int(0.1 * X.shape[1])
    if eps is not None and eps < 0:
        raise ValueError("Epsilon cannot be negative")
    if eps is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if eps is None and n_nonzero_coefs > X.shape[1]:
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    if precompute_gram == 'auto':
        precompute_gram = X.shape[0] > X.shape[1]
    if precompute_gram:
        G = np.dot(X.T, X)
        G = np.asfortranarray(G)
        Xy = np.dot(X.T, y)
        if eps is not None:
            norms_squared = np.sum((y ** 2), axis=0)
        else:
            norms_squared = None
        return orthogonal_mp_gram(G, Xy, n_nonzero_coefs, eps, norms_squared,
                                  overwrite_gram=overwrite_X,
                                  overwrite_Xy=True)

    coef = np.zeros((X.shape[1], y.shape[1]))
    for k in xrange(y.shape[1]):
        x, idx = _cholesky_omp(X, y[:, k], n_nonzero_coefs, eps,
                               overwrite_X=overwrite_X)
        coef[idx, k] = x
    return np.squeeze(coef)


def orthogonal_mp_gram(Gram, Xy, n_nonzero_coefs=None, eps=None,
                       norms_squared=None, overwrite_gram=False,
                       overwrite_Xy=False):
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

    eps: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    norms_squared: array-like, shape = (n_targets,)
        Squared L2 norms of the lines of y. Required if eps is not None.

    overwrite_gram: bool,
        Whether the gram matrix can be overwritten by the algorithm. This
        is only helpful if it is already Fortran-ordered, otherwise a copy
        is made anyway.

    overwrite_Xy: bool,
        Whether the covariance vector Xy can be overwritten by the algorithm.

    Returns
    -------
    coef: array, shape = (n_features,) or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    OrthogonalMatchingPursuit
    orthogonal_mp
    lars_path

    Notes
    -----
    Orthogonal matching pursuit was introduced in G. Mallat, Z. Zhang,
    Matching pursuits with time-frequency dictionaries, IEEE Transactions on
    Signal Processing, Vol. 41, No. 12. (December 1993), pp. 3397-3415.
    (http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf)

    This implementation is based on Rubinstein, R., Zibulevsky, M. and Elad,
    M., Efficient Implementation of the K-SVD Algorithm using Batch Orthogonal
    Matching Pursuit Technical Report - CS Technion, April 2008.
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVX-OMP-v2.pdf

    """
    Gram, Xy = map(np.asanyarray, (Gram, Xy))
    if Xy.ndim == 1:
        Xy = Xy[:, np.newaxis]
        if eps is not None:
            norms_squared = [norms_squared]

    if n_nonzero_coefs == None and eps is None:
        n_nonzero_coefs = int(0.1 * len(Gram))
    if eps is not None and norms_squared == None:
        raise ValueError('Gram OMP needs the precomputed norms in order \
                          to evaluate the error sum of squares.')
    if eps is not None and eps < 0:
        raise ValueError("Epsilon cennot be negative")
    if eps is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if eps is None and n_nonzero_coefs > len(Gram):
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    coef = np.zeros((len(Gram), Xy.shape[1]))
    for k in range(Xy.shape[1]):
        x, idx = _gram_omp(Gram, Xy[:, k], n_nonzero_coefs,
                           norms_squared[k] if eps is not None else None, eps,
                           overwrite_gram=overwrite_gram,
                           overwrite_Xy=overwrite_Xy)
        coef[idx, k] = x
    return np.squeeze(coef)


class OrthogonalMatchingPursuit(LinearModel):
    """Orthogonal Mathching Pursuit model (OMP)

    Parameters
    ----------
    n_nonzero_coefs: int, optional
        Desired number of non-zero entries in the solution. If None (by
        default) this value is set to 10% of n_features.

    eps: float, optional
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    fit_intercept: boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize: boolean, optional
        If False, the regressors X are assumed to be already normalized.

    precompute_gram: {True, False, 'auto'},
        Whether to use a precomputed Gram and Xy matrix to speed up
        calculations. Improves performance when `n_targets` or `n_samples` is
        very large. Note that if you already have such matrices, you can pass
        them directly to the fit method.

    overwrite_X: bool,
        Whether the design matrix X can be overwritten by the algorithm.
        This is only helpful if X is already Fortran-ordered, otherwise a
        copy is made anyway.

    overwrite_gram: bool,
        Whether the gram matrix can be overwritten by the algorithm. This
        is only helpful if it is already Fortran-ordered, otherwise a copy
        is made anyway.

    overwrite_Xy: bool,
        Whether the covariance vector Xy can be overwritten by the
        algorithm.


    Attributes
    ----------
    coef_: array, shape = (n_features,) or (n_features, n_targets)
        parameter vector (w in the fomulation formula)

    intercept_: float or array, shape =(n_targets,)
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
    http://www.cs.technion.ac.il/~ronrubin/Publications/KSVX-OMP-v2.pdf

    See also
    --------
    orthogonal_mp
    orthogonal_mp_gram
    lars_path
    Lars
    LassoLars

    """

    def __init__(self, overwrite_X=False, overwrite_gram=False,
            overwrite_Xy=False, n_nonzero_coefs=None, eps=None,
            fit_intercept=True, normalize=True, precompute_gram=False):
        self.n_nonzero_coefs = n_nonzero_coefs
        self.eps = eps
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute_gram = precompute_gram
        self.overwrite_gram = overwrite_gram
        self.overwrite_Xy = overwrite_Xy
        self.overwrite_X = overwrite_X
        self.eps = eps

    def fit(self, X, y, Gram=None, Xy=None, overwrite_x=False,
            overwrite_gram=False, overwrite_xy=False):
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
        X = np.atleast_2d(X)
        y = np.atleast_1d(y)
        n_features = X.shape[1]
        X = as_float_array(X, self.overwrite_X)

        X, y, X_mean, y_mean, X_std = self._center_data(X, y, self.fit_intercept,
                self.normalize)

        if self.n_nonzero_coefs == None and self.eps is None:
            self.n_nonzero_coefs = int(0.1 * n_features)

        if Gram is not None:
            Gram = np.atleast_2d(Gram)

            if not self.overwrite_gram:
                overwrite_gram = True
                Gram = Gram.copy('F')
            else:
                Gram = np.asfortranarray(Gram)

            if y.shape[1] > 1:  # subsequent targets will be affected
                overwrite_gram = False

            if Xy is None:
                Xy = np.dot(X.T, y)
            else:
                if not self.overwrite_Xy:
                    Xy = Xy.copy()
                if self.normalize:
                    if len(Xy.shape) == 1:
                        Xy /= X_std
                    else:
                        Xy /= X_std[:, np.newaxis]

            if self.normalize:
                Gram /= X_std
                Gram /= X_std[:, np.newaxis]

            norms_sq = np.sum((y ** 2), axis=0) if self.eps is not None else None
            self.coef_ = orthogonal_mp_gram(Gram, Xy, self.n_nonzero_coefs,
                                            self.eps, norms_sq,
                                            overwrite_gram, True).T
        else:
            precompute_gram = self.precompute_gram
            if precompute_gram == 'auto':
                precompute_gram = X.shape[0] > X.shape[1]
            self.coef_ = orthogonal_mp(X, y, self.n_nonzero_coefs, self.eps,
                                       precompute_gram=precompute_gram,
                                       overwrite_X=self.overwrite_X).T

        self._set_intercept(X_mean, y_mean, X_std)
        return self
