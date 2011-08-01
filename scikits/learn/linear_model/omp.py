"""
Orthogonal matching pursuit algorithms based on
http://www.cs.technion.ac.il/~ronrubin/Publications/KSVX-OMP-v2.pdf

OMP introduced in S. G. Mallat, Z. Zhang, Matching pursuits with time-frequency
dictionaries, IEEE Transactions on Signal Processing, Vol. 41, No. 12.
(December 1993), pp. 3397-3415.
http://blanche.polytechnique.fr/~mallat/papiers/MallatPursuit93.pdf
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

premature = """ Orthogonal matching pursuit ended prematurely due to linear
dependence in the dictionary. The requested precision might not have been met.
"""


def _cholesky_omp(X, y, n_nonzero_coefs, eps=None):
    """
    Solves a single Orthogonal Matching Pursuit problem using
    the Cholesky decomposition.

    Parameters:
    -----------
    X: array, shape (n_samples, n_features)
        Input dictionary. Columns are assumed to have unit norm.

    y: array, shape: n_samples
        Input targets

    n_nonzero_coefs: int
        Targeted number of non-zero elements

    eps: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    Returns:
    --------
    gamma: array, shape n_nonzero_coefs
        Non-zero elements of the solution

    idx: array, shape n_nonzero_coefs
        Indices of the positions of the elements in gamma within the solution
        vector

    """

    X = X.copy('F')

    min_float = np.finfo(X.dtype).eps
    potrs, = get_lapack_funcs(('potrs',), (X,))

    alpha = np.dot(X.T, y)
    residual = y
    n_active = 0
    idx = []

    max_features = X.shape[1] if eps is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=X.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.abs(np.dot(X.T, residual)).argmax()
        if lam in idx or alpha[lam] ** 2 < min_float:
            # atom already selected or inner product too small
            warn(premature)
            break
        if n_active > 0:
            # Updates the Cholesky decomposition of X' X
            L[n_active, :n_active] = np.dot(X[:, idx].T, X[:, lam])
            v = np.dot(L[n_active, :n_active], L[n_active, :n_active])
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            if 1 - v <= min_float:  # selected atoms are dependent
                warn(premature)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        idx.append(lam)
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], alpha[idx], lower=True,
                         overwrite_b=False)

        residual = y - np.dot(X[:, idx], gamma)
        if eps is not None and np.dot(residual.T, residual) <= eps:
            break
        elif n_active == max_features:
            break

    return gamma, idx


def _gram_omp(G, Xy, n_nonzero_coefs, eps_0=None, eps=None):
    """
    Solves a single Orthogonal Matching Pursuit problem using
    the Cholesky decomposition, based on the Gram matrix and more
    precomputations.

    Parameters:
    -----------
    G: array, shape (n_features, n_features)
        Gram matrix of the input data matrix

    Xy: array, shape: n_features
        Input targets

    n_nonzero_coefs: int
        Targeted number of non-zero elements

    eps_0: float
        Squared norm of y, required if eps is not None.

    eps: float
        Targeted squared error, if not None overrides n_nonzero_coefs.

    Returns:
    --------
    gamma: array, shape n_nonzero_coefs
        Non-zero elements of the solution

    idx: array, shape n_nonzero_coefs
        Indices of the positions of the elements in gamma within the solution
        vector

    """

    min_float = np.finfo(G.dtype).eps
    potrs, = get_lapack_funcs(('potrs',), (G,))

    idx = []
    alpha = Xy
    eps_curr = eps_0
    delta = 0
    n_active = 0

    max_features = len(G) if eps is not None else n_nonzero_coefs
    L = np.empty((max_features, max_features), dtype=G.dtype)
    L[0, 0] = 1.

    while True:
        lam = np.abs(alpha).argmax()
        if lam in idx or alpha[lam] ** 2 < min_float:
            # selected same atom twice, or inner product too small
            warn(premature)
            break
        if n_active > 0:
            L[n_active, :n_active] = G[lam, idx]
            v = np.dot(L[n_active, :n_active], L[n_active, :n_active])
            solve_triangular(L[:n_active, :n_active], L[n_active, :n_active])
            if 1 - v <= min_float:  # selected atoms are dependent
                warn(premature)
                break
            L[n_active, n_active] = np.sqrt(1 - v)
        idx.append(lam)
        n_active += 1
        # solves LL'x = y as a composition of two triangular systems
        gamma, _ = potrs(L[:n_active, :n_active], Xy[idx], lower=True,
                         overwrite_b=False)

        beta = np.dot(gamma, G[idx])
        alpha = Xy - beta
        if eps is not None:
            eps_curr += delta
            delta = np.inner(gamma, beta[idx])
            eps_curr -= delta
            if eps_curr <= eps:
                break
        elif n_active == max_features:
            break

    return gamma, idx


def orthogonal_mp(X, y, n_nonzero_coefs=None, eps=None, compute_gram=False):
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
    X: array of shape (n_samples, n_features)
        Input dictionary. Columns are assumed to have unit norm.

    y: array, shape: n_samples or (n_samples, n_targets)
        Input targets

    n_nonzero_coefs: int
        Desired number of non-zero entries in the solution

    eps: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    compute_gram, bool:
        Whether to perform precomputations. Improves performance when n_targets
        or n_samples is very large.

    Returns
    -------
    coef: array of shape: n_features or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    orthogonal_mp_gram
    lars_path

    """
    X = np.asanyarray(X)
    y = np.asanyarray(y)
    if y.ndim == 1:
        y = y[:, np.newaxis]
    if n_nonzero_coefs == None and eps == None:
        n_nonzero_coefs = int(0.1 * X.shape[1])
    if eps is not None and eps < 0:
        raise ValueError("Epsilon cannot be negative")
    if eps is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if eps is None and n_nonzero_coefs > X.shape[1]:
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    if compute_gram:
        G = np.dot(X.T, X)
        Xy = np.dot(X.T, y)
        if eps is not None:
            norms_squared = np.sum((y ** 2), axis=0)
        else:
            norms_squared = None
        return orthogonal_mp_gram(G, Xy, n_nonzero_coefs, eps, norms_squared)

    coef = np.zeros((X.shape[1], y.shape[1]))
    for k in xrange(y.shape[1]):
        x, idx = _cholesky_omp(X, y[:, k], n_nonzero_coefs, eps)
        coef[idx, k] = x
    return np.squeeze(coef)


def orthogonal_mp_gram(G, Xy, n_nonzero_coefs=None, eps=None,
                       norms_squared=None):
    """Gram Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems using only
    the Gram matrix X.T * X and the product X * y.

    Parameters
    ----------
    G: array of shape (n_features, n_features)
        Gram matrix of the input data: X.T * X

    Xy: array, shape: n_features or (n_features, n_targets)
        Input targets multiplied by X: X * y

    n_nonzero_coefs: int
        Desired number of non-zero entries in the solution

    eps: float
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    norms_squared: array, shape: n_targets
        Squared L2 norms of the lines of y. Required if eps is not None.

    Returns
    -------
    coef: array of shape: n_features or (n_features, n_targets)
        Coefficients of the OMP solution

    See also
    --------
    orthogonal_mp
    lars_path

    """
    G = np.asanyarray(G)
    Xy = np.asanyarray(Xy)
    if Xy.ndim == 1:
        Xy = Xy[:, np.newaxis]
        if eps is not None:
            norms_squared = [norms_squared]

    if n_nonzero_coefs == None and eps is None:
        n_nonzero_coefs = int(0.1 * len(G))
    if eps is not None and norms_squared == None:
        raise ValueError('Gram OMP needs the precomputed norms in order \
                          to evaluate the error sum of squares.')
    if eps is not None and eps < 0:
        raise ValueError("Epsilon cennot be negative")
    if eps is None and n_nonzero_coefs <= 0:
        raise ValueError("The number of atoms must be positive")
    if eps is None and n_nonzero_coefs > len(G):
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    coef = np.zeros((len(G), Xy.shape[1]))
    for k in range(Xy.shape[1]):
        x, idx = _gram_omp(G, Xy[:, k], n_nonzero_coefs,
                           norms_squared[k] if eps is not None else None, eps)
        coef[idx, k] = x
    return np.squeeze(coef)


class OrthogonalMatchingPursuit(LinearModel):
    """Orthogonal Mathching Pursuit model (OMP)

    Parameters
    ----------
    n_nonzero_coefs: int, optional
        number of selected active features

    eps: float, optional
        Maximum norm of the residual. If not None, overrides n_nonzero_coefs.

    fit_intercept: boolean, optional
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize: boolean, optional
        If False, the regressors X are assumed to be already normalized.

    precompute : True | False | array-like, optional
        Whether to use a precomputed Gram matrix to speed up
        calculations. The Gram matrix can also be passed as argument.

    Attributes
    ----------
    `coef_` : array, shape = [n_features]
        parameter vector (w in the fomulation formula)

    `intercept_` : float
        independent term in decision function.

    References
    ----------
    TODO
    """

    def __init__(self, n_nonzero_coefs=None, eps=None, fit_intercept=True,
                 normalize=True, precompute=False):
        self.n_nonzero_coefs = n_nonzero_coefs
        self.eps = eps
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.precompute = precompute

    def fit(self, X, y, **params):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Training data.

        y : array-like, shape = (n_samples,) or (n_samples, n_targets)
            Target values.

        Returns
        -------
        self : object
            returns an instance of self.
        """
        self._set_params(**params)

        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)
        if self.normalize:
            norms = np.sqrt(np.sum(X ** 2, axis=0))
            nonzeros = np.flatnonzero(norms)
            X[:, nonzeros] /= norms[nonzeros]

        self.coef_ = orthogonal_mp(X, y, self.n_nonzero_coefs, self.eps,
                                   self.precompute).T

        if self.normalize:
            self.coef_ /= norms
        self._set_intercept(Xmean, ymean)
        return self
