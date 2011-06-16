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

import numpy as np
from scipy import linalg


def _cholesky_omp(X, y, n_atoms, eps=None):
    """
    Solves a single Orthogonal Matching Pursuit problem using
    the Cholesky decomposition.

    Parameters:
    -----------
    X: array, shape (n_samples, n_features)
        Input dictionary

    y: array, shape: n_samples
        Input targets

    n_atoms: int
        Targeted number of non-zero elements

    eps: float
        Targeted squared error, if not None overrides n_atoms.

    Returns:
    --------
    gamma: array, shape n_atoms
        Non-zero elements of the solution

    idx: array, shape n_atoms
        Indices of the positions of the elements in gamma within the solution
        vector

    """

    alpha = np.dot(X.T, y)
    residual = y
    it = 0
    idx = []
    L = np.ones((1, 1))

    for it in xrange(X.shape[1]):
        lam = np.abs(np.dot(X.T, residual)).argmax()
        if len(idx) > 0:
            w = linalg.solve_triangular(L, np.dot(X[:, idx].T, X[:, lam]),
                                        lower=True, unit_diagonal=True)
            # Updates the Cholesky decomposition of X' X
            L = np.r_[np.c_[L, np.zeros(len(L))],
                      np.atleast_2d(np.append(w, np.sqrt(1 - np.dot(w.T, w))))]
        idx.append(lam)
        # solves LL'x = y as a composition of two triangular systems
        Lc = linalg.solve_triangular(L, alpha[idx], lower=True)
        gamma = linalg.solve_triangular(L, Lc, trans=1, lower=True)
        residual = y - np.dot(X[:, idx], gamma)
        if eps != None and np.dot(residual.T, residual) <= eps:
            break
        elif n_atoms != None and it == n_atoms - 1:
            break

    return gamma, idx


def _gram_omp(G, Xy, n_atoms, eps_0=None, eps=None):
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

    n_atoms: int
        Targeted number of non-zero elements

    eps_0: float
        Squared norm of y, required if eps is not None.

    eps: float
        Targeted squared error, if not None overrides n_atoms.

    Returns:
    --------
    gamma: array, shape n_atoms
        Non-zero elements of the solution

    idx: array, shape n_atoms
        Indices of the positions of the elements in gamma within the solution
        vector

    """

    idx = []
    L = np.ones((1, 1))
    alpha = Xy
    eps_curr = eps_0
    delta = 0
    it = 0

    for it in xrange(len(G)):
        lam = np.abs(alpha).argmax()
        if len(idx) > 0:
            w = linalg.solve_triangular(L, G[idx, lam],
                                        lower=True, unit_diagonal=True)
            L = np.r_[np.c_[L, np.zeros(len(L))],
                      np.atleast_2d(np.append(w, np.sqrt(1 - np.inner(w, w))))]
        idx.append(lam)
        Lc = linalg.solve_triangular(L, Xy[idx], lower=True)
        gamma = linalg.solve_triangular(L, Lc, trans=1, lower=True)
        beta = np.dot(G[:, idx], gamma)
        alpha = Xy - beta
        if eps != None:
            eps_curr += delta
            delta = np.inner(gamma, beta[idx])
            eps_curr -= delta
            if eps_curr <= eps:
                break
        elif n_atoms != None and it == n_atoms - 1:
            break
    return gamma, idx


def orthogonal_mp(X, y, n_atoms=None, eps=None, compute_gram=False):
    """Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems.

    Parameters:
    -----------
    X: array of shape (n_samples, n_features)
        Input dictionary

    y: array, shape: n_samples or (n_samples, n_targets)
        Input targets

    n_atoms: int
        Desired number of non-zero entries in the solution

    eps: float
        Maximum norm of the residual. If not None, overrides n_atoms.

    compute_gram, bool:
        Whether to perform precomputations. Improves performance when n_targets
        or n_samples is very large.

    Returns:
    --------
    coef: array of shape: n_features or (n_features, n_targets)
        Coefficients of the OMP solution
    """
    X = np.asanyarray(X)
    y = np.asanyarray(y)
    if y.ndim == 1:
        y = y[:, np.newaxis]
    if n_atoms == None and eps == None:
        raise ValueError('OMP needs either a target number of atoms (n_atoms) \
                         or a target residual error (eps)')
    if eps != None and eps < 0:
        raise ValueError("Epsilon must be positive")
    if eps == None and n_atoms < 0:
        raise ValueError("The number of atoms must be positive")
    if eps == None and n_atoms > X.shape[1]:
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    if compute_gram:
        G = np.dot(X.T, X)
        Xy = np.dot(X.T, y)
        if eps != None:
            norms_squared = np.sum((y ** 2), axis=0)
        else:
            norms_squared = None
        return orthogonal_mp_gram(G, Xy, n_atoms, eps, norms_squared)

    coef = np.zeros((X.shape[1], y.shape[1]))
    for k in xrange(y.shape[1]):
        x, idx = _cholesky_omp(X, y[:, k], n_atoms, eps)
        coef[idx, k] = x
    return np.squeeze(coef)


def orthogonal_mp_gram(G, Xy, n_atoms=None, eps=None, norms_squared=None):
    """Gram Orthogonal Matching Pursuit (OMP)

    Solves n_targets Orthogonal Matching Pursuit problems using only
    the Gram matrix X.T * X and the product X * y.

    Parameters:
    -----------

    G: array of shape (n_features, n_features)
        Gram matrix of the input data: X.T * X

    Xy: array, shape: n_features or (n_features, n_targets)
        Input targets multiplied by X: X * y

    n_atoms: int
        Desired number of non-zero entries in the solution

    eps: float
        Maximum norm of the residual. If not None, overrides n_atoms.

    norms_squared: array, shape: n_targets
        Squared L2 norms of the lines of y. Required if eps is not None.

    Returns:
    --------
    coef: array of shape: n_features or (n_features, n_targets)
        Coefficients of the OMP solution

    """
    G = np.asanyarray(G)
    Xy = np.asanyarray(Xy)
    if Xy.ndim == 1:
        Xy = Xy[:, np.newaxis]

    if n_atoms == None and eps == None:
        raise ValueError('OMP needs either a target number of atoms (n_atoms) \
                         or a target residual error (eps)')
    if eps != None and norms_squared == None:
        raise ValueError('Gram OMP needs the precomputed norms in order \
                          to evaluate the error sum of squares.')
    if eps != None and eps < 0:
        raise ValueError("Epsilon must be positive")
    if eps == None and n_atoms < 0:
        raise ValueError("The number of atoms must be positive")
    if eps == None and n_atoms > len(G):
        raise ValueError("The number of atoms cannot be more than the number \
                          of features")
    coef = np.zeros((len(G), Xy.shape[1]))
    for k in range(Xy.shape[1]):
        x, idx = _gram_omp(G, Xy[:, k], n_atoms,
                           norms_squared[k] if eps else None, eps)
        coef[idx, k] = x
    return np.squeeze(coef)
