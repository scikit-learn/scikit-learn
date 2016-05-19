import numpy as np
from ..utils.extmath import randomized_svd
from numpy.linalg import norm


def _soft_thresh(X, threshold):
    "Apply soft thresholding to an array"

    sign = np.sign(X)
    return np.multiply(sign, np.maximum(np.abs(X) - threshold, 0))


def norm_1(X):
    return np.sum(np.abs(X))


def _sv_thresh(X, threshold, num_svalue):
    """
    Perform singular value thresholding.

    Parameters
    ---------
    X : array of shape [n_samples, n_features]
        The input array.
    threshold : float
        The threshold for the singualar values.
    num_svalue : int
        The number of singular values to compute.

    Returns
    -------
    X_thresh : array of shape [n_samples, n_features]
        The output after performing singular value thresholding.
    grater_sv : int
        The number of singular values of `X` which were greater than
        `threshold`
    (U, s, V): tuple
        The singular value decomposition
    """
    m, n = X.shape
    U, s, V = randomized_svd(X, num_svalue)
    greater_sv = np.count_nonzero(s > threshold)
    s = _soft_thresh(s, threshold)
    S = np.diag(s)
    X_thresh = np.dot(U, np.dot(S, V))
    return X_thresh, greater_sv, (U, s, V)


def rpca(M, lam=None, mu=None, max_iter=1000, eps_primal=1e-7, eps_dual=1e-5,
         rho=1.6, initial_sv=10, max_mu=1e6, verbose=False):
    # See http://arxiv.org/pdf/1009.5055v3.pdf

    # This implementation follows Algorithm 5 from the paper with minor
    # modifications

    if lam is None:
        lam = 1.0/np.sqrt(max(M.shape))

    d = min(M.shape)

    # See "Choosing Parameters" paragraph in section 4
    mu = 1.25/norm(M, 2)

    # The sparse matrix
    S = np.zeros_like(M)

    # The low rank matrix
    L = np.zeros_like(M)

    # See equation 10
    J = min(norm(M, 2), np.max(np.abs(M)))
    Y = M/J

    M_fro_norm = norm(M, 'fro')

    # This variable tried to predict how many singular values will be required.
    sv = initial_sv

    for iter_ in range(max_iter):
        # See Section 4, paragraph "Order of Updating A and E" to see why
        # `S` iterate is computed before `L` ierate.
        S_old = S
        S = _soft_thresh(M - L + (Y/mu), lam/mu)
        L, svp, (U, s, V) = _sv_thresh(M - S + (Y/mu), 1/mu, sv)
        Y = Y + mu*(M - L - S)

        mu_old = mu
        mu = rho*mu
        mu = min(mu, max_mu)

        # See Equation 18
        if svp < sv:
            sv = svp + 1
        else:
            sv = svp + int(round(0.05*d))

        sv = min(sv, M.shape[0], M.shape[1])

        primal_error = norm(M - L - S, 'fro')/M_fro_norm
        dual_error = mu_old*norm(S - S_old, 'fro')/M_fro_norm

        if verbose:
            print('rpca: Iteration %d - Primal Error = %e Dual Error = %e' %
                  (iter_, primal_error, dual_error))
        if primal_error < eps_primal and dual_error < eps_dual:
            break

    if iter_ >= max_iter:
        if verbose:
            print('rpca: Failed to converge within %d iterations' % max_iter)

    return L, S, (U, s, V)
