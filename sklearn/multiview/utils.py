""" Utilities to be used in this package.

There are mainly four functions to help other classes:
Hbeta, x2p and whiten for the tsne algorithm."""

___version___ = '1.0'
___author___ = 'Maria Araceli Burgueño Caballero'
___email___ = "mburgueno@uoc.edu"
___status___ = "Pre-Production"


import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils import check_array


def Hbeta(D, beta):
    """Compute H beta for matrix D, given beta.

    Parameters
    ----------

    D : ndarray or nscalar.
        Input value.
    beta : numeric.
        Number to operate with in the exponential distribution.

    Returns
    -------

    (H, P) : tuple
        Return a tuple with H and P values. H is a float value and P
        can be either an scalar or a ndarray.

    Examples
    --------

    >>> import numpy as np
    >>> matrix = np.array([[1, 2, 3], [2, 3, 4], [5, 6, 2]])
    >>> hbeta = Hbeta(matrix, 3)
    >>> print("H -> %g\nP-> %a" % hbeta)
        H -> 21.6814
        P-> array([[  8.66214422e-01,   4.31262766e-02,   2.14713088e-03],
                   [  4.31262766e-02,   2.14713088e-03,   1.06899352e-04],
                   [  5.32220535e-06,   2.64977002e-07,   4.31262766e-02]])
    """
    P = np.exp(-D * beta)
    sumP = np.sum(P)
    if sumP == 0:
        H = 0
        P = D * 0
    else:
        H = np.log(sumP) + beta * np.sum(np.dot(D, P)) / sumP
        P /= sumP
    return (H, P)


def x2p(X, perplexity=15, tol=1e-5):
    """ Compute a pair-wise conditional matrix given input matrix.

    Parameters
    ----------
    X : ndarray.
        Input data. Matrix with n rows.
    perplexity : numeric.
        Target perplexity. Loosely translates into the expected number of
        behibours per point.
    tol : float
        Tolerance to be used in the computation. It must be a small value.

    Returns
    -------
    (P, beta) : tuple
        P: Pair-wise conditional probability matrix.
        beta: array used during computation.

    Examples
    --------

    >>> x = np.array(([2, 4, 3], [1, 5, 7], [8, 6, 9])).T
    >>> x2p(x)
        (array([[ 0. ,  0.5,  0.5],
       [ 0.5,  0. ,  0.5],
       [ 0.5,  0.5,  0. ]]), array([  8.88178420e-16,
                                      8.88178420e-16,   8.88178420e-16]))
    """
    X = check_array(X)
    n = X.shape[0]
    D = euclidean_distances(X)
    P = np.zeros((n, n))
    beta = np.array([1] * n, dtype="float64")
    logU = np.log(perplexity)
    for i in np.arange(n):
        betamin = -float("inf")
        betamax = float("inf")
        Di = np.delete(D[i], i)
        H, this_P = Hbeta(Di, beta[i])
        H_diff = H - logU
        tries = 0
        while abs(H_diff) > tol and tries < 50:
            if H_diff > 0:
                betamin = beta[i]
                if betamax == float("inf"):
                    beta[i] = beta[i] * 2
                else:
                    beta[i] = (beta[i] + betamax) / 2
            else:
                betamax = beta[i]
                if betamin == -float("inf"):
                    beta[i] = beta[i] / 2
                else:
                    beta[i] = (beta[i] + betamin) / 2
            H, this_P = Hbeta(Di, beta[i])
            H_diff = H - logU
            tries += 1
        P[i, np.arange(n) != i] = this_P
    return (P, beta)


def whiten(X, row_norm=False, verbose=0, n_comp=-1):
    """Whitening of matrix X, and return that new matrix whitened.

    Parameters
    ----------

    X : ndarray
        Input data (2D).
    row_norm : int, default 0
        If row_norm is True, then input data is scaled as well as centered.
        If it is False, then no center or scale is computed.
    verbose : int, default 0
        Verbosity mode.
    n_comp : int, default -1
        Number of rows of output matrix. If n_comp is -1, the number
        of columns of input data is set.

    Returns
    -------

    X : ndarray
        Whitened matrix.

    Examples
    --------

    >>> x = np.array(([2, 4, 3], [1, 5, 7], [8, 6, 9])).T
    >>> whiten(x, row_norm=True)
        [[  1.22474487e+00   7.07106781e-01  -7.26265414e-08]
         [ -1.22474487e+00   7.07106781e-01  -3.69804654e-08]
         [  2.10705234e-16  -1.41421356e+00  -5.30203132e-09]]
    >>> y = np.array(([1, 2, 3], [4, 5, 6], [7, 8, 9])).T
    >>> whiten(y)
        [[  1.22474487e+00  -2.98023224e-08   0.00000000e+00]
         [  0.00000000e+00   0.00000000e+00   0.00000000e+00]
         [ -1.22474487e+00   2.98023224e-08   0.00000000e+00]]
    """
    if n_comp == -1:
        n_comp = X.shape[1]
    if verbose:
        print("Centering")
    n = X.shape[0]
    p = X.shape[1]
    # Centering matrix columns
    mean = X.mean(axis=0)
    sd = np.std(X, axis=0, ddof=1)
    X -= mean
    if row_norm:
        X = (X / sd).T
    else:
        X = X.T
    if verbose:
        print("Whitening")
    V = np.dot(X, X.T) / n
    u, s, v = np.linalg.svd(V)
    D = (np.diag(1 / s**(1 / 2)))
    K = np.dot(D, u.T)
    K = np.matrix(K[:n_comp, :].reshape((n_comp, p)))
    X = np.dot(K, X).T
    return X
