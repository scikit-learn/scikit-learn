"""Utilities to evaluate pairwise distances or affinity of sets of samples"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np
from scipy.sparse import issparse
from ..utils import safe_asanyarray, atleast2d_or_csr
from ..utils.extmath import safe_sparse_dot


def euclidean_distances(X, Y, Y_norm_squared=None, squared=False):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    Parameters
    ----------
    X: array-like, shape = [n_samples_1, n_features]

    Y: array-like, shape = [n_samples_2, n_features]

    Y_norm_squared: array-like, shape = [n_samples_2], optional
        Pre-computed (Y**2).sum(axis=1)

    squared: boolean, optional
        Return squared Euclidean distances.

    Returns
    -------
    distances: array, shape = [n_samples_1, n_samples_2]

    Examples
    --------
    >>> from scikits.learn.metrics.pairwise import euclidean_distances
    >>> X = [[0, 1], [1, 1]]
    >>> # distance between rows of X
    >>> euclidean_distances(X, X)
    array([[ 0.,  1.],
           [ 1.,  0.]])
    >>> # get distance to origin
    >>> euclidean_distances(X, [[0, 0]])
    array([[ 1.        ],
           [ 1.41421356]])
    """
    # should not need X_norm_squared because if you could precompute that as
    # well as Y, then you should just pre-compute the output and not even
    # call this function.
    if X is Y:
        X = Y = safe_asanyarray(X)
    else:
        X = safe_asanyarray(X)
        Y = safe_asanyarray(Y)

    if X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices")

    if issparse(X):
        XX = X.multiply(X).sum(axis=1)
    else:
        XX = np.sum(X * X, axis=1)[:, np.newaxis]

    if X is Y:  # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    elif Y_norm_squared is None:
        if issparse(Y):
            YY = Y.copy()
            YY.data **= 2
            YY = YY.sum(axis=1).T
        else:
            YY = np.sum(Y ** 2, axis=1)[np.newaxis, :]
    else:
        YY = atleast2d_or_csr(Y_norm_squared)
        if YY.shape != (1, Y.shape[0]):
            raise ValueError(
                        "Incompatible dimensions for Y and Y_norm_squared")

    # TODO:
    # a faster cython implementation would do the dot product first,
    # and then add XX, add YY, and do the clipping of negative values in
    # a single pass over the output matrix.
    distances = XX + YY  # Using broadcasting
    distances -= 2 * safe_sparse_dot(X, Y.T)
    distances = np.maximum(distances, 0)
    return distances if squared else np.sqrt(distances)

euclidian_distances = euclidean_distances  # both spelling for backward compat


def linear_kernel(X, Y):
    """
    Compute the linear kernel between X and Y.

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    return np.dot(X, Y.T)


def polynomial_kernel(X, Y, degree=3, gamma=0, coef0=1):
    """
    Compute the polynomial kernel between X and Y.

    K(X, Y) = (gamma <X, Y> + coef0)^degree

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    degree: int

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = linear_kernel(X, Y)
    K *= gamma
    K += coef0
    K **= degree
    return K


def sigmoid_kernel(X, Y, gamma=0, coef0=1):
    """
    Compute the sigmoid kernel between X and Y.

    K(X, Y) = tanh(gamma <X, Y> + coef0)

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    degree: int

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = linear_kernel(X, Y)
    K *= gamma
    K += coef0
    np.tanh(K, K)   # compute tanh in-place
    return K


def rbf_kernel(X, Y, gamma=0):
    """
    Compute the rbf (gaussian) kernel between X and Y.

    K(X, Y) = exp(-gamma ||X-Y||^2)

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    gamma: float

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = euclidean_distances(X, Y, squared=True)
    K *= -gamma
    np.exp(K, K)    # exponentiate K in-place
    return K
