"""Utilities to evaluate pairwise distances or metrics between 2
sets of points.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np

from ..utils.extmath import safe_sparse_dot

def euclidean_distances(X, Y, Y_norm_squared=None, squared=False):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    Y_norm_squared: array [n_samples_2], optional
        pre-computed (Y**2).sum(axis=1)

    squared: boolean, optional
        This routine will return squared Euclidean distances instead.

    Returns
    -------
    distances: array of shape (n_samples_1, n_samples_2)

    Examples
    --------
    >>> from scikits.learn.metrics.pairwise import euclidean_distances
    >>> X = [[0, 1], [1, 1]]
    >>> # distrance between rows of X
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
        X = Y = np.asanyarray(X)
    else:
        X = np.asanyarray(X)
        Y = np.asanyarray(Y)

    if X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices")

    XX = np.sum(X * X, axis=1)[:, np.newaxis]
    if X is Y: # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    elif Y_norm_squared is None:
        YY = Y.copy()
        YY **= 2
        YY = np.sum(YY, axis=1)[np.newaxis, :]
    else:
        YY = np.asanyarray(Y_norm_squared)
        if YY.shape != (Y.shape[0],):
            raise ValueError("Incompatible dimension for Y and Y_norm_squared")
        YY = YY[np.newaxis, :]

    # TODO:
    # a faster cython implementation would do the dot product first,
    # and then add XX, add YY, and do the clipping of negative values in
    # a single pass over the output matrix.
    distances = XX + YY # Using broadcasting
    distances -= 2 * np.dot(X, Y.T)
    distances = np.maximum(distances, 0)
    if squared:
        return distances
    else:
        return np.sqrt(distances)

euclidian_distances = euclidean_distances # both spelling for backward compat


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


def polynomial_kernel(X, Y, degree=3):
    """
    Compute the polynomial kernel between X and Y.

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    degree: int

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    K = linear_kernel(X, Y)
    K += 1
    K **= degree
    return K


def rbf_kernel(X, Y, sigma=1.0):
    """
    Compute the rbf (gaussian) kernel between X and Y.

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    sigma: float

    Returns
    -------
    Gram matrix: array of shape (n_samples_1, n_samples_2)
    """
    K = -euclidean_distances(X, Y, squared=True)
    K /= (2 * (sigma ** 2))
    np.exp(K, K)
    return K


def cosine_similarity(X, Y):
    """Compute pairwise cosine similarities between rows in X and Y

    Cosine similarity is a normalized linear kernel with value ranging
      - -1: similar vectors with opposite signs
      - 0: completely dissimilar (orthogonal) vectors
      - 1: similar vectors (same sign)

    In practice, cosine similarity is often used to measure the
    relatedness of text documents represented by sparse vectors of word
    counts, frequencies or TF-IDF weights. In this cases all features
    are non negative and the similarities range from 0 to 1 instead.

    Cosine similarity can be used as an affinity matrix for spectral
    and power iteration clustering algorithms.

    Parameters
    ----------

    X: array or sparse matrix of shape (n_samples_1, n_features)

    Y: array or sparse matrix of shape (n_samples_2, n_features)

    Returns
    -------

    array or sparse matrix of shape (n_samples_1, n_samples_2)

    Examples
    --------

    >>> from scikits.learn.metrics.pairwise import cosine_similarity
    >>> X = np.asarray([[0, 1], [1, 1], [0, -1], [0, 0]], dtype=np.float64)
    >>> cosine_similarity(X, X).round(decimals=2)
    array([[ 1.  ,  0.71, -1.  ,  0.  ],
           [ 0.71,  1.  , -0.71,  0.  ],
           [-1.  , -0.71,  1.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ]])

    >>> from scipy.sparse import csr_matrix
    >>> X_sparse = csr_matrix(X)
    >>> cosine_similarity(X_sparse, X_sparse).toarray().round(decimals=2)

    It is possible to use the cosine similarity to perform similarity
    queries:

    >>> query = [[0.5, 0.9]]
    >>> cosine_similarity(X, query)

    """
    similarities = safe_sparse_dot(X, Y.T)

    if hasattr(X, 'multiply'):
        norms_x = np.sqrt((X.multiply(X)).sum(axis=1))
    else:
        norms_x = np.sqrt((X ** 2).sum(axis=1))
    nnzeros_x = np.where(norms_x > 0)
    similarities[nnzeros_x, :] /= norms_x[nnzeros_x].reshape(
        (nnzeros_x[0].sum(), -1))

    if hasattr(Y, 'multiply'):
        norms_y = np.sqrt((Y.multiply(Y)).sum(axis=1))
    else:
        norms_y = np.sqrt((Y ** 2).sum(axis=1))
    nnzeros_y = np.where(norms_y > 0)
    similarities[:, nnzeros_y,] /= norms_y[nnzeros_y].reshape(
        (-1, nnzeros_y[0].sum()))

    return similarities

