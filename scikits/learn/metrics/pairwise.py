"""Utilities to evaluate pairwise distances or affinity of sets of samples"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np
from scipy.spatial import distance
from scipy.sparse import csr_matrix, issparse
from ..utils import safe_asanyarray, atleast2d_or_csr, deprecated
from ..utils.extmath import safe_sparse_dot


# Utility Functions
def check_pairwise_arrays(XA, XB):
    """ Sets XA and XB appropriately and checks inputs

    If XB is None, it is set as a pointer to XA (i.e. not a copy).
    If XB is given, this does not happen.
    All distance metrics should use this function first to assert that the
    given parameters are correct and safe to use.

    Specifically, this function first ensures that both XA and XB are arrays,
    then checkes that they are at least two dimensional. Finally, the function
    checks that the size of the second dimension of the two arrays is equal.

    Parameters
    ----------
    XA: {array-like, sparse matrix}, shape = [n_samples_a, n_features]

    XB: {array-like, sparse matrix}, shape = [n_samples_b, n_features]

    Returns
    -------
    safe_XA: {array-like, sparse matrix}, shape = [n_samples_a, n_features]
        An array equal to XA, guarenteed to be a numpy array.

    safe_XB: {array-like, sparse matrix}, shape = [n_samples_b, n_features]
        An array equal to XB if XB was not None, guarenteed to be a numpy array.
        If XB was None, safe_XB will be a pointer to XA.

    """
    if XB is XA or XB is None:
        XA = XB = safe_asanyarray(XA)
    else:
        XA = safe_asanyarray(XA)
        XB = safe_asanyarray(XB)
    if len(XA.shape) < 2:
        raise ValueError("XA is required to be at least two dimensional.")
    if len(XB.shape) < 2:
        raise ValueError("XB is required to be at least two dimensional.")
    if XA.shape[1] != XB.shape[1]:
        raise ValueError("Incompatible dimension for XA and XB matrices")
    return XA, XB


# Distances
def euclidean_distances(XA, XB=None, Y_norm_squared=None, squared=False):
    """
    Considering the rows of XA (and XB=XA) as vectors, compute the
    distance matrix between each pair of vectors.

    Parameters
    ----------
    XA: {array-like, sparse matrix}, shape = [n_samples_1, n_features]

    XB: {array-like, sparse matrix}, shape = [n_samples_2, n_features]

    Y_norm_squared: array-like, shape = [n_samples_2], optional
        Pre-computed (XB**2).sum(axis=1)

    squared: boolean, optional
        Return squared Euclidean distances.

    Returns
    -------
    distances: {array, sparse matrix}, shape = [n_samples_1, n_samples_2]

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
    XA, XB = check_pairwise_arrays(XA, XB)
    if issparse(XA):
        XX = XA.multiply(XA).sum(axis=1)
    else:
        XX = np.sum(XA * XA, axis=1)[:, np.newaxis]

    if XA is XB:  # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    elif Y_norm_squared is None:
        if issparse(XB):
            # scipy.sparse matrices don't have element-wise scalar
            # exponentiation, and tocsr has a copy kwarg only on CSR matrices.
            YY = XB.copy() if isinstance(XB, csr_matrix) else XB.tocsr()
            YY.data **= 2
            YY = np.asarray(YY.sum(axis=1)).T
        else:
            YY = np.sum(XB ** 2, axis=1)[np.newaxis, :]
    else:
        YY = atleast2d_or_csr(Y_norm_squared)
        if YY.shape != (1, XB.shape[0]):
            raise ValueError(
                        "Incompatible dimensions for XB and Y_norm_squared")

    # TODO: a faster Cython implementation would do the clipping of negative
    # values in a single pass over the output matrix.
    distances = safe_sparse_dot(XA, XB.T, dense_output=True)
    distances *= -2
    distances += XX
    distances += YY
    distances = np.maximum(distances, 0)
    return distances if squared else np.sqrt(distances)


@deprecated("use euclidean_distances instead")
def euclidian_distances(*args, **kwargs):
    return euclidean_distances(*args, **kwargs)


def l1_distances(X, Y=None):
    """
    Computes the componentwise L1 pairwise-distances between the vectors
    in X and Y.

    Parameters
    ----------
    X: array_like
        An array with shape (n_samples_X, n_features)

    Y: array_like, optional
        An array with shape (n_samples_Y, n_features).

    Returns
    -------
    D: array with shape (n_samples_X * n_samples_Y, n_features)
        The array of componentwise L1 pairwise-distances.

    Examples
    --------
    >>> from scikits.learn.metrics.pairwise import l1_distances
    >>> l1_distances(3, 3)
    array([[0]])
    >>> l1_distances(3, 2)
    array([[1]])
    >>> l1_distances(2, 3)
    array([[1]])
    >>> import numpy as np
    >>> X = np.ones((1, 2))
    >>> y = 2*np.ones((2, 2))
    >>> l1_distances(X, y)
    array([[ 1.,  1.],
           [ 1.,  1.]])
    """
    X, Y = check_pairwise_arrays(X, Y)
    n_samples_X, n_features_X = X.shape
    n_samples_Y, n_features_Y = Y.shape
    if n_features_X != n_features_Y:
        raise Exception("X and Y should have the same number of features!")
    else:
        n_features = n_features_X
    D = np.abs(X[:, np.newaxis, :] - Y[np.newaxis, :, :])
    D = D.reshape((n_samples_X * n_samples_Y, n_features))

    return D


# Kernels
def linear_kernel(X, Y=None):
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
    X, Y = check_pairwise_arrays(X, Y)
    return safe_sparse_dot(X, Y.T, dense_output=True)


def polynomial_kernel(X, Y=None, degree=3, gamma=0, coef0=1):
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
    X, Y = check_pairwise_arrays(X, Y)
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = linear_kernel(X, Y)
    K *= gamma
    K += coef0
    K **= degree
    return K


def sigmoid_kernel(X, Y=None, gamma=0, coef0=1):
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
    X, Y = check_pairwise_arrays(X, Y)
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = linear_kernel(X, Y)
    K *= gamma
    K += coef0
    np.tanh(K, K)   # compute tanh in-place
    return K


def rbf_kernel(X, Y=None, gamma=0):
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
    X, Y = check_pairwise_arrays(X, Y)
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = euclidean_distances(X, Y, squared=True)
    K *= -gamma
    np.exp(K, K)    # exponentiate K in-place
    return K


# Helper functions - distance
pairwise_function_map = {}
pairwise_function_map['euclidean'] = euclidean_distances
pairwise_function_map['l1'] = l1_distances


def pairwise_distances(XA, XB=None, metric="euclidean", **kwds):
    """ Calculates the distance matrix from a vector matrix XA and optional XB.

    This method takes either a vector array or a distance matrix, and returns
    a distance matrix. If the input is a vector array, the distances are
    computed. If the input is a distances matrix, it is returned instead.

    This method provides a safe way to take a distance matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    If XB is given (default is None), then the returned matrix is the pairwise
    distance between the arrays from both XA and XB.

    Parameters
    ----------
    XA: array [n_samples_a, n_samples_a] if metric == "precomputed", or,
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    XB: array [n_samples_b, n_features]
        A second feature array only if XA has shape [n_samples_a, n_features].

    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        If metric is "precomputed", XA is assumed to be a distance matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from XA as input and return a value indicating
        the distance between them.

    **kwds: optional keyword parameters
        Any further parameters are passed directly to the distance metric.

    Returns
    -------
    D: array [n_samples_a, n_samples_a] or [n_samples_a, n_samples_b]
        A distance matrix D such that D_{i, j} is the distance between the
        ith and jth vectors of the given matrix XA, if XB is None.
        If XB is not None, then D_{i, j} is the distance between the ith array
        from XA and the jth array from XB.

    """
    if metric == "precomputed":
        if XA.shape[0] != XA.shape[1]:
            raise ValueError("X is not square!")
        return XA
    elif metric in pairwise_function_map:
        return pairwise_function_map[metric](XA, XB, **kwds)
    else:
        # FIXME: the distance module doesn't support sparse matrices!
        if XB is None:
            return distance.squareform(distance.pdist(XA, metric=metric),
                                       **kwds)
        else:
            return distance.cdist(XA, XB, metric=metric, **kwds)
