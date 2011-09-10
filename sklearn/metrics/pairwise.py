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
def check_pairwise_arrays(X, Y):
    """ Set X and Y appropriately and checks inputs

    If Y is None, it is set as a pointer to X (i.e. not a copy).
    If Y is given, this does not happen.
    All distance metrics should use this function first to assert that the
    given parameters are correct and safe to use.

    Specifically, this function first ensures that both X and Y are arrays,
    then checkes that they are at least two dimensional. Finally, the function
    checks that the size of the second dimension of the two arrays is equal.

    Parameters
    ----------
    X: {array-like, sparse matrix}, shape = [n_samples_a, n_features]

    Y: {array-like, sparse matrix}, shape = [n_samples_b, n_features]

    Returns
    -------
    safe_X: {array-like, sparse matrix}, shape = [n_samples_a, n_features]
        An array equal to X, guarenteed to be a numpy array.

    safe_Y: {array-like, sparse matrix}, shape = [n_samples_b, n_features]
        An array equal to Y if Y was not None, guarenteed to be a numpy array.
        If Y was None, safe_Y will be a pointer to X.

    """
    if Y is X or Y is None:
        X = Y = safe_asanyarray(X)
    else:
        X = safe_asanyarray(X)
        Y = safe_asanyarray(Y)
    X = atleast2d_or_csr(X)
    Y = atleast2d_or_csr(Y)
    if len(X.shape) < 2:
        raise ValueError("X is required to be at least two dimensional.")
    if len(Y.shape) < 2:
        raise ValueError("Y is required to be at least two dimensional.")
    if X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices")
    return X, Y


# Distances
def euclidean_distances(X, Y=None, Y_norm_squared=None, squared=False):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    Parameters
    ----------
    X: {array-like, sparse matrix}, shape = [n_samples_1, n_features]

    Y: {array-like, sparse matrix}, shape = [n_samples_2, n_features]

    Y_norm_squared: array-like, shape = [n_samples_2], optional
        Pre-computed (Y**2).sum(axis=1)

    squared: boolean, optional
        Return squared Euclidean distances.

    Returns
    -------
    distances: {array, sparse matrix}, shape = [n_samples_1, n_samples_2]

    Examples
    --------
    >>> from sklearn.metrics.pairwise import euclidean_distances
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
    X, Y = check_pairwise_arrays(X, Y)
    if issparse(X):
        XX = X.multiply(X).sum(axis=1)
    else:
        XX = np.sum(X * X, axis=1)[:, np.newaxis]

    if X is Y:  # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    elif Y_norm_squared is None:
        if issparse(Y):
            # scipy.sparse matrices don't have element-wise scalar
            # exponentiation, and tocsr has a copy kwarg only on CSR matrices.
            YY = Y.copy() if isinstance(Y, csr_matrix) else Y.tocsr()
            YY.data **= 2
            YY = np.asarray(YY.sum(axis=1)).T
        else:
            YY = np.sum(Y ** 2, axis=1)[np.newaxis, :]
    else:
        YY = atleast2d_or_csr(Y_norm_squared)
        if YY.shape != (1, Y.shape[0]):
            raise ValueError(
                        "Incompatible dimensions for Y and Y_norm_squared")

    # TODO: a faster Cython implementation would do the clipping of negative
    # values in a single pass over the output matrix.
    distances = safe_sparse_dot(X, Y.T, dense_output=True)
    distances *= -2
    distances += XX
    distances += YY
    distances = np.maximum(distances, 0)
    return distances if squared else np.sqrt(distances)


@deprecated("use euclidean_distances instead")
def euclidian_distances(*args, **kwargs):
    return euclidean_distances(*args, **kwargs)


def manhattan_distances(X, Y=None, sum_over_features=True):
    """ Compute the L1 distances between the vectors in X and Y.

    With sum_over_features equal to False it returns the componentwise
    distances.

    Parameters
    ----------
    X: array_like
        An array with shape (n_samples_X, n_features).

    Y: array_like, optional
        An array with shape (n_samples_Y, n_features).

    sum_over_features: bool, default=True
        If True the function returns the pairwise distance matrix
        else it returns the componentwise L1 pairwise-distances.

    Returns
    -------
    D: array
        If sum_over_features is False shape is
        (n_samples_X, n_samples_Y, n_features) and D contains the
        componentwise L1 pairwise-distances (ie. absolute difference),
        else shape is (n_samples_X, n_samples_Y) and D contains
        the pairwise l1 distances.

    Examples
    --------
    >>> from sklearn.metrics.pairwise import manhattan_distances
    >>> manhattan_distances(3, 3)
    array([[0]])
    >>> manhattan_distances(3, 2)
    array([[1]])
    >>> manhattan_distances(2, 3)
    array([[1]])
    >>> import numpy as np
    >>> X = np.ones((1, 2))
    >>> y = 2 * np.ones((2, 2))
    >>> manhattan_distances(X, y, sum_over_features=False)
    array([[[ 1.,  1.],
            [ 1.,  1.]]])
    """
    X, Y = check_pairwise_arrays(X, Y)
    n_samples_X, n_features_X = X.shape
    n_samples_Y, n_features_Y = Y.shape
    if n_features_X != n_features_Y:
        raise Exception("X and Y should have the same number of features!")
    D = np.abs(X[:, np.newaxis, :] - Y[np.newaxis, :, :])
    if sum_over_features:
        D = np.sum(D, axis=2)
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
pairwise_distance_functions = {
    'euclidean': euclidean_distances,
    'l2': euclidean_distances,
    'l1': manhattan_distances,
    'manhattan': manhattan_distances,
    'cityblock': manhattan_distances
    }


def pairwise_distances(X, Y=None, metric="euclidean", **kwds):
    """ Compute the distance matrix from a vector array X and optional Y.

    This method takes either a vector array or a distance matrix, and returns
    a distance matrix. If the input is a vector array, the distances are
    computed. If the input is a distances matrix, it is returned instead.

    This method provides a safe way to take a distance matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    If Y is given (default is None), then the returned matrix is the pairwise
    distance between the arrays from both X and Y.

    Please note that support for sparse matrices is currently limited to those
    metrics listed in pairwise.pairwise_distance_functions.

    Valid values for metric are:
    - from scikits.learn: ['euclidean', 'l2', 'l1', 'manhattan', 'cityblock']
    - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
      'correlation', 'cosine', 'dice', 'hamming', 'jaccard', 'kulsinski',
      'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
      'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeucludean', 'yule']
      See the documentation for scipy.spatial.distance for details on these
      metrics.
    Note in the case of 'euclidean' and 'cityblock' (which are valid
    scipy.spatial.distance metrics), the values will use the scikits.learn
    implementation, which is faster and has support for sparse matrices.

    Parameters
    ----------
    X: array [n_samples_a, n_samples_a] if metric == "precomputed", or,
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    Y: array [n_samples_b, n_features]
        A second feature array only if X has shape [n_samples_a, n_features].

    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter, or
        a metric listed in pairwise.pairwise_distance_functions.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    **kwds: optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a scipy.spatial.distance metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    D: array [n_samples_a, n_samples_a] or [n_samples_a, n_samples_b]
        A distance matrix D such that D_{i, j} is the distance between the
        ith and jth vectors of the given matrix X, if Y is None.
        If Y is not None, then D_{i, j} is the distance between the ith array
        from X and the jth array from Y.

    """
    if metric == "precomputed":
        if X.shape[0] != X.shape[1]:
            raise ValueError("X is not square!")
        return X
    elif metric in pairwise_distance_functions:
        return pairwise_distance_functions[metric](X, Y, **kwds)
    elif callable(metric):
        # Check matrices first (this is usually done by the metric).
        X, Y = check_pairwise_arrays(X, Y)
        n_x, n_y = X.shape[0], Y.shape[0]
        # Calculate distance for each element in X and Y.
        D = np.zeros((n_x, n_y), dtype='float')
        for i in range(n_x):
            start = 0
            if X is Y:
                start = i
            for j in range(start, n_y):
                # Kernel assumed to be symmetric.
                D[i][j] = metric(X[i], Y[j], **kwds)
                if X is Y:
                    D[j][i] = D[i][j]
        return D
    else:
        # Note that if the metric is callable
        # FIXME: the distance module doesn't support sparse matrices!
        if Y is None:
            return distance.squareform(distance.pdist(X, metric=metric,
                                                      **kwds))
        else:
            return distance.cdist(X, Y, metric=metric, **kwds)


# Helper functions - distance
pairwise_kernel_functions = {
    'rbf': rbf_kernel,
    'sigmoid': sigmoid_kernel,
    'polynomial': polynomial_kernel,
    'poly': polynomial_kernel,
    'linear': linear_kernel
    }


def pairwise_kernels(X, Y=None, metric="linear", **kwds):
    """ Compute the kernel between arrays X and optional array Y.

    This method takes either a vector array or a kernel matrix, and returns
    a kernel matrix. If the input is a vector array, the kernels are
    computed. If the input is a kernel matrix, it is returned instead.

    This method provides a safe way to take a kernel matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    If Y is given (default is None), then the returned matrix is the pairwise
    kernel between the arrays from both X and Y.

    Valid values for metric are:
    ['rbf', 'sigmoid', 'polynomial', 'poly', 'linear']

    Parameters
    ----------
    X: array [n_samples_a, n_samples_a] if metric == "precomputed", or,
             [n_samples_a, n_features] otherwise
        Array of pairwise kernels between samples, or a feature array.

    Y: array [n_samples_b, n_features]
        A second feature array only if X has shape [n_samples_a, n_features].

    metric: string, or callable
        The metric to use when calculating kernel between instances in a
        feature array. If metric is a string, it must be one of the metrics
        in pairwise.pairwise_kernel_functions.
        If metric is "precomputed", X is assumed to be a kernel matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    **kwds: optional keyword parameters
        Any further parameters are passed directly to the kernel function.

    Returns
    -------
    K: array [n_samples_a, n_samples_a] or [n_samples_a, n_samples_b]
        A kernel matrix K such that K_{i, j} is the kernel between the
        ith and jth vectors of the given matrix X, if Y is None.
        If Y is not None, then K_{i, j} is the kernel between the ith array
        from X and the jth array from Y.

    """
    if metric == "precomputed":
        if X.shape[0] != X.shape[1]:
            raise ValueError("X is not square!")
        return X
    elif metric in pairwise_kernel_functions:
        return pairwise_kernel_functions[metric](X, Y, **kwds)
    elif callable(metric):
        # Check matrices first (this is usually done by the metric).
        X, Y = check_pairwise_arrays(X, Y)
        n_x, n_y = X.shape[0], Y.shape[0]
        # Calculate kernel for each element in X and Y.
        K = np.zeros((n_x, n_y), dtype='float')
        for i in range(n_x):
            start = 0
            if X is Y:
                start = i
            for j in range(start, n_y):
                # Kernel assumed to be symmetric.
                K[i][j] = metric(X[i], Y[j], **kwds)
                if X is Y:
                    K[j][i] = K[i][j]
        return K
    else:
        raise AttributeError("Unknown metric %s" % metric)
