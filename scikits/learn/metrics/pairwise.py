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
def return_self_if_square(X):
    """ Returns the argument iff it is a square numpy array

    If the argument is a square numpy array, it is returned.
    Otherwise, raise a ValueError instead.

    Examples
    --------
    >>> from scikits.learn.metrics.pairwise import return_self_if_square
    >>> import numpy as np
    >>> X = np.array([[0., 1.], [1., 1.]])
    >>> return_self_if_square(X)
    array([[ 0.,  1.],
           [ 1.,  1.]])
    >>> X = np.array([[0, 1], [1, 1], [2, 2]])
    >>> return_self_if_square(X)
    Traceback (most recent call last):
        ...
    ValueError: X is not square!
    """
    if X.shape[0] != X.shape[1]:
            raise ValueError("X is not square!")
    return X

def check_set_Y(X, Y):
    """ Sets Y appropriately and checks inputs

    If Y is None, it is set as a pointer to X (i.e. not a copy).
    If Y is given, this does not happen.
    All distance metrics should use this function first to assert that the
    given parameters are correct and safe to use.

    Parameters
    ----------
    X: {array-like, sparse matrix}, shape = [n_samples_1, n_features]

    Y: {array-like, sparse matrix}, shape = [n_samples_2, n_features]

    Returns
    -------
    safe_X: {array-like, sparse matrix}, shape = [n_samples_1, n_features]
        An array equal to X, guarenteed to be a numpy array.
        
    safe_Y: {array-like, sparse matrix}, shape = [n_samples_2, n_features]
        An array equal to Y if Y was not None, guarenteed to be a numpy array.
        If Y was None, safe_Y will be a pointer to X.
        
    """
    if Y is X or Y is None:
        X = Y = safe_asanyarray(X)
    else:
        X = safe_asanyarray(X)
        Y = safe_asanyarray(Y)
    X, Y = np.atleast_2d(X), np.atleast_2d(Y)
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
    X, Y = check_set_Y(X, Y)
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
    X, Y = check_set_Y(X, Y)
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
    X, Y = check_set_Y(X, Y)
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
    X, Y = check_set_Y(X, Y)
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
    X, Y = check_set_Y(X, Y)
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
    X, Y = check_set_Y(X, Y)
    if gamma == 0:
        gamma = 1.0 / X.shape[1]

    K = euclidean_distances(X, Y, squared=True)
    K *= -gamma
    np.exp(K, K)    # exponentiate K in-place
    return K

