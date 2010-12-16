"""Utilities to evaluate pairwise distances or metrics between 2
sets of points.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np


def euclidian_distances(X, Y):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    Parameters
    ----------
    X: array of shape (n_samples_1, n_features)

    Y: array of shape (n_samples_2, n_features)

    Returns
    -------
    distances: array of shape (n_samples_1, n_samples_2)

    Examples
    --------
    >>> X = [[0, 1], [1, 1]]
    >>> # distrance between rows of X
    >>> euclidian_distances(X, X)
    array([[ 0.,  1.],
           [ 1.,  0.]])
    >>> # get distance to origin
    >>> euclidian_distances(X, [[0, 0]])
    array([[ 1.        ],
           [ 1.41421356]])
     """
    # shortcut in the common case euclidean_distances(X, X)
    compute_Y = X is not Y

    X = np.asanyarray(X)
    Y = np.asanyarray(Y)

    if X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices")

    XX = np.sum(X * X, axis=1)[:, np.newaxis]
    if compute_Y:
        YY = np.sum(Y * Y, axis=1)[np.newaxis, :]
    else:
        YY = XX.T

    distances = XX + YY # Using broadcasting
    distances -= 2 * np.dot(X, Y.T)
    distances = np.maximum(distances, 0)
    return np.sqrt(distances)
