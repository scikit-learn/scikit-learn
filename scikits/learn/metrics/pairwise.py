"""Utilities to evaluate pairwise distances or metrics between 2
sets of points.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np


def euclidian_distances(X, Y=None):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vector

    Parameters
    ----------
    X, array of shape (n_samples_1, n_features)

    Y, array of shape (n_samples_2, n_features), default None
            if Y is None, then Y=X is used instead

    Returns
    -------
    distances, array of shape (n_samples_1, n_samples_2)
    """
    X = np.asanyarray(X)
    Y = np.asanyarray(Y)
    if Y is None:
        Y = X
    if X.shape[1] != Y.shape[1]:
        raise ValueError, "incompatible dimension for X and Y matrices"

    XX = np.sum(X * X, axis=1)[:, np.newaxis]
    if Y is None:
        YY = XX.T
    else:
        YY = np.sum(Y * Y, axis=1)[np.newaxis, :]
    distances = XX + YY # Using broadcasting
    distances -= 2 * np.dot(X, Y.T)
    distances = np.maximum(distances, 0)
    distances = np.sqrt(distances)
    return distances
