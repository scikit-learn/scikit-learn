"""Utilities to evaluate pairwise distances or metrics between 2
sets of points.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np

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
