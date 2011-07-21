"""
Clustering algorithms
"""
import numpy as np
from scipy.spatial import distance

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .hierarchical import ward_tree, Ward, WardAgglomeration
from .k_means_ import k_means, KMeans, MiniBatchKMeans


def calculate_similarity(X, metric="euclidean"):
    """ Calculates the similarity matrix from a vector matrix X.

    This method takes either a vector array or a similarity matrix, and returns
    a similarity matrix. If the input is a vector array, the similarities are
    computed. If the input is a similarity matrix, it is returned instead.

    This method provides a safe way to take a similarity matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    Parameters
    ----------
    X: array [n_points, n_points] if metric == "precomputed", or,
             [n_points, n_features] otherwise
        Array of similarities between points, or a feature array.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        If metric is "precomputed", X is assumed to be a similarity matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.

    Returns
    -------
    S: array [n_points, n_points]
        A similarity matrix S such that S_{i, j} is the distance between the
        ith and jth vectors of the given matrix X.

    """
    if metric == "precomputed":
        if X.shape[0] != X.shape[1]:
            raise ValueError("X is not square!")
        return X

    # In all other cases, the array is to be considered as a feature array.
    D = distance.squareform(distance.pdist(X, metric=metric))
    # Convert distance to similarity (FIXME: use heat kernel?)
    S = 1. - (D / np.max(D))
    return S

# Must be after calculate_similarity, as that method is required in dbscan_.
from .dbscan_ import dbscan, DBSCAN
