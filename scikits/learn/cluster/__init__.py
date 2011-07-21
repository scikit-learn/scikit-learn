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

def calculate_similarity(X, metric=None, is_similarity=None):
    """ Calculates the similarity matrix from a vector matrix X.

    This method takes either a vector array or a similarity matrix, and returns
    a similarity matrix. If the input is a vector array, the similarities are
    computed. If the input is a similarity matrix, it is returned instead.

    This method provides a safe way to take a similarity matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    Parameters
    ----------
    X: array [n_points, n_points] or [n_points, n_features]
        Array of similarities between points, or a feature array.
        If the array is square, it is treated as a similarity array,
        otherwise it is treated as a feature array. Use is_similarity to
        override this pattern.
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded.
    is_similarity: boolean, optional (default=None)
        Overrides the behaviour of the array handling of S.
        If is_similarity is None, any square array is handled as a similarity
        array and any non-square array is a feature array.
        If is_similarity is True, any array is handled as a similarity array,
        and the procedure will raise a ValueError if the array is not square.
        If is_similarity is False, any array will be handled as a feature
        array, including square matrices.

    Returns
    -------
    S: array [n_points, n_points]
        A similarity matrix S such that S_{i, j} is the distance between the
        ith and jth vectors of the given matrix X.

    """
    n, d = X.shape
    # If the array looks square, it may be a similarity array.
    if n == d:
        if is_similarity in (None, True):
            return X
    elif is_similarity:
        # Array is not square, so it cannot be a similarity array.
        raise ValueError("Array not square, cannot be a similarity array."
                         " Shape = %s" % repr((n, d)))
    # In all other cases, the array is to be considered as a feature array.
    D = distance.squareform(distance.pdist(X, metric=metric))
    S = 1. - (D / np.max(D))
    return S

# Must be after calculate_similarity, as that method is required in dbscan_.
from .dbscan_ import dbscan, DBSCAN
