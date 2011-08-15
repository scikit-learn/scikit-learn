"""
Metrics module with score functions, performance metrics and
pairwise metrics or distances computation
"""
from scipy.spatial import distance

from .metrics import confusion_matrix, roc_curve, auc, precision_score, \
                recall_score, fbeta_score, f1_score, zero_one_score, \
                precision_recall_fscore_support, classification_report, \
                precision_recall_curve, explained_variance_score, r2_score, \
                zero_one, mean_square_error, hinge_loss

from .cluster import homogeneity_completeness_v_measure
from .cluster import homogeneity_score
from .cluster import completeness_score
from .cluster import v_measure_score
from . import pairwise
from .pairwise import euclidean_distances

pairwise_function_map = {}
pairwise_function_map['euclidean'] = pairwise.euclidean_distances
pairwise_function_map['precomputed'] = pairwise.return_self_if_square
pairwise_function_map['l1'] = pairwise.l1_distances


def pairwise_distances(X, Y=None, metric="euclidean"):
    """ Calculates the distance matrix from a vector matrix X.

    This method takes either a vector array or a distance matrix, and returns
    a distance matrix. If the input is a vector array, the distances are
    computed. If the input is a distances matrix, it is returned instead.

    This method provides a safe way to take a distance matrix as input, while
    preserving compatability with many other algorithms that take a vector
    array.

    Parameters
    ----------
    X: array [n_samples, n_samples] if metric == "precomputed", or,
             [n_samples, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    X: array [n_samples, n_features]
        A second feature array only if X has shape [n_samples, n_features].

    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    Returns
    -------
    D: array [n_samples, n_samples]
        A distance matrix D such that D_{i, j} is the distance between the
        ith and jth vectors of the given matrix X.

    """
    if metric in pairwise_function_map:
        return pairwise_function_map[metric](X, Y)
    else:
        # FIXME: the distance module doesn't support sparse matrices!
        if Y is None:
            return distance.squareform(distance.pdist(X, metric=metric))
        else:
            return distance.cdist(X, Y, metric=metric)
