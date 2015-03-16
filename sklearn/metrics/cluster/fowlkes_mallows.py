from __future__ import division

from .adjacency_matrix import adjacency_matrix

from scipy.spatial.distance import cosine


def fowlkes_mallows_index(clustering_1, clustering_2):
    """
    Mesure the similarity of two clusterings of a set of points.
    Let:
    - TP be the number of pair of points (x_i, x_j) that belongs
        in the same clusters in both clustering_1 and clustering_2
    - FP be the number of pair of points (x_i, x_j) that belongs
        in the same clusters in clustering_1 and not in clustering_2
    - FN be the number of pair of points (x_i, x_j) that belongs
        in the same clusters in clustering_2 and not in clustering_1

    The Fowlkes-Mallows index has the following formula:

        fowlkes_mallows_index = TP / sqrt((TP + FP) * (TP + FN))

    Parameter
    ---------
    clustering_1: list of int.
        "clustering_1[i] = c" means that point i is assigned to cluster c
    clustering_2: list of int.
        "clustering_2[i] = c" means that point i is assigned to cluster c

    Return
    ------
    fowlkes_mallows_index: float between 0 and 1. 1 means that both
        clusterings perfectly match, 0 means that they totally disconnect
    """
    adj_mat_1 = adjacency_matrix(clustering_1)
    adj_mat_2 = adjacency_matrix(clustering_2)
    return 1 - cosine(adj_mat_1.flatten(), adj_mat_2.flatten())
