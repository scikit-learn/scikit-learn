import numpy as np


def adjacency_matrix(cluster_assignement):
    """
    Parameter
    ---------
    cluster_assignement: vector (n_samples) of int i, 0 <= i < k

    Return
    ------
    adj_matrix: matrix (n_samples, n_samples)
        adji_matrix[i, j] = cluster_assignement[i] == cluster_assignement[j]
    """
    n_samples = len(cluster_assignement)
    adj_matrix = np.zeros((n_samples, n_samples))
    for i, val in enumerate(cluster_assignement):
        for j in range(i, n_samples):
            linked = val == cluster_assignement[j]
            adj_matrix[i, j] = linked
            adj_matrix[j, i] = linked
    return adj_matrix
