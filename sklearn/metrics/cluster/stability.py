from __future__ import division

from math import sqrt

import numpy as np
from .adjacency_matrix import adjacency_matrix
from .fowlkes_mallows import fowlkes_mallows_index
from scipy.spatial.distance import pdist

from ...utils import check_random_state


def stability(X, cluster_estimator, k_max=None, nb_draw=100, prop_subset=.8,
              random_state=None, p=None, distance='fowlkes-mallows',
              verbose=False):
    """Stability algorithm.
    For k from 2 to k_max, compute stability of cluster estimator to produce k
    clusters. Stability measures if the estimator produces the same clusters
    given small variations in the input data. It draws two overlapping subsets
    A and B of input data. For points in the two subsets, we compute the
    clustering C_A and C_B done on subsets A and B. We then compute the
    similarity of those clustering. We can use the opposite of a distance
    as a similarity

    The stability of cluster_estimator with k cluster is the expectation of
    similarity(C_A, C_B)

    Ref: Ben-Hur, Elisseeff, Guyon: a stability based method for discovering
    structure in clusterd data, 2002
    Overview of stability: Luxburg: clustering stability: an overview

    Parameters
    ----------
    X : array-like or sparse matrix, shape (n_samples, n_features)
        The observations to cluster.
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    k_max: int: maximum number of clusters (default = n_samples / 2)
    nb_draw: number of draws to estimate expectation of expectation of
        similarity(C_A, C_B)
    prop_subset: 0 < float < 1: proportion of input data taken in each subset
    distance: a string naming a distance or a cluster similarity. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski', 'fowlkes-mallows'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)
    Return
    ------
    k: int
    """
    rng = check_random_state(random_state)
    cluster_similarity = function_cluster_similarity(distance, p)

    n_samples, n_features = X.shape
    if not k_max:
        k_max = n_samples // 2

    best_stab, best_k = 0, 0
    for k in range(2, k_max + 1):
        cluster_estimator.set_params(n_clusters=k)
        this_score = sum(
            _one_stability_measure(cluster_estimator, X, prop_subset,
                                   cluster_similarity)
            for _ in range(nb_draw)) / nb_draw
        if verbose:
            print('for %d cluster, stability is %f' % (k, this_score))

        if this_score >= best_stab:
            best_stab = this_score
            best_k = k

    return best_k


def _one_stability_measure(cluster_estimator, X, prop_sample,
                           cluster_similarity, random_state=None):
    """
    Draws two subsets A and B from X, compute C_A, clustering on subset
    A, and C_B, clustering on subset B, then returns

    similarity(C_A, C_B)

    Parameter
    ---------
    X: array of size n_samples, n_features
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    prop_sample: 0 < float < 1, proportion of X taken in each subset
    cluster_similarity: function (list, list) -> float
    """
    rng = check_random_state(random_state)

    n_sample = X.shape[0]
    set_1 = rng.uniform(size=n_sample) < prop_sample
    set_2 = rng.uniform(size=n_sample) < prop_sample
    nb_points_1, nb_points_2 = 0, 0
    points_1, points_2 = [], []
    common_points_1, common_points_2 = [], []
    for i, (is_1, is_2) in enumerate(zip(set_1, set_2)):
        if is_1 and is_2:
            common_points_1.append(nb_points_1)
            common_points_2.append(nb_points_2)
        if is_1:
            points_1.append(i)
            nb_points_1 += 1
        if is_2:
            points_2.append(i)
            nb_points_2 += 1

    assi_1 = cluster_estimator.fit_predict(X[np.ix_(points_1)])
    assi_2 = cluster_estimator.fit_predict(X[np.ix_(points_2)])

    clustering_1 = [assi_1[c] for c in common_points_1]
    clustering_2 = [assi_2[c] for c in common_points_2]
    return cluster_similarity(clustering_1, clustering_2)


def function_cluster_similarity(metric='fowlkes-mallows', p=None):
    """
    Given the name of a distance, return function to estimate
    two clusterings  similarity

    Parameter
    --------
    metric: a string naming a distance or a cluster similarity. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski', 'fowlkes-mallows'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)

    Return:
    function (clustering_1, clustering_2) -> similarity; with:
        clustering_k: a list. clustering_k[i] = c means that
           point x_i belongs to cluster c in clustering k
        similarity: float
    """
    if metric == 'fowlkes-mallows':
        return fowlkes_mallows_index
    if metric == 'l2':
        # Translate to something understood by scipy
        metric = 'euclidean'
    elif metric in ('l1', 'manhattan'):
        metric = 'cityblock'

    def cluster_dist(clustering_1, clustering_2):
        adj_mat_1 = adjacency_matrix(clustering_1).flatten()
        adj_mat_2 = adjacency_matrix(clustering_2).flatten()
        return -pdist([adj_mat_1, adj_mat_2], metric=metric, p=p)
    return cluster_dist
