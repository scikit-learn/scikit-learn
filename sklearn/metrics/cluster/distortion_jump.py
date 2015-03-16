from math import pow

from .distortion import distortion

import numpy as np


def distortion_jump(X, cluster_estimator, k_max=None,
                    distortion_meth='sqeuclidean', p=2):
    """
    Find the number of clusters that maximizes efficiency while minimizing
    error by information theoretic standards (wikipedia). For each number of
    cluster, it calculates the distortion reduction. Roughly, it selects k such
    as the difference between distortion with k clusters minus distortion with
    k-1 clusters is maximal.

    More precisely, let d(k) equals distortion with k clusters.
    Let Y=nb_feature/2, let D[k] = d(k)^{-Y}
    k^* = argmax(D[k] - D[k-1])

    Parameters
    ----------
    X: numpy array of shape (nb_date, nb_features)
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    k_max: int: maximum number of clusters
    distortion_meth: can be a function X, labels -> float,
        can be a string naming a scipy.spatial distance. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)

    Return
    ------
    k_star: int: optimal number of cluster
    """
    nb_data, nb_feature = X.shape
    # if no maximum number of clusters set, take datasize divided by 2
    if not k_max:
        k_max = nb_data // 2

    Y = - nb_feature / 2
    info_gain = 0
    old_dist = pow(
        distortion(X, np.zeros(nb_data), distortion_meth, p) / nb_feature, Y)
    for k in range(2, k_max + 1):
        cluster_estimator.set_params(n_clusters=k)
        labs = cluster_estimator.fit_predict(X)
        new_dist = pow(
            distortion(X, labs, distortion_meth, p) / nb_feature, Y)
        if new_dist - old_dist >= info_gain:
            k_star = k
            info_gain = new_dist - old_dist
        old_dist = new_dist
    return k_star
