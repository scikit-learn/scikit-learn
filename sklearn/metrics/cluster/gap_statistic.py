from __future__ import division
from collections import defaultdict

from math import sqrt, log

import numpy as np

from scipy.spatial.distance import cdist

from ...utils import check_random_state


def inertia(X, labels, metric='sqeuclidean', p=2):
    """
    Given data and their cluster assigment, compute the sum of
    within-cluster mean distance to cluster's mean

    Parameter
    ---------
    X: numpy array of shape (nb_data, nb_feature)
    labels: list of int of length nb_data
    metric: a string naming a scipy.spatial distance. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)

    Return
    ------
    distortion: float
    """
    if metric == 'l2':
        # Translate to something understood by scipy
        metric = 'euclidean'
    elif metric in ('l1', 'manhattan'):
        metric = 'cityblock'

    assi = defaultdict(list)
    for i, l in enumerate(labels):
        assi[l].append(i)

    inertia = .0
    nb_feature = X.shape[1]
    for points in assi.values():
        clu_points = X[points, :]
        clu_center = np.mean(clu_points, axis=0).reshape(1, nb_feature)
        inertia += (np.sum(cdist(clu_points, clu_center, metric=metric, p=p)) /
                    (2 * len(clu_points)))
    return inertia


def normal_inertia(X, cluster_estimator, nb_draw=100,
                   metric='sqeuclidean', p=2, random_state=None,
                   mu=None, sigma=None):
    """
    Draw multivariate normal data of size data_shape = (nb_data, nb_feature),
    with same mean and covariance as X.
    Clusterize data using cluster_estimator and compute inertia

    Parameter
    ---------
    X numpy array of size (nb_data, nb_feature)
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    nb_draw: number of samples to calculate expected_inertia
    metric: a string naming a scipy.spatial distance. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)
    mu: mean of drawn data
    sigma: covariance matrix of drawn data

    Return
    ------
    dist: list of inertias (float) obtained on random dataset
    """
    rng = check_random_state(random_state)
    nb_data, nb_feature = X.shape

    if mu is None:
        # data mean has no influence on distortion
        mu = np.zeros(nb_feature)
    if sigma is None:
        sigma = np.cov(X.transpose())

    dist = []
    for i in range(nb_draw):
        X_rand = rng.multivariate_normal(mu, sigma, size=nb_data)
        dist.append(inertia(
            X_rand, cluster_estimator.fit_predict(X_rand),
            metric, p))

    return dist


def uniform_inertia(X, cluster_estimator, nb_draw=100, val_min=None,
                    val_max=None, metric='sqeuclidean', p=2,
                    random_state=None):
    """
    Uniformly draw data of size data_shape = (nb_data, nb_feature)
    in the smallest hyperrectangle containing real data X.
    Clusterize data using cluster_estimator and compute inertia

    Parameter
    ---------
    X: numpy array of shape (nb_data, nb_feature)
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    nb_draw: number of samples to calculate expected_inertia
    val_min: minimum values of each dimension of input data
        array of length nb_feature
    val_max: maximum values of each dimension of input data
        array of length nb_feature
    metric: a string naming a scipy.spatial distance. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)

    Return
    ------
    dist: list of distortions (float) obtained on random dataset
    """
    rng = check_random_state(random_state)
    if val_min is None:
        val_min = np.min(X, axis=0)
    if val_max is None:
        val_max = np.max(X, axis=0)

    dist = []
    for i in range(nb_draw):
        X_rand = rng.uniform(size=X.shape) * (val_max - val_min) + val_min
        dist.append(inertia(X_rand, cluster_estimator.fit_predict(X_rand),
                            metric, p))

    return dist


def gap_statistic(X, cluster_estimator, k_max=None, nb_draw=100,
                  random_state=None, draw_model='uniform',
                  metric='sqeuclidean', p=2):
    """
    Estimating optimal number of cluster for data X with cluster_estimator by
    comparing inertia of clustered real data with inertia of clustered
    random data. Let W_rand(k) be the inertia of random data in k clusters,
    W_real(k) inertia of real data in k clusters, statistic gap is defined
    as

    Gap(k) = E(log(W_rand(k))) - log(W_real(k))

    We draw nb_draw random data "shapened-like X" (shape depend on draw_model)
    We select the smallest k such as the gap between inertia of k clusters
    of random data and k clusters of real data is superior to the gap with
    k + 1 clusters minus a "standard-error" safety. Precisely:

    k_star = min_k k
         s.t. Gap(k) >= Gap(k + 1) - s(k + 1)
              s(k) = stdev(log(W_rand)) * sqrt(1 + 1 / nb_draw)

    From R.Tibshirani, G. Walther and T.Hastie, Estimating the number of
    clusters in a dataset via the Gap statistic, Journal of the Royal
    Statistical Socciety: Seris (B) (Statistical Methodology), 63(2), 411-423

    Parameter
    ---------
    X: data. array nb_data * nb_feature
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
    nb_draw: int: number of random data of shape (nb_data, nb_feature) drawn
        to estimate E(log(D_rand(k)))
    draw_model: under which i.i.d data are draw. default: uniform data
        (following Tibshirani et al.)
        can be 'uniform', 'normal' (Gaussian distribution)
    metric: a string naming a scipy.spatial distance. can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'canberra', 'wminkowski'])
    p : double
        The p-norm to apply (for Minkowski, weighted and unweighted)

    Return
    ------
    k: int: number of cluster that maximizes the gap statistic
    """
    rng = check_random_state(random_state)

    # if no maximum number of clusters set, take datasize divided by 2
    if not k_max:
        k_max = X.shape[0] // 2
    if draw_model == 'uniform':
        val_min = np.min(X, axis=0)
        val_max = np.max(X, axis=0)
    elif draw_model == 'normal':
        mu = np.mean(X, axis=0)
        sigma = np.cov(X.transpose())

    old_gap = - float("inf")
    for k in range(2, k_max + 2):
        cluster_estimator.set_params(n_clusters=k)
        real_dist = inertia(X, cluster_estimator.fit_predict(X),
                            metric, p)
        # expected distortion
        if draw_model == 'uniform':
            rand_dist = uniform_inertia(X, cluster_estimator, nb_draw,
                                        val_min, val_max, metric,
                                        p)
        elif draw_model == 'normal':
            rand_dist = normal_inertia(X, cluster_estimator, nb_draw,
                                       metric=metric,
                                       p=p, mu=mu, sigma=sigma)
        else:
            raise ValueError(
                "For gap statistic, model for random data is unknown")
        rand_dist = np.log(rand_dist)
        exp_dist = np.mean(rand_dist)
        std_dist = np.std(rand_dist)
        gap = exp_dist - log(real_dist)
        safety = std_dist * sqrt(1 + 1 / nb_draw)
        if k > 2 and old_gap >= gap - safety:
            return k - 1
        old_gap = gap
    # if k was found, the function would have returned
    # no clusters were found -> only 1 cluster
    return 1
