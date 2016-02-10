from collections import defaultdict

import numpy as np


def calinski_harabaz_index(X, labels):
    """
    Compute the Calinski and Harabaz (1974). It a ratio between the
    within-cluster dispersion and the between-cluster dispersion

    CH(k) = trace(B_k) / (k -1) * (n - k) / trace(W_k)

    With B_k the between group dispersion matrix, W_k the within-cluster
    dispersion matrix

    B_k = \sum_q n_q (c_q - c) (c_q -c)^T
    W_k = \sum_q \sum_{x \in C_q} (x - c_q) (x - c_q)^T

    Ref: R.B.Calinsky, J.Harabasz: A dendrite method for cluster analysis 1974

    Parameter
    ---------
    X: numpy array of size (nb_data, nb_feature)
    labels: list of int of length nb_data: labels[i] is the cluster
        assigned to X[i, :]

    Return
    ------
    res: float: mean silhouette of this clustering
    """
    assi = defaultdict(list)
    for i, l in enumerate(labels):
        assi[l].append(i)

    nb_data, nb_feature = X.shape
    disp_intra = np.zeros((nb_feature, nb_feature))
    disp_extra = np.zeros((nb_feature, nb_feature))
    center = np.mean(X, axis=0)

    for points in assi.values():
        clu_points = X[points, :]
        # unbiaised estimate of variace is \sum (x - mean_x)^2 / (n - 1)
        # so, if I want sum of dispersion, I need
        # W_k = cov(X) * (n - 1)
        nb_point = clu_points.shape[0]
        disp_intra += np.cov(clu_points, rowvar=0) * (nb_point - 1)
        extra_var = (np.mean(clu_points, axis=0) - center).reshape(
            (nb_feature, 1))
        disp_extra += np.multiply(extra_var, extra_var.transpose()) * nb_point
    return (disp_extra.trace() * (nb_data - len(assi)) /
            (disp_intra.trace() * (len(assi) - 1)))


def calc_calinski_harabaz(X, cluster_estimator, n_clusters):
    """
    Compute calinski harabaz for clusters made by cluster estimator

    Parameter
    ---------
    X numpy array of size (nb_data, nb_feature)
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    n_clusters: number of clusters

    """
    cluster_estimator.set_params(n_clusters=n_clusters)
    return calinski_harabaz_index(X, cluster_estimator.fit_predict(X))


def max_CH_index(X, cluster_estimator, k_max=None):
    """
    Select number of cluster maximizing the Calinski and Harabasz (1974).
    It a ratio between the within-cluster dispersion and the between-cluster
    dispersion

    Ref: R.B.Calinsky, J.Harabasz: A dendrite method for cluster analysis 1974

    Parameters
    ----------
    X: numpy array of shape (nb_date, nb_features)
    cluster_estimator: ClusterMixing estimator object.
        need parameter n_clusters
        need method fit_predict: X -> labels
    k_max: int: maximum number of clusters

    Return
    ------
    k_star: int: optimal number of cluster
    """
    # if no maximum number of clusters set, take datasize divided by 2
    if not k_max:
        k_max = X.shape[0] // 2

    return max((k for k in range(2, k_max + 1)),
               key=lambda k: calc_calinski_harabaz(X, cluster_estimator, k))
