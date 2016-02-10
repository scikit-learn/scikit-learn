from collections import defaultdict

import numpy as np
from scipy.spatial.distance import cdist


def distortion(X, labels, distortion_meth='sqeuclidean', p=2):
    """
    Given data and their cluster assigment, compute the distortion D

    Parameter
    ---------
    X: numpy array of shape (nb_data, nb_feature)
    labels: list of int of length nb_data
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
    distortion: float
    """
    if isinstance(distortion_meth, str):
        return distortion_metrics(X, labels, distortion_meth, p)
    else:
        return distortion_meth(X, labels)


def distortion_metrics(X, labels, metric='sqeuclidean', p=2):
    """
    Given data and their cluster assigment, compute the distortion D

    D = \sum_{x \in X} distance(x, c_x)

    With c_x the center of the cluster containing x, distance is the distance
    defined by metrics

    Parameter
    ---------
    X: numpy array of shape (nb_data, nb_feature)
    labels: list of int of length nb_data
    metric: string naming a scipy.spatial distance. metric can be in
        ['euclidian', 'minkowski', 'seuclidiean', 'sqeuclidean', 'chebyshev'
         'cityblock', 'cosine', 'correlation', 'hamming', 'jaccard',
         'Bray-Curtis', 'mahalanobis', 'yule', 'matching', 'dice', 'kulsinski',
         'rogerstanimoto', 'russellrao', 'sokalmichener', 'sokalsneath',
         'wminkowski', 'canberra']
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

    distance_sum = .0
    nb_feature = X.shape[1]
    for points in assi.values():
        clu_points = X[points, :]
        clu_center = np.mean(clu_points, axis=0).reshape(1, nb_feature)
        distance_sum += np.sum(cdist(
            clu_points, clu_center, metric=metric, p=p))

    return distance_sum / X.shape[1]
