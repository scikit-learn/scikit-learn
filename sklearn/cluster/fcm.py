# -*- coding: utf-8 -*-
"""Fuzzy C-Means clustering"""

# Authors: Junyi Li (lijy263@mail2.sysu.edu.cn)
# License: BSD 3 clause
# The import part this file is copy from k_means_.py,
# some functions are not used for now.
import numpy as np

from ..base import BaseEstimator, ClusterMixin
from ..utils import check_array
from ..utils import check_random_state
from ..preprocessing import normalize


def reconstruct_label(labels):
    if len(labels) <= 0:
        return labels
    tmp = []
    for i, xi in enumerate(labels):
        if xi in tmp:
            labels[i] = tmp.index(xi)
        else:
            tmp.append(xi)
            labels[i] = len(tmp) - 1
    return labels


def fcm(X, n_clusters, eps=3, m=2, random_state=None, max_iter=300,
        sample_weight=None):
    """Fuzzy CMeans clustering

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)

    n_clusters : integer
        The number of seeds to choose

    eps : float, optional (default = 3)
        If the sum of the abs of the (new_cluster_probability_matrix \
        - old_cluster_probability_matrix) is smaller than eps, \
         then the algorithm will stop.

    m : float, optional (default = 2)
        The fuzzy number.

    random_state : int, RandomState instance or None (default)
        The generator used to initialize the centers. Use an int to make the
        randomness deterministic.

    max_iter : int, optional (default = 300)
        The maximum iteration time for FCM clustering algorthm.

    sample_weight : array, shape (n_samples,), optional
        Weight of each sample, such that a sample with a weight of at least
        ``min_samples`` is by itself a core sample; a sample with negative
        weight may inhibit its eps-neighbor from being core.
        Note that weights are absolute, and default to 1.

    Returns
    -------
    centroid : float ndarray with shape (n_clusters, n_features)
        Centroids found at the last iteration of fuzzy c-means.

    labels : array [n_samples]
        Cluster labels for each point.

    """
    if m <= 1:
        raise ValueError("Invalid number of m."
                         " m=%d must be bigger than 1." % m)
    if eps <= 0:
        raise ValueError("Invalid number of eps."
                         " eps=%d must be bigger than 0." % eps)
    random_state = check_random_state(random_state)

    if max_iter <= 0:
        raise ValueError('Number of iterations should be a positive number,'
                         ' got %d instead' % max_iter)

    if len(X) < n_clusters:
        n_clusters = len(X)

    X = check_array(X, accept_sparse='csr')
    X = normalize(X)
    membership_mat = np.random.random((len(X), n_clusters))
    membership_mat = membership_mat / np.sum(membership_mat,
                                             axis=1)[:, np.newaxis]

    for iter_time in range(max_iter):
        working_membership_mat = membership_mat ** m
        Centroids = np.dot(working_membership_mat.T,
                           X) / np.sum(working_membership_mat.T,
                                       axis=1)[:, np.newaxis]

        n_c_distance_mat = np.zeros((len(X), n_clusters))
        for i, x in enumerate(X):
            for j, c in enumerate(Centroids):
                n_c_distance_mat[i][j] = np.linalg.norm(x - c, 2)

        new_membership_mat = np.zeros((len(X), n_clusters))

        for i, x in enumerate(X):
            for j, c in enumerate(Centroids):
                new_membership_mat[i][j] = 1. / np.sum(
                    (
                            n_c_distance_mat[i][j] /
                            n_c_distance_mat[i]) ** (2 / (m - 1)
                                                     )
                )
        if np.sum(abs(new_membership_mat - membership_mat)) < eps:
            break
        membership_mat = new_membership_mat
    labels = np.argmax(new_membership_mat, axis=1)
    reconstruct_label(labels)
    return Centroids, labels, iter_time


class FCM(BaseEstimator, ClusterMixin):
    """Perform fuzzy c-means clustering

    Parameters
    ----------
    X : array or sparse (CSR) matrix of shape (n_samples, n_features), or \
            array of shape (n_samples, n_samples)

    n_clusters : integer, optional (default = 3)
        The number of seeds to choose

    eps : float, optional (default = 3)
        If the sum of the abs of the (new_cluster_probability_matrix \
        - old_cluster_probability_matrix) is smaller than eps, \
         then the algorithm will stop.

    m : float, optional (default = 2)
        The fuzzy number.

    random_state : int, RandomState instance or None (default)
        The generator used to initialize the centers. Use an int to make the
        randomness deterministic.

    max_iter : int, optional (default = 300)
        The maximum iteration time for FCM clustering algorthm.

    sample_weight : array, shape (n_samples,), optional
        Weight of each sample, such that a sample with a weight of at least
        ``min_samples`` is by itself a core sample; a sample with negative
        weight may inhibit its eps-neighbor from being core.
        Note that weights are absolute, and default to 1.

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Coordinates of cluster centers.

    labels_ :
        Labels of each point

    Notes
    -----
    Now, something remains implementing: sample weighted \
    section and parallel run the model.
    """
    def __init__(self, n_clusters=3, eps=3, m=2, init='random',
                 random_state=None, max_iter=300):
        self.n_clusters = n_clusters
        self.init = init
        self.m = m
        self.eps = eps
        self.max_iter = max_iter
        self.random_state = random_state

    def _check_test_data(self, X):
        X = check_array(X, accept_sparse='csr')
        n_samples, n_features = X.shape
        expected_n_features = self.cluster_centers_.shape[1]
        if not n_features == expected_n_features:
            raise ValueError("Incorrect number of features. "
                             "Got %d features, expected %d" % (
                                 n_features, expected_n_features))

        return X

    def fit(self, X, y=None, sample_weight=None):
        """Compute Fuzzy c-means clustering.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Training instances to cluster.

        y : Ignored
            not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None)

        """

        random_state = check_random_state(self.random_state)
        self.cluster_centers_, self.labels_, self.n_iter_ = \
            fcm(
                X,
                n_clusters=self.n_clusters,
                m=self.m,
                eps=self.eps,
                random_state=random_state,
                max_iter=self.max_iter,
                sample_weight=sample_weight
                )
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        """Compute cluster centers and predict cluster index for each sample.

        Convenience method; equivalent to calling fit(X) followed by
        predict(X).

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to transform.

        y : Ignored
            not used, present here for API consistency by convention.

        sample_weight : array-like, shape (n_samples,), optional
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None)

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """
        return self.fit(X, sample_weight=sample_weight).labels_
