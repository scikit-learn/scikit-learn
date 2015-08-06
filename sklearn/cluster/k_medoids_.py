# -*- coding: utf-8 -*-
"""K-medoids clustering"""

# Authors: Timo Erkkilä <timo.erkkila@gmail.com>
#          Antti Lehmussola <antti.lehmussola@gmail.com>
# License: BSD 3 clause

import numpy as np
import warnings

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from ..utils import check_array

# With Python 2.x ValueError needs to be imported,
# otherwise that's not needed
try:
    from exceptions import ValueError
except ImportError:
    pass


class KMedoids(BaseEstimator, ClusterMixin, TransformerMixin):
    """
    k-medoids class.

    Parameters
    ----------
    n_clusters: int, optional, default: 8
      How many medoids. Must be positive.

    distance_metric: string, optional, default: 'euclidean'
      What distance metric to use.

    clustering: {'pam'}, optional, default: 'pam'
      What clustering mode to use.

    init: {'random', 'heuristic'}, optional, default: 'heuristic'
      Specify medoid initialization.
    """

    # Supported clustering methods
    CLUSTERING_METHODS = ['pam']

    # Supported initialization methods
    INIT_METHODS = ['random', 'heuristic']

    def __init__(self, n_clusters=8, distance_metric='euclidean',
                 clustering_method='pam', init='heuristic'):

        # Check n_clusters
        if (n_clusters is None or 
            n_clusters <= 0 or 
            not isinstance(n_clusters, int)):
            raise ValueError("n_clusters has to be nonnegative integer")

        # Check distance_metric
        if callable(distance_metric):
            self.distance_func = distance_metric
        elif distance_metric in PAIRWISE_DISTANCE_FUNCTIONS:
            self.distance_func = \
                PAIRWISE_DISTANCE_FUNCTIONS[distance_metric]
        else:
            raise ValueError("distance_metric needs to be " +
                             "callable or one of the " +
                             "following strings: " +
                             "{}".format(PAIRWISE_DISTANCE_FUNCTIONS.keys()) +
                             ". Instead, '{}' ".format(distance_metric) +
                             "was given.")

        # Check clustering_method
        if clustering_method not in self.CLUSTERING_METHODS:
            raise ValueError("clustering must be one of the following: " +
                             "{}".format(self.CLUSTERING_METHODS))

        # Check init
        if init not in self.INIT_METHODS:
            raise ValueError("init needs to be one of " +
                             "the following: " +
                             "{}".format(self.INIT_METHODS))

        self.n_clusters = n_clusters

        self.distance_metric = distance_metric

        self.init = init

        self.clustering_method = clustering_method

    def fit(self, X, y=None):
        """Fit K-Medoids to the provided data.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)

        Returns
        -------
        self
        """

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = check_array(X)

        # Apply distance metric to get the distance matrix
        D = self.distance_func(X)

        # Check that the number of clusters is less than or equal to
        # the number of clusters
        if self.n_clusters > D.shape[0]:
            raise ValueError("The number of medoids " +
                             "({}) ".format(self.n_clusters) +
                             "must be larger than the sides " +
                             "of the distance matrix ({})".format(D.shape[0]))

        medoid_ics = self.__get_initial_medoid_indices(D, self.n_clusters)

        # Old medoids will be stored here for reference
        old_medoid_ics = np.zeros((self.n_clusters,))

        # Continue the algorithm as long as
        # the medoids keep changing
        while not np.all(old_medoid_ics == medoid_ics):

            # Keep a copy of the old medoid assignments
            old_medoid_ics = np.copy(medoid_ics)

            # Get cluster indices
            cluster_ics = self.__get_cluster_ics(D, medoid_ics)

            # Update medoids with the new cluster indices
            self.__update_medoid_ics_in_place(D, cluster_ics, medoid_ics)

        # Expose labels_ which are the assignments of
        # the training data to clusters
        self.labels_ = cluster_ics

        # Expose cluster centers, i.e. medoids
        self.cluster_centers_ = X.take(medoid_ics, axis=0)

        # Return self to enable method chaining
        return self

    def __get_cluster_ics(self, D, medoid_ics):
        """Returns cluster indices for D and current medoid indices"""

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        cluster_ics = np.argmin(D[medoid_ics, :], axis=0)

        return cluster_ics

    def __update_medoid_ics_in_place(self, D, cluster_ics, medoid_ics):
        """In-place update of the medoid indices"""

        # Update the medoids for each cluster
        for cluster_idx in range(self.n_clusters):

            if sum(cluster_ics == cluster_idx) == 0:
                warnings.warn("Cluster {} is empty!".format(cluster_idx))
                continue

            # Find current cost that is associated with cluster_idx.
            # Cost is the sum of the distance from the cluster
            # members to the medoid.
            curr_cost = np.sum(D[medoid_ics[cluster_idx],
                                 cluster_ics == cluster_idx])

            # Extract the distance matrix between the data points
            # inside the cluster_idx
            D_in = D[cluster_ics == cluster_idx, :]
            D_in = D_in[:, cluster_ics == cluster_idx]

            # Calculate all costs there exists between all
            # the data points in the cluster_idx
            all_costs = np.sum(D_in, axis=1)

            # Find the index for the smallest cost in cluster_idx
            min_cost_idx = np.argmin(all_costs)

            # find the value of the minimum cost in cluster_idx
            min_cost = all_costs[min_cost_idx]

            # If the minimum cost is smaller than that
            # exhibited by the currently used medoid,
            # we switch to using the new medoid in cluster_idx
            if min_cost < curr_cost:

                # Find data points that belong to cluster_idx,
                # and assign the newly found medoid as the medoid
                # for cluster c
                medoid_ics[cluster_idx] = \
                    np.where(cluster_ics == cluster_idx)[0][min_cost_idx]

    def transform(self, X):
        """Transforms X to cluster-distance space.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
          Data to transform.

        Returns
        -------
        X_new : array, shape=(n_samples, n_clusters)
          X transformed in the new space.
        """

        # Apply distance metric wrt. cluster centers (medoids),
        # and return these distances
        return self.distance_func(X, Y=self.cluster_centers_)

    def predict(self, X):

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = check_array(X)

        # Apply distance metric wrt. cluster centers (medoids)
        D = self.distance_func(X, Y=self.cluster_centers_)

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        labels = np.argmin(D, axis=1)

        return labels

    def inertia(self, X):

        # Map the original X to the distance-space
        Xt = self.transform(X)

        # Define inertia as the sum of the sample-distances
        # to closest cluster centers
        inertia = np.sum(np.min(Xt, axis=1))

        return inertia

    def __get_initial_medoid_indices(self, D, n_clusters):

        if self.init == 'random':  # Random initialization

            # Pick random k medoids as the initial ones.
            medoids = np.random.permutation(D.shape[0])[:n_clusters]

        elif self.init == 'heuristic':  # Initialization by heuristic

            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = list(map(lambda x: x[0],
                               sorted(enumerate(np.sum(D, axis=1)),
                                      key=lambda x: x[1])))[:n_clusters]

        else:

            raise ValueError("Initialization not implemented for method: " +
                             "'{}'".format(self.init))

        return medoids
