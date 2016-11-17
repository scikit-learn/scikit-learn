# -*- coding: utf-8 -*-
"""K-medoids clustering"""

# Authors: Timo Erkkilä <timo.erkkila@gmail.com>
#          Antti Lehmussola <antti.lehmussola@gmail.com>
#          Kornel Kiełczewski <kornel.k@plusnet.pl>
# License: BSD 3 clause

import warnings

import numpy as np

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import pairwise_distances
from ..utils import check_array, check_random_state
from ..utils.validation import check_is_fitted


class KMedoids(BaseEstimator, ClusterMixin, TransformerMixin):
    """
    k-medoids class.

    Read more in the :ref:`User Guide <k_medoids>`.

    Parameters
    ----------
    n_clusters : int, optional, default: 8
        How many medoids. Must be positive.

    distance_metric : string, optional, default: 'euclidean'
        What distance metric to use.

    init : {'random', 'heuristic'}, optional, default: 'heuristic'
        Specify medoid initialization method. Random randomly selects
        n_clusters elements from the dataset, while heuristic picks n_clusters
        first data points that have the smallest sum distance to every other
        point.

    max_iter : int, optional, default : 300
        Specify the maximum number of iterations when fitting.

    random_state : int, optional, default: None
        Specify random state for the random number generator.

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Cluster centers, i.e. medoids (elements from the original dataset)

    labels_ :
        Labels of each point

    inertia_ : float
        Sum of distances of samples to their closest cluster center.

    Examples
    --------

    >>> from sklearn.cluster import KMedoids
    >>> import numpy as np

    >>> X = np.asarray([[1, 2], [1, 4], [1, 0],
    ...                 [4, 2], [4, 4], [4, 0]])
    >>> kmedoids = KMedoids(n_clusters=2, random_state=0).fit(X)
    >>> kmedoids.labels_
    array([0, 0, 0, 1, 1, 1])
    >>> kmedoids.predict([[0,0], [4,4]])
    array([0, 1])
    >>> kmedoids.cluster_centers_
    array([[1, 2],
           [4, 2]])
    >>> kmedoids.inertia_
    8.0

    >>> kmedoids = KMedoids(n_clusters=2, random_state=0, distance_metric='manhattan').fit(X)
    >>> kmedoids.labels_
    array([0, 0, 0, 1, 1, 1])
    >>> kmedoids.predict([[0,0], [4,4]])
    array([0, 1])
    >>> kmedoids.cluster_centers_
    array([[1, 2],
           [4, 2]])
    >>> kmedoids.inertia_
    8.0
    """

    def __init__(self, n_clusters=8, distance_metric='euclidean',
                 init='heuristic', max_iter=300, random_state=None):
        self.n_clusters = n_clusters
        self.distance_metric = distance_metric
        self.init = init
        self.max_iter = max_iter
        self.random_state = random_state

    def _check_init_args(self):
        """Validates the input arguments. """

        # Check n_clusters
        if (self.n_clusters is None or self.n_clusters <= 0 or
                not isinstance(self.n_clusters, int)):
            raise ValueError("n_clusters should be a nonnegative integer. "
                             "%s was given" % self.n_clusters)

        # Check init
        init_methods = ['random', 'heuristic']
        if self.init not in init_methods:
            raise ValueError("init needs to be one of " +
                             "the following: " +
                             "%s" % init_methods)

        # Check random state
        self.random_state_ = check_random_state(self.random_state)

    def fit(self, X, y=None):
        """Fit K-Medoids to the provided data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape=(n_samples, n_features).
        Dataset to cluster.

        Returns
        -------
        self
        """

        self._check_init_args()

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = self._check_input(X)

        # Apply distance metric to get the distance matrix
        distances = pairwise_distances(X, metric=self.distance_metric)

        medoid_ics = self._get_initial_medoid_indices(distances,
                                                      self.n_clusters)
        cluster_ics = self._get_cluster_ics(distances, medoid_ics)

        # Old medoids will be stored here for reference
        old_medoid_ics = np.zeros((self.n_clusters,))

        # Continue the algorithm as long as
        # the medoids keep changing and the maximum number
        # of iterations is not exceeded
        self.n_iter_ = 0
        while (not np.all(old_medoid_ics == medoid_ics) and
                self.n_iter_ < self.max_iter):
            self.n_iter_ += 1

            # Keep a copy of the old medoid assignments
            old_medoid_ics = np.copy(medoid_ics)

            # Get cluster indices
            cluster_ics = self._get_cluster_ics(distances, medoid_ics)

            # Update medoids with the new cluster indices
            self._update_medoid_ics_in_place(distances, cluster_ics,
                                             medoid_ics)

        # Expose labels_ which are the assignments of
        # the training data to clusters
        self.labels_ = cluster_ics

        # Expose cluster centers, i.e. medoids
        self.cluster_centers_ = X[medoid_ics]

        # Expose intertia
        self.inertia_ = self._compute_inertia(X)

        # Return self to enable method chaining
        return self

    def _check_input(self, X):

        X = check_array(X, accept_sparse=['csr', 'csc'])
        # Check that the number of clusters is less than or equal to
        # the number of samples
        if self.n_clusters > X.shape[0]:
            raise ValueError("The number of medoids %d must be larger "
                             "than the number of samples %d."
                             % (self.n_clusters, X.shape[0]))

        return X

    def _get_cluster_ics(self, distances, medoid_ics):
        """Returns cluster indices for distances and current medoid indices"""

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        cluster_ics = np.argmin(distances[medoid_ics, :], axis=0)

        return cluster_ics

    def _update_medoid_ics_in_place(self, distances, cluster_ics, medoid_ics):
        """In-place update of the medoid indices"""

        # Update the medoids for each cluster
        for cluster_idx in range(self.n_clusters):

            if sum(cluster_ics == cluster_idx) == 0:
                warnings.warn("Cluster %d is empty!" % cluster_idx)
                continue

            # Find current cost that is associated with cluster_idx.
            # Cost is the sum of the distance from the cluster
            # members to the medoid.
            curr_cost = np.sum(distances[medoid_ics[cluster_idx],
                                         cluster_ics == cluster_idx])

            # Extract the distance matrix between the data points
            # inside the cluster_idx
            in_cluster_distances = distances[cluster_ics == cluster_idx, :]
            in_cluster_distances = (
                in_cluster_distances[:, cluster_ics == cluster_idx])

            # Calculate all costs there exists between all
            # the data points in the cluster_idx
            all_costs = np.sum(in_cluster_distances, axis=1)

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
                medoid_ics[cluster_idx] = np.where(
                    cluster_ics == cluster_idx)[0][min_cost_idx]

    def transform(self, X):
        """Transforms X to cluster-distance space.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape=(n_samples, n_features)
            Data to transform.

        Returns
        -------
        X_new : {array-like, sparse matrix}, shape=(n_samples, n_clusters)
            X transformed in the new space.
        """
        X = check_array(X, accept_sparse=['csr', 'csc'])

        check_is_fitted(self, "cluster_centers_")

        # Apply distance metric wrt. cluster centers (medoids),
        # and return these distances
        return pairwise_distances(X, Y=self.cluster_centers_,
                                  metric=self.distance_metric)

    def predict(self, X):
        """Predict the closest cluster each sample in X belongs to.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            New data to predict.

        Returns
        -------
        labels : array, shape [n_samples,]
            Index of the cluster each sample belongs to.
        """

        check_is_fitted(self, "cluster_centers_")

        # Check that the array is good and attempt to convert it to
        # Numpy array if possible
        X = check_array(X, accept_sparse=['csr', 'csc'])

        # Apply distance metric wrt. cluster centers (medoids)
        distances = pairwise_distances(X, Y=self.cluster_centers_,
                                       metric=self.distance_metric)

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        labels = np.argmin(distances, axis=1)

        return labels

    def _compute_inertia(self, X):
        """Compute inertia of new samples. Inertia is defined as the sum of the
        sample distances to closest cluster centers.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape=(n_samples, n_features)
            Samples to compute inertia for.

        Returns
        -------
        Sum of sample distances to closest cluster centers.

        """

        # Map the original X to the distance-space
        distances = self.transform(X)

        # Define inertia as the sum of the sample-distances
        # to closest cluster centers
        inertia = np.sum(np.min(distances, axis=1))

        return inertia

    def _get_initial_medoid_indices(self, distances, n_clusters):
        """Initialize medoid indices using either a random or
        heuristic approach.
        """

        if self.init == 'random':  # Random initialization
            # Pick random k medoids as the initial ones.
            medoids = self.random_state_.permutation(
                distances.shape[0])[:n_clusters]
        elif self.init == 'heuristic':  # Initialization by heuristic
            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = list(np.argsort(np.sum(distances, axis=1))[:n_clusters])
        else:
            raise ValueError("Initialization not implemented for method: '%s'"
                             % self.init)

        return medoids
