# -*- coding: utf-8 -*-
"""K-medoids clustering"""

# Authors: Timo Erkkilä <timo.erkkila@gmail.com>
#          Antti Lehmussola <antti.lehmussola@gmail.com>
#          Kornel Kiełczewski <kornel.k@plusnet.pl>
# License: BSD 3 clause

import warnings

import numpy as np

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import pairwise_distances, pairwise_distances_argmin
from ..utils import check_array, check_random_state
from ..utils.validation import check_is_fitted
from ..exceptions import ConvergenceWarning


class KMedoids(BaseEstimator, ClusterMixin, TransformerMixin):
    """k-medoids clustering.

    Read more in the :ref:`User Guide <k_medoids>`.

    Parameters
    ----------
    n_clusters : int, optional, default: 8
        The number of clusters to form as well as the number of medoids to
        generate.


    metric : string, or callable, optional, default: 'euclidean'
        What distance metric to use. See :func:metrics.pairwise_distances

    init : {'random', 'heuristic'}, optional, default: 'heuristic'
        Specify medoid initialization method. Random selects n_clusters
        elements from the dataset, while heuristic picks the n_clusters points
        with the smallest sum distance to every other point.

    max_iter : int, optional, default : 300
        Specify the maximum number of iterations when fitting.

    random_state : int, optional, default: None
        Specify random state for the random number generator. Used to
        initialise medoids when init='random'.

    Attributes
    ----------
    cluster_centers_ : array, shape = (n_clusters, n_features)
        Cluster centers, i.e. medoids (elements from the original dataset)

    labels_ : array, shape = (n_samples,)
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

    >>> kmedoids = KMedoids(n_clusters=2, random_state=0,
    ...                     metric='manhattan').fit(X)
    >>> kmedoids.labels_
    array([0, 0, 0, 1, 1, 1])
    >>> kmedoids.predict([[0,0], [4,4]])
    array([0, 1])
    >>> kmedoids.cluster_centers_
    array([[1, 2],
           [4, 2]])
    >>> kmedoids.inertia_
    8.0

    References
    ----------
    Kaufman, L. and Rousseeuw, P.J., Statistical Data Analysis Based on
    the L1–Norm and Related Methods, edited by Y. Dodge, North-Holland,
    405–416. 1987

    See also
    --------

    KMeans
        The KMeans algorithm minimizes the within-cluster sum-of-squares
        criterion. It scales well to large number of samples.


    Notes
    -----
        Since all pairwise distances are calculated and stored in memory for
        the duration of fit, the space complexity is O(n_samples ** 2).
    """

    def __init__(self, n_clusters=8, metric='euclidean',
                 init='heuristic', max_iter=300, random_state=None):
        self.n_clusters = n_clusters
        self.metric = metric
        self.init = init
        self.max_iter = max_iter
        self.random_state = random_state

    def _check_nonnegative_int(self, value, desc):
        """Validates if value is a valid integer > 0"""

        if (value is None or value <= 0 or
                not isinstance(value, (int, np.integer))):
            raise ValueError("%s should be a nonnegative integer. "
                             "%s was given" % (desc, value))

    def _check_init_args(self):
        """Validates the input arguments. """

        # Check n_clusters and max_iter
        self._check_nonnegative_int(self.n_clusters, "n_clusters")
        self._check_nonnegative_int(self.max_iter, "max_iter")

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
        X : {array-like, sparse matrix}, shape = (n_samples, n_features).
            Dataset to cluster.

        Returns
        -------
        self
        """

        self._check_init_args()
        X = check_array(X, accept_sparse=['csr', 'csc'])
        if self.n_clusters > X.shape[0]:
            raise ValueError("The number of medoids (%d) must be less "
                             "than the number of samples %d."
                             % (self.n_clusters, X.shape[0]))

        D = pairwise_distances(X, metric=self.metric)
        medoid_idxs = self._initialize_medoids(D, self.n_clusters)
        labels = None

        # Continue the algorithm as long as
        # the medoids keep changing and the maximum number
        # of iterations is not exceeded
        self.n_iter_ = 0
        for self.n_iter_ in range(0, self.max_iter):
            old_medoid_idxs = np.copy(medoid_idxs)
            labels = np.argmin(D[medoid_idxs, :], axis=0)

            # Update medoids with the new cluster indices
            self._update_medoid_idxs_in_place(D, labels, medoid_idxs)
            if np.all(old_medoid_idxs == medoid_idxs):
                break
            elif self.n_iter_ == self.max_iter - 1:
                warnings.warn("Maximum number of iteration reached before "
                              "convergence. Consider increasing max_iter to "
                              "improve the fit.",
                              ConvergenceWarning)

        # Expose labels_ which are the assignments of
        # the training data to clusters
        self.labels_ = labels
        self.cluster_centers_ = X[medoid_idxs]
        self.inertia_ = self._compute_inertia(self.transform(X))

        # Return self to enable method chaining
        return self

    def _update_medoid_idxs_in_place(self, D, labels, medoid_idxs):
        """In-place update of the medoid indices"""

        # Update the medoids for each cluster
        for k in range(self.n_clusters):
            # Extract the distance matrix between the data points
            # inside the cluster k
            cluster_k_idxs = np.where(labels == k)[0]

            if len(cluster_k_idxs) == 0:
                warnings.warn("Cluster %d is empty!" % k)
                continue

            in_cluster_distances = D[cluster_k_idxs, cluster_k_idxs[:, None]]

            # Calculate all costs from each point to all others in the cluster
            in_cluster_all_costs = np.sum(in_cluster_distances, axis=1)

            min_cost_idx = np.argmin(in_cluster_all_costs)
            min_cost = in_cluster_all_costs[min_cost_idx]
            curr_cost = in_cluster_all_costs[
                np.argmax(cluster_k_idxs == medoid_idxs[k])]

            # Adopt a new medoid if its distance is smaller then the current
            if min_cost < curr_cost:
                medoid_idxs[k] = cluster_k_idxs[min_cost_idx]

    def transform(self, X):
        """Transforms X to cluster-distance space.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Data to transform.

        Returns
        -------
        X_new : {array-like, sparse matrix}, shape=(n_samples, n_clusters)
            X transformed in the new space of distances to cluster centers.
        """
        X = check_array(X, accept_sparse=['csr', 'csc'])
        check_is_fitted(self, "cluster_centers_")

        return pairwise_distances(X, Y=self.cluster_centers_,
                                  metric=self.metric)

    def predict(self, X):
        """Predict the closest cluster for each sample in X

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            New data to predict.

        Returns
        -------
        labels : array, shape = (n_samples,)
            Index of the cluster each sample belongs to.
        """
        check_is_fitted(self, "cluster_centers_")
        X = check_array(X, accept_sparse=['csr', 'csc'])

        # Return data points to clusters based on which cluster assignment
        # yields the smallest distance
        return pairwise_distances_argmin(X, Y=self.cluster_centers_,
                                         metric=self.metric)

    def _compute_inertia(self, distances):
        """Compute inertia of new samples. Inertia is defined as the sum of the
        sample distances to closest cluster centers.

        Parameters
        ----------
        distances : {array-like, sparse matrix}, shape=(n_samples, n_clusters)
            Distances to cluster centers.

        Returns
        -------
        Sum of sample distances to closest cluster centers.
        """

        # Define inertia as the sum of the sample-distances
        # to closest cluster centers
        inertia = np.sum(np.min(distances, axis=1))

        return inertia

    def _initialize_medoids(self, D, n_clusters):
        """Select initial mediods randomly or heuristically"""

        if self.init == 'random':  # Random initialization
            # Pick random k medoids as the initial ones.
            medoids = self.random_state_.choice(len(D), n_clusters)
        else:  # Initialization by heuristic
            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = np.argpartition(np.sum(D, axis=1),
                                      range(0, n_clusters))[:n_clusters]

        return medoids
