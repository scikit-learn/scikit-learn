
import numpy as np
import warnings

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS

from exceptions import ValueError

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

    init: {'random', 'heuristic', or ndarray}, optional, default: 'heuristic'
      Specify medoid initialization.
    """

    # Supported clustering methods
    CLUSTERING_METHODS = ['pam']

    # Supported initialization methods
    INIT_METHODS = ['random','heuristic']
    
    def __init__(self, n_clusters = 8, distance_metric='euclidean', 
                 clustering_method='pam', init='heuristic'):

        # Check n_clusters
        if n_clusters <= 0 or not isinstance(n_clusters, int):
            raise ValueError("n_clusters has to be nonnegative integer")
        
        # Check distance_metric
        if callable(distance_metric):
            self.__distance_metric = distance_metric
        elif distance_metric in PAIRWISE_DISTANCE_FUNCTIONS:
            self.__distance_metric = PAIRWISE_DISTANCE_FUNCTIONS[distance_metric]
        else:
            raise ValueError("distance_metric needs to be callable or one of the following strings: " + 
                             "{}".format(PAIRWISE_DISTANCE_FUNCTIONS.keys()))

        # Check clustering_method
        if clustering_method not in self.CLUSTERING_METHODS:
            raise ValueError("clustering must be one of the following: {}".format(self.CLUSTERING_METHODS))

        # Check init
        if init not in self.INIT_METHODS:
            raise ValueError("init needs to be one of the following: {}".format(self.INIT_METHODS))

        self.n_clusters = n_clusters
        
        self.init = init

        self.clustering_method = clustering_method
        

    def fit(self, X):

        D = self.__distance_metric(X)

        if self.n_clusters > D.shape[0]:
            raise ValueError("The number of medoids ({}) ".format(self.n_clusters) +
                             "must be larger than the sides " +
                             "of the distance matrix ({})".format(D.shape[0]))

        medoidIcs = self.__get_initial_medoid_indices(D,self.n_clusters)

        # Old medoids will be stored here for reference
        oldMedoidIcs = np.zeros((self.n_clusters,))

        # Continue the algorithm as long as
        # the medoids keep changing
        while not np.all(oldMedoidIcs == medoidIcs):

            # Keep a copy of the old medoid assignments
            oldMedoidIcs = np.copy(medoidIcs)
            
            # Assign data points to clusters based on
            # which cluster assignment yields
            # the smallest distance
            clusters = np.argmin(D[medoidIcs, :], axis=0)
            
            # Update the medoids for each cluster
            for c in xrange(self.n_clusters):
                
                if sum(clusters == c) == 0:
                    warnings.warn("Cluster {} is empty!".format(c))
                    continue

                # Find current cost that is associated with cluster c.
                # Cost is the sum of the distance from the cluster
                # members to the medoid.
                currCost = np.sum(D[medoidIcs[c], clusters == c])
                
                # Extract the distance matrix between the data points
                # inside the cluster c
                Din = D[clusters == c, :]
                Din = Din[:, clusters == c]
                
                # Calculate all costs there exists between all
                # the data points in the cluster c
                allCosts = np.sum(Din, axis=1)
                
                # Find the index for the smallest cost in cluster c
                minCostIdx = np.argmin(allCosts)
                
                # find the value of the minimum cost in cluster c
                minCost = allCosts[minCostIdx]

                # If the minimum cost is smaller than that
                # exhibited by the currently used medoid,
                # we switch to using the new medoid in cluster c
                if minCost < currCost:

                    # Find data points that belong to cluster c,
                    # and assign the newly found medoid as the medoid
                    # for cluster c
                    medoidIcs[c] = np.where(clusters == c)[0][minCostIdx]

            # Expose labels_ which are the assignments of
            # the training data to clusters
            self.labels_ = clusters

            # Expose cluster centers, i.e. medoids
            self.medoids_ = X.take(medoidIcs,axis=0)

        # Return self to enable method chaining
        return self

    
    def transform(self, X):

        print self.medoids_.shape
        
        return self.__distance_metric(X, Y=self.medoids_)


    def predict(self, X):
        
        D = self.__distance_metric(X, Y=self.medoids_)

        # Assign data points to clusters based on
        # which cluster assignment yields
        # the smallest distance
        labels = np.argmin(D, axis=1)

        return labels


    def __get_initial_medoid_indices(self, D, n_clusters):

        if self.init == 'random':

            # Pick random k medoids as the initial ones.
            medoids = np.random.permutation(D.shape[0])[:n_clusters]

        elif self.init == 'heuristic':

            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = map(lambda x: x[0],
                          sorted(enumerate(np.sum(D, axis=1)),
                                 key=lambda x: x[1]))[:n_clusters]

        else:
            
            raise ValueError("Initialization not implemented for method: '{}'".format(self.init))

        return medoids
    
    
