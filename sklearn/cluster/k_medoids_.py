
import numpy as np
import warnings

from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS

from exceptions import ValueError

class KMedoids(object):
    """
    k-medoids class.

    Parameters
    ----------
    n_clusters: int, optional, default: 8
      How many medoids. Must be positive.

    distance_metric: string, optional, default: 'euclidean'
      What distance metric to use.

    clustering: string, optional, default: 'pam'
      What clustering mode to use.

    random_init: bool, optional, default: False
      Use random medoid initialization.
    """

    def __init__(self, n_clusters = 8, distance_metric='euclidean', 
                 clustering='pam', random_init=False):

        if distance_metric not in PAIRWISE_DISTANCE_FUNCTIONS.keys():
            raise ValueError("distance_metric needs to be one of the following: " + 
                             "{}".format(PAIRWISE_DISTANCE_FUNCTIONS.keys()))

        self.dist_func = PAIRWISE_DISTANCE_FUNCTIONS[distance_metric]
        
        if clustering != 'pam':
            raise ValueError("clustering must be 'pam'")

        if n_clusters <= 0:
            raise ValueError("n_clusters has to be nonnegative integer")

        self.__n_clusters = n_clusters
        
        self.__random_init = random_init
        

    def fit(self, X):

        D = self.dist_func(X)

        if self.__n_clusters >= D.shape[0]:
            raise ValueError("The number of medoids ({}) ".format(self.__n_clusters) +
                             "must be larger than the sides " +
                             "of the distance matrix ({})".format(D.shape[0]))

        medoids = self.__get_initial_medoids(D,self.__n_clusters)

        # Old medoids will be stored here for reference
        oldMedoids = np.zeros((self.__n_clusters,))

        # Continue the algorithm as long as
        # the medoids keep changing
        while not np.all(oldMedoids == medoids):

            # Keep a copy of the old medoid assignments
            oldMedoids = np.copy(medoids)
            
            # Assign data points to clusters based on
            # which cluster assignment yields
            # the smallest distance
            clusters = np.argmin(D[medoids, :], axis=0)
            
            # Update the medoids for each cluster
            for c in xrange(k):
                
                if sum(clusters == c) == 0:
                    warnings.warn("Cluster {} is empty!".format(c))
                    continue

                # Find current cost that is associated with cluster c.
                # Cost is the sum of the distance from the cluster
                # members to the medoid.
                currCost = np.sum(D[medoids[c], clusters == c])
                
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
                    medoids[c] = np.where(clusters == c)[0][minCostIdx]

        return clusters, medoids

    def __get_initial_medoids(self,D,n_clusters):

        if self.__random_init:

            # Pick random k medoids as the initial ones.
            medoids = np.random.permutation(D.shape[0])[:n_clusters]

        else:

            # Pick K first data points that have the smallest sum distance
            # to every other point. These are the initial medoids.
            medoids = map(lambda x: x[0],
                          sorted(enumerate(np.sum(D, axis=1)),
                                 key=lambda x: x[1]))[:n_clusters]
            
        return medoids
            
