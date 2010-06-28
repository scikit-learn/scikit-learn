"""Mean-shift algorithm for clustering

Author: Alexandre Gramfort alexandre.gramfort@inria.fr
"""

import numpy as np

def meanshift(X, bandwidth):
    """Perform MeanShift Clustering of data using a flat kernel

    Arguments:
    =========

    X : array [n_samples, n_features]
        Input points

    bandwidth : float
        kernel bandwidth

    Output:
    ======

    cluster_centers: array [n_clusters, n_features]

    labels : array [n_samples]
        cluster labels for each point

    References:
    ===========

    K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
    Density Function, with Applications in Pattern Recognition"

    Notes:
    =====
    See examples/plot_meanshift.py for an example.
    
    """

    n_points, n_features = X.shape

    n_clusters = 0
    bandwidth_squared = bandwidth**2
    points_idx_init = np.arange(n_points)
    stop_thresh     = 1e-3*bandwidth # when mean has converged
    cluster_centers = [] # center of clusters
    been_visited_flag = np.zeros(n_points, dtype=np.bool) # track if a points been seen already
    n_points_init = n_points # number of points to possibly use as initilization points
    cluster_votes = [] # used to resolve conflicts on cluster membership

    np.random.seed(0)

    while n_points_init:

        tmp_index   = np.random.randint(n_points_init) # pick a random seed point
        start_idx   = points_idx_init[tmp_index] # use this point as start of mean
        my_mean     = X[start_idx,:] # intilize mean to this points location
        my_members  = np.zeros(n_points, dtype=np.bool) # points that will get added to this cluster
        this_cluster_votes = np.zeros(n_points, dtype=np.uint16) # used to resolve conflicts on cluster membership

        while True: # loop until convergence

            # dist squared from mean to all points still active
            sqrt_dist_to_all = np.sum((my_mean - X)**2, axis=1)

            in_idx = sqrt_dist_to_all < bandwidth_squared # points within bandwidth
            this_cluster_votes[in_idx] += 1 # add a vote for all the in points belonging to this cluster

            my_old_mean  = my_mean # save the old mean
            my_mean      = np.mean(X[in_idx,:], axis=0) # compute the new mean
            my_members   = np.logical_or(my_members, in_idx) # add any point within bandwidth to the cluster
            been_visited_flag[my_members] = True # mark that these points have been visited

            if np.linalg.norm(my_mean-my_old_mean) < stop_thresh:

                # check for merge possibilities
                merge_with = -1
                for c in range(n_clusters):
                     # distance from possible new clust max to old clust max
                    dist_to_other = np.linalg.norm(my_mean - cluster_centers[c])
                    if dist_to_other < bandwidth/2: # if its within bandwidth/2 merge new and old
                        merge_with = c
                        break

                if merge_with >= 0: # something to merge
                    # record the max as the mean of the two merged (I know biased twoards new ones)
                    cluster_centers[merge_with] = 0.5 * (my_mean+cluster_centers[merge_with])
                    cluster_votes[merge_with] += this_cluster_votes # add these votes to the merged cluster
                else: # its a new cluster
                    n_clusters += 1 # increment clusters
                    cluster_centers.append(my_mean) # record the mean
                    cluster_votes.append(this_cluster_votes)

                break

        # we can initialize with any of the points not yet visited
        points_idx_init = np.where(been_visited_flag == False)[0]
        n_points_init   = points_idx_init.size # number of active points in set

    labels = np.argmax(cluster_votes, axis=0) # a point belongs to the cluster with the most votes

    return cluster_centers, labels

class MeanShift(object):
    """MeanShift clustering"""
    def __init__(self, bandwidth):
        super(MeanShift, self).__init__()
        self.bandwidth = bandwidth
    
    def fit(self, X):
        """compute MeanShift"""
        self.cluster_centers, self.labels = meanshift(X, self.bandwidth)
        return self
