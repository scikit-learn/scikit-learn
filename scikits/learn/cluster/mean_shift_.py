""" Algorithms for clustering : Meanshift,  Affinity propagation and spectral
clustering.

Author: Alexandre Gramfort alexandre.gramfort@inria.fr
        Gael Varoquaux gael.varoquaux@normalesup.org
"""

from math import floor
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
from itertools import izip

from ..base import BaseEstimator
from ..metrics.pairwise import euclidean_distances


def estimate_bandwidth(X, quantile=0.3):
    """Estimate the bandwith ie the radius to use with an RBF kernel
    in the MeanShift algorithm

    X: array [n_samples, n_features]
        Input points

    quantile: float, default 0.3
        should be between [0, 1]
        0.5 means that the median is all pairwise distances is used
    """
    distances = euclidean_distances(X, X)
    distances = np.triu(distances, 1)
    distances_sorted = np.sort(distances[distances > 0])
    bandwidth = distances_sorted[floor(quantile * len(distances_sorted))]
    return bandwidth


def mean_shift(X, bandwidth=None, seeds=None, cluster_all_points=True, max_iterations=300):
    """Perform MeanShift Clustering of data using a flat kernel

    Seed using a bucketing/binning/discretizing technique for scalability.
    
    Parameters
    ----------

    X : array [n_samples, n_features]
        Input points

    bandwidth : float, optional
        kernel bandwidth
        If bandwidth is not defined, it is set using
        a heuristic given by the median of all pairwise distances

    bin_size: float, optional
        The algorithm is seeded using points in a bucketed (discretized)
        version of the points. bin_size controls the coarseness of the
        discretization. If not defined, it is set to one half of the
        bandwidth (this is rather arbitrary).

    min_bin_freq: int, optional
       To speed up the algorithm, accept only those bins with at least
       min_bin_freq points as seeds. If not defined, set to 1.
       
    Returns
    -------

    cluster_centers : array [n_clusters, n_features]
        Coordinates of cluster centers

    labels : array [n_samples]
        cluster labels for each point

    Notes
    -----
    See examples/plot_meanshift.py for an example.

    Dorin Comaniciu and Peter Meer, "Mean Shift: A robust approach toward
    feature space analysis". IEEE Transactions on Pattern Analysis and
    Machine Intelligence. 2002. pp. 603-619.

    """
    if bandwidth is None:
        bandwidth = estimate_bandwidth(X)
    if seeds is None:
        seeds = get_bucket_seeds(X, bandwidth)
        
    n_points, n_features = X.shape
    stop_thresh = 0.1 * bandwidth  # when mean has converged
    cluster_center_intensity_dict = {}

    # used to efficiently look up nearest neighbors (effective in lower dimensions)
    kd_tree = cKDTree(X)

    # For each seed, climb gradient until convergence
    for my_mean in seeds:
        completed_iterations = 0
        while True:  
            # Find mean of points within bandwidth
            points_within = [X[idx] for idx in get_points_within_range(kd_tree, my_mean, bandwidth)]
            if completed_iterations == 0 and len(points_within) == 0:
                break # Depending on seeding strategy, this condition may occur
            my_old_mean = my_mean  # save the old mean
            my_mean = np.mean(points_within, axis=0)

            # If converged or at max_iterations, add the cluster ---
            # we deal with duplicates below
            if np.linalg.norm(my_mean - my_old_mean) < stop_thresh or \
                   completed_iterations == max_iterations:
                # record the point and intensity, duplicates implicitly ignored
                cluster_center_intensity_dict[tuple(my_mean)] = len(points_within)
                break
            completed_iterations += 1
            if completed_iterations == max_iterations:
                print "Max iterations reached"

    # POST PROCESSING: remove near duplicate points
    # If the distance between two kernels is less than the bandwidth,
    # then we have to remove one because it is a duplicate. Remove the
    # one with fewer points.
    sorted_by_intensity = sorted(cluster_center_intensity_dict.items() , key=lambda tup: tup[1], reverse=True)
    sorted_centers = [center for center, intensity in sorted_by_intensity]
    is_unique = np.ones(len(sorted_centers), dtype=np.bool)
    cluster_center_kd_tree = cKDTree(sorted_centers)
    for i, center in enumerate(sorted_centers):
        if is_unique[i]:
            neighbor_idxs = get_points_within_range(cluster_center_kd_tree, center, bandwidth)
            for neighbor_idx in neighbor_idxs[1:]: # skip nearest point because it is the current point
                is_unique[neighbor_idx] = 0
    cluster_centers = [center for center, unique in izip(sorted_centers, is_unique) if unique]

    # ASSIGN LABELS: a point belongs to the cluster that it is closest to
    centers_tree = cKDTree(cluster_centers)
    # Every point is assigned a label, so keep these small using 4byte ints if possible
    if len(cluster_centers) < 65535:  
        labels = np.zeros(n_points, dtype=np.uint16)
    else:
        labels = np.zeros(n_points, dtype=np.uint32)
    for point_idx in xrange(len(X)):
        if cluster_all_points:
            distance, idx = centers_tree.query(X[point_idx], 1)
            labels[point_idx] = idx
        # If not forced to cluster all points, put those that are not within a kernel in cluster -1
        else:
            distance, idx = centers_tree.query(X[point_idx], 1, distance_upper_bound=bandwidth)
            if distance <= bandwidth:
                labels[point_idx] = idx
            else:
                labels[point_idx] = -1
    return cluster_centers, labels

def get_bucket_seeds(X, bin_size, min_bin_freq=1):
    """
    Finds seeds for clustering.mean_shift by first bucketing/discretizing
    data onto a grid whose lines are spaced bin_size apart, and then
    choosing those buckets with at least min_bin_freq points.
    Parameters
    ----------

    X : array [n_samples, n_features]
        Input points, the same points that will be used in mean_shift

    bin_size: float
        Controls the coarseness of the discretization. Smaller values lead
        to more seeding (which is computationally more expensive). If you're
        not sure how to set this, set it to the value of the bandwidth used
        in clustering.mean_shift

    min_bin_freq: integer, default 1
        Only bins with at least min_bin_freq will be selected as seeds.
        Raising this value decreases the number of seeds found, which
        makes mean_shift computationally cheaper.
    Returns
    -------

    bin_seeds : array [n_samples, n_features]
        points used as initial kernel posistions in clustering.mean_shift
    """
    
    # Discretize (i.e., quantize, bin) points to bins
    bin_sizes = defaultdict(int)
    discretized_points = X.copy() / bin_size
    discretized_points = np.cast[np.int32](discretized_points)
    for discretized_point in discretized_points:
        bin_sizes[tuple(discretized_point)] += 1
        
    # Select only those bins as seeds which have enough members
    bin_seeds = np.array([point for point, freq in bin_sizes.iteritems() if \
                          freq >= min_bin_freq], dtype=np.float32)
    bin_seeds = bin_seeds * bin_size
    return bin_seeds

def get_points_within_range(kd_tree, query_point, max_distance):
    """
    This function adds a feature to the cKDTree that is missing:
    efficiently query all points within a given range of a coordinate.
    Instead of doing this, cKDTree efficiently queries the k-nearest
    points.
    
    Parameters
    ----------

    kd_tree : a cKDTree object from scipy.spatial, required
        Should contain all n-dimensional points of interest

    query_point : n-dimensional point (e.g., tuple, list, np array)
        The coordinate of interest. All points within max_distance
        of this coordinate will be returned.

    max_distance: float
        All points within max_distance of the query_point will
        be returned
       
    Returns
    -------

    cluster_centers : array [n_matchin_points]
        Contains the indices of the points used to create
        kd_tree which are within max_distance of query_point
    """
    max_neighbors = min(kd_tree.n, 100)
    # Because we don't want the kd_tree to return a list with all distances,
    # we start wtih a low value of max_neighbors and increase it until it
    # is sufficient to contain all points within max_distance
    while True:
        distances, indices = kd_tree.query(query_point, max_neighbors,
                                           distance_upper_bound=max_distance)
        max_distance_idx = distances.argmax()
        # If max distance is infinite, then we've found all neighbors within the
        # distance upper bound. Otherwise we haven't gone out far enough and need
        # to requery with more max results
        if distances[max_distance_idx] == np.inf or max_neighbors == kd_tree.n:
            return indices[:max_distance_idx]
        max_neighbors = min(max_neighbors*5, kd_tree.n)
        
##############################################################################

class MeanShift(BaseEstimator):
    """MeanShift clustering


    Parameters
    ----------

    bandwidth: float, optional
        Bandwith used in the RBF kernel
        If not set, the bandwidth is estimated.
        See clustering.estimate_bandwidth

    seeds: array [n_samples, n_features], optional
        Seeds used to initialize kernels. If not set,
        the seeds are calculated by clustering.get_bucket_seeds
        with bandwidth as the grid size and default values for
        other parameters.

    cluster_all_points: boolean, default True
        If true, then all points are clustered, even those orphans that are
        not within any kernel. Orphans are assigned to the nearest kernel.
        If false, then orphans are given cluster label -1.
    Methods
    -------

    fit(X):
        Compute MeanShift clustering

    Attributes
    ----------

    cluster_centers_: array, [n_clusters, n_features]
        Coordinates of cluster centers

    labels_:
        Labels of each point

    Notes
    -----

    Reference:

    K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a
    Density Function, with Applications in Pattern Recognition"

    The algorithmic complexity of the mean shift algorithm is O(T n^2)
    with n the number of samples and T the number of iterations. It is
    not adviced for a large number of samples.
    """

    def __init__(self, bandwidth=None, seeds=None, cluster_all_points=True):
        self.bandwidth = bandwidth
        self.seeds = seeds
        self.cluster_all_points = cluster_all_points
        
    def fit(self, X, **params):
        """ Compute MeanShift

            Parameters
            -----------
            X : array [n_samples, n_features]
                Input points

        """
        self._set_params(**params)
        self.cluster_centers_, self.labels_ = mean_shift(X,
                                                         seeds = self.seeds,
                                                         bandwidth = self.bandwidth,
                                                         cluster_all_points = self.cluster_all_points)
        return self
