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


def mean_shift(X, bandwidth=None, bin_size=None, min_bin_freq=1, cluster_all_points=True, max_per_kernel=15000, max_iterations=300):
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
    if bin_size is None:
        bin_size  = bandwidth / 2.
        
    n_points, n_features = X.shape
    n_clusters = 0
    stop_thresh = 0.1 * bandwidth  # when mean has converged
    cluster_centers = []  # center of clusters
    cluster_intensities = [] # intensities of clusters used in case of merging

    # used to efficiently look up nearest neighbors (effective in lower dimensions)
    kd_tree = cKDTree(X)
    
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

    # For each seed, climb gradient until convergence
    for my_mean in bin_seeds:
        completed_iterations = 0
        while True:  
            # Find mean of points within bandwidth
            distances, indices = kd_tree.query(my_mean, max_per_kernel,
                                               distance_upper_bound=bandwidth+0.000001)
            points_within = np.array(list(finite_point_generator(distances, indices, X, max_per_kernel)))
            if completed_iterations == 0 and len(points_within) == 0:
                break
            my_old_mean = my_mean  # save the old mean
            my_mean = np.mean(points_within, axis=0)

            # If converged or at max_iterations, check if cluster is or better
            if np.linalg.norm(my_mean - my_old_mean) < stop_thresh or \
                   completed_iterations == max_iterations:
                # check for nearby (duplicate) clusters
                merge_with = -1
                closest_distance = np.inf
                for c in range(n_clusters):
                    dist_to_other = np.linalg.norm(my_mean -
                                                        cluster_centers[c])
                    # if its within the bandwidth mark it
                    if dist_to_other <= bandwidth and dist_to_other < closest_distance:
                        merge_with = c
                        closest_distance = dist_to_other
                if merge_with >= 0:  # In case of duplicates
                    # use cluster which has higher kernel value
                    if len(points_within) > cluster_intensities[merge_with]:
                        cluster_centers[merge_with] = my_mean
                        cluster_intensities[merge_with] = len(points_within)
                else:  # If there's no duplicate than accepted it
                    n_clusters += 1  # increment clusters
                    cluster_centers.append(my_mean)  # record the mean
                    cluster_intensities.append(len(points_within))
                break
            completed_iterations += 1
            if completed_iterations == max_iterations:
                print "Max iterations reached"
    # a point belongs to the cluster that it is closest to
    centers_tree = cKDTree(cluster_centers)
    if n_clusters < 65535:
        labels = np.zeros(n_points, dtype=np.uint16)
    else:
        labels = np.zeros(n_points, dtype=np.uint32)
    for point_idx in xrange(len(X)):
        if cluster_all_points:
            distance, idx = centers_tree.query(X[point_idx], 1)
            labels[point_idx] = idx
        else:
            distance, idx = centers_tree.query(X[point_idx], 1, distance_upper_bound=bandwidth)
            if distance <= bandwidth: # make sure distance is not inf
                labels[point_idx] = idx
            else:
                labels[point_idx] = -1
    return cluster_centers, labels

def finite_point_generator(distances, indices, X, max_per_kernel):
    i = 0
    for distance, index in izip(distances, indices):
        if distance == np.inf:
            break
        yield X[index]
        i += 1
        if i == max_per_kernel:
            print "\t\tmore points should have been returned but max_per_kernel kicked in"
##############################################################################

class MeanShift(BaseEstimator):
    """MeanShift clustering


    Parameters
    ----------

    bandwidth: float, optional
        Bandwith used in the RBF kernel
        If not set, the bandwidth is estimated.
        See clustering.estimate_bandwidth

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

    def __init__(self, bandwidth=None, bin_size=None, min_bin_freq=None, cluster_all_points=True):
        self.bandwidth = bandwidth
        self.bin_size = bin_size
        self.min_bin_freq = min_bin_freq
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
                                                         bandwidth = self.bandwidth,
                                                         bin_size = self.bin_size,
                                                         min_bin_freq = self.min_bin_freq,
                                                         cluster_all_points = self.cluster_all_points)
        return self


class MeanShiftGrid(BaseEstimator):
    def __init__(self, bin_size, min_bin_freq, bandwidth=None):
        self.bandwidth = bandwidth
        self.bin_size = bin_size
        self.min_bin_freq = min_bin_freq
    def fit(self, X, **params):
        """ Compute MeanShift

            Parameters
            -----------
            X : array [n_samples, n_features]
                Input points

        """
        self._set_params(**params)
        self.cluster_centers_, self.labels_ = mean_shift_grid(X, self.bin_size,
                                                              self.min_bin_freq,
                                                              self.bandwidth)
        return self
