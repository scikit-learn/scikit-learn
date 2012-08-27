"""Meanshift clustering.

Authors: Conrad Lee conradlee@gmail.com
         Alexandre Gramfort alexandre.gramfort@inria.fr
         Gael Varoquaux gael.varoquaux@normalesup.org
"""

from collections import defaultdict
import numpy as np

from ..utils import extmath, check_random_state
from ..base import BaseEstimator, ClusterMixin
from ..neighbors import NearestNeighbors


def estimate_bandwidth(X, quantile=0.3, n_samples=None, random_state=0):
    """Estimate the bandwith to use with MeanShift algorithm

    Parameters
    ----------
    X: array [n_samples, n_features]
        Input points

    quantile: float, default 0.3
        should be between [0, 1]
        0.5 means that the median is all pairwise distances is used

    n_samples: int
        The number of samples to use. If None, all samples are used.

    random_state: int or RandomState
        Pseudo number generator state used for random sampling.

    Returns
    -------
    bandwidth: float
        The bandwidth parameter
    """
    random_state = check_random_state(random_state)
    if n_samples is not None:
        idx = random_state.permutation(X.shape[0])[:n_samples]
        X = X[idx]
    nbrs = NearestNeighbors(n_neighbors=int(X.shape[0] * quantile))
    nbrs.fit(X)
    d, _ = nbrs.kneighbors(X, return_distance=True)

    bandwidth = np.mean(np.max(d, axis=1))
    return bandwidth


def mean_shift(X, bandwidth=None, seeds=None, bin_seeding=False,
               cluster_all=True, max_iterations=300):
    """Perform MeanShift Clustering of data using a flat kernel

    Seed using a binning technique for scalability.

    Parameters
    ----------

    X : array [n_samples, n_features]
        Input points

    bandwidth : float, optional
        kernel bandwidth
        If bandwidth is not defined, it is set using
        a heuristic given by the median of all pairwise distances

    seeds: array [n_seeds, n_features]
        point used as initial kernel locations

    bin_seeding: boolean
        If true, initial kernel locations are not locations of all
        points, but rather the location of the discretized version of
        points, where points are binned onto a grid whose coarseness
        corresponds to the bandwidth. Setting this option to True will speed
        up the algorithm because fewer seeds will be initialized.
        default value: False
        Ignored if seeds argument is not None

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

    """
    if bandwidth is None:
        bandwidth = estimate_bandwidth(X)
    if seeds is None:
        if bin_seeding:
            seeds = get_bin_seeds(X, bandwidth)
        else:
            seeds = X
    n_samples, n_features = X.shape
    stop_thresh = 1e-3 * bandwidth  # when mean has converged
    center_intensity_dict = {}
    nbrs = NearestNeighbors(radius=bandwidth).fit(X)

    # For each seed, climb gradient until convergence or max_iterations
    for my_mean in seeds:
        completed_iterations = 0
        while True:
            # Find mean of points within bandwidth
            i_nbrs = nbrs.radius_neighbors([my_mean], bandwidth,
                                           return_distance=False)[0]
            points_within = X[i_nbrs]
            if len(points_within) == 0:
                break  # Depending on seeding strategy this condition may occur
            my_old_mean = my_mean  # save the old mean
            my_mean = np.mean(points_within, axis=0)
            # If converged or at max_iterations, addS the cluster
            if extmath.norm(my_mean - my_old_mean) < stop_thresh or \
                   completed_iterations == max_iterations:
                center_intensity_dict[tuple(my_mean)] = len(points_within)
                break
            completed_iterations += 1

    # POST PROCESSING: remove near duplicate points
    # If the distance between two kernels is less than the bandwidth,
    # then we have to remove one because it is a duplicate. Remove the
    # one with fewer points.
    sorted_by_intensity = sorted(center_intensity_dict.items(),
                                 key=lambda tup: tup[1], reverse=True)
    sorted_centers = np.array([tup[0] for tup in sorted_by_intensity])
    unique = np.ones(len(sorted_centers), dtype=np.bool)
    nbrs = NearestNeighbors(radius=bandwidth).fit(sorted_centers)
    for i, center in enumerate(sorted_centers):
        if unique[i]:
            neighbor_idxs = nbrs.radius_neighbors([center],
                                                  return_distance=False)[0]
            unique[neighbor_idxs] = 0
            unique[i] = 1  # leave the current point as unique
    cluster_centers = sorted_centers[unique]

    # ASSIGN LABELS: a point belongs to the cluster that it is closest to
    nbrs = NearestNeighbors(n_neighbors=1).fit(cluster_centers)
    labels = np.zeros(n_samples, dtype=np.int)
    distances, idxs = nbrs.kneighbors(X)
    if cluster_all:
        labels = idxs.flatten()
    else:
        labels.fill(-1)
        bool_selector = distances.flatten() <= bandwidth
        labels[bool_selector] = idxs.flatten()[bool_selector]
    return cluster_centers, labels


def get_bin_seeds(X, bin_size, min_bin_freq=1):
    """Finds seeds for mean_shift

    Finds seeds by first binning data onto a grid whose lines are
    spaced bin_size apart, and then choosing those bins with at least
    min_bin_freq points.

    Parameters
    ----------

    X : array [n_samples, n_features]
        Input points, the same points that will be used in mean_shift

    bin_size: float
        Controls the coarseness of the binning. Smaller values lead
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

    # Bin points
    bin_sizes = defaultdict(int)
    for point in X:
        binned_point = np.cast[np.int32](point / bin_size)
        bin_sizes[tuple(binned_point)] += 1

    # Select only those bins as seeds which have enough members
    bin_seeds = np.array([point for point, freq in bin_sizes.iteritems() if \
                          freq >= min_bin_freq], dtype=np.float32)
    bin_seeds = bin_seeds * bin_size
    return bin_seeds


class MeanShift(BaseEstimator, ClusterMixin):
    """MeanShift clustering

    Parameters
    ----------
    bandwidth: float, optional
        Bandwith used in the RBF kernel
        If not set, the bandwidth is estimated.
        See clustering.estimate_bandwidth

    seeds: array [n_samples, n_features], optional
        Seeds used to initialize kernels. If not set,
        the seeds are calculated by clustering.get_bin_seeds
        with bandwidth as the grid size and default values for
        other parameters.

    cluster_all: boolean, default True
        If true, then all points are clustered, even those orphans that are
        not within any kernel. Orphans are assigned to the nearest kernel.
        If false, then orphans are given cluster label -1.

    Attributes
    ----------
    `cluster_centers_` : array, [n_clusters, n_features]
        Coordinates of cluster centers

    `labels_` :
        Labels of each point

    Notes
    -----

    Scalability:

    Because this implementation uses a flat kernel and
    a Ball Tree to look up members of each kernel, the complexity will is
    to O(T*n*log(n)) in lower dimensions, with n the number of samples
    and T the number of points. In higher dimensions the complexity will
    tend towards O(T*n^2).

    Scalability can be boosted by using fewer seeds, for examply by using
    a higher value of min_bin_freq in the get_bin_seeds function.

    Note that the estimate_bandwidth function is much less scalable than
    the mean shift algorithm and will be the bottleneck if it is used.

    References
    ----------

    Dorin Comaniciu and Peter Meer, "Mean Shift: A robust approach toward
    feature space analysis". IEEE Transactions on Pattern Analysis and
    Machine Intelligence. 2002. pp. 603-619.

    """
    def __init__(self, bandwidth=None, seeds=None, bin_seeding=False,
                 cluster_all=True):
        self.bandwidth = bandwidth
        self.seeds = seeds
        self.bin_seeding = bin_seeding
        self.cluster_all = cluster_all
        self.cluster_centers_ = None
        self.labels_ = None

    def fit(self, X):
        """ Compute MeanShift

        Parameters
        -----------
        X : array [n_samples, n_features]
            Input points
        """
        self.cluster_centers_, self.labels_ = \
                               mean_shift(X,
                                          bandwidth=self.bandwidth,
                                          seeds=self.seeds,
                                          bin_seeding=self.bin_seeding,
                                          cluster_all=self.cluster_all)
        return self
