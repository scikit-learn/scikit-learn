# -*- coding: utf-8 -*-

#
# Author:   Shane Grigsby         #
# Email:    refuge@rocktalus.com  #
# Date:     May 2013              #
# Updated:  Nov 2014              #


# Imports #

import scipy
from ..utils import check_array, check_consistent_length
from sklearn.neighbors import BallTree
from sklearn.base import BaseEstimator, ClusterMixin


class setOfObjects(BallTree):

    """Build balltree data structure with processing index from given data
    in preparation for OPTICS Algorithm

    Parameters
    ----------
    data_points: array [n_samples, n_features]"""

    def __init__(self, data_points, **kwargs):

        super(setOfObjects, self).__init__(data_points, **kwargs)

        self._n = len(self.data)
        # Start all points as 'unprocessed' ##
        self._processed = scipy.zeros((self._n, 1), dtype=bool)
        self._reachability = scipy.ones(self._n) * scipy.inf
        self._core_dist = scipy.ones(self._n) * scipy.nan
        self._index = scipy.array(range(self._n))
        self._nneighbors = scipy.ones(self._n, dtype=int)
        # Start all points as noise ##
        self._cluster_id = -scipy.ones(self._n, dtype=int)
        self._is_core = scipy.ones(self._n, dtype=bool)
        # Ordering is important below... ###
        self._ordered_list = []


def _prep_optics(self, epsilon, MinPts):
    """Prep data set for main OPTICS loop

    Parameters
    ----------
    epsilon: float or int
        Determines maximum object size that can be extracted.
        Smaller epsilons reduce run time
    MinPts: int
        The minimum number of samples in a neighborhood to be
        considered a core point

    Returns
    -------
    Modified setOfObjects tree structure"""

    self._nneighbors = self.query_radius(self.data, r=epsilon,
                                         count_only=True)

    # Only core distance lookups for points capable of being core
    mask_idx = self._nneighbors >= MinPts 
    core_query = self.get_arrays()[0][mask_idx]
    # Check to see if that there is at least one cluster
    if len(core_query) >= 1 :
        core_dist = self.query(core_query,k=MinPts)[0][:,-1]
        self._core_dist[mask_idx] = core_dist


# Main OPTICS loop #


def _build_optics(SetOfObjects, epsilon, MinPts):
    """Builds OPTICS ordered list of clustering structure

    Parameters
    ----------
    SetofObjects: Instantiated and prepped instance of 'setOfObjects' class
    epsilon: float or int
        Determines maximum object size that can be extracted. Smaller
        epsilons reduce run time. This should be equal to epsilon
        in 'prep_optics'
    MinPts: int
        The minimum number of samples in a neighborhood to be considered a
        core point. Must be equal to MinPts used in 'prep_optics'
    Output_file_name: string
        Valid path where write access is available.
        Stores cluster structure"""

    for point in SetOfObjects._index:
        if not SetOfObjects._processed[point]:
            _expandClusterOrder(SetOfObjects, point, epsilon,
                                MinPts)

# OPTICS helper functions; these should not be public #

# Not parallelizable. The order that entries are written to
# the '_ordered_list' is important!


def _expandClusterOrder(SetOfObjects, point, epsilon, MinPts):
    if SetOfObjects._core_dist[point] <= epsilon:
        while not SetOfObjects._processed[point]:
            SetOfObjects._processed[point] = True
            SetOfObjects._ordered_list.append(point)
            point = _set_reach_dist(SetOfObjects, point, epsilon)
    else:
        SetOfObjects._processed[point] = True    # Probably not needed... #


# As above, not parallelizable. Parallelizing would allow items in
# the 'unprocessed' list to switch to 'processed'
def _set_reach_dist(SetOfObjects, point_index, epsilon):

    # Assumes that the query returns ordered (smallest distance first)
    # entries. This is the case for the balltree query...

    distances, indices = SetOfObjects.query(SetOfObjects.data[point_index],
                                            SetOfObjects._nneighbors[point_index])
    # Above line is not pep8, but is more readable than shortening it...

    # Checks to see if there more than one member in the neighborhood ##
    if scipy.iterable(distances):

        # Masking processed values ##
        unprocessed = indices[(SetOfObjects._processed[indices] < 1)[0].T]
        rdistances = scipy.maximum(
            distances[(SetOfObjects._processed[indices] < 1)[0].T],
            SetOfObjects._core_dist[point_index])
        SetOfObjects._reachability[
            unprocessed] = scipy.minimum(
                SetOfObjects._reachability[
                    unprocessed],
                rdistances)

        # Checks to see if everything is already processed;
        # if so, return control to main loop ##
        if unprocessed.size > 0:
            # Define return order based on reachability distance ###
            return unprocessed[scipy.argmin(SetOfObjects._reachability[unprocessed])]
        else:
            return point_index

# End Algorithm #

# class OPTICS(object):


class OPTICS(BaseEstimator, ClusterMixin):

    """Estimate clustering structure from vector array

    OPTICS: Ordering Points To Identify the Clustering Structure
    Equivalent to DBSCAN, finds core sample of high density and expands
    clusters from them. Unlike DBSCAN, keeps cluster hierarchy for varying
    epsilon values. Optimized for usage on large point datasets.

    Parameters
    ----------
    eps : float, optional
    The maximum distance between two samples for them to be considered
    as in the same neighborhood. This is also the largest object size
    expected within the dataset. Lower eps values can be used after
    OPTICS is run the first time, with fast returns of labels.
    min_samples : int, optional
    The number of samples in a neighborhood for a point to be considered
    as a core point.

    Attributes
    ----------
    `core_samples` : array, shape = [n_core_samples]
        Indices of core samples.
    `labels_` : array, shape = [n_samples]
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    References
    ----------
    Ankerst, Mihael, Markus M. Breunig, Hans-Peter Kriegel, and Jörg Sander.
    "OPTICS: ordering points to identify the clustering structure." ACM SIGMOD
    Record 28, no. 2 (1999): 49-60.
    """

    def __init__(self, eps=0.5, min_samples=5):  # , metric='euclidean'):
        self.eps = eps
        self.min_samples = min_samples
        # self.metric = metric
        self.processed = False

    def fit(self, X, y=None):
        """Perform OPTICS clustering

        Extracts an ordered list of points and reachability distances, and
        performs initial clustering using 'eps' distance specified at OPTICS
        object instantiation.

        Parameters
        ----------
        X : array [n_samples, n_features]"""

        #  Checks for sparse matrices
        X = check_array(X)

        tree = setOfObjects(X)  # ,self.metric)
        _prep_optics(tree, self.eps * 10.0, self.min_samples)
        _build_optics(tree, self.eps * 10.0, self.min_samples)
        self._index = tree._index[:]
        self._reachability = tree._reachability[:]
        self._core_dist = tree._core_dist[:]
        self._cluster_id = tree._cluster_id[:]
        self._is_core = tree._is_core[:]
        self._ordered_list = tree._ordered_list[:]
        _ExtractDBSCAN(self, self.eps)  # extraction needs to be < eps
        self.labels_ = self._cluster_id[:]
        self.core_samples = self._index[self._is_core[:] == True]
        self.n_clusters = max(self._cluster_id)
        self.processed = True
        return self  # self.core_samples, self.labels_

    def extract(self, epsPrime):
        """Performs DBSCAN equivalent extraction for arbitrary epsilon.
        Can be run multiple times.

        Parameters
        ----------
        epsilon_prime: float or int
        Must be less than or equal to what was used for prep and build steps

        Returns
        -------
        New core_samples and labels_ arrays. Modifies OPTICS object and stores
        core_samples and lables_ as attributes."""

        if self.processed is True:
            if epsPrime > self.eps * 10.0:
                print('Specify an epsilon smaller than ' + str(self.eps * 10))
            else:
                self.eps_prime = epsPrime
                _ExtractDBSCAN(self, epsPrime)
                self.labels_ = self._cluster_id[:]
                self.core_samples = self._index[self._is_core[:] == True]
                self.labels_ = self._cluster_id[:]
                self.n_clusters = max(self._cluster_id)
                if epsPrime > (self.eps * 1.1):
                    print("Warning, eps is close to epsPrime:")
                    print("Output may be unstable")
                return self.core_samples, self.labels_
        else:
            print("Run fit method first")

# Extract DBSCAN Equivalent cluster structure ##

# Important: Epsilon prime should be less than epsilon used in OPTICS #


def _ExtractDBSCAN(SetOfObjects, epsilon_prime):

    # Start Cluster_id at zero, incremented to '1' for first cluster
    cluster_id = 0
    for entry in SetOfObjects._ordered_list:
        if SetOfObjects._reachability[entry] > epsilon_prime:
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                cluster_id += 1
                SetOfObjects._cluster_id[entry] = cluster_id
            else:
                # This is only needed for compatibility for repeated scans.
                # -1 is Noise points
                SetOfObjects._cluster_id[entry] = -1
                SetOfObjects._is_core[entry] = 0
        else:
            SetOfObjects._cluster_id[entry] = cluster_id
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                # One (i.e., 'True') for core points #
                SetOfObjects._is_core[entry] = 1
            else:
                # Zero (i.e., 'False') for non-core, non-noise points #
                SetOfObjects._is_core[entry] = 0
