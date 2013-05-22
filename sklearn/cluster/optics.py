# -*- coding: utf-8 -*-

#
# Written by Shane Grigsby     ##
# Email: refuge@rocktalus.com  ##
# Date:  May 2013              ##
#


# Imports #

import sys
import scipy
from sklearn.neighbors import BallTree

# Main Class #


class setOfObjects(BallTree):

    """Build balltree data structure with processing index from given data
    in preparation for OPTICS Algorithm

    Parameters
    ----------
    data_points: array [n_samples, n_features]"""

    def __init__(self, data_points):

        super(setOfObjects, self).__init__(data_points)

        self._n = len(self.data)
        # Start all points as 'unprocessed' ##
        self._processed = scipy.zeros((self._n, 1), dtype=bool)
        self._reachability = scipy.ones(self._n) * scipy.inf
        self._core_dist = scipy.ones(self._n) * scipy.nan
        # Might be faster to use a list below? ##
        self._index = scipy.array(range(self._n))
        self._nneighbors = scipy.ones(self._n, dtype=int)
        # Start all points as noise ##
        self._cluster_id = -scipy.ones(self._n, dtype=int)
        self._is_core = scipy.ones(self._n, dtype=bool)
        # Ordering is important below... ###
        self._ordered_list = []

    # Used in prep step #
    def _set_neighborhood(self, point, epsilon):
        self._nneighbors[
            point] = self.query_radius(
                self.data[
                    point],
                epsilon,
                count_only=1)[
                    0]

    # Used in prep step #
    def _set_core_dist(self, point, MinPts):
        self._core_dist[point] = self.query(
            self.data[point], MinPts)[0][0][-1]

# Prep Method #

# Paralizeable! #


def prep_optics(SetofObjects, epsilon, MinPts):
    """Prep data set for main OPTICS loop

    Parameters
    ----------
    SetofObjects: Instantiated instance of 'setOfObjects' class
    epsilon: float or int
        Determines maximum object size that can be extracted.
        Smaller epsilons reduce run time
    MinPts: int
        The minimum number of samples in a neighborhood to be
        considered a core point

    Returns
    -------
    Modified setOfObjects tree structure"""

    for i in SetofObjects._index:
        SetofObjects._set_neighborhood(i, epsilon)
    for j in SetofObjects._index:
        if SetofObjects._nneighbors[j] >= MinPts:
            SetofObjects._set_core_dist(j, MinPts)
    print(
        'Core distances and neighborhoods prepped for ' + str(
        SetofObjects._n) + ' points.')

# Main OPTICS loop #


def build_optics(SetOfObjects, epsilon, MinPts, Output_file_name):
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
        if SetOfObjects._processed[point] == False:
            expandClusterOrder(SetOfObjects, point, epsilon,
                               MinPts, Output_file_name)

# OPTICS helper functions; these should not be public #

# NOT Paralizeable! The order that entries are written to
# the '_ordered_list' is important!


def expandClusterOrder(SetOfObjects, point, epsilon, MinPts, Output_file_name):
    if SetOfObjects._core_dist[point] <= epsilon:
        while not SetOfObjects._processed[point]:
            SetOfObjects._processed[point] = True
            SetOfObjects._ordered_list.append(point)
            # Comment following two lines to not write to a text file ##
            with open(Output_file_name, 'a') as file:
                file.write((str(point) + ', ' + str(
                    SetOfObjects._reachability[point]) + '\n'))
                # Keep following line! ##
                point = set_reach_dist(SetOfObjects, point, epsilon)
        print('Object Found!')
    else:
        SetOfObjects._processed[point] = True    # Probably not needed... #


# As above, NOT paralizable! Paralizing would allow items in
# 'unprocessed' list to switch to 'processed' ###
def set_reach_dist(SetOfObjects, point_index, epsilon):

    # Assumes that the query returns ordered (smallest distance first)
    # entries. This is the case for the balltree query...

    distances, indices = SetOfObjects.query(SetOfObjects.data[point_index],
                                            SetOfObjects._nneighbors[point_index])

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
            return sorted(zip(SetOfObjects._reachability[unprocessed], unprocessed), key=lambda reachability: reachability[0])[0][1]
        else:
            return point_index
    else:  # Not sure if this else statement is actaully needed... ##
        return point_index

# Extract DBSCAN Equivalent cluster structure ##

# Important: Epsilon prime should be less than epsilon used in OPTICS #


def ExtractDBSCAN(SetOfObjects, epsilon_prime):
    """Performs DBSCAN equivalent extraction for arbitrary epsilon.
    Can be run multiple times.

    Parameters
    ----------
    SetOfObjects: Prepped and build instance of setOfObjects
    epsilon_prime: float or int
        Must be less than or equal to what was used for prep and build steps

    Returns
    -------
    Modified setOfObjects with cluster_id and is_core attributes."""

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
        else:
            SetOfObjects._cluster_id[entry] = cluster_id
            if SetOfObjects._core_dist[entry] <= epsilon_prime:
                # One (i.e., 'True') for core points #
                SetOfObjects._is_core[entry] = 1
            else:
                # Zero (i.e., 'False') for non-core, non-noise points #
                SetOfObjects._is_core[entry] = 0


# End Algorithm #
