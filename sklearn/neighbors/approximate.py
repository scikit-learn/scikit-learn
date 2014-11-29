"""Approximate nearest neighbor search"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
import warnings
from .base import KNeighborsMixin, RadiusNeighborsMixin
from ..base import BaseEstimator
from ..utils.validation import check_array
from ..utils import check_random_state
from ..metrics.pairwise import pairwise_distances

from ..random_projection import GaussianRandomProjection

__all__ = ["LSHForest"]

HASH_DTYPE = '>u4'
MAX_HASH_SIZE = np.dtype(HASH_DTYPE).itemsize * 8


def _find_matching_indices(sorted_array, item, left_mask, right_mask):
    """Finds indices in sorted array of integers.

    Most significant h bits in the binary representations of the
    integers are matched with the items' most significant h bits.
    """
    left_index = np.searchsorted(sorted_array, item & left_mask)
    right_index = np.searchsorted(sorted_array, item | right_mask,
                                  side='right')
    return np.arange(left_index, right_index)


def _find_longest_prefix_match(bit_string_array, query, hash_size,
                               left_masks, right_masks):
    """Private function to find the longest prefix match for query.

    Most significant bits are considered as the prefix.
    """
    hi = hash_size
    lo = 0

    if _find_matching_indices(bit_string_array, query, left_masks[hi],
                              right_masks[hi]).shape[0] > 0:
        return hi

    while lo < hi:
        mid = (lo+hi) // 2

        k = _find_matching_indices(bit_string_array, query,
                                   left_masks[mid],
                                   right_masks[mid]).shape[0]
        if k > 0:
            lo = mid + 1
            res = mid
        else:
            hi = mid

    return res


class ProjectionToHashMixin(object):
    """Turn a transformed real-valued array into a hash"""
    @staticmethod
    def _to_hash(projected):
        if projected.shape[1] % 8 != 0:
            raise ValueError('Require reduced dimensionality to be a multiple '
                             'of 8 for hashing')
        # XXX: perhaps non-copying operation better
        out = np.packbits((projected > 0).astype(int)).view(dtype=HASH_DTYPE)
        return out.reshape(projected.shape[0], -1)

    def fit_transform(self, X, y=None):
        self.fit(X)
        return self.transform(X)

    def transform(self, X, y=None):
        return self._to_hash(super(ProjectionToHashMixin, self).transform(X))


class GaussianRandomProjectionHash(ProjectionToHashMixin,
                                   GaussianRandomProjection):
    """Use GaussianRandomProjection to produce a cosine LSH fingerprint"""
    def __init__(self,
                 n_components=8,
                 random_state=None):
        super(GaussianRandomProjectionHash, self).__init__(
            n_components=n_components,
            random_state=random_state)


def _array_of_arrays(list_of_arrays):
    """Creates an array of array from list of arrays."""
    out = np.empty(len(list_of_arrays), dtype=object)
    out[:] = list_of_arrays
    return out


class LSHForest(BaseEstimator, KNeighborsMixin, RadiusNeighborsMixin):
    """Performs approximate nearest neighbor search using LSH forest.

    LSH Forest: Locality Sensitive Hashing forest [1] is an alternative
    method for vanilla approximate nearest neighbor search methods.
    LSH forest data structure has been implemented using sorted
    arrays and binary search and 32 bit fixed-length hashes.
    Random projection is used as the hash family which approximates
    cosine distance.

    Parameters
    ----------

    n_estimators : int (default = 10)
        Number of trees in the LSH Forest.

    min_hash_match : int (default = 4)
        lowest hash length to be searched when candidate selection is
        performed for nearest neighbors.

    n_candidates : int (default = 10)
        Minimum number of candidates evaluated per estimator, assuming enough
        items meet the `min_hash_match` constraint.

    n_neighbors : int (default = 5)
        Number of neighbors to be returned from query function when
        it is not provided to the :meth:`kneighbors` method.

    radius : float, optinal (default = 1.0)
        Radius from the data point to its neighbors. This is the parameter
        space to use by default for the :meth`radius_neighbors` queries.

    radius_cutoff_ratio : float, optional (defualt = 0.9)
        A value ranges from 0 to 1. Radius neighbors will be searched until
        the ratio between total neighbors within the radius and the total
        candidates becomes less than this value unless it is terminated by
        hash length reaching `min_hash_match`.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------

    hash_functions_ : list of GaussianRandomProjectionHash objects
        Hash function g(p,x) for a tree is an array of 32 randomly generated
        float arrays with the same dimenstion as the data set. This array is
        stored in GaussianRandomProjectionHash object and can be obtained
        from ``components_`` attribute.

    trees_ : array, shape (n_estimators, n_samples)
        Each tree (corresponding to a hash function) contains an array of
        sorted hashed values. The array representation may change in future
        versions.

    original_indices_ : array, shape (n_estimators, n_samples)
        Original indices of sorted hashed values in the fitted index.

    References
    ----------

    .. [1] M. Bawa, T. Condie and P. Ganesan, "LSH Forest: Self-Tuning
           Indexes for Similarity Search", WWW '05 Proceedings of the
           14th international conference on World Wide Web,  651-660,
           2005.

    Examples
    --------
      >>> import numpy as np
      >>> from sklearn.neighbors import LSHForest

      >>> X_train = [[5, 5, 2], [21, 5, 5], [1, 1, 1], [8, 9, 1], [6, 10, 2]]
      >>> X_test = [[9, 1, 6], [3, 1, 10], [7, 10, 3]]
      >>> lshf = LSHForest()
      >>> lshf.fit(X_train)
      LSHForest(min_hash_match=4, n_candidates=50, n_estimators=10, n_neighbors=5,
           radius=1.0, radius_cutoff_ratio=0.9, random_state=None)
      >>> distances, indices = lshf.kneighbors(X_test, n_neighbors=2)
      >>> distances                                        # doctest: +ELLIPSIS
      array([[ 0.069...,  0.149...],
             [ 0.229...,  0.481...],
             [ 0.004...,  0.014...]])
      >>> indices
      array([[1, 2],
             [2, 0],
             [4, 0]])

    """

    def __init__(self, n_estimators=10, radius=1.0, n_candidates=50,
                 n_neighbors=5, min_hash_match=4, radius_cutoff_ratio=.9,
                 random_state=None):
        self.n_estimators = n_estimators
        self.radius = radius
        self.random_state = random_state
        self.n_candidates = n_candidates
        self.n_neighbors = n_neighbors
        self.min_hash_match = min_hash_match
        self.radius_cutoff_ratio = radius_cutoff_ratio

    def _compute_distances(self, query, candidates):
        """Computes the cosine distance.

        Distance is from the query to points in the candidates array.
        Returns argsort of distances in the candidates
        array and sorted distances.
        """
        distances = pairwise_distances(query, self._fit_X[candidates],
                                       metric='cosine')[0]
        distance_positions = np.argsort(distances)
        return distance_positions, distances[distance_positions]

    def _generate_masks(self):
        """Creates left and right masks for all hash lengths."""
        tri_size = MAX_HASH_SIZE + 1
        # Called once on fitting, output is independent of hashes
        left_mask = np.tril(np.ones((tri_size, tri_size), dtype=int))[:, 1:]
        right_mask = left_mask[::-1, ::-1]

        self._left_mask = np.packbits(left_mask).view(dtype=HASH_DTYPE)
        self._right_mask = np.packbits(right_mask).view(dtype=HASH_DTYPE)

    def _get_candidates(self, query, max_depth, bin_queries, n_neighbors):
        """Performs the Synchronous ascending phase.

        Returns an array of candidates, their distance ranks and
        distances.
        """
        index_size = self._fit_X.shape[0]
        candidates = []
        min_candidates = self.n_candidates * self.n_estimators
        while max_depth > self.min_hash_match and (len(candidates)
                                                   < min_candidates or
                                                   len(set(candidates))
                                                   < n_neighbors):

            for i in range(self.n_estimators):
                candidates.extend(
                    self.original_indices_[i][_find_matching_indices(
                        self.trees_[i],
                        bin_queries[i],
                        self._left_mask[max_depth],
                        self._right_mask[max_depth])].tolist())
            max_depth = max_depth - 1

        candidates = np.unique(candidates)
        # For insufficient candidates, candidates are filled.
        # Candidates are filled from unselected indices uniformly.
        if candidates.shape[0] < n_neighbors:
            warnings.warn(
                "Number of candidates is not sufficient to retrieve"
                " %i neighbors with"
                " min_hash_match = %i. Candidates are filled up"
                " uniformly from unselected"
                " indices." % (n_neighbors, self.min_hash_match))
            remaining = np.setdiff1d(np.arange(0, index_size), candidates)
            to_fill = n_neighbors - candidates.shape[0]
            candidates = np.concatenate((candidates, remaining[:to_fill]))

        ranks, distances = self._compute_distances(query,
                                                   candidates.astype(int))

        return (candidates[ranks[:n_neighbors]],
                distances[:n_neighbors])

    def _get_radius_neighbors(self, query, max_depth, bin_queries, radius):
        """Finds radius neighbors from the candidates obtained.

        Their distances from query are smaller than radius.
        Returns radius neighbors and distances.
        """
        ratio_within_radius = 1
        threshold = 1 - self.radius_cutoff_ratio
        total_candidates = np.array([], dtype=int)
        total_neighbors = np.array([], dtype=int)
        total_distances = np.array([], dtype=float)

        while max_depth > self.min_hash_match and (ratio_within_radius
                                                   > threshold):
            candidates = []
            for i in range(self.n_estimators):
                candidates.extend(
                    self.original_indices_[i][_find_matching_indices(
                        self.trees_[i],
                        bin_queries[i],
                        self._left_mask[max_depth],
                        self._right_mask[max_depth])].tolist())
            candidates = np.setdiff1d(candidates, total_candidates)
            total_candidates = np.append(total_candidates, candidates)
            ranks, distances = self._compute_distances(query, candidates)
            m = np.searchsorted(distances, radius, side='right')
            positions = np.searchsorted(total_distances, distances[:m])
            total_neighbors = np.insert(total_neighbors, positions,
                                        candidates[ranks[:m]])
            total_distances = np.insert(total_distances, positions,
                                        distances[:m])
            ratio_within_radius = (total_neighbors.shape[0] /
                                   float(total_candidates.shape[0]))
            max_depth = max_depth - 1
        return total_neighbors, total_distances

    def fit(self, X):
        """Fit the LSH forest on the data.

        This creates binary hashes of input data points by getting the
        dot product of input points and hash_function then
        transforming the projection into a binary string array based
        on the sign (positive/negative) of the projection.
        A sorted array of binary hashes is created.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        self : object
            Returns self.
        """

        self._fit_X = check_array(X)

        # Creates a g(p,x) for each tree
        self.hash_functions_ = []
        self.trees_ = []
        self.original_indices_ = []

        rng = check_random_state(self.random_state)
        int_max = np.iinfo(np.int32).max

        for i in range(self.n_estimators):
            # This is g(p,x) for a particular tree.
            # Builds a single tree. Hashing is done on an array of data points.
            # `GaussianRandomProjection` is used for hashing.
            # `n_components=hash size and n_features=n_dim.
            hasher = GaussianRandomProjectionHash(MAX_HASH_SIZE,
                                                  rng.randint(0, int_max))
            hashes = hasher.fit_transform(self._fit_X)[:, 0]
            original_index = np.argsort(hashes)
            bin_hashes = hashes[original_index]
            self.original_indices_.append(original_index)
            self.trees_.append(bin_hashes)
            self.hash_functions_.append(hasher)

        self._generate_masks()

        return self

    def _query(self, query, n_neighbors=None, radius=None, is_radius=False):
        """Finds neighbors and distances.

        Returns the neighbors whose distances from the query is less
        than radius if is_radius is True.
        Otherwise returns m number of neighbors and the distances
        for a given query.
        """
        bin_queries = []
        query = query.reshape((1, query.shape[0]))

        # descend phase
        max_depth = 0
        for i in range(self.n_estimators):
            bin_query = self.hash_functions_[i].transform(query)[0][0]
            k = _find_longest_prefix_match(self.trees_[i], bin_query,
                                           MAX_HASH_SIZE,
                                           self._left_mask,
                                           self._right_mask)
            if k > max_depth:
                max_depth = k
            bin_queries.append(bin_query)

        if is_radius:
            return self._get_radius_neighbors(query, max_depth,
                                              bin_queries, radius)

        else:
            return self._get_candidates(query, max_depth,
                                        bin_queries, n_neighbors)

    def kneighbors(self, X, n_neighbors=None, return_distance=True):
        """
        Returns the n_number of approximated nearest neighbors

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single query.

        n_neighbors : int, opitonal (default = None)
            Number of neighbors required. If not provided, this will
            return the number specified at the initialization.

        return_distance : boolean, optional (default = False)
            Returns the distances of neighbors if set to True.

        Returns
        -------
        dist : array, shape (n_samples, n_neighbors)
            Array representing the cosine distances to each point,
            only present if return_distance=True.

        ind : array, shape (n_samples, n_neighbors)
            Indices of the approximate nearest points in the population
            matrix.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        X = check_array(X)

        neighbors, distances = [], []
        for i in range(X.shape[0]):
            neighs, dists = self._query(X[i], n_neighbors)
            neighbors.append(neighs)
            distances.append(dists)

        if return_distance:
            return np.array(distances), np.array(neighbors)
        else:
            return np.array(neighbors)

    def radius_neighbors(self, X, radius=None, return_distance=True):
        """
        Returns the approximated nearest neighbors within the radius

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single query.

        radius : float
            Limiting distance of neighbors to return.
            (default is the value passed to the constructor).

        return_distance : boolean, optional (default = False)
            Returns the distances of neighbors if set to True.

        Returns
        -------
        dist : array, shape (n_samples,) of arrays
            Array representing the cosine distances to each point,
            only present if return_distance=True.

        ind : array, shape (n_samples,) of arrays
            An array of arrays of indices of the approximated nearest points
            with in the `radius` to the query in the population matrix.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if radius is None:
            radius = self.radius

        X = check_array(X)

        neighbors, distances = [], []
        for i in range(X.shape[0]):
            neighs, dists = self._query(X[i], radius=radius,
                                        is_radius=True)
            neighbors.append(neighs)
            distances.append(dists)

        if return_distance:
            return _array_of_arrays(distances), _array_of_arrays(neighbors)
        else:
            return _array_of_arrays(neighbors)

    def partial_fit(self, X):
        """
        Inserts new data into the already fitted LSH Forest.
        Cost is proportional to new total size, so additions
        should be batched.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            New data point to be inserted into the LSH Forest.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted before"
                             " inserting.")

        X = check_array(X)

        if X.shape[1] != self._fit_X.shape[1]:
            raise ValueError("Number of features in X and"
                             " fitted array does not match.")
        n_samples = X.shape[0]
        input_array_size = self._fit_X.shape[0]

        for i in range(self.n_estimators):
            bin_X = self.hash_functions_[i].transform(X)[:, 0]
            # gets the position to be added in the tree.
            positions = self.trees_[i].searchsorted(bin_X)
            # adds the hashed value into the tree.
            self.trees_[i] = np.insert(self.trees_[i],
                                       positions, bin_X)
            # add the entry into the original_indices_.
            self.original_indices_[i] = np.insert(self.original_indices_[i],
                                                  positions,
                                                  np.arange(input_array_size,
                                                            input_array_size +
                                                            n_samples))

        # adds the entry into the input_array.
        self._fit_X = np.row_stack((self._fit_X, X))
