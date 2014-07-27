"""
Locality Sensitive Hashing Forest for Approximate Nearest Neighbor Search
-------------------------------------------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
import itertools
from ..base import BaseEstimator
from ..utils.validation import safe_asarray
from ..utils import check_random_state

from ..random_projection import GaussianRandomProjection

__all__ = ["LSHForest"]


def _find_matching_indices(sorted_array, item, left_mask, right_mask):
    """
    Finds indices in sorted array of strings where their first
    h elements match the items' first h elements
    """
    left_index = np.searchsorted(sorted_array,
                                 item & left_mask)
    right_index = np.searchsorted(sorted_array,
                                  item | right_mask,
                                  side='right')
    return np.arange(left_index, right_index)


def _find_longest_prefix_match(bit_string_array, query, hash_size,
                               left_masks, right_masks):
    """
    Private function to find the longest prefix match for query
    in the bit_string_array
    """
    hi = hash_size
    lo = 0

    if _find_matching_indices(bit_string_array, query, left_masks[hi],
                              right_masks[hi]).shape[0] > 0:
        return hi

    while lo < hi:
        mid = (lo+hi)//2

        k = _find_matching_indices(bit_string_array, query,
                                   left_masks[mid],
                                   right_masks[mid]).shape[0]
        if k > 0:
            lo = mid + 1
            res = mid
        else:
            hi = mid

    return res


def _simple_euclidean_distance(query, candidates):
    """
    Private function to calculate Euclidean distances between each
    point in candidates and query
    """
    distances = np.zeros(candidates.shape[0])
    for i in range(candidates.shape[0]):
        distances[i] = np.linalg.norm(candidates[i] - query)
    return distances


class LSHForest(BaseEstimator):

    """
    Performs approximate nearest neighbor search using LSH forest.

    LSH Forest: Locality Sensitive Hashing forest [1] is an alternative
    method for vanilla approximate nearest neighbor search methods.
    LSH forest data structure has been implemented using sorted
    arrays and binary search.

    Parameters
    ----------

    max_label_length: int, optional(default = 32)
        Maximum length of a binary hash

    n_trees: int, optional (default = 10)
        Number of trees in the LSH Forest.

    c: int, optional(default = 10)
        Threshold value to select candidates for nearest neighbors.
        Number of candidates is often greater than c*n_trees(unless
        restricted by lower_bound)

    n_neighbors: int, optional(default = 1)
        Number of neighbors to be returned from query funcitond when
        it is not provided with the query.

    lower_bound: int, optional(defualt = 4)
        lowerest hash length to be searched when candidate selection is
        performed for nearest neighbors.

    radius : float, optinal(default = 1.0)
        Range of parameter space to use by default for :meth`radius_neighbors`
        queries.

    radius_cutoff_ratio: float, optional(defualt = 0.9)
        Cut off ratio of radius neighbors to candidates at the radius
        neighbor search

    random_state: float, optional(default = 1)
        A random value to initialize random number generator.

    Attributes
    ----------

    `hash_functions`: list of arrays
        Randomly generated LS hash function g(p,x) for each tree.


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

      >>> X = np.logspace(0, 3, num=5000)
      >>> X = X.reshape((100,50))
      >>> lshf = LSHForest()
      >>> lshf.fit(X)
      LSHForest(c=50, lower_bound=4, max_label_length=32, n_neighbors=1, n_trees=10,
           radius=1.0, radius_cutoff_ratio=0.9, random_state=None)

      >>> lshf.kneighbors(X[:5], n_neighbors=3, return_distance=True)
      (array([[0, 1, 2],
             [1, 0, 2],
             [2, 1, 3],
             [3, 2, 4],
             [4, 3, 5]]), array([[ 0.        ,  0.52344831,  1.08434102],
             [ 0.        ,  0.52344831,  0.56089272],
             [ 0.        ,  0.56089272,  0.60101568],
             [ 0.        ,  0.60101568,  0.6440088 ],
             [ 0.        ,  0.6440088 ,  0.6900774 ]]))

    """

    def __init__(self, max_label_length=32, n_trees=10,
                 radius=1.0, c=50, n_neighbors=1,
                 lower_bound=4, radius_cutoff_ratio=.9,
                 random_state=None):
        self.max_label_length = int(max_label_length/2*2)
        self.n_trees = n_trees
        self.radius = radius
        self.random_state = random_state
        self.c = c
        self.n_neighbors = n_neighbors
        self.lower_bound = lower_bound
        self.radius_cutoff_ratio = radius_cutoff_ratio

    def _generate_hash_function(self):
        """
        Fits a `GaussianRandomProjections` with `n_components=hash_size
        and n_features=n_dim.
        """
        random_state = check_random_state(self.random_state)
        grp = GaussianRandomProjection(n_components=self.max_label_length,
                                       random_state=random_state.randint(0,
                                                                         10))
        X = np.zeros((2, self._n_dim), dtype=float)
        grp.fit(X)
        return grp

    def _do_hash(self, input_array=None):
        """
        Does hashing on an array of data points.
        This creates a binary hash by getting the dot product of
        input_point and hash_function then transforming the projection
        into a binary string array based on the sign(positive/negative)
        of the projection.

        Parameters
        ----------

        input_array: array_like, shape (n_samples, n_features)
            A matrix of dimensions (n_samples, n_features), which is being
            hashed.
        """
        if input_array is None:
            raise ValueError("input_array cannot be None.")

        grp = self._generate_hash_function()
        res = np.array(grp.transform(input_array) > 0, dtype=int)

        return res, grp.components_

    def _create_tree(self):
        """
        Builds a single tree (in this case creates a sorted array of
        binary hashes).
        """
        hashes, hash_function = self._do_hash(self._input_array)
        binary_hashes = []
        for i in range(hashes.shape[0]):
            xx = tuple(hashes[i])
            binary_hashes.append(self.cache[xx[:self.cache_N]] * self.k
                                 + self.cache[xx[self.cache_N:]])

        return np.argsort(binary_hashes), np.sort(binary_hashes), hash_function

    def _compute_distances(self, query, candidates):
        """
        Computes the Euclidean distance from the query
        to points in the candidates array.
        Returns argsort of distances in the candidates
        array and sorted distances.
        """
        distances = _simple_euclidean_distance(
            query, self._input_array[candidates])
        return np.argsort(distances), np.sort(distances)

    def _generate_masks(self):
        """
        Creates left and right masks for all hash lengths
        """
        self._left_mask, self._right_mask = [], []

        for length in range(self.max_label_length+1):
            left_mask = int("".join(['1' for i in range(length)])
                            + "".join(['0' for i in
                                       range(self.max_label_length-length)]),
                            2)
            self._left_mask.append(left_mask)
            right_mask = int("".join(['0' for i in range(length)])
                             + "".join(['1' for i in
                                        range(self.max_label_length-length)]),
                             2)
            self._right_mask.append(right_mask)

        self._left_mask = np.array(self._left_mask)
        self._right_mask = np.array(self._right_mask)

    def _get_candidates(self, query, max_depth, bin_queries, m):
        """
        Performs the Synchronous ascending phase in the LSH Forest
        paper.
        Returns an array of candidates, their distance rancks and
        distances.
        """
        candidates = []
        n_candidates = self.c * self.n_trees
        while max_depth > self.lower_bound and (len(candidates) < n_candidates
                                                or len(set(candidates)) < m):
            for i in range(self.n_trees):
                candidates.extend(
                    self._original_indices[i, _find_matching_indices(
                        self._trees[i],
                        bin_queries[i],
                        self._left_mask[max_depth],
                        self._right_mask[max_depth])].tolist())
            max_depth = max_depth - 1
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)

        return candidates, ranks, distances

    def _get_radius_neighbors(self, query, max_depth, bin_queries, radius):
        """
        Finds neighbors of which the distances from query are smaller than
        radius, from the candidates obtained.
        Returns radius neighbors and distances.
        """
        ratio_within_radius = 1
        threshold = 1 - self.radius_cutoff_ratio
        total_candidates = np.array([], dtype=int)
        total_neighbors = np.array([], dtype=int)
        total_distances = np.array([], dtype=float)

        while max_depth > self.lower_bound and ratio_within_radius > threshold:
            candidates = []
            for i in range(self.n_trees):
                candidates.extend(
                    self._original_indices[i, _find_matching_indices(
                        self._trees[i],
                        bin_queries[i],
                        self._left_mask[max_depth],
                        self._right_mask[max_depth])].tolist())
            candidates = np.setdiff1d(candidates, total_candidates)
            total_candidates = np.append(total_candidates, candidates)
            ranks, distances = self._compute_distances(query, candidates)
            m = np.searchsorted(distances, radius, side='right')
            total_neighbors = np.append(total_neighbors,
                                        candidates[ranks[:m]])
            total_distances = np.append(total_distances, distances[:m])
            ratio_within_radius = (total_neighbors.shape[0] /
                                   float(total_candidates.shape[0]))
            max_depth = max_depth - 1
        return total_neighbors, total_distances

    def fit(self, X=None):
        """
        Fit the LSH forest on the data.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.
       """
        if X is None:
            raise ValueError("X cannot be None")

        self._input_array = safe_asarray(X)
        self._n_dim = self._input_array.shape[1]

        digits = ['0', '1']
        # Creates a g(p,x) for each tree
        self.hash_functions_ = []
        self._trees = []
        self._original_indices = []

        self.cache_N = int(self.max_label_length/2)
        hashes = [x for x in itertools.product((0, 1),
                                               repeat=self.cache_N)]

        self.cache = {}
        for item in hashes:
            self.cache[tuple(item)] = int("".join([digits[y] for y in item]),
                                          2)

        self.k = 2 ** self.cache_N

        for i in range(self.n_trees):
            # This is g(p,x) for a particular tree.
            original_index, bin_hashes, hash_function = self._create_tree()
            self._original_indices.append(original_index)
            self._trees.append(bin_hashes)
            self.hash_functions_.append(hash_function)

        self.hash_functions_ = np.array(self.hash_functions_)
        self._trees = np.array(self._trees)
        self._original_indices = np.array(self._original_indices)
        self._generate_masks()

        return self

    def _query(self, query, m=None, radius=None, is_radius=False):
        """
        Returns the neighbors whose distances from the query is less
        than radius if is_radius is True.
        Otherwise returns m number of neighbors and the distances
        for a given query.
        """
        bin_queries = []

        # descend phase
        max_depth = 0
        for i in range(self.n_trees):
            projections = np.array(np.dot(self.hash_functions_[i],
                                          query) > 0, dtype=int)
            xx = tuple(projections)
            bin_query = (self.cache[xx[:self.cache_N]] * self.k +
                         self.cache[xx[self.cache_N:]])
            k = _find_longest_prefix_match(self._trees[i], bin_query,
                                           self.max_label_length,
                                           self._left_mask,
                                           self._right_mask)
            if k > max_depth:
                max_depth = k
            bin_queries.append(bin_query)

        if not is_radius:
            candidates, ranks, distances = self._get_candidates(query,
                                                                max_depth,
                                                                bin_queries,
                                                                m)

            return candidates[ranks[:m]], distances[:m]
        else:
            return self._get_radius_neighbors(query, max_depth,
                                              bin_queries, radius)

    def kneighbors(self, X, n_neighbors=None, return_distance=False):
        """
        Returns the n_number of approximated nearest neighbors

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single query.

        n_neighbors: int, opitonal (default = None)
            Number of neighbors required. If not provided, this will
            return the number specified at the initialization.

        return_distance: boolean, optional (default = False)
            Returns the distances of neighbors if set to True.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if X is None:
            raise ValueError("X cannot be None.")

        if n_neighbors is not None:
            self.n_neighbors = n_neighbors

        X = safe_asarray(X)
        x_dim = X.ndim

        if x_dim == 1:
            neighbors, distances = self._query(X, self.n_neighbors)
            if return_distance:
                return np.array([neighbors]), np.array([distances])
            else:
                return np.array([neighbors])
        else:
            neighbors, distances = [], []
            for i in range(X.shape[0]):
                neighs, dists = self._query(X[i], self.n_neighbors)
                neighbors.append(neighs)
                distances.append(dists)

            if return_distance:
                return np.array(neighbors), np.array(distances)
            else:
                return np.array(neighbors)

    def radius_neighbors(self, X, radius=None, return_distance=False):
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

        return_distance: boolean, optional (default = False)
            Returns the distances of neighbors if set to True.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if X is None:
            raise ValueError("X cannot be None.")

        if radius is not None:
            self.radius = radius

        X = safe_asarray(X)
        x_dim = X.ndim

        if x_dim == 1:
            neighbors, distances = self._query(X, radius=self.radius,
                                               is_radius=True)
            if return_distance:
                return np.array([neighbors]), np.array([distances])
            else:
                return np.array([neighbors])
        else:
            neighbors, distances = [], []
            for i in range(X.shape[0]):
                neighs, dists = self._query(X[i], radius=self.radius,
                                            is_radius=True)
                neighbors.append(neighs)
                distances.append(dists)

            if return_distance:
                return np.array(neighbors), np.array(distances)
            else:
                return np.array(neighbors)

    def insert(self, item):
        """
        Inserts a new data point into the LSH Forest.

        Parameters
        ----------
        item: array_like, shape (n_features, )
            New data point to be inserted into the LSH Forest.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted before"
                             " inserting.")

        item = safe_asarray(item)

        if item.ndim != 1:
            raise ValueError("item shoud be a 1-D vector.")
        if item.shape[0] != self._input_array.shape[1]:
            raise ValueError("Number of features in item and"
                             " fitted array does not match.")

        input_array_shape = self._input_array.shape[0]
        trees = np.empty((self.n_trees, input_array_shape + 1),
                         dtype=int)
        original_incides = np.empty((self.n_trees,
                                     input_array_shape + 1),
                                    dtype=int)
        for i in range(self.n_trees):
            projections = np.array(np.dot(self.hash_functions_[i],
                                          item) > 0, dtype=int)
            xx = tuple(projections)
            bin_query = (self.cache[xx[:self.cache_N]] * self.k +
                         self.cache[xx[self.cache_N:]])
            # gets the position to be added in the tree.
            position = self._trees[i].searchsorted(bin_query)
            # adds the hashed value into the tree.
            trees[i] = np.insert(self._trees[i],
                                 position, bin_query)
            # add the entry into the original_indices.
            original_incides[i] = np.insert(self._original_indices[i],
                                            position,
                                            input_array_shape)
        self._trees = trees
        self._original_indices = original_incides

        # adds the entry into the input_array.
        self._input_array = np.row_stack((self._input_array, item))
