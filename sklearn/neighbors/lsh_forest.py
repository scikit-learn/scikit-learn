"""
Locality Sensitive Hashing Forest for Approximate Nearest Neighbor Search
-------------------------------------------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from ..base import BaseEstimator
from ..utils.validation import safe_asarray
from ..feature_extraction.lshashing import RandomProjections

__all__ = ["LSHForest"]


def _bisect_left(a, x):
    """Private function to perform bisect left operation"""
    return np.searchsorted(a, x)


def _bisect_right(a, x):
    """Private function to perform bisect right operation"""
    return np.searchsorted(np.array([item[:len(x)] for item in a]),
                           x, side='right')


def _find_matching_indices(sorted_array, item, h):
    """
    Finds indices in sorted array of strings where their first
    h elements match the items' first h elements
    """
    left_index = _bisect_left(sorted_array, item[:h])
    right_index = _bisect_right(sorted_array, item[:h])
    return np.arange(left_index, right_index)


def _find_longest_prefix_match(bit_string_array, query):
    """
    Private function to find the longest prefix match for query
    in the bit_string_array
    """
    hi = len(query)
    lo = 0

    if _find_matching_indices(bit_string_array, query, hi).shape[0] > 0:
        return hi

    while lo < hi:
        mid = (lo + hi) // 2
        k = _find_matching_indices(bit_string_array, query, mid).shape[0]
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

    max_label_length: int, optional (default = 32)
        Maximum length of a binary hash

    n_trees: int, optional (default = 10)
        Number of trees in the LSH Forest.

    hashing_algorithm: {'random_projections'},
        optional (default = 'random_projections')
        Algorithm of LSH family by which the hashing is performed on
        the data vectors.

        -'random_projections': hash using :class:`RandomProjections`

    c: int, optional (default = 10)
        Threshold value to select candidates for nearest neighbors.
        Number of candidates is often greater than c*n_trees(unless
        restricted by lower_bound)

    n_neighbors: int, optional(default = 1)
        Number of neighbors to be returned from query funcitond when
        it is not provided with the query.

    lower_bound: int, optional(defualt = 4)
        lowerest hash length to be searched when candidate selection is
        performed for nearest neighbors.

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

      >>> X = np.logspace(0, 3, num=50)
      >>> X = X.reshape((10,5))
      >>> lshf = LSHForest()
      >>> lshf.fit(X)
      LSHForest(c=50, hashing_algorithm='random_projections', lower_bound=4,
           max_label_length=32, n_neighbors=None, n_trees=10, seed=None)

      >>> lshf.kneighbors(X[:5], n_neighbors=3, return_distance=True)
      (array([[0, 1, 2],
             [1, 0, 2],
             [2, 1, 0],
             [3, 2, 1],
             [4, 3, 2]]), array([[  0.        ,   3.15525015,   9.54018168],
             [  0.        ,   3.15525015,   6.38493153],
             [  0.        ,   6.38493153,   9.54018168],
             [  0.        ,  12.92048135,  19.30541288],
             [  0.        ,  26.1457523 ,  39.06623365]]))

    """

    def __init__(self, max_label_length=32, n_trees=10,
                 hashing_algorithm='random_projections',
                 c=50, n_neighbors=1, lower_bound=4, seed=1):
        self.max_label_length = max_label_length
        self.n_trees = n_trees
        self.hashing_algorithm = hashing_algorithm
        self.random_state = np.random.RandomState(seed)
        self.c = c
        self.m = n_neighbors
        self.lower_bound = lower_bound

    def _select_hashing_algorithm(self, n_dim, hash_size):
        """ Selectes the LSH algorithm """
        if n_dim is None or hash_size is None:
            raise ValueError("n_dim or hash_size cannot be None.")

        if self.hashing_algorithm == 'random_projections':
            return RandomProjections(n_dim=n_dim, hash_size=hash_size)
        else:
            raise ValueError("Unknown hashing algorithm: %s"
                             % (self.hashing_algorithm))

    def _create_tree(self, hash_function=None):
        """
        Builds a single tree (in this case creates a sorted array of
        binary hashes).
        """
        if hash_function is None:
            raise ValueError("hash_funciton cannot be None.")

        n_points = self._input_array.shape[0]
        binary_hashes = []
        for i in range(n_points):
            binary_hashes.append(
                self.hash_generator.do_hash(self._input_array[i],
                                            hash_function))

        return np.argsort(binary_hashes), np.sort(binary_hashes)

    def _compute_distances(self, query, candidates):
        distances = _simple_euclidean_distance(
            query, self._input_array[candidates])
        return np.argsort(distances), np.sort(distances)

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
        n_dim = self._input_array.shape[1]

        self.hash_generator = self._select_hashing_algorithm(
            n_dim, self.max_label_length)

        # Creates a g(p,x) for each tree
        self.hash_functions_ = []
        self._trees = []
        self._original_indices = []
        for i in range(self.n_trees):
            # This is g(p,x) for a particular tree.
            hash_function = self.hash_generator.generate_hash_function()
            original_index, bin_hashes = self._create_tree(hash_function)
            self._original_indices.append(original_index)
            self._trees.append(bin_hashes)
            self.hash_functions_.append(hash_function)

        self.hash_functions_ = np.array(self.hash_functions_)
        self._trees = np.array(self._trees)
        self._original_indices = np.array(self._original_indices)

        return self

    def _query(self, query, m):
        """
        returns self.m number of neighbors and the distances
        for a given query.
        """
        if query is None:
            raise ValueError("query cannot be None.")
        query = np.array(query)

        bin_queries = []

        # descend phase
        max_depth = 0
        for i in range(self.n_trees):
            bin_query = self.hash_generator.do_hash(
                query, self.hash_functions_[i])
            k = _find_longest_prefix_match(self._trees[i], bin_query)
            if k > max_depth:
                max_depth = k
            bin_queries.append(bin_query)

        # Synchronous ascend phase
        candidates = []
        n_candidates = self.c * self.n_trees
        while max_depth > self.lower_bound and (len(candidates) < n_candidates
                                                or len(set(candidates)) < m):
            for i in range(self.n_trees):
                candidates.extend(
                    self._original_indices[i, _find_matching_indices(
                        self._trees[i],
                        bin_queries[i],
                        max_depth)].tolist())
            max_depth = max_depth - 1
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)

        return candidates[ranks[:m]], distances[:m]

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
        if n_neighbors is not None:
            self.m = n_neighbors
        X = safe_asarray(X)
        x_dim = X.ndim

        if x_dim == 1:
            neighbors, distances = self._query(X, self.m)
            if return_distance:
                return np.array([neighbors]), np.array([distances])
            else:
                return np.array([neighbors])
        else:
            neighbors, distances = [], []
            for i in range(X.shape[0]):
                neighs, dists = self._query(X[i], self.m)
                neighbors.append(neighs)
                distances.append(dists)
            return np.array(neighbors), np.array(distances)
