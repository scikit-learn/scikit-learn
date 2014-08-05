"""
Locality Sensitive Hashing Forest for Approximate Nearest Neighbor Search
-------------------------------------------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from bisect import bisect_left, bisect_right
from ..base import BaseEstimator
from ..utils.validation import check_array
from ..utils import check_random_state

from ..random_projection import GaussianRandomProjection

__all__ = ["LSHForest"]


def _find_matching_indices(sorted_array, item, left_mask, right_mask):
    """Finds indices in sorted array of integers.

    Most significant h bits in the binary representations of the
    integers are matched with the items' most significant h bits.
    """
    left_index = bisect_left(sorted_array, item & left_mask)
    right_index = bisect_right(sorted_array, item | right_mask)
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
    """Private function to calculate Euclidean distances.

    Distance is calculated between each point in candidates
    and query.
    """
    distances = np.zeros(candidates.shape[0])
    for i in range(candidates.shape[0]):
        distances[i] = np.linalg.norm(candidates[i] - query)
    return distances


class LSHForest(BaseEstimator):

    """Performs approximate nearest neighbor search using LSH forest.

    LSH Forest: Locality Sensitive Hashing forest [1] is an alternative
    method for vanilla approximate nearest neighbor search methods.
    LSH forest data structure has been implemented using sorted
    arrays and binary search. 32 bit fixed length hashes are used in
    this implementation.

    Parameters
    ----------

    n_trees: int (default = 10)
        Number of trees in the LSH Forest.

    c: int (default = 10)
        Value to restrict candidates selection for nearest neighbors.
        Number of candidates is often greater than c*n_trees(unless
        restricted by lower_bound)

    n_neighbors: int (default = 1)
        Number of neighbors to be returned from query function when
        it is not provided to :meth:`k_neighbors`

    lower_bound: int (defualt = 4)
        lowest hash length to be searched when candidate selection is
        performed for nearest neighbors.

    radius : float, optinal (default = 1.0)
        Radius from the data point to its neighbors. This is the parameter
        space to use by default for :meth`radius_neighbors` queries.

    radius_cutoff_ratio: float, optional (defualt = 0.9)
        Cut off ratio of radius neighbors to candidates at the radius
        neighbor search

    random_state: numpy.RandomState, optional
        The generator used to initialize random projections.
        Defaults to numpy.random.

    Attributes
    ----------

    `hash_functions_`: list of arrays
        Hash function g(p,x) for a tree is an array of 32 randomly generated
        float arrays with the same dimenstion as the data set.


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
      LSHForest(c=50, lower_bound=4, n_neighbors=1, n_trees=10, radius=1.0,
           radius_cutoff_ratio=0.9, random_state=None)

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

    def __init__(self, n_trees=10, radius=1.0, c=50, n_neighbors=1,
                 lower_bound=4, radius_cutoff_ratio=.9,
                 random_state=None):
        self.n_trees = n_trees
        self.radius = radius
        self.random_state = random_state
        self.c = c
        self.n_neighbors = n_neighbors
        self.lower_bound = lower_bound
        self.radius_cutoff_ratio = radius_cutoff_ratio

    def _generate_hash_function(self):
        """Fits a `GaussianRandomProjection`.

        `n_components=hash_size and n_features=n_dim.
        """
        random_state = check_random_state(self.random_state)
        grp = GaussianRandomProjection(n_components=self.max_label_length,
                                       random_state=random_state.randint(0,
                                                                         10))
        X = np.zeros((1, self._n_dim), dtype=float)
        grp.fit(X)
        return grp

    def _create_tree(self):
        """Builds a single tree.

        Here, it creates a sorted array of binary hashes.
        Hashing is done on an array of data points.
        This creates a binary hashes by getting the dot product of
        input points and hash_function then transforming the projection
        into a binary string array based on the sign (positive/negative)
        of the projection.
        """
        grp = self._generate_hash_function()
        hashes = np.array(grp.transform(self._input_array) > 0, dtype=int)
        hash_function = grp.components_

        binary_hashes = np.packbits(hashes).view(dtype='>u4')

        return np.argsort(binary_hashes), np.sort(binary_hashes), hash_function

    def _compute_distances(self, query, candidates):
        """Computes the Euclidean distance.

        Distance is from the queryto points in the candidates array.
        Returns argsort of distances in the candidates
        array and sorted distances.
        """
        distances = _simple_euclidean_distance(
            query, self._input_array[candidates])
        return np.argsort(distances), np.sort(distances)

    def _generate_masks(self):
        """Creates left and right masks for all hash lengths."""
        tri_size = self.max_label_length + 1
        left_mask = np.tril(np.ones((tri_size, tri_size), dtype=int))[:, 1:]
        right_mask = np.triu(np.ones((tri_size, tri_size), dtype=int))[:, :-1]

        self._left_mask = np.packbits(left_mask).view(dtype='>u4')
        self._right_mask = np.packbits(right_mask).view(dtype='>u4')

    def _get_candidates(self, query, max_depth, bin_queries, m):
        """Performs the Synchronous ascending phase.

        Returns an array of candidates, their distance ranks and
        distances.
        """
        candidates = []
        n_candidates = self.c * self.n_trees
        while max_depth > self.lower_bound and (len(candidates) < n_candidates
                                                or len(set(candidates)) < m):
            for i in range(self.n_trees):
                candidates.extend(
                    self._original_indices[i][_find_matching_indices(
                        self._trees[i],
                        bin_queries[i],
                        self._left_mask[max_depth],
                        self._right_mask[max_depth])].tolist())
            max_depth = max_depth - 1
        candidates = np.unique(candidates)
        ranks, distances = self._compute_distances(query, candidates)

        return candidates, ranks, distances

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

        while max_depth > self.lower_bound and ratio_within_radius > threshold:
            candidates = []
            for i in range(self.n_trees):
                candidates.extend(
                    self._original_indices[i][_find_matching_indices(
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

    def _convert_to_hash(self, item, tree_n):
        """Converts item(a date point) into an integer.

        Value of the integer is the value represented by the
        binary hashed value.
        """
        projections = np.array(np.dot(self.hash_functions_[tree_n],
                                      item) > 0, dtype=int)

        return np.packbits(projections).view(dtype='>u4')[0]

    def fit(self, X):
        """Fit the LSH forest on the data.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.
        """

        self._input_array = check_array(X)
        self._n_dim = self._input_array.shape[1]

        self.max_label_length = 32

        # Creates a g(p,x) for each tree
        self.hash_functions_ = []
        self._trees = []
        self._original_indices = []

        for i in range(self.n_trees):
            # This is g(p,x) for a particular tree.
            original_index, bin_hashes, hash_function = self._create_tree()
            self._original_indices.append(original_index)
            self._trees.append(bin_hashes)
            self.hash_functions_.append(hash_function)

        self.hash_functions_ = np.array(self.hash_functions_)
        self._generate_masks()

        return self

    def _query(self, query, m=None, radius=None, is_radius=False):
        """Finds neighbors and distances.
        Returns the neighbors whose distances from the query is less
        than radius if is_radius is True.
        Otherwise returns m number of neighbors and the distances
        for a given query.
        """
        bin_queries = []

        # descend phase
        max_depth = 0
        for i in range(self.n_trees):
            bin_query = self._convert_to_hash(query, i)
            k = _find_longest_prefix_match(self._trees[i], bin_query,
                                           self.max_label_length,
                                           self._left_mask,
                                           self._right_mask)
            if k > max_depth:
                max_depth = k
            bin_queries.append(bin_query)

        if is_radius:
            return self._get_radius_neighbors(query, max_depth,
                                              bin_queries, radius)

        else:
            candidates, ranks, distances = self._get_candidates(query,
                                                                max_depth,
                                                                bin_queries,
                                                                m)

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
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        X = check_array(X)
        x_dim = X.ndim

        if x_dim == 1:
            neighbors, distances = self._query(X, n_neighbors)
            if return_distance:
                return np.array([neighbors]), np.array([distances])
            else:
                return np.array([neighbors])
        else:
            neighbors, distances = [], []
            for i in range(X.shape[0]):
                neighs, dists = self._query(X[i], n_neighbors)
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

        if radius is None:
            radius = self.radius

        X = check_array(X)
        x_dim = X.ndim

        if x_dim == 1:
            neighbors, distances = self._query(X, radius=radius,
                                               is_radius=True)
            if return_distance:
                return np.array([neighbors]), np.array([distances])
            else:
                return np.array([neighbors])
        else:
            neighbors, distances = [], []
            for i in range(X.shape[0]):
                neighs, dists = self._query(X[i], radius=radius,
                                            is_radius=True)
                neighbors.append(neighs)
                distances.append(dists)

            if return_distance:
                return np.array(neighbors), np.array(distances)
            else:
                return np.array(neighbors)

    def insert(self, X):
        """
        Inserts new data into the LSH Forest. Cost is proportional
        to new total size, so additions should be batched.

        Parameters
        ----------
        X: array_like, shape (n_samples, n_features)
            New data point to be inserted into the LSH Forest.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted before"
                             " inserting.")

        X = check_array(X)

        if X.shape[1] != self._input_array.shape[1]:
            raise ValueError("Number of features in X and"
                             " fitted array does not match.")
        n_samples = X.shape[0]
        input_array_size = self._input_array.shape[0]

        for i in range(self.n_trees):
            bin_X = [self._convert_to_hash(X[j], i) for j in range(n_samples)]
            # gets the position to be added in the tree.
            positions = self._trees[i].searchsorted(bin_X)
            # adds the hashed value into the tree.
            self._trees[i] = np.insert(self._trees[i],
                                       positions, bin_X)
            # add the entry into the original_indices.
            self._original_indices[i] = np.insert(self._original_indices[i],
                                                  positions,
                                                  np.arange(input_array_size,
                                                            input_array_size +
                                                            n_samples))

        # adds the entry into the input_array.
        self._input_array = np.row_stack((self._input_array, X))
