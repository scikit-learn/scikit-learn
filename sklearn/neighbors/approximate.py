"""Approximate nearest neighbor search"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from ..base import BaseEstimator
from ..utils.validation import check_array
from ..utils import check_random_state
from ..utils.extmath import row_norms

from ..random_projection import GaussianRandomProjection

__all__ = ["LSHForest"]


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


class LSHForest(BaseEstimator):
    """Performs approximate nearest neighbor search using LSH forest.

    LSH Forest: Locality Sensitive Hashing forest [1] is an alternative
    method for vanilla approximate nearest neighbor search methods.
    LSH forest data structure has been implemented using sorted
    arrays and binary search. 32 bit fixed length hashes are used in
    this implementation.

    Parameters
    ----------

    n_estimators : int (default = 10)
        Number of trees in the LSH Forest.

    n_candidates : int (default = 10)
        Value to restrict candidates selected from a single esitimator(tree)
        for nearest neighbors. Number of total candidates is often greater
        than n_candidates*n_estimators(unless restricted by min_hash_length)

    n_neighbors : int (default = 5)
        Number of neighbors to be returned from query function when
        it is not provided to :meth:`k_neighbors`

    min_hash_length : int (default = 4)
        lowest hash length to be searched when candidate selection is
        performed for nearest neighbors.

    radius : float, optinal (default = 1.0)
        Radius from the data point to its neighbors. This is the parameter
        space to use by default for :meth`radius_neighbors` queries.

    radius_cutoff_ratio : float, optional (defualt = 0.9)
        Cut off ratio of radius neighbors to candidates at the radius
        neighbor search

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------

    `hash_functions_` : list of arrays
        Hash function g(p,x) for a tree is an array of 32 randomly generated
        float arrays with the same dimension as the data set.


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
      LSHForest(min_hash_length=4, n_candidates=50, n_estimators=10, n_neighbors=5,
           radius=1.0, radius_cutoff_ratio=0.9, random_state=None)
      >>> distances, indices = lshf.kneighbors(X[:5], n_neighbors=3)
      >>> distances
      array([[ 1.08434102,  0.52344831,  0.        ],
             [ 0.56089272,  0.52344831,  0.        ],
             [ 0.60101568,  0.56089272,  0.        ],
             [ 0.6440088 ,  0.60101568,  0.        ],
             [ 0.6900774 ,  0.6440088 ,  0.        ]])
      >>> indices
      array([[2, 1, 0],
             [2, 0, 1],
             [3, 1, 2],
             [4, 2, 3],
             [5, 3, 4]])

    """

    def __init__(self, n_estimators=10, radius=1.0, n_candidates=50,
                 n_neighbors=5, min_hash_length=4, radius_cutoff_ratio=.9,
                 random_state=None):
        self.n_estimators = n_estimators
        self.radius = radius
        self.random_state = random_state
        self.n_candidates = n_candidates
        self.n_neighbors = n_neighbors
        self.min_hash_length = min_hash_length
        self.radius_cutoff_ratio = radius_cutoff_ratio

    def _create_tree(self, seed):
        """Builds a single tree.

        Here, it creates a sorted array of binary hashes.
        Hashing is done on an array of data points.
        `GaussianRandomProjection` is used for hashing.
        `n_components=hash_size and n_features=n_dim.

        This creates binary hashes by getting the dot product of
        input points and hash_function then transforming the projection
        into a binary string array based on the sign (positive/negative)
        of the projection.
        """
        grp = GaussianRandomProjection(n_components=self.max_label_length,
                                       random_state=seed)
        X = np.zeros((1, self._fit_X.shape[1]), dtype=float)
        grp.fit(X)

        hashes = (grp.transform(self._fit_X) > 0).astype(int)
        hash_function = grp.components_

        binary_hashes = np.packbits(hashes).view(dtype='>u4')
        original_indices = np.argsort(binary_hashes)

        return original_indices, binary_hashes[original_indices], hash_function

        return original_indices, binary_hashes[original_indices], hash_function

    def _compute_distances(self, query, candidates):
        """Computes the Euclidean distance.

        Distance is from the queryto points in the candidates array.
        Returns argsort of distances in the candidates
        array and sorted distances.
        """
        distances = row_norms(self._fit_X[candidates] - query)
        distance_positions = np.argsort(distances)
        return distance_positions, distances[distance_positions]

    def _generate_masks(self):
        """Creates left and right masks for all hash lengths."""
        tri_size = self.max_label_length + 1
        left_mask = np.tril(np.ones((tri_size, tri_size), dtype=int))[:, 1:]
        right_mask = np.triu(np.ones((tri_size, tri_size), dtype=int))[:, :-1]

        self._left_mask = np.packbits(left_mask).view(dtype='>u4')
        self._right_mask = np.packbits(right_mask).view(dtype='>u4')

    def _get_candidates(self, query, max_depth, bin_queries, n_neighbors):
        """Performs the Synchronous ascending phase.

        Returns an array of candidates, their distance ranks and
        distances.
        """
        candidates = []
        max_candidates = self.n_candidates * self.n_estimators
        while max_depth > self.min_hash_length and (len(candidates)
                                                    < max_candidates or
                                                    len(set(candidates))
                                                    < n_neighbors):

            for i in range(self.n_estimators):
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

        while max_depth > self.min_hash_length and (ratio_within_radius
                                                    > threshold):
            candidates = []
            for i in range(self.n_estimators):
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
            positions = np.searchsorted(total_distances, distances[:m])
            total_neighbors = np.insert(total_neighbors, positions,
                                        candidates[ranks[:m]])
            total_distances = np.insert(total_distances, positions,
                                        distances[:m])
            ratio_within_radius = (total_neighbors.shape[0] /
                                   float(total_candidates.shape[0]))
            max_depth = max_depth - 1
        return total_neighbors[::-1], total_distances[::-1]

    def _convert_to_hash(self, y, tree_n):
        """Converts item(a data point) into an integer.

        Value of the integer is the value represented by the
        binary hashed value.
        """
        projections = (np.dot(self.hash_functions_[tree_n],
                              y) > 0).astype(int)

        return np.packbits(projections).view(dtype='>u4')[0]

    def fit(self, X):
        """Fit the LSH forest on the data.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.
        """

        self._fit_X = check_array(X)
        self.max_label_length = 32

        # Creates a g(p,x) for each tree
        self.hash_functions_ = []
        self._trees = []
        self._original_indices = []

        random_state = check_random_state(self.random_state)

        for i in range(self.n_estimators):
            # This is g(p,x) for a particular tree.
            original_index, bin_hashes, hash_function = self._create_tree(
                random_state.randint(0, np.iinfo(np.int32).max))
            self._original_indices.append(original_index)
            self._trees.append(bin_hashes)
            self.hash_functions_.append(hash_function)

        self.hash_functions_ = np.array(self.hash_functions_)
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

        # descend phase
        max_depth = 0
        for i in range(self.n_estimators):
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
                                                                n_neighbors)

            return (candidates[ranks[:n_neighbors]][::-1],
                    distances[:n_neighbors][::-1])

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
            return np.array(distances), np.array(neighbors)
        else:
            return np.array(neighbors)

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
        self._fit_X = np.row_stack((self._fit_X, X))
