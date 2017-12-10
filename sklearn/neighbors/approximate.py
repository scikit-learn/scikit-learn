"""Approximate nearest neighbor search"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
#         Joel Nothman <joel.nothman@gmail.com>

import numpy as np
import warnings

from scipy import sparse

from .base import KNeighborsMixin, RadiusNeighborsMixin
from ..base import BaseEstimator
from ..utils.validation import check_array
from ..utils import check_random_state
from ..metrics.pairwise import pairwise_distances

from ..random_projection import GaussianRandomProjection

__all__ = ["LSHForest"]

HASH_DTYPE = '>u4'
MAX_HASH_SIZE = np.dtype(HASH_DTYPE).itemsize * 8


def _find_matching_indices(tree, bin_X, left_mask, right_mask):
    """Finds indices in sorted array of integers.

    Most significant h bits in the binary representations of the
    integers are matched with the items' most significant h bits.
    """
    left_index = np.searchsorted(tree, bin_X & left_mask)
    right_index = np.searchsorted(tree, bin_X | right_mask,
                                  side='right')
    return left_index, right_index


def _find_longest_prefix_match(tree, bin_X, hash_size,
                               left_masks, right_masks):
    """Find the longest prefix match in tree for each query in bin_X

    Most significant bits are considered as the prefix.
    """
    hi = np.empty_like(bin_X, dtype=np.intp)
    hi.fill(hash_size)
    lo = np.zeros_like(bin_X, dtype=np.intp)
    res = np.empty_like(bin_X, dtype=np.intp)

    left_idx, right_idx = _find_matching_indices(tree, bin_X,
                                                 left_masks[hi],
                                                 right_masks[hi])
    found = right_idx > left_idx
    res[found] = lo[found] = hash_size

    r = np.arange(bin_X.shape[0])
    kept = r[lo < hi]  # indices remaining in bin_X mask
    while kept.shape[0]:
        mid = (lo.take(kept) + hi.take(kept)) // 2

        left_idx, right_idx = _find_matching_indices(tree,
                                                     bin_X.take(kept),
                                                     left_masks[mid],
                                                     right_masks[mid])
        found = right_idx > left_idx
        mid_found = mid[found]
        lo[kept[found]] = mid_found + 1
        res[kept[found]] = mid_found
        hi[kept[~found]] = mid[~found]

        kept = r[lo < hi]

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

    def transform(self, X):
        return self._to_hash(super(ProjectionToHashMixin, self).transform(X))


class GaussianRandomProjectionHash(ProjectionToHashMixin,
                                   GaussianRandomProjection):
    """Use GaussianRandomProjection to produce a cosine LSH fingerprint"""
    def __init__(self,
                 n_components=32,
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

    The cosine distance is defined as ``1 - cosine_similarity``: the lowest
    value is 0 (identical point) but it is bounded above by 2 for the farthest
    points. Its value does not depend on the norm of the vector points but
    only on their relative angles.

    Parameters
    ----------

    n_estimators : int (default = 10)
        Number of trees in the LSH Forest.

    radius : float, optinal (default = 1.0)
        Radius from the data point to its neighbors. This is the parameter
        space to use by default for the :meth:`radius_neighbors` queries.

    n_candidates : int (default = 50)
        Minimum number of candidates evaluated per estimator, assuming enough
        items meet the `min_hash_match` constraint.

    n_neighbors : int (default = 5)
        Number of neighbors to be returned from query function when
        it is not provided to the :meth:`kneighbors` method.

    min_hash_match : int (default = 4)
        lowest hash length to be searched when candidate selection is
        performed for nearest neighbors.

    radius_cutoff_ratio : float, optional (default = 0.9)
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
        float arrays with the same dimension as the data set. This array is
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
      >>> from sklearn.neighbors import LSHForest

      >>> X_train = [[5, 5, 2], [21, 5, 5], [1, 1, 1], [8, 9, 1], [6, 10, 2]]
      >>> X_test = [[9, 1, 6], [3, 1, 10], [7, 10, 3]]
      >>> lshf = LSHForest(random_state=42)
      >>> lshf.fit(X_train)  # doctest: +NORMALIZE_WHITESPACE
      LSHForest(min_hash_match=4, n_candidates=50, n_estimators=10,
                n_neighbors=5, radius=1.0, radius_cutoff_ratio=0.9,
                random_state=42)
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

        warnings.warn("LSHForest has poor performance and has been deprecated "
                      "in 0.19. It will be removed in version 0.21.",
                      DeprecationWarning)

    def _compute_distances(self, query, candidates):
        """Computes the cosine distance.

        Distance is from the query to points in the candidates array.
        Returns argsort of distances in the candidates
        array and sorted distances.
        """
        if candidates.shape == (0,):
            # needed since _fit_X[np.array([])] doesn't work if _fit_X sparse
            return np.empty(0, dtype=np.int), np.empty(0, dtype=float)

        if sparse.issparse(self._fit_X):
            candidate_X = self._fit_X[candidates]
        else:
            candidate_X = self._fit_X.take(candidates, axis=0, mode='clip')
        distances = pairwise_distances(query, candidate_X,
                                       metric='cosine')[0]
        distance_positions = np.argsort(distances)
        distances = distances.take(distance_positions, mode='clip', axis=0)
        return distance_positions, distances

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
        # Number of candidates considered including duplicates
        # XXX: not sure whether this is being calculated correctly wrt
        #      duplicates from different iterations through a single tree
        n_candidates = 0
        candidate_set = set()
        min_candidates = self.n_candidates * self.n_estimators
        while (max_depth > self.min_hash_match and
               (n_candidates < min_candidates or
                len(candidate_set) < n_neighbors)):

            left_mask = self._left_mask[max_depth]
            right_mask = self._right_mask[max_depth]
            for i in range(self.n_estimators):
                start, stop = _find_matching_indices(self.trees_[i],
                                                     bin_queries[i],
                                                     left_mask, right_mask)
                n_candidates += stop - start
                candidate_set.update(
                    self.original_indices_[i][start:stop].tolist())
            max_depth -= 1

        candidates = np.fromiter(candidate_set, count=len(candidate_set),
                                 dtype=np.intp)
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

        while (max_depth > self.min_hash_match and
               ratio_within_radius > threshold):
            left_mask = self._left_mask[max_depth]
            right_mask = self._right_mask[max_depth]
            candidates = []
            for i in range(self.n_estimators):
                start, stop = _find_matching_indices(self.trees_[i],
                                                     bin_queries[i],
                                                     left_mask, right_mask)
                candidates.extend(
                    self.original_indices_[i][start:stop].tolist())
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

    def fit(self, X, y=None):
        """Fit the LSH forest on the data.

        This creates binary hashes of input data points by getting the
        dot product of input points and hash_function then
        transforming the projection into a binary string array based
        on the sign (positive/negative) of the projection.
        A sorted array of binary hashes is created.

        Parameters
        ----------
        X : array_like or sparse (CSR) matrix, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        self : object
            Returns self.
        """

        self._fit_X = check_array(X, accept_sparse='csr')

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

    def _query(self, X):
        """Performs descending phase to find maximum depth."""
        # Calculate hashes of shape (n_samples, n_estimators, [hash_size])
        bin_queries = np.asarray([hasher.transform(X)[:, 0]
                                  for hasher in self.hash_functions_])
        bin_queries = np.rollaxis(bin_queries, 1)

        # descend phase
        depths = [_find_longest_prefix_match(tree, tree_queries, MAX_HASH_SIZE,
                                             self._left_mask, self._right_mask)
                  for tree, tree_queries in zip(self.trees_,
                                                np.rollaxis(bin_queries, 1))]

        return bin_queries, np.max(depths, axis=0)

    def kneighbors(self, X, n_neighbors=None, return_distance=True):
        """Returns n_neighbors of approximate nearest neighbors.

        Parameters
        ----------
        X : array_like or sparse (CSR) matrix, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single query.

        n_neighbors : int, optional (default = None)
            Number of neighbors required. If not provided, this will
            return the number specified at the initialization.

        return_distance : boolean, optional (default = True)
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

        X = check_array(X, accept_sparse='csr')

        neighbors, distances = [], []
        bin_queries, max_depth = self._query(X)
        for i in range(X.shape[0]):

            neighs, dists = self._get_candidates(X[[i]], max_depth[i],
                                                 bin_queries[i],
                                                 n_neighbors)
            neighbors.append(neighs)
            distances.append(dists)

        if return_distance:
            return np.array(distances), np.array(neighbors)
        else:
            return np.array(neighbors)

    def radius_neighbors(self, X, radius=None, return_distance=True):
        """Finds the neighbors within a given radius of a point or points.

        Return the indices and distances of some points from the dataset
        lying in a ball with size ``radius`` around the points of the query
        array. Points lying on the boundary are included in the results.

        The result points are *not* necessarily sorted by distance to their
        query point.

        LSH Forest being an approximate method, some true neighbors from the
        indexed dataset might be missing from the results.

        Parameters
        ----------
        X : array_like or sparse (CSR) matrix, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single query.

        radius : float
            Limiting distance of neighbors to return.
            (default is the value passed to the constructor).

        return_distance : boolean, optional (default = False)
            Returns the distances of neighbors if set to True.

        Returns
        -------
        dist : array, shape (n_samples,) of arrays
            Each element is an array representing the cosine distances
            to some points found within ``radius`` of the respective query.
            Only present if ``return_distance=True``.

        ind : array, shape (n_samples,) of arrays
            Each element is an array of indices for neighbors within ``radius``
            of the respective query.
        """
        if not hasattr(self, 'hash_functions_'):
            raise ValueError("estimator should be fitted.")

        if radius is None:
            radius = self.radius

        X = check_array(X, accept_sparse='csr')

        neighbors, distances = [], []
        bin_queries, max_depth = self._query(X)
        for i in range(X.shape[0]):

            neighs, dists = self._get_radius_neighbors(X[[i]], max_depth[i],
                                                       bin_queries[i], radius)
            neighbors.append(neighs)
            distances.append(dists)

        if return_distance:
            return _array_of_arrays(distances), _array_of_arrays(neighbors)
        else:
            return _array_of_arrays(neighbors)

    def partial_fit(self, X, y=None):
        """
        Inserts new data into the already fitted LSH Forest.
        Cost is proportional to new total size, so additions
        should be batched.

        Parameters
        ----------
        X : array_like or sparse (CSR) matrix, shape (n_samples, n_features)
            New data point to be inserted into the LSH Forest.
        """
        X = check_array(X, accept_sparse='csr')
        if not hasattr(self, 'hash_functions_'):
            return self.fit(X)

        if X.shape[1] != self._fit_X.shape[1]:
            raise ValueError("Number of features in X and"
                             " fitted array does not match.")
        n_samples = X.shape[0]
        n_indexed = self._fit_X.shape[0]

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
                                                  np.arange(n_indexed,
                                                            n_indexed +
                                                            n_samples))

        # adds the entry into the input_array.
        if sparse.issparse(X) or sparse.issparse(self._fit_X):
            self._fit_X = sparse.vstack((self._fit_X, X))
        else:
            self._fit_X = np.row_stack((self._fit_X, X))

        return self
