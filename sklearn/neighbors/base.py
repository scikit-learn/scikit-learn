"""Base and mixin classes for nearest neighbors"""
# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam
from functools import partial
from distutils.version import LooseVersion

import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import csr_matrix, issparse

from .ball_tree import BallTree
from .kd_tree import KDTree
from ..base import BaseEstimator
from ..metrics import pairwise_distances_chunked
from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from ..utils import check_X_y, check_array, gen_even_slices
from ..utils.multiclass import check_classification_targets
from ..utils.validation import check_is_fitted
from ..utils.validation import check_non_negative
from ..externals import six
from ..utils import Parallel, delayed, effective_n_jobs
from ..utils._joblib import __version__ as joblib_version
from ..exceptions import DataConversionWarning, EfficiencyWarning

VALID_METRICS = dict(ball_tree=BallTree.valid_metrics,
                     kd_tree=KDTree.valid_metrics,
                     # The following list comes from the
                     # sklearn.metrics.pairwise doc string
                     brute=(list(PAIRWISE_DISTANCE_FUNCTIONS.keys()) +
                            ['braycurtis', 'canberra', 'chebyshev',
                             'correlation', 'cosine', 'dice', 'hamming',
                             'jaccard', 'kulsinski', 'mahalanobis',
                             'matching', 'minkowski', 'rogerstanimoto',
                             'russellrao', 'seuclidean', 'sokalmichener',
                             'sokalsneath', 'sqeuclidean',
                             'yule', 'wminkowski']))


VALID_METRICS_SPARSE = dict(ball_tree=[],
                            kd_tree=[],
                            brute=PAIRWISE_DISTANCE_FUNCTIONS.keys())


def _check_weights(weights):
    """Check to make sure weights are valid"""
    if weights in (None, 'uniform', 'distance'):
        return weights
    elif callable(weights):
        return weights
    else:
        raise ValueError("weights not recognized: should be 'uniform', "
                         "'distance', or a callable function")


def _get_weights(dist, weights):
    """Get the weights from an array of distances and a parameter ``weights``

    Parameters
    ===========
    dist : ndarray
        The input distances
    weights : {'uniform', 'distance' or a callable}
        The kind of weighting used

    Returns
    ========
    weights_arr : array of the same shape as ``dist``
        if ``weights == 'uniform'``, then returns None
    """
    if weights in (None, 'uniform'):
        return None
    elif weights == 'distance':
        # if user attempts to classify a point that was zero distance from one
        # or more training points, those training points are weighted as 1.0
        # and the other points as 0.0
        if dist.dtype is np.dtype(object):
            for point_dist_i, point_dist in enumerate(dist):
                # check if point_dist is iterable
                # (ex: RadiusNeighborClassifier.predict may set an element of
                # dist to 1e-6 to represent an 'outlier')
                if hasattr(point_dist, '__contains__') and 0. in point_dist:
                    dist[point_dist_i] = point_dist == 0.
                else:
                    dist[point_dist_i] = 1. / point_dist
        else:
            with np.errstate(divide='ignore'):
                dist = 1. / dist
            inf_mask = np.isinf(dist)
            inf_row = np.any(inf_mask, axis=1)
            dist[inf_row] = inf_mask[inf_row]
        return dist
    elif callable(weights):
        return weights(dist)
    else:
        raise ValueError("weights not recognized: should be 'uniform', "
                         "'distance', or a callable function")


def _is_sorted_by_data(graph):
    """Returns whether the graph's non-zero entries are sorted by data

    The non-zero entries are stored in graph.data and graph.indices.
    For each row (or sample), the non-zero entries can be either:
        - sorted by indices, as after graph.sort_indices()
        - sorted by data, as after _check_precomputed(graph)
        - not sorted.

    Parameters
    ----------
    graph : CSR sparse matrix, shape (n_samples, n_samples)
        Neighbors graph as given by kneighbors_graph or radius_neighbors_graph

    Returns
    -------
    res : boolean
        Whether input graph is sorted by data
    """
    assert graph.format == 'csr'
    out_of_order = graph.data[:-1] > graph.data[1:]
    line_change = np.unique(graph.indptr[1:-1] - 1)
    line_change = line_change[line_change < out_of_order.shape[0]]
    return (out_of_order.sum() == out_of_order[line_change].sum())


def _check_precomputed(X):
    """Check precomputed distance matrix

    If the precomputed distance matrix is sparse, it checks that the non-zero
    entries are sorted by distances. If not, the matrix is copied and sorted.

    Parameters
    ----------
    X : {sparse matrix, array-like}, (n_samples, n_samples)
        Distance matrix to other samples. X may be a sparse matrix, in which
        case only non-zero elements may be considered neighbors.

    Returns
    -------
    X : {sparse matrix, array-like}, (n_samples, n_samples)
        Distance matrix to other samples. X may be a sparse matrix, in which
        case only non-zero elements may be considered neighbors.
    """
    if not issparse(X):
        X = check_array(X)
        check_non_negative(X, whom="precomputed distance matrix.")
        return X
    else:
        graph = X

    if graph.format not in ('csr', 'csc', 'coo', 'lil'):
        raise TypeError('Sparse matrix in {!r} format is not supported due to '
                        'its handling of explicit zeros'.format(graph.format))
    copied = graph.format != 'csr'
    graph = check_array(graph, accept_sparse='csr')
    check_non_negative(graph, whom="precomputed distance matrix.")

    if not _is_sorted_by_data(graph):
        warnings.warn('Precomputed sparse input was not sorted by data.',
                      EfficiencyWarning)
        if not copied:
            graph = graph.copy()

        # if each sample has the same number of provided neighbors
        row_nnz = np.diff(graph.indptr)
        if row_nnz.max() == row_nnz.min():
            n_samples = graph.shape[0]
            distances = graph.data.reshape(n_samples, -1)

            order = np.argsort(distances, kind='mergesort')
            order += np.arange(n_samples)[:, None] * row_nnz[0]
            order = order.ravel()
            graph.data = graph.data[order]
            graph.indices = graph.indices[order]

        else:
            for start, stop in zip(graph.indptr, graph.indptr[1:]):
                order = np.argsort(graph.data[start:stop], kind='mergesort')
                graph.data[start:stop] = graph.data[start:stop][order]
                graph.indices[start:stop] = graph.indices[start:stop][order]
    return graph


def _kneighbors_from_graph(graph, n_neighbors, return_distance):
    """Decompose a nearest neighbors sparse graph into distances and indices

    Parameters
    ----------
    graph : CSR sparse matrix, shape (n_samples, n_samples)
        Neighbors graph as given by kneighbors_graph or radius_neighbors_graph

    n_neighbors : int
        Number of neighbors required for each sample.

    return_distance : boolean
        If False, distances will not be returned

    Returns
    -------
    neigh_dist : array, shape (n_samples, n_neighbors)
        Distances to nearest neighbors. Only present if return_distance=True.

    neigh_ind : array, shape (n_samples, n_neighbors)
        Indices of nearest neighbors.
    """
    n_samples = graph.shape[0]
    assert graph.format == 'csr'

    # number of neighbors by samples
    row_nnz = np.diff(graph.indptr)
    row_nnz_min = row_nnz.min()
    if n_neighbors is not None and row_nnz_min < n_neighbors:
        raise ValueError(
            '%d neighbors per samples are required, but some samples have only'
            ' %d neighbors in precomputed graph matrix. Decrease number of '
            'neighbors used or recompute the graph with more neighbors.'
            % (n_neighbors, row_nnz_min))

    # if each sample has the same number of provided neighbors
    if row_nnz.max() == row_nnz_min:

        def extract(a):
            return a.reshape(n_samples, -1)[:, :n_neighbors]

    else:
        idx = np.tile(np.arange(n_neighbors), n_samples)
        idx = idx.reshape(n_samples, n_neighbors)
        idx += graph.indptr[:-1, np.newaxis]

        def extract(a):
            return a.take(idx, mode='clip').reshape(n_samples, n_neighbors)

    if return_distance:
        return extract(graph.data), extract(graph.indices)
    else:
        return extract(graph.indices)


def _radius_neighbors_from_graph(graph, radius, return_distance):
    """Decompose a nearest neighbors sparse graph into distances and indices

    Parameters
    ----------
    graph : CSR sparse matrix, shape (n_samples, n_samples)
        Neighbors graph as given by kneighbors_graph or radius_neighbors_graph

    radius : float > 0
        Radius of neighborhoods.

    return_distance : boolean
        If False, distances will not be returned

    Returns
    -------
    neigh_dist : array, shape (n_samples,) of arrays
        Distances to nearest neighbors. Only present if return_distance=True.

    neigh_ind :array, shape (n_samples,) of arrays
        Indices of nearest neighbors.
    """
    assert graph.format == 'csr'

    if graph.data.max() <= radius:
        data, indices, indptr = graph.data, graph.indices, graph.indptr
    else:
        mask = graph.data <= radius
        if return_distance:
            data = np.compress(mask, graph.data)
        indices = np.compress(mask, graph.indices)
        indptr = np.concatenate(([0], np.cumsum(mask)))[graph.indptr]

    if return_distance:
        neigh_dist = np.array(np.split(data, indptr[1:-1]))
    neigh_ind = np.array(np.split(indices, indptr[1:-1]))

    if return_distance:
        return neigh_dist, neigh_ind
    else:
        return neigh_ind


class NeighborsBase(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for nearest neighbors estimators."""

    @abstractmethod
    def __init__(self, n_neighbors=None, radius=None,
                 algorithm='auto', leaf_size=30, metric='minkowski',
                 p=2, metric_params=None, n_jobs=None):

        self.n_neighbors = n_neighbors
        self.radius = radius
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.metric = metric
        self.metric_params = metric_params
        self.p = p
        self.n_jobs = n_jobs
        self._check_algorithm_metric()

    def _check_algorithm_metric(self):
        if self.algorithm not in ['auto', 'brute',
                                  'kd_tree', 'ball_tree']:
            raise ValueError("unrecognized algorithm: '%s'" % self.algorithm)

        if self.algorithm == 'auto':
            if self.metric == 'precomputed':
                alg_check = 'brute'
            elif (callable(self.metric) or
                  self.metric in VALID_METRICS['ball_tree']):
                alg_check = 'ball_tree'
            else:
                alg_check = 'brute'
        else:
            alg_check = self.algorithm

        if callable(self.metric):
            if self.algorithm == 'kd_tree':
                # callable metric is only valid for brute force and ball_tree
                raise ValueError(
                    "kd_tree algorithm does not support callable metric '%s'"
                    % self.metric)
        elif self.metric not in VALID_METRICS[alg_check]:
            raise ValueError("Metric '%s' not valid for algorithm '%s'"
                             % (self.metric, self.algorithm))

        if self.metric_params is not None and 'p' in self.metric_params:
            warnings.warn("Parameter p is found in metric_params. "
                          "The corresponding parameter from __init__ "
                          "is ignored.", SyntaxWarning, stacklevel=3)
            effective_p = self.metric_params['p']
        else:
            effective_p = self.p

        if self.metric in ['wminkowski', 'minkowski'] and effective_p < 1:
            raise ValueError("p must be greater than one for minkowski metric")

    def _fit(self, X):
        self._check_algorithm_metric()
        if self.metric_params is None:
            self.effective_metric_params_ = {}
        else:
            self.effective_metric_params_ = self.metric_params.copy()

        effective_p = self.effective_metric_params_.get('p', self.p)
        if self.metric in ['wminkowski', 'minkowski']:
            self.effective_metric_params_['p'] = effective_p

        self.effective_metric_ = self.metric
        # For minkowski distance, use more efficient methods where available
        if self.metric == 'minkowski':
            p = self.effective_metric_params_.pop('p', 2)
            if p < 1:
                raise ValueError("p must be greater than one "
                                 "for minkowski metric")
            elif p == 1:
                self.effective_metric_ = 'manhattan'
            elif p == 2:
                self.effective_metric_ = 'euclidean'
            elif p == np.inf:
                self.effective_metric_ = 'chebyshev'
            else:
                self.effective_metric_params_['p'] = p

        if isinstance(X, NeighborsBase):
            self._fit_X = X._fit_X
            self._tree = X._tree
            self._fit_method = X._fit_method
            return self

        elif isinstance(X, BallTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'ball_tree'
            return self

        elif isinstance(X, KDTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'kd_tree'
            return self

        if self.effective_metric_ == 'precomputed':
            X = _check_precomputed(X)
        else:
            X = check_array(X, accept_sparse='csr')

        n_samples = X.shape[0]
        if n_samples == 0:
            raise ValueError("n_samples must be greater than 0")

        if issparse(X):
            if self.algorithm not in ('auto', 'brute'):
                warnings.warn("cannot use tree with sparse input: "
                              "using brute force")
            if self.effective_metric_ not in VALID_METRICS_SPARSE['brute'] \
                    and not callable(self.effective_metric_):

                raise ValueError("metric '%s' not valid for sparse input"
                                 % self.effective_metric_)
            self._fit_X = X.copy()
            self._tree = None
            self._fit_method = 'brute'
            return self

        self._fit_method = self.algorithm
        self._fit_X = X

        if self._fit_method == 'auto':
            # A tree approach is better for small number of neighbors,
            # and KDTree is generally faster when available
            if ((self.n_neighbors is None or
                 self.n_neighbors < self._fit_X.shape[0] // 2) and
                    self.metric != 'precomputed'):
                if self.effective_metric_ in VALID_METRICS['kd_tree']:
                    self._fit_method = 'kd_tree'
                elif (callable(self.effective_metric_) or
                        self.effective_metric_ in VALID_METRICS['ball_tree']):
                    self._fit_method = 'ball_tree'
                else:
                    self._fit_method = 'brute'
            else:
                self._fit_method = 'brute'

        if self._fit_method == 'ball_tree':
            self._tree = BallTree(X, self.leaf_size,
                                  metric=self.effective_metric_,
                                  **self.effective_metric_params_)
        elif self._fit_method == 'kd_tree':
            self._tree = KDTree(X, self.leaf_size,
                                metric=self.effective_metric_,
                                **self.effective_metric_params_)
        elif self._fit_method == 'brute':
            self._tree = None
        else:
            raise ValueError("algorithm = '%s' not recognized"
                             % self.algorithm)

        if self.n_neighbors is not None:
            if self.n_neighbors <= 0:
                raise ValueError(
                    "Expected n_neighbors > 0. Got %d" %
                    self.n_neighbors
                )
            else:
                if not np.issubdtype(type(self.n_neighbors), np.integer):
                    raise TypeError(
                        "n_neighbors does not take %s value, "
                        "enter integer value" %
                        type(self.n_neighbors))

        return self

    @property
    def _pairwise(self):
        # For cross-validation routines to split data correctly
        return self.metric == 'precomputed'


class KNeighborsMixin(object):
    """Mixin for k-neighbors searches"""

    def _kneighbors_reduce_func(self, dist, start,
                                n_neighbors, return_distance):
        """Reduce a chunk of distances to the nearest neighbors

        Callback to :func:`sklearn.metrics.pairwise.pairwise_distances_chunked`

        Parameters
        ----------
        dist : array of shape (n_samples_chunk, n_samples)
        start : int
            The index in X which the first row of dist corresponds to.
        n_neighbors : int
        return_distance : bool

        Returns
        -------
        dist : array of shape (n_samples_chunk, n_neighbors), optional
            Returned only if return_distance
        neigh : array of shape (n_samples_chunk, n_neighbors)
        """
        sample_range = np.arange(dist.shape[0])[:, None]
        neigh_ind = np.argpartition(dist, n_neighbors - 1, axis=1)
        neigh_ind = neigh_ind[:, :n_neighbors]
        # argpartition doesn't guarantee sorted order, so we sort again
        neigh_ind = neigh_ind[
            sample_range, np.argsort(dist[sample_range, neigh_ind])]
        if return_distance:
            if self.effective_metric_ == 'euclidean':
                result = np.sqrt(dist[sample_range, neigh_ind]), neigh_ind
            else:
                result = dist[sample_range, neigh_ind], neigh_ind
        else:
            result = neigh_ind
        return result

    def kneighbors(self, X=None, n_neighbors=None, return_distance=True):
        """Finds the K-neighbors of a point.
        Returns indices of and distances to the neighbors of each point.

        Parameters
        ----------
        X : array-like, shape (n_query, n_features), \
                or (n_query, n_indexed) if metric == 'precomputed'
            The query point or points.
            If not provided, neighbors of each indexed point are returned.
            In this case, the query point is not considered its own neighbor.

        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned

        Returns
        -------
        neigh_dist : array, shape (n_query, n_neighbors)
            Array representing the lengths to points, only present if
            return_distance=True

        neigh_ind : array, shape (n_query, n_neighbors)
            Indices of the nearest points in the population matrix.

        Examples
        --------
        In the following example, we construct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(n_neighbors=1)
        >>> neigh.fit(samples) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> print(neigh.kneighbors([[1., 1., 1.]])) # doctest: +ELLIPSIS
        (array([[0.5]]), array([[2]]))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False) # doctest: +ELLIPSIS
        array([[1],
               [2]]...)

        """
        check_is_fitted(self, "_fit_method")

        if n_neighbors is None:
            n_neighbors = self.n_neighbors
        elif n_neighbors <= 0:
            raise ValueError(
                "Expected n_neighbors > 0. Got %d" %
                n_neighbors
            )
        else:
            if not np.issubdtype(type(n_neighbors), np.integer):
                raise TypeError(
                    "n_neighbors does not take %s value, "
                    "enter integer value" %
                    type(n_neighbors))

        if X is not None:
            query_is_train = False
            if self.effective_metric_ == 'precomputed':
                X = _check_precomputed(X)
            else:
                X = check_array(X, accept_sparse='csr')
        else:
            query_is_train = True
            X = self._fit_X
            # Include an extra neighbor to account for the sample itself being
            # returned, which is removed later
            n_neighbors += 1

        train_size = self._fit_X.shape[0]
        if n_neighbors > train_size:
            raise ValueError(
                "Expected n_neighbors <= n_samples, "
                " but n_samples = %d, n_neighbors = %d" %
                (train_size, n_neighbors)
            )

        n_jobs = effective_n_jobs(self.n_jobs)
        chunked_results = None
        if self._fit_method == 'brute':

            if self.effective_metric_ == 'precomputed' and issparse(X):
                results = _kneighbors_from_graph(
                    X, n_neighbors=n_neighbors,
                    return_distance=return_distance)
            else:
                reduce_func = partial(self._kneighbors_reduce_func,
                                      n_neighbors=n_neighbors,
                                      return_distance=return_distance)

                # for efficiency, use squared euclidean distances
                if self.effective_metric_ == 'euclidean':
                    kwds = {'squared': True}
                else:
                    kwds = self.effective_metric_params_

                chunked_results = pairwise_distances_chunked(
                    X, self._fit_X, reduce_func=reduce_func,
                    metric=self.effective_metric_, n_jobs=n_jobs,
                    **kwds)

        elif self._fit_method in ['ball_tree', 'kd_tree']:
            if issparse(X):
                raise ValueError(
                    "%s does not work with sparse matrices. Densify the data, "
                    "or set algorithm='brute'" % self._fit_method)
            if LooseVersion(joblib_version) < LooseVersion('0.12'):
                # Deal with change of API in joblib
                delayed_query = delayed(self._tree.query,
                                        check_pickle=False)
                parallel_kwargs = {"backend": "threading"}
            else:
                delayed_query = delayed(self._tree.query)
                parallel_kwargs = {"prefer": "threads"}
            chunked_results = Parallel(n_jobs, **parallel_kwargs)(
                delayed_query(
                    X[s], n_neighbors, return_distance)
                for s in gen_even_slices(X.shape[0], n_jobs)
            )
        else:
            raise ValueError("internal: _fit_method not recognized")

        if chunked_results is not None:
            if return_distance:
                neigh_dist, neigh_ind = zip(*chunked_results)
                results = np.vstack(neigh_dist), np.vstack(neigh_ind)
            else:
                results = np.vstack(chunked_results)

        if not query_is_train:
            return results
        else:
            # If the query data is the same as the indexed data, we would like
            # to ignore the first nearest neighbor of every sample, i.e
            # the sample itself.
            if return_distance:
                neigh_dist, neigh_ind = results
            else:
                neigh_ind = results

            n_samples, _ = X.shape
            sample_range = np.arange(n_samples)[:, None]
            sample_mask = neigh_ind != sample_range

            # Corner case: When the number of duplicates are more
            # than the number of neighbors, the first NN will not
            # be the sample, but a duplicate.
            # In that case mask the first duplicate.
            dup_gr_nbrs = np.all(sample_mask, axis=1)
            sample_mask[:, 0][dup_gr_nbrs] = False
            neigh_ind = np.reshape(
                neigh_ind[sample_mask], (n_samples, n_neighbors - 1))

            if return_distance:
                neigh_dist = np.reshape(
                    neigh_dist[sample_mask], (n_samples, n_neighbors - 1))
                return neigh_dist, neigh_ind
            return neigh_ind

    def kneighbors_graph(self, X=None, n_neighbors=None,
                         mode='connectivity'):
        """Computes the (weighted) graph of k-Neighbors for points in X

        Parameters
        ----------
        X : array-like, shape (n_query, n_features), \
                or (n_query, n_indexed) if metric == 'precomputed'
            The query point or points.
            If not provided, neighbors of each indexed point are returned.
            In this case, the query point is not considered its own neighbor.

        n_neighbors : int
            Number of neighbors for each sample.
            (default is value passed to the constructor).

        mode : {'connectivity', 'distance'}, optional
            Type of returned matrix: 'connectivity' will return the
            connectivity matrix with ones and zeros, in 'distance' the
            edges are Euclidean distance between points.

        Returns
        -------
        A : sparse matrix in CSR format, shape = [n_query, n_samples_fit]
            n_samples_fit is the number of samples in the fitted data
            A[i, j] is assigned the weight of edge that connects i to j.

        Examples
        --------
        >>> X = [[0], [3], [1]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(n_neighbors=2)
        >>> neigh.fit(X) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> A = neigh.kneighbors_graph(X)
        >>> A.toarray()
        array([[1., 0., 1.],
               [0., 1., 1.],
               [1., 0., 1.]])

        See also
        --------
        NearestNeighbors.radius_neighbors_graph
        """
        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        # check the input only in self.kneighbors

        # construct CSR matrix representation of the k-NN graph
        if mode == 'connectivity':
            A_ind = self.kneighbors(X, n_neighbors, return_distance=False)
            n_samples1 = A_ind.shape[0]
            A_data = np.ones(n_samples1 * n_neighbors)

        elif mode == 'distance':
            A_data, A_ind = self.kneighbors(
                X, n_neighbors, return_distance=True)
            A_data = np.ravel(A_data)

        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity" '
                'or "distance" but got "%s" instead' % mode)

        n_samples1 = A_ind.shape[0]
        n_samples2 = self._fit_X.shape[0]
        n_nonzero = n_samples1 * n_neighbors
        A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

        kneighbors_graph = csr_matrix((A_data, A_ind.ravel(), A_indptr),
                                      shape=(n_samples1, n_samples2))

        return kneighbors_graph


class RadiusNeighborsMixin(object):
    """Mixin for radius-based neighbors searches"""

    def _radius_neighbors_reduce_func(self, dist, start,
                                      radius, return_distance):
        """Reduce a chunk of distances to the nearest neighbors

        Callback to :func:`sklearn.metrics.pairwise.pairwise_distances_chunked`

        Parameters
        ----------
        dist : array of shape (n_samples_chunk, n_samples)
        start : int
            The index in X which the first row of dist corresponds to.
        radius : float
        return_distance : bool

        Returns
        -------
        dist : list of n_samples_chunk 1d arrays, optional
            Returned only if return_distance
        neigh : list of n_samples_chunk 1d arrays
        """
        neigh_ind = [np.where(d <= radius)[0] for d in dist]

        if return_distance:
            if self.effective_metric_ == 'euclidean':
                dist = [np.sqrt(d[neigh_ind[i]])
                        for i, d in enumerate(dist)]
            else:
                dist = [d[neigh_ind[i]]
                        for i, d in enumerate(dist)]
            results = dist, neigh_ind
        else:
            results = neigh_ind
        return results

    def radius_neighbors(self, X=None, radius=None, return_distance=True,
                         sort_results=False):
        """Finds the neighbors within a given radius of a point or points.

        Return the indices and distances of each point from the dataset
        lying in a ball with size ``radius`` around the points of the query
        array. Points lying on the boundary are included in the results.

        The result points are *not* necessarily sorted by distance to their
        query point.

        Parameters
        ----------
        X : array-like, (n_samples, n_features), optional
            The query point or points.
            If not provided, neighbors of each indexed point are returned.
            In this case, the query point is not considered its own neighbor.

        radius : float
            Limiting distance of neighbors to return.
            (default is the value passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned.

        sort_results : boolean, optional. Defaults to False.
            If True, the distances and indices will be sorted before being
            returned. If False, the results will not be sorted. If
            return_distance == False, setting sort_results = True will
            result in an error.

        Returns
        -------
        neigh_dist : array, shape (n_samples,) of arrays
            Array representing the distances to each point, only present if
            return_distance=True. The distance values are computed according
            to the ``metric`` constructor parameter.

        neigh_ind : array, shape (n_samples,) of arrays
            An array of arrays of indices of the approximate nearest points
            from the population matrix that lie within a ball of size
            ``radius`` around the query points.

        Examples
        --------
        In the following example, we construct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1, 1, 1]:

        >>> import numpy as np
        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(radius=1.6)
        >>> neigh.fit(samples) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> rng = neigh.radius_neighbors([[1., 1., 1.]])
        >>> print(np.asarray(rng[0][0])) # doctest: +ELLIPSIS
        [1.5 0.5]
        >>> print(np.asarray(rng[1][0])) # doctest: +ELLIPSIS
        [1 2]

        The first array returned contains the distances to all points which
        are closer than 1.6, while the second array returned contains their
        indices.  In general, multiple points can be queried at the same time.

        Notes
        -----
        Because the number of neighbors of each point is not necessarily
        equal, the results for multiple query points cannot be fit in a
        standard data array.
        For efficiency, `radius_neighbors` returns arrays of objects, where
        each object is a 1D array of indices or distances.
        """
        check_is_fitted(self, "_fit_method")

        if X is not None:
            query_is_train = False
            if self.effective_metric_ == 'precomputed':
                X = _check_precomputed(X)
            else:
                X = check_array(X, accept_sparse='csr')
        else:
            query_is_train = True
            X = self._fit_X

        if radius is None:
            radius = self.radius

        if self._fit_method == 'brute':

            if self.effective_metric_ == 'precomputed' and issparse(X):
                results = _radius_neighbors_from_graph(
                    X, radius=radius, return_distance=return_distance)

            else:
                # for efficiency, use squared euclidean distances
                if self.effective_metric_ == 'euclidean':
                    radius *= radius
                    kwds = {'squared': True}
                else:
                    kwds = self.effective_metric_params_

                reduce_func = partial(self._radius_neighbors_reduce_func,
                                      radius=radius,
                                      return_distance=return_distance)

                chunked_results = pairwise_distances_chunked(
                    X, self._fit_X, reduce_func=reduce_func,
                    metric=self.effective_metric_, n_jobs=self.n_jobs,
                    **kwds)
                if return_distance:
                    neigh_dist_chunks, neigh_ind_chunks = zip(*chunked_results)
                    neigh_dist_list = sum(neigh_dist_chunks, [])
                    neigh_ind_list = sum(neigh_ind_chunks, [])
                    # See https://github.com/numpy/numpy/issues/5456
                    # to understand why this is initialized this way.
                    neigh_dist = np.empty(len(neigh_dist_list), dtype='object')
                    neigh_dist[:] = neigh_dist_list
                    neigh_ind = np.empty(len(neigh_ind_list), dtype='object')
                    neigh_ind[:] = neigh_ind_list
                    results = neigh_dist, neigh_ind
                else:
                    neigh_ind_list = sum(chunked_results, [])
                    results = np.empty(len(neigh_ind_list), dtype='object')
                    results[:] = neigh_ind_list

        elif self._fit_method in ['ball_tree', 'kd_tree']:
            if issparse(X):
                raise ValueError(
                    "%s does not work with sparse matrices. Densify the data, "
                    "or set algorithm='brute'" % self._fit_method)

            n_jobs = effective_n_jobs(self.n_jobs)
            if LooseVersion(joblib_version) < LooseVersion('0.12'):
                # Deal with change of API in joblib
                delayed_query = delayed(self._tree.query_radius,
                                        check_pickle=False)
                parallel_kwargs = {"backend": "threading"}
            else:
                delayed_query = delayed(self._tree.query_radius)
                parallel_kwargs = {"prefer": "threads"}
            chunked_results = Parallel(n_jobs, **parallel_kwargs)(
                delayed_query(X[s], radius, return_distance,
                              sort_results=sort_results)
                for s in gen_even_slices(X.shape[0], n_jobs)
            )
            if return_distance:
                neigh_ind, neigh_dist = tuple(zip(*chunked_results))
                results = np.hstack(neigh_dist), np.hstack(neigh_ind)
            else:
                results = np.hstack(chunked_results)
        else:
            raise ValueError("internal: _fit_method not recognized")

        if not query_is_train:
            return results
        else:
            # If the query data is the same as the indexed data, we would like
            # to ignore the first nearest neighbor of every sample, i.e
            # the sample itself.
            if return_distance:
                neigh_dist, neigh_ind = results
            else:
                neigh_ind = results

            for ind, ind_neighbor in enumerate(neigh_ind):
                mask = ind_neighbor != ind

                neigh_ind[ind] = ind_neighbor[mask]
                if return_distance:
                    neigh_dist[ind] = neigh_dist[ind][mask]

            if return_distance:
                return neigh_dist, neigh_ind
            return neigh_ind

    def radius_neighbors_graph(self, X=None, radius=None, mode='connectivity',
                               sort_results=False):
        """Computes the (weighted) graph of Neighbors for points in X

        Neighborhoods are restricted the points at a distance lower than
        radius.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features], optional
            The query point or points.
            If not provided, neighbors of each indexed point are returned.
            In this case, the query point is not considered its own neighbor.

        radius : float
            Radius of neighborhoods.
            (default is the value passed to the constructor).

        mode : {'connectivity', 'distance'}, optional
            Type of returned matrix: 'connectivity' will return the
            connectivity matrix with ones and zeros, in 'distance' the
            edges are Euclidean distance between points.

        sort_results : boolean, optional. Defaults to False.
            If True, the distances and indices will be sorted before being
            returned. If False, the results will not be sorted.
            Only used with mode='distance'.

        Returns
        -------
        A : sparse matrix in CSR format, shape = [n_samples, n_samples]
            A[i, j] is assigned the weight of edge that connects i to j.

        Examples
        --------
        >>> X = [[0], [3], [1]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(radius=1.5)
        >>> neigh.fit(X) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> A = neigh.radius_neighbors_graph(X)
        >>> A.toarray()
        array([[1., 0., 1.],
               [0., 1., 0.],
               [1., 0., 1.]])

        See also
        --------
        kneighbors_graph
        """
        # check the input only in self.radius_neighbors

        if radius is None:
            radius = self.radius

        # construct CSR matrix representation of the NN graph
        if mode == 'connectivity':
            A_ind = self.radius_neighbors(X, radius,
                                          return_distance=False)
            A_data = None
        elif mode == 'distance':
            dist, A_ind = self.radius_neighbors(X, radius,
                                                return_distance=True,
                                                sort_results=sort_results)
            A_data = np.concatenate(list(dist))
        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity", '
                'or "distance" but got %s instead' % mode)

        n_samples1 = A_ind.shape[0]
        n_samples2 = self._fit_X.shape[0]
        n_neighbors = np.array([len(a) for a in A_ind])
        A_ind = np.concatenate(list(A_ind))
        if A_data is None:
            A_data = np.ones(len(A_ind))
        A_indptr = np.concatenate((np.zeros(1, dtype=int),
                                   np.cumsum(n_neighbors)))

        return csr_matrix((A_data, A_ind, A_indptr),
                          shape=(n_samples1, n_samples2))


class SupervisedFloatMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape [n_samples, n_features],
            or [n_samples, n_samples] if metric='precomputed'.

        y : {array-like, sparse matrix}
            Target values, array of float values, shape = [n_samples]
             or [n_samples, n_outputs]
        """
        if not isinstance(X, (KDTree, BallTree)):
            X, y = check_X_y(X, y, "csr", multi_output=True)
        self._y = y
        return self._fit(X)


class SupervisedIntegerMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape [n_samples, n_features],
            or [n_samples, n_samples] if metric='precomputed'.

        y : {array-like, sparse matrix}
            Target values of shape = [n_samples] or [n_samples, n_outputs]

        """
        if not isinstance(X, (KDTree, BallTree)):
            X, y = check_X_y(X, y, "csr", multi_output=True)

        if y.ndim == 1 or y.ndim == 2 and y.shape[1] == 1:
            if y.ndim != 1:
                warnings.warn("A column-vector y was passed when a 1d array "
                              "was expected. Please change the shape of y to "
                              "(n_samples, ), for example using ravel().",
                              DataConversionWarning, stacklevel=2)

            self.outputs_2d_ = False
            y = y.reshape((-1, 1))
        else:
            self.outputs_2d_ = True

        check_classification_targets(y)
        self.classes_ = []
        self._y = np.empty(y.shape, dtype=np.int)
        for k in range(self._y.shape[1]):
            classes, self._y[:, k] = np.unique(y[:, k], return_inverse=True)
            self.classes_.append(classes)

        if not self.outputs_2d_:
            self.classes_ = self.classes_[0]
            self._y = self._y.ravel()

        return self._fit(X)


class UnsupervisedMixin(object):
    def fit(self, X, y=None):
        """Fit the model using X as training data

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape [n_samples, n_features],
            or [n_samples, n_samples] if metric='precomputed'.
        """
        return self._fit(X)
