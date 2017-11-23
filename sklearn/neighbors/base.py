"""Base and mixin classes for nearest neighbors"""
# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam
import warnings
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import csr_matrix, issparse

from .ball_tree import BallTree
from .kd_tree import KDTree
from ..base import BaseEstimator
from ..metrics import pairwise_distances
from ..metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from ..utils import check_X_y, check_array, _get_n_jobs, gen_even_slices
from ..utils.multiclass import check_classification_targets
from ..externals import six
from ..externals.joblib import Parallel, delayed
from ..exceptions import NotFittedError
from ..exceptions import DataConversionWarning

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


class NeighborsBase(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for nearest neighbors estimators."""

    @abstractmethod
    def __init__(self):
        pass

    def _init_params(self, n_neighbors=None, radius=None,
                     algorithm='auto', leaf_size=30, metric='minkowski',
                     p=2, metric_params=None, n_jobs=1):

        self.n_neighbors = n_neighbors
        self.radius = radius
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.metric = metric
        self.metric_params = metric_params
        self.p = p
        self.n_jobs = n_jobs

        if algorithm not in ['auto', 'brute',
                             'kd_tree', 'ball_tree']:
            raise ValueError("unrecognized algorithm: '%s'" % algorithm)

        if algorithm == 'auto':
            if metric == 'precomputed':
                alg_check = 'brute'
            elif callable(metric) or metric in VALID_METRICS['ball_tree']:
                alg_check = 'ball_tree'
            else:
                alg_check = 'brute'
        else:
            alg_check = algorithm

        if callable(metric):
            if algorithm == 'kd_tree':
                # callable metric is only valid for brute force and ball_tree
                raise ValueError(
                    "kd_tree algorithm does not support callable metric '%s'"
                    % metric)
        elif metric not in VALID_METRICS[alg_check]:
            raise ValueError("Metric '%s' not valid for algorithm '%s'"
                             % (metric, algorithm))

        if self.metric_params is not None and 'p' in self.metric_params:
            warnings.warn("Parameter p is found in metric_params. "
                          "The corresponding parameter from __init__ "
                          "is ignored.", SyntaxWarning, stacklevel=3)
            effective_p = metric_params['p']
        else:
            effective_p = self.p

        if self.metric in ['wminkowski', 'minkowski'] and effective_p < 1:
            raise ValueError("p must be greater than one for minkowski metric")

        self._fit_X = None
        self._tree = None
        self._fit_method = None

    def _fit(self, X):
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

        X = check_array(X, accept_sparse='csr')

        n_samples = X.shape[0]
        if n_samples == 0:
            raise ValueError("n_samples must be greater than 0")

        if issparse(X):
            if self.algorithm not in ('auto', 'brute'):
                warnings.warn("cannot use tree with sparse input: "
                              "using brute force")
            if self.effective_metric_ not in VALID_METRICS_SPARSE['brute']:
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

        return self

    @property
    def _pairwise(self):
        # For cross-validation routines to split data correctly
        return self.metric == 'precomputed'


class KNeighborsMixin(object):
    """Mixin for k-neighbors searches"""

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
        dist : array
            Array representing the lengths to points, only present if
            return_distance=True

        ind : array
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
        (array([[ 0.5]]), array([[2]]...))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False) # doctest: +ELLIPSIS
        array([[1],
               [2]]...)

        """
        if self._fit_method is None:
            raise NotFittedError("Must fit neighbors before querying.")

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        if X is not None:
            query_is_train = False
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
        n_samples, _ = X.shape
        sample_range = np.arange(n_samples)[:, None]

        n_jobs = _get_n_jobs(self.n_jobs)
        if self._fit_method == 'brute':
            # for efficiency, use squared euclidean distances
            if self.effective_metric_ == 'euclidean':
                dist = pairwise_distances(X, self._fit_X, 'euclidean',
                                          n_jobs=n_jobs, squared=True)
            else:
                dist = pairwise_distances(
                    X, self._fit_X, self.effective_metric_, n_jobs=n_jobs,
                    **self.effective_metric_params_)

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

        elif self._fit_method in ['ball_tree', 'kd_tree']:
            if issparse(X):
                raise ValueError(
                    "%s does not work with sparse matrices. Densify the data, "
                    "or set algorithm='brute'" % self._fit_method)
            result = Parallel(n_jobs, backend='threading')(
                delayed(self._tree.query, check_pickle=False)(
                    X[s], n_neighbors, return_distance)
                for s in gen_even_slices(X.shape[0], n_jobs)
            )
            if return_distance:
                dist, neigh_ind = tuple(zip(*result))
                result = np.vstack(dist), np.vstack(neigh_ind)
            else:
                result = np.vstack(result)
        else:
            raise ValueError("internal: _fit_method not recognized")

        if not query_is_train:
            return result
        else:
            # If the query data is the same as the indexed data, we would like
            # to ignore the first nearest neighbor of every sample, i.e
            # the sample itself.
            if return_distance:
                dist, neigh_ind = result
            else:
                neigh_ind = result

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
                dist = np.reshape(
                    dist[sample_mask], (n_samples, n_neighbors - 1))
                return dist, neigh_ind
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
        A : sparse matrix in CSR format, shape = [n_samples, n_samples_fit]
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
        array([[ 1.,  0.,  1.],
               [ 0.,  1.,  1.],
               [ 1.,  0.,  1.]])

        See also
        --------
        NearestNeighbors.radius_neighbors_graph
        """
        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        # kneighbors does the None handling.
        if X is not None:
            X = check_array(X, accept_sparse='csr')
            n_samples1 = X.shape[0]
        else:
            n_samples1 = self._fit_X.shape[0]

        n_samples2 = self._fit_X.shape[0]
        n_nonzero = n_samples1 * n_neighbors
        A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

        # construct CSR matrix representation of the k-NN graph
        if mode == 'connectivity':
            A_data = np.ones(n_samples1 * n_neighbors)
            A_ind = self.kneighbors(X, n_neighbors, return_distance=False)

        elif mode == 'distance':
            A_data, A_ind = self.kneighbors(
                X, n_neighbors, return_distance=True)
            A_data = np.ravel(A_data)

        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity" '
                'or "distance" but got "%s" instead' % mode)

        kneighbors_graph = csr_matrix((A_data, A_ind.ravel(), A_indptr),
                                      shape=(n_samples1, n_samples2))

        return kneighbors_graph


class RadiusNeighborsMixin(object):
    """Mixin for radius-based neighbors searches"""

    def radius_neighbors(self, X=None, radius=None, return_distance=True):
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
            If False, distances will not be returned

        Returns
        -------
        dist : array, shape (n_samples,) of arrays
            Array representing the distances to each point, only present if
            return_distance=True. The distance values are computed according
            to the ``metric`` constructor parameter.

        ind : array, shape (n_samples,) of arrays
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
        [ 1.5  0.5]
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
        if self._fit_method is None:
            raise NotFittedError("Must fit neighbors before querying.")

        if X is not None:
            query_is_train = False
            X = check_array(X, accept_sparse='csr')
        else:
            query_is_train = True
            X = self._fit_X

        if radius is None:
            radius = self.radius

        n_samples = X.shape[0]
        if self._fit_method == 'brute':
            # for efficiency, use squared euclidean distances
            if self.effective_metric_ == 'euclidean':
                dist = pairwise_distances(X, self._fit_X, 'euclidean',
                                          n_jobs=self.n_jobs, squared=True)
                radius *= radius
            else:
                dist = pairwise_distances(X, self._fit_X,
                                          self.effective_metric_,
                                          n_jobs=self.n_jobs,
                                          **self.effective_metric_params_)

            neigh_ind_list = [np.where(d <= radius)[0] for d in dist]

            # See https://github.com/numpy/numpy/issues/5456
            # if you want to understand why this is initialized this way.
            neigh_ind = np.empty(n_samples, dtype='object')
            neigh_ind[:] = neigh_ind_list

            if return_distance:
                dist_array = np.empty(n_samples, dtype='object')
                if self.effective_metric_ == 'euclidean':
                    dist_list = [np.sqrt(d[neigh_ind[i]])
                                 for i, d in enumerate(dist)]
                else:
                    dist_list = [d[neigh_ind[i]]
                                 for i, d in enumerate(dist)]
                dist_array[:] = dist_list

                results = dist_array, neigh_ind
            else:
                results = neigh_ind

        elif self._fit_method in ['ball_tree', 'kd_tree']:
            if issparse(X):
                raise ValueError(
                    "%s does not work with sparse matrices. Densify the data, "
                    "or set algorithm='brute'" % self._fit_method)
            results = self._tree.query_radius(X, radius,
                                              return_distance=return_distance)
            if return_distance:
                results = results[::-1]
        else:
            raise ValueError("internal: _fit_method not recognized")

        if not query_is_train:
            return results
        else:
            # If the query data is the same as the indexed data, we would like
            # to ignore the first nearest neighbor of every sample, i.e
            # the sample itself.
            if return_distance:
                dist, neigh_ind = results
            else:
                neigh_ind = results

            for ind, ind_neighbor in enumerate(neigh_ind):
                mask = ind_neighbor != ind

                neigh_ind[ind] = ind_neighbor[mask]
                if return_distance:
                    dist[ind] = dist[ind][mask]

            if return_distance:
                return dist, neigh_ind
            return neigh_ind

    def radius_neighbors_graph(self, X=None, radius=None, mode='connectivity'):
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
        array([[ 1.,  0.,  1.],
               [ 0.,  1.,  0.],
               [ 1.,  0.,  1.]])

        See also
        --------
        kneighbors_graph
        """
        if X is not None:
            X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])

        n_samples2 = self._fit_X.shape[0]
        if radius is None:
            radius = self.radius

        # construct CSR matrix representation of the NN graph
        if mode == 'connectivity':
            A_ind = self.radius_neighbors(X, radius,
                                          return_distance=False)
            A_data = None
        elif mode == 'distance':
            dist, A_ind = self.radius_neighbors(X, radius,
                                                return_distance=True)
            A_data = np.concatenate(list(dist))
        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity", '
                'or "distance" but got %s instead' % mode)

        n_samples1 = A_ind.shape[0]
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
