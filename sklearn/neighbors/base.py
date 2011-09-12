"""Base and mixin classes for nearest neighbors"""
# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam
import warnings

import re

import numpy as np
from scipy.sparse import csr_matrix, issparse
from scipy.spatial.ckdtree import cKDTree

from .ball_tree import BallTree
from ..base import BaseEstimator
from ..metrics import euclidean_distances
from ..utils import safe_asanyarray, atleast2d_or_csr

ALGORITHM_DESCR = \
    ("algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional\n"
     "    Algorithm used to compute the nearest neighbors:\n"
     "\n"
     "    - 'ball_tree' will use :class:`BallTree`\n"
     "    - 'kd_tree' will use :class:`scipy.spatial.cKDtree`\n"
     "    - 'brute' will use a brute-force search.\n"
     "    - 'auto' will attempt to decide the most appropriate algorithm\n"
     "      based on the values passed to :meth:`fit` method.\n"
     "\n"
     "    Note: fitting on sparse input will override the setting of\n"
     "    this parameter, using brute force.")

LEAF_SIZE_DESCR = \
    ("leaf_size : int, optional (default = 30)\n"
     "    Leaf size passed to BallTree or cKDTree.  This can affect the\n"
     "    speed of the construction and query, as well as the memory\n"
     "    required to store the tree.  The optimal value depends on the\n"
     "    nature of the problem.\n")

NOTES_DESCR = \
    ("Notes\n"
     "-----\n"
     "See :ref:`Nearest Neighbors <neighbors>` in the online documentation\n"
     "for a discussion of the choice of ``algorithm`` and ``leaf_size``.\n"
     "\n"
     "References\n"
     "----------\n"
     "http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm\n")

N_NEIGHBORS_DESCR = \
    ("n_neighbors : int, optional (default = 5)\n"
     "    Number of neighbors to use by default for "
     ":meth:`k_neighbors` queries.\n")

RADIUS_DESCR = \
    ("radius : float, optional (default = 1.0)\n"
     "    Range of parameter space to use by default "
     "for :meth`radius_neighbors`\n"
     "    queries.\n")

WEIGHTS_DESCR = \
    ("weights : str or callable\n"
     "    weight function used in prediction.  Possible values:\n"
     "\n"
     "    - 'uniform' : uniform weights.  All points in each neighborhood\n"
     "      are weighted equally.\n"
     "    - 'distance' : weight points by the inverse of their distance.\n"
     "      in this case, closer neighbors of a query point will have a\n"
     "      greater influence than neighbors which are further away.\n"
     "    - [callable] : a user-defined function which accepts an\n"
     "      array of distances, and returns an array of the same shape\n"
     "      containing the weights.\n"
     "\n"
     "    Uniform weights are used by default.")


NEIGHBORS_DOC_DICT = {'algorithm': ALGORITHM_DESCR,
                      'leaf_size': LEAF_SIZE_DESCR,
                      'notes': NOTES_DESCR,
                      'n_neighbors': N_NEIGHBORS_DESCR,
                      'radius': RADIUS_DESCR,
                      'weights': WEIGHTS_DESCR}


def construct_docstring(template,
                        flag='@INCLUDE',
                        doc_dict=NEIGHBORS_DOC_DICT,
                        indent=4):
    """
    Construct a documentation string using information from mixins and
    subclasses.

    Parameters
    ----------
    template : string
        This string gives the template for the docstring.  Anywhere the
        characters in `flag` are present, the rest of the line will be
        replaced with the argument of `flag`
    flag : string
        The character string indicating when a replacement should be made
    doc_dict : dict object
        dictionary of replacement strings for the documentation
    indent : int
        The indentation level of the documentation.  Default is 4.

    Returns
    -------
    docstr : string
    """
    # find all instances of the flag in the template
    indices = [match.start() for match in re.finditer(flag, template)]

    # build the final docstring from the template
    docstr = ""
    i_last = 0
    for i in indices:
        docstr += template[i_last:i]

        current = template[i:]

        i_endline = current.find('\n')

        key = current[:i_endline].strip(flag).strip()

        docstr += doc_dict[key].strip().replace('\n', '\n' + indent * ' ')
        docstr += '\n'

        i_last = i + i_endline + 1

    docstr += template[i_last:]

    return docstr


def _check_weights(weights):
    """
    check to make sure weights are valid
    """
    if weights in (None, 'uniform', 'distance'):
        return weights
    elif callable(weights):
        return weights
    else:
        raise ValueError("weights not recognized: should be 'uniform', "
                         "'distance', or a callable function")


def _get_weights(dist, weights):
    """
    get the weights from an array of distances and a parameter ``weights``,
    which can be either a string or an executable

    returns ``weights_arr``, an array of the same size as ``dist``
    if ``weights == 'uniform'``, then returns None
    """
    if dist.dtype == np.dtype(object):
        if weights in (None, 'uniform'):
            return None
        elif weights == 'distance':
            return [1. / d for d in dist]
        elif callable(weights):
            return [weights(d) for d in dist]
        else:
            raise ValueError("weights not recognized: should be 'uniform', "
                             "'distance', or a callable function")
    else:
        if weights in (None, 'uniform'):
            return None
        elif weights == 'distance':
            return 1. / dist
        elif callable(weights):
            return weights(dist)
        else:
            raise ValueError("weights not recognized: should be 'uniform', "
                             "'distance', or a callable function")


class NeighborsBase(BaseEstimator):
    """Base class for nearest neighbors estimators."""
    #FIXME: include float parameter p for using different distance metrics.
    # this can be passed directly to BallTree and cKDTree.  Brute-force will
    # rely on soon-to-be-updated functionality in the pairwise module.
    def _init_params(self, n_neighbors=None, radius=None,
                     algorithm='auto', leaf_size=30):
        self.n_neighbors = n_neighbors
        self.radius = radius
        self.algorithm = algorithm
        self.leaf_size = leaf_size

        if algorithm not in ['auto', 'brute', 'kd_tree', 'ball_tree']:
            raise ValueError("unrecognized algorithm")

        self._fit_X = None
        self._tree = None
        self._fit_method = None

    def _fit(self, X):
        if isinstance(X, BallTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'ball_tree'
            return self

        if isinstance(X, cKDTree):
            self._fit_X = X.data
            self._tree = X
            self._fit_method = 'kd_tree'
            return self

        X = safe_asanyarray(X)

        if issparse(X):
            if self.algorithm not in ('auto', 'brute'):
                warnings.warn("cannot use tree with sparse input: "
                              "using brute force")
            self._fit_X = X.tocsr()
            self._tree = None
            self._fit_method = 'brute'
            return self

        self._fit_method = self.algorithm
        self._fit_X = X

        if self._fit_method == 'auto':
            # BallTree outperforms the others in nearly any circumstance.
            self._fit_method = 'ball_tree'

        if self._fit_method == 'kd_tree':
            self._tree = cKDTree(X, self.leaf_size)
        elif self._fit_method == 'ball_tree':
            self._tree = BallTree(X, self.leaf_size)
        elif self._fit_method == 'brute':
            self._tree = None
        else:
            raise ValueError("algorithm = '%s' not recognized"
                             % self.algorithm)
        return self


class KNeighborsMixin(object):
    """Mixin for k-neighbors searches"""

    def kneighbors(self, X, n_neighbors=None, return_distance=True):
        """Finds the K-neighbors of a point.

        Returns distance

        Parameters
        ----------
        X : array-like, last dimension same as that of fit data
            The new point.

        n_neighbors : int
            Number of neighbors to get (default is the value
            passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned

        Returns
        -------
        dist : array
            Array representing the lengths to point, only present if
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
        >>> print neigh.kneighbors([1., 1., 1.]) # doctest: +ELLIPSIS
        (array([[ 0.5]]), array([[2]]...))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False) # doctest: +ELLIPSIS
        array([[1],
               [2]]...)

        """
        if self._fit_method == None:
            raise ValueError("must fit neighbors before querying")

        X = atleast2d_or_csr(X)

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        if self._fit_method == 'brute':
            dist = euclidean_distances(X, self._fit_X, squared=True)
            # XXX: should be implemented with a partial sort
            neigh_ind = dist.argsort(axis=1)[:, :n_neighbors]
            if return_distance:
                j = np.arange(neigh_ind.shape[0])[:, None]
                return np.sqrt(dist[j, neigh_ind]), neigh_ind
            else:
                return neigh_ind
        elif self._fit_method == 'ball_tree':
            return self._tree.query(X, n_neighbors,
                                    return_distance=return_distance)
        elif self._fit_method == 'kd_tree':
            dist, ind = self._tree.query(X, n_neighbors)
            # kd_tree returns a 1D array for n_neighbors = 1
            if n_neighbors == 1:
                dist = dist[:, None]
                ind = ind[:, None]
            if return_distance:
                return dist, ind
            else:
                return ind
        else:
            raise ValueError("internal: _fit_method not recognized")

    def kneighbors_graph(self, X, n_neighbors=None,
                         mode='connectivity'):
        """Computes the (weighted) graph of k-Neighbors for points in X

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Sample data

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
        >>> A.todense()
        matrix([[ 1.,  0.,  1.],
                [ 0.,  1.,  1.],
                [ 1.,  0.,  1.]])

        See also
        --------
        NearestNeighbors.radius_neighbors_graph
        """
        X = np.asanyarray(X)

        if n_neighbors is None:
            n_neighbors = self.n_neighbors

        n_samples1 = X.shape[0]
        n_samples2 = self._fit_X.shape[0]
        n_nonzero = n_samples1 * n_neighbors
        A_indptr = np.arange(0, n_nonzero + 1, n_neighbors)

        # construct CSR matrix representation of the k-NN graph
        if mode == 'connectivity':
            A_data = np.ones((n_samples1, n_neighbors))
            A_ind = self.kneighbors(X, n_neighbors, return_distance=False)

        elif mode == 'distance':
            data, ind = self.kneighbors(X, n_neighbors + 1,
                                        return_distance=True)
            A_data, A_ind = data[:, 1:], ind[:, 1:]

        else:
            raise ValueError(
                'Unsupported mode, must be one of "connectivity" '
                'or "distance" but got "%s" instead' % mode)

        return csr_matrix((A_data.ravel(), A_ind.ravel(), A_indptr),
                          shape=(n_samples1, n_samples2))


class RadiusNeighborsMixin(object):
    """Mixin for radius-based neighbors searches"""

    def radius_neighbors(self, X, radius=None, return_distance=True):
        """Finds the neighbors of a point within a given radius.

        Returns distance

        Parameters
        ----------
        X : array-like, last dimension same as that of fit data
            The new point.

        radius : float
            Limiting distance of neighbors to return.
            (default is the value passed to the constructor).

        return_distance : boolean, optional. Defaults to True.
            If False, distances will not be returned

        Returns
        -------
        dist : array
            Array representing the lengths to point, only present if
            return_distance=True

        ind : array
            Indices of the nearest points in the population matrix.

        Examples
        --------
        In the following example, we construnct a NeighborsClassifier
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(radius=1.6)
        >>> neigh.fit(samples) # doctest: +ELLIPSIS
        NearestNeighbors(algorithm='auto', leaf_size=30, ...)
        >>> print neigh.radius_neighbors([1., 1., 1.]) # doctest: +ELLIPSIS
        (array([[ 1.5  0.5]]...), array([[1 2]]...)

        The first array returned contains the distances to all points which
        are closer than 1.6, while the second array returned contains their
        indices.  In general, multiple points can be queried at the same time.
        Because the number of neighbors of each point is not necessarily
        equal, `radius_neighbors` returns an array of objects, where each
        object is a 1D array of indices.

        """
        if self._fit_method == None:
            raise ValueError("must fit neighbors before querying")

        X = atleast2d_or_csr(X)

        if radius is None:
            radius = self.radius

        if self._fit_method == 'brute':
            dist = euclidean_distances(X, self._fit_X, squared=True)
            rad2 = radius ** 2

            neigh_ind = [np.where(d < rad2)[0] for d in dist]

            # if there are the same number of neighbors for each point,
            # we can do a normal array.  Otherwise, we return an object
            # array with elements that are numpy arrays
            try:
                neigh_ind = np.asarray(neigh_ind, dtype=int)
                dtype_F = float
            except ValueError:
                neigh_ind = np.asarray(neigh_ind, dtype='object')
                dtype_F = object

            if return_distance:
                dist = np.array([np.sqrt(d[neigh_ind[i]]) \
                                     for i, d in enumerate(dist)],
                                dtype=dtype_F)
                return dist, neigh_ind
            else:
                return neigh_ind
        elif self._fit_method == 'ball_tree':
            if return_distance:
                ind, dist = self._tree.query_radius(X, radius,
                                                    return_distance=True)
                return dist, ind
            else:
                ind = self._tree.query_radius(X, radius,
                                              return_distance=False)
                return ind
        elif self._fit_method == 'kd_tree':
            Npts = self._fit_X.shape[0]
            dist, ind = self._tree.query(X, Npts,
                                         distance_upper_bound=radius)

            ind = [ind_i[:ind_i.searchsorted(Npts)] for ind_i in ind]

            # if there are the same number of neighbors for each point,
            # we can do a normal array.  Otherwise, we return an object
            # array with elements that are numpy arrays
            try:
                ind = np.asarray(ind, dtype=int)
                dtype_F = float
            except ValueError:
                ind = np.asarray(ind, dtype='object')
                dtype_F = object

            if return_distance:
                dist = np.array([dist_i[:len(ind[i])]
                                 for i, dist_i in enumerate(dist)],
                                dtype=dtype_F)
                return dist, ind
            else:
                return ind
        else:
            raise ValueError("internal: _fit_method not recognized")

    def radius_neighbors_graph(self, X, radius=None, mode='connectivity'):
        """Computes the (weighted) graph of Neighbors for points in X

        Neighborhoods are restricted the points at a distance lower than
        radius.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Sample data

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
        >>> A.todense()
        matrix([[ 1.,  0.,  1.],
                [ 0.,  1.,  0.],
                [ 1.,  0.,  1.]])

        See also
        --------
        kneighbors_graph
        """
        X = np.asanyarray(X)

        if radius is None:
            radius = self.radius

        n_samples1 = X.shape[0]
        n_samples2 = self._fit_X.shape[0]

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

        n_neighbors = np.array([len(a) for a in A_ind])
        n_nonzero = np.sum(n_neighbors)
        if A_data is None:
            A_data = np.ones(n_nonzero)
        A_ind = np.concatenate(list(A_ind))
        A_indptr = np.concatenate((np.zeros(1, dtype=int),
                                   np.cumsum(n_neighbors)))

        return csr_matrix((A_data, A_ind, A_indptr),
                          shape=(n_samples1, n_samples2))


class SupervisedFloatMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, cKDTree}
            Training data. If array or matrix, then the shape
            is [n_samples, n_features]

        y : {array-like, sparse matrix}, shape = [n_samples]
            Target values, array of float values.
        """
        self._y = np.asanyarray(y)
        return self._fit(X)


class SupervisedIntegerMixin(object):
    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, cKDTree}
            Training data. If array or matrix, then the shape
            is [n_samples, n_features]

        y : {array-like, sparse matrix}, shape = [n_samples]
            Target values, array of integer values.
        """
        self._y = np.asanyarray(y)
        return self._fit(X)


class UnsupervisedMixin(object):
    def fit(self, X, y=None):
        """Fit the model using X as training data

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, cKDTree}
            Training data. If array or matrix, shape = [n_samples, n_features]
        """
        return self._fit(X)
