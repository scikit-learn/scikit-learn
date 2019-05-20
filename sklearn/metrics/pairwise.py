# -*- coding: utf-8 -*-

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Robert Layton <robertlayton@gmail.com>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
#          Philippe Gervais <philippe.gervais@inria.fr>
#          Lars Buitinck
#          Joel Nothman <joel.nothman@gmail.com>
# License: BSD 3 clause

import itertools
from functools import partial
import warnings

import numpy as np
from scipy.spatial import distance
from scipy.sparse import csr_matrix
from scipy.sparse import issparse

from ..utils.validation import _num_samples
from ..utils.validation import check_non_negative
from ..utils import check_array
from ..utils import gen_even_slices
from ..utils import gen_batches, get_chunk_n_rows
from ..utils.extmath import row_norms, safe_sparse_dot
from ..preprocessing import normalize
from ..utils._joblib import Parallel
from ..utils._joblib import delayed
from ..utils._joblib import effective_n_jobs

from .pairwise_fast import _chi2_kernel_fast, _sparse_manhattan
from ..exceptions import DataConversionWarning


# Utility Functions
def _return_float_dtype(X, Y):
    """
    1. If dtype of X and Y is float32, then dtype float32 is returned.
    2. Else dtype float is returned.
    """
    if not issparse(X) and not isinstance(X, np.ndarray):
        X = np.asarray(X)

    if Y is None:
        Y_dtype = X.dtype
    elif not issparse(Y) and not isinstance(Y, np.ndarray):
        Y = np.asarray(Y)
        Y_dtype = Y.dtype
    else:
        Y_dtype = Y.dtype

    if X.dtype == Y_dtype == np.float32:
        dtype = np.float32
    else:
        dtype = np.float

    return X, Y, dtype


def check_pairwise_arrays(X, Y, precomputed=False, dtype=None):
    """ Set X and Y appropriately and checks inputs

    If Y is None, it is set as a pointer to X (i.e. not a copy).
    If Y is given, this does not happen.
    All distance metrics should use this function first to assert that the
    given parameters are correct and safe to use.

    Specifically, this function first ensures that both X and Y are arrays,
    then checks that they are at least two dimensional while ensuring that
    their elements are floats (or dtype if provided). Finally, the function
    checks that the size of the second dimension of the two arrays is equal, or
    the equivalent check for a precomputed distance matrix.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples_a, n_features)

    Y : {array-like, sparse matrix}, shape (n_samples_b, n_features)

    precomputed : bool
        True if X is to be treated as precomputed distances to the samples in
        Y.

    dtype : string, type, list of types or None (default=None)
        Data type required for X and Y. If None, the dtype will be an
        appropriate float type selected by _return_float_dtype.

        .. versionadded:: 0.18

    Returns
    -------
    safe_X : {array-like, sparse matrix}, shape (n_samples_a, n_features)
        An array equal to X, guaranteed to be a numpy array.

    safe_Y : {array-like, sparse matrix}, shape (n_samples_b, n_features)
        An array equal to Y if Y was not None, guaranteed to be a numpy array.
        If Y was None, safe_Y will be a pointer to X.

    """
    X, Y, dtype_float = _return_float_dtype(X, Y)

    estimator = 'check_pairwise_arrays'
    if dtype is None:
        dtype = dtype_float

    if Y is X or Y is None:
        X = Y = check_array(X, accept_sparse='csr', dtype=dtype,
                            estimator=estimator)
    else:
        X = check_array(X, accept_sparse='csr', dtype=dtype,
                        estimator=estimator)
        Y = check_array(Y, accept_sparse='csr', dtype=dtype,
                        estimator=estimator)

    if precomputed:
        if X.shape[1] != Y.shape[0]:
            raise ValueError("Precomputed metric requires shape "
                             "(n_queries, n_indexed). Got (%d, %d) "
                             "for %d indexed." %
                             (X.shape[0], X.shape[1], Y.shape[0]))
    elif X.shape[1] != Y.shape[1]:
        raise ValueError("Incompatible dimension for X and Y matrices: "
                         "X.shape[1] == %d while Y.shape[1] == %d" % (
                             X.shape[1], Y.shape[1]))

    return X, Y


def check_paired_arrays(X, Y):
    """ Set X and Y appropriately and checks inputs for paired distances

    All paired distance metrics should use this function first to assert that
    the given parameters are correct and safe to use.

    Specifically, this function first ensures that both X and Y are arrays,
    then checks that they are at least two dimensional while ensuring that
    their elements are floats. Finally, the function checks that the size
    of the dimensions of the two arrays are equal.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples_a, n_features)

    Y : {array-like, sparse matrix}, shape (n_samples_b, n_features)

    Returns
    -------
    safe_X : {array-like, sparse matrix}, shape (n_samples_a, n_features)
        An array equal to X, guaranteed to be a numpy array.

    safe_Y : {array-like, sparse matrix}, shape (n_samples_b, n_features)
        An array equal to Y if Y was not None, guaranteed to be a numpy array.
        If Y was None, safe_Y will be a pointer to X.

    """
    X, Y = check_pairwise_arrays(X, Y)
    if X.shape != Y.shape:
        raise ValueError("X and Y should be of same shape. They were "
                         "respectively %r and %r long." % (X.shape, Y.shape))
    return X, Y


# Pairwise distances
def euclidean_distances(X, Y=None, Y_norm_squared=None, squared=False,
                        X_norm_squared=None):
    """
    Considering the rows of X (and Y=X) as vectors, compute the
    distance matrix between each pair of vectors.

    For efficiency reasons, the euclidean distance between a pair of row
    vector x and y is computed as::

        dist(x, y) = sqrt(dot(x, x) - 2 * dot(x, y) + dot(y, y))

    This formulation has two advantages over other ways of computing distances.
    First, it is computationally efficient when dealing with sparse data.
    Second, if one argument varies but the other remains unchanged, then
    `dot(x, x)` and/or `dot(y, y)` can be pre-computed.

    However, this is not the most precise way of doing this computation, and
    the distance matrix returned by this function may not be exactly
    symmetric as required by, e.g., ``scipy.spatial.distance`` functions.

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples_1, n_features)

    Y : {array-like, sparse matrix}, shape (n_samples_2, n_features)

    Y_norm_squared : array-like, shape (n_samples_2, ), optional
        Pre-computed dot-products of vectors in Y (e.g.,
        ``(Y**2).sum(axis=1)``)
        May be ignored in some cases, see the note below.

    squared : boolean, optional
        Return squared Euclidean distances.

    X_norm_squared : array-like, shape = [n_samples_1], optional
        Pre-computed dot-products of vectors in X (e.g.,
        ``(X**2).sum(axis=1)``)
        May be ignored in some cases, see the note below.

    Notes
    -----
    To achieve better accuracy, `X_norm_squared` and `Y_norm_squared` may be
    unused if they are passed as ``float32``.

    Returns
    -------
    distances : array, shape (n_samples_1, n_samples_2)

    Examples
    --------
    >>> from sklearn.metrics.pairwise import euclidean_distances
    >>> X = [[0, 1], [1, 1]]
    >>> # distance between rows of X
    >>> euclidean_distances(X, X)
    array([[0., 1.],
           [1., 0.]])
    >>> # get distance to origin
    >>> euclidean_distances(X, [[0, 0]])
    array([[1.        ],
           [1.41421356]])

    See also
    --------
    paired_distances : distances betweens pairs of elements of X and Y.
    """
    X, Y = check_pairwise_arrays(X, Y)

    # If norms are passed as float32, they are unused. If arrays are passed as
    # float32, norms needs to be recomputed on upcast chunks.
    # TODO: use a float64 accumulator in row_norms to avoid the latter.
    if X_norm_squared is not None:
        XX = check_array(X_norm_squared)
        if XX.shape == (1, X.shape[0]):
            XX = XX.T
        elif XX.shape != (X.shape[0], 1):
            raise ValueError(
                "Incompatible dimensions for X and X_norm_squared")
        if XX.dtype == np.float32:
            XX = None
    elif X.dtype == np.float32:
        XX = None
    else:
        XX = row_norms(X, squared=True)[:, np.newaxis]

    if X is Y and XX is not None:
        # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    elif Y_norm_squared is not None:
        YY = np.atleast_2d(Y_norm_squared)

        if YY.shape != (1, Y.shape[0]):
            raise ValueError(
                "Incompatible dimensions for Y and Y_norm_squared")
        if YY.dtype == np.float32:
            YY = None
    elif Y.dtype == np.float32:
        YY = None
    else:
        YY = row_norms(Y, squared=True)[np.newaxis, :]

    if X.dtype == np.float32:
        # To minimize precision issues with float32, we compute the distance
        # matrix on chunks of X and Y upcast to float64
        distances = _euclidean_distances_upcast(X, XX, Y, YY)
    else:
        # if dtype is already float64, no need to chunk and upcast
        distances = - 2 * safe_sparse_dot(X, Y.T, dense_output=True)
        distances += XX
        distances += YY
    np.maximum(distances, 0, out=distances)

    # Ensure that distances between vectors and themselves are set to 0.0.
    # This may not be the case due to floating point rounding errors.
    if X is Y:
        np.fill_diagonal(distances, 0)

    return distances if squared else np.sqrt(distances, out=distances)


def _euclidean_distances_upcast(X, XX=None, Y=None, YY=None, batch_size=None):
    """Euclidean distances between X and Y

    Assumes X and Y have float32 dtype.
    Assumes XX and YY have float64 dtype or are None.

    X and Y are upcast to float64 by chunks, which size is chosen to limit
    memory increase by approximately 10% (at least 10MiB).
    """
    n_samples_X = X.shape[0]
    n_samples_Y = Y.shape[0]
    n_features = X.shape[1]

    distances = np.empty((n_samples_X, n_samples_Y), dtype=np.float32)

    if batch_size is None:
        x_density = X.nnz / np.prod(X.shape) if issparse(X) else 1
        y_density = Y.nnz / np.prod(Y.shape) if issparse(Y) else 1

        # Allow 10% more memory than X, Y and the distance matrix take (at
        # least 10MiB)
        maxmem = max(
            ((x_density * n_samples_X + y_density * n_samples_Y) * n_features
             + (x_density * n_samples_X * y_density * n_samples_Y)) / 10,
            10 * 2 ** 17)

        # The increase amount of memory in 8-byte blocks is:
        # - x_density * batch_size * n_features (copy of chunk of X)
        # - y_density * batch_size * n_features (copy of chunk of Y)
        # - batch_size * batch_size (chunk of distance matrix)
        # Hence x² + (xd+yd)kx = M, where x=batch_size, k=n_features, M=maxmem
        #                                 xd=x_density and yd=y_density
        tmp = (x_density + y_density) * n_features
        batch_size = (-tmp + np.sqrt(tmp ** 2 + 4 * maxmem)) / 2
        batch_size = max(int(batch_size), 1)

    x_batches = gen_batches(n_samples_X, batch_size)

    for i, x_slice in enumerate(x_batches):
        X_chunk = X[x_slice].astype(np.float64)
        if XX is None:
            XX_chunk = row_norms(X_chunk, squared=True)[:, np.newaxis]
        else:
            XX_chunk = XX[x_slice]

        y_batches = gen_batches(n_samples_Y, batch_size)

        for j, y_slice in enumerate(y_batches):
            if X is Y and j < i:
                # when X is Y the distance matrix is symmetric so we only need
                # to compute half of it.
                d = distances[y_slice, x_slice].T

            else:
                Y_chunk = Y[y_slice].astype(np.float64)
                if YY is None:
                    YY_chunk = row_norms(Y_chunk, squared=True)[np.newaxis, :]
                else:
                    YY_chunk = YY[:, y_slice]

                d = -2 * safe_sparse_dot(X_chunk, Y_chunk.T, dense_output=True)
                d += XX_chunk
                d += YY_chunk

            distances[x_slice, y_slice] = d.astype(np.float32, copy=False)

    return distances


def _argmin_min_reduce(dist, start):
    indices = dist.argmin(axis=1)
    values = dist[np.arange(dist.shape[0]), indices]
    return indices, values


def pairwise_distances_argmin_min(X, Y, axis=1, metric="euclidean",
                                  batch_size=None, metric_kwargs=None):
    """Compute minimum distances between one point and a set of points.

    This function computes for each row in X, the index of the row of Y which
    is closest (according to the specified distance). The minimal distances are
    also returned.

    This is mostly equivalent to calling:

        (pairwise_distances(X, Y=Y, metric=metric).argmin(axis=axis),
         pairwise_distances(X, Y=Y, metric=metric).min(axis=axis))

    but uses much less memory, and is faster for large arrays.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples1, n_features)
        Array containing points.

    Y : {array-like, sparse matrix}, shape (n_samples2, n_features)
        Arrays containing points.

    axis : int, optional, default 1
        Axis along which the argmin and distances are to be computed.

    metric : string or callable, default 'euclidean'
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    batch_size : integer
        .. deprecated:: 0.20
            Deprecated for removal in 0.22.
            Use sklearn.set_config(working_memory=...) instead.

    metric_kwargs : dict, optional
        Keyword arguments to pass to specified metric function.

    Returns
    -------
    argmin : numpy.ndarray
        Y[argmin[i], :] is the row in Y that is closest to X[i, :].

    distances : numpy.ndarray
        distances[i] is the distance between the i-th row in X and the
        argmin[i]-th row in Y.

    See also
    --------
    sklearn.metrics.pairwise_distances
    sklearn.metrics.pairwise_distances_argmin
    """
    if batch_size is not None:
        warnings.warn("'batch_size' is ignored. It was deprecated in version "
                      "0.20 and will be removed in version 0.22. "
                      "Use sklearn.set_config(working_memory=...) instead.",
                      DeprecationWarning)
    X, Y = check_pairwise_arrays(X, Y)

    if metric_kwargs is None:
        metric_kwargs = {}

    if axis == 0:
        X, Y = Y, X

    indices, values = zip(*pairwise_distances_chunked(
        X, Y, reduce_func=_argmin_min_reduce, metric=metric,
        **metric_kwargs))
    indices = np.concatenate(indices)
    values = np.concatenate(values)

    return indices, values


def pairwise_distances_argmin(X, Y, axis=1, metric="euclidean",
                              batch_size=None, metric_kwargs=None):
    """Compute minimum distances between one point and a set of points.

    This function computes for each row in X, the index of the row of Y which
    is closest (according to the specified distance).

    This is mostly equivalent to calling:

        pairwise_distances(X, Y=Y, metric=metric).argmin(axis=axis)

    but uses much less memory, and is faster for large arrays.

    This function works with dense 2D arrays only.

    Parameters
    ----------
    X : array-like
        Arrays containing points. Respective shapes (n_samples1, n_features)
        and (n_samples2, n_features)

    Y : array-like
        Arrays containing points. Respective shapes (n_samples1, n_features)
        and (n_samples2, n_features)

    axis : int, optional, default 1
        Axis along which the argmin and distances are to be computed.

    metric : string or callable
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    batch_size : integer
        .. deprecated:: 0.20
            Deprecated for removal in 0.22.
            Use sklearn.set_config(working_memory=...) instead.

    metric_kwargs : dict
        keyword arguments to pass to specified metric function.

    Returns
    -------
    argmin : numpy.ndarray
        Y[argmin[i], :] is the row in Y that is closest to X[i, :].

    See also
    --------
    sklearn.metrics.pairwise_distances
    sklearn.metrics.pairwise_distances_argmin_min
    """
    if metric_kwargs is None:
        metric_kwargs = {}

    return pairwise_distances_argmin_min(X, Y, axis, metric,
                                         metric_kwargs=metric_kwargs,
                                         batch_size=batch_size)[0]


def haversine_distances(X, Y=None):
    """Compute the Haversine distance between samples in X and Y

    The Haversine (or great circle) distance is the angular distance between
    two points on the surface of a sphere. The first distance of each point is
    assumed to be the latitude, the second is the longitude, given in radians.
    The dimension of the data must be 2.

    .. math::
       D(x, y) = 2\\arcsin[\\sqrt{\\sin^2((x1 - y1) / 2)
                                + \\cos(x1)\\cos(y1)\\sin^2((x2 - y2) / 2)}]

    Parameters
    ----------
    X : array_like, shape (n_samples_1, 2)

    Y : array_like, shape (n_samples_2, 2), optional

    Returns
    -------
    distance : {array}, shape (n_samples_1, n_samples_2)

    Notes
    -----
    As the Earth is nearly spherical, the haversine formula provides a good
    approximation of the distance between two points of the Earth surface, with
    a less than 1% error on average.

    Examples
    --------
    We want to calculate the distance between the Ezeiza Airport
    (Buenos Aires, Argentina) and the Charles de Gaulle Airport (Paris, France)

    >>> from sklearn.metrics.pairwise import haversine_distances
    >>> bsas = [-34.83333, -58.5166646]
    >>> paris = [49.0083899664, 2.53844117956]
    >>> result = haversine_distances([bsas, paris])
    >>> result * 6371000/1000  # multiply by Earth radius to get kilometers
    array([[    0.        , 11279.45379464],
           [11279.45379464,     0.        ]])
    """
    from sklearn.neighbors import DistanceMetric
    return DistanceMetric.get_metric('haversine').pairwise(X, Y)


def manhattan_distances(X, Y=None, sum_over_features=True):
    """ Compute the L1 distances between the vectors in X and Y.

    With sum_over_features equal to False it returns the componentwise
    distances.

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array_like
        An array with shape (n_samples_X, n_features).

    Y : array_like, optional
        An array with shape (n_samples_Y, n_features).

    sum_over_features : bool, default=True
        If True the function returns the pairwise distance matrix
        else it returns the componentwise L1 pairwise-distances.
        Not supported for sparse matrix inputs.

    Returns
    -------
    D : array
        If sum_over_features is False shape is
        (n_samples_X * n_samples_Y, n_features) and D contains the
        componentwise L1 pairwise-distances (ie. absolute difference),
        else shape is (n_samples_X, n_samples_Y) and D contains
        the pairwise L1 distances.

    Examples
    --------
    >>> from sklearn.metrics.pairwise import manhattan_distances
    >>> manhattan_distances([[3]], [[3]])#doctest:+ELLIPSIS
    array([[0.]])
    >>> manhattan_distances([[3]], [[2]])#doctest:+ELLIPSIS
    array([[1.]])
    >>> manhattan_distances([[2]], [[3]])#doctest:+ELLIPSIS
    array([[1.]])
    >>> manhattan_distances([[1, 2], [3, 4]],\
         [[1, 2], [0, 3]])#doctest:+ELLIPSIS
    array([[0., 2.],
           [4., 4.]])
    >>> import numpy as np
    >>> X = np.ones((1, 2))
    >>> y = np.full((2, 2), 2.)
    >>> manhattan_distances(X, y, sum_over_features=False)#doctest:+ELLIPSIS
    array([[1., 1.],
           [1., 1.]])
    """
    X, Y = check_pairwise_arrays(X, Y)

    if issparse(X) or issparse(Y):
        if not sum_over_features:
            raise TypeError("sum_over_features=%r not supported"
                            " for sparse matrices" % sum_over_features)

        X = csr_matrix(X, copy=False)
        Y = csr_matrix(Y, copy=False)
        D = np.zeros((X.shape[0], Y.shape[0]))
        _sparse_manhattan(X.data, X.indices, X.indptr,
                          Y.data, Y.indices, Y.indptr,
                          X.shape[1], D)
        return D

    if sum_over_features:
        return distance.cdist(X, Y, 'cityblock')

    D = X[:, np.newaxis, :] - Y[np.newaxis, :, :]
    D = np.abs(D, D)
    return D.reshape((-1, X.shape[1]))


def cosine_distances(X, Y=None):
    """Compute cosine distance between samples in X and Y.

    Cosine distance is defined as 1.0 minus the cosine similarity.

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array_like, sparse matrix
        with shape (n_samples_X, n_features).

    Y : array_like, sparse matrix (optional)
        with shape (n_samples_Y, n_features).

    Returns
    -------
    distance matrix : array
        An array with shape (n_samples_X, n_samples_Y).

    See also
    --------
    sklearn.metrics.pairwise.cosine_similarity
    scipy.spatial.distance.cosine : dense matrices only
    """
    # 1.0 - cosine_similarity(X, Y) without copy
    S = cosine_similarity(X, Y)
    S *= -1
    S += 1
    np.clip(S, 0, 2, out=S)
    if X is Y or Y is None:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        S[np.diag_indices_from(S)] = 0.0
    return S


# Paired distances
def paired_euclidean_distances(X, Y):
    """
    Computes the paired euclidean distances between X and Y

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    Y : array-like, shape (n_samples, n_features)

    Returns
    -------
    distances : ndarray (n_samples, )
    """
    X, Y = check_paired_arrays(X, Y)
    return row_norms(X - Y)


def paired_manhattan_distances(X, Y):
    """Compute the L1 distances between the vectors in X and Y.

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    Y : array-like, shape (n_samples, n_features)

    Returns
    -------
    distances : ndarray (n_samples, )
    """
    X, Y = check_paired_arrays(X, Y)
    diff = X - Y
    if issparse(diff):
        diff.data = np.abs(diff.data)
        return np.squeeze(np.array(diff.sum(axis=1)))
    else:
        return np.abs(diff).sum(axis=-1)


def paired_cosine_distances(X, Y):
    """
    Computes the paired cosine distances between X and Y

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    Y : array-like, shape (n_samples, n_features)

    Returns
    -------
    distances : ndarray, shape (n_samples, )

    Notes
    -----
    The cosine distance is equivalent to the half the squared
    euclidean distance if each sample is normalized to unit norm
    """
    X, Y = check_paired_arrays(X, Y)
    return .5 * row_norms(normalize(X) - normalize(Y), squared=True)


PAIRED_DISTANCES = {
    'cosine': paired_cosine_distances,
    'euclidean': paired_euclidean_distances,
    'l2': paired_euclidean_distances,
    'l1': paired_manhattan_distances,
    'manhattan': paired_manhattan_distances,
    'cityblock': paired_manhattan_distances}


def paired_distances(X, Y, metric="euclidean", **kwds):
    """
    Computes the paired distances between X and Y.

    Computes the distances between (X[0], Y[0]), (X[1], Y[1]), etc...

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : ndarray (n_samples, n_features)
        Array 1 for distance computation.

    Y : ndarray (n_samples, n_features)
        Array 2 for distance computation.

    metric : string or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        specified in PAIRED_DISTANCES, including "euclidean",
        "manhattan", or "cosine".
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    Returns
    -------
    distances : ndarray (n_samples, )

    Examples
    --------
    >>> from sklearn.metrics.pairwise import paired_distances
    >>> X = [[0, 1], [1, 1]]
    >>> Y = [[0, 1], [2, 1]]
    >>> paired_distances(X, Y)
    array([0., 1.])

    See also
    --------
    pairwise_distances : Computes the distance between every pair of samples
    """

    if metric in PAIRED_DISTANCES:
        func = PAIRED_DISTANCES[metric]
        return func(X, Y)
    elif callable(metric):
        # Check the matrix first (it is usually done by the metric)
        X, Y = check_paired_arrays(X, Y)
        distances = np.zeros(len(X))
        for i in range(len(X)):
            distances[i] = metric(X[i], Y[i])
        return distances
    else:
        raise ValueError('Unknown distance %s' % metric)


# Kernels
def linear_kernel(X, Y=None, dense_output=True):
    """
    Compute the linear kernel between X and Y.

    Read more in the :ref:`User Guide <linear_kernel>`.

    Parameters
    ----------
    X : array of shape (n_samples_1, n_features)

    Y : array of shape (n_samples_2, n_features)

    dense_output : boolean (optional), default True
        Whether to return dense output even when the input is sparse. If
        ``False``, the output is sparse if both input arrays are sparse.

        .. versionadded:: 0.20

    Returns
    -------
    Gram matrix : array of shape (n_samples_1, n_samples_2)
    """
    X, Y = check_pairwise_arrays(X, Y)
    return safe_sparse_dot(X, Y.T, dense_output=dense_output)


def polynomial_kernel(X, Y=None, degree=3, gamma=None, coef0=1):
    """
    Compute the polynomial kernel between X and Y::

        K(X, Y) = (gamma <X, Y> + coef0)^degree

    Read more in the :ref:`User Guide <polynomial_kernel>`.

    Parameters
    ----------
    X : ndarray of shape (n_samples_1, n_features)

    Y : ndarray of shape (n_samples_2, n_features)

    degree : int, default 3

    gamma : float, default None
        if None, defaults to 1.0 / n_features

    coef0 : float, default 1

    Returns
    -------
    Gram matrix : array of shape (n_samples_1, n_samples_2)
    """
    X, Y = check_pairwise_arrays(X, Y)
    if gamma is None:
        gamma = 1.0 / X.shape[1]

    K = safe_sparse_dot(X, Y.T, dense_output=True)
    K *= gamma
    K += coef0
    K **= degree
    return K


def sigmoid_kernel(X, Y=None, gamma=None, coef0=1):
    """
    Compute the sigmoid kernel between X and Y::

        K(X, Y) = tanh(gamma <X, Y> + coef0)

    Read more in the :ref:`User Guide <sigmoid_kernel>`.

    Parameters
    ----------
    X : ndarray of shape (n_samples_1, n_features)

    Y : ndarray of shape (n_samples_2, n_features)

    gamma : float, default None
        If None, defaults to 1.0 / n_features

    coef0 : float, default 1

    Returns
    -------
    Gram matrix : array of shape (n_samples_1, n_samples_2)
    """
    X, Y = check_pairwise_arrays(X, Y)
    if gamma is None:
        gamma = 1.0 / X.shape[1]

    K = safe_sparse_dot(X, Y.T, dense_output=True)
    K *= gamma
    K += coef0
    np.tanh(K, K)  # compute tanh in-place
    return K


def rbf_kernel(X, Y=None, gamma=None):
    """
    Compute the rbf (gaussian) kernel between X and Y::

        K(x, y) = exp(-gamma ||x-y||^2)

    for each pair of rows x in X and y in Y.

    Read more in the :ref:`User Guide <rbf_kernel>`.

    Parameters
    ----------
    X : array of shape (n_samples_X, n_features)

    Y : array of shape (n_samples_Y, n_features)

    gamma : float, default None
        If None, defaults to 1.0 / n_features

    Returns
    -------
    kernel_matrix : array of shape (n_samples_X, n_samples_Y)
    """
    X, Y = check_pairwise_arrays(X, Y)
    if gamma is None:
        gamma = 1.0 / X.shape[1]

    K = euclidean_distances(X, Y, squared=True)
    K *= -gamma
    np.exp(K, K)  # exponentiate K in-place
    return K


def laplacian_kernel(X, Y=None, gamma=None):
    """Compute the laplacian kernel between X and Y.

    The laplacian kernel is defined as::

        K(x, y) = exp(-gamma ||x-y||_1)

    for each pair of rows x in X and y in Y.
    Read more in the :ref:`User Guide <laplacian_kernel>`.

    .. versionadded:: 0.17

    Parameters
    ----------
    X : array of shape (n_samples_X, n_features)

    Y : array of shape (n_samples_Y, n_features)

    gamma : float, default None
        If None, defaults to 1.0 / n_features

    Returns
    -------
    kernel_matrix : array of shape (n_samples_X, n_samples_Y)
    """
    X, Y = check_pairwise_arrays(X, Y)
    if gamma is None:
        gamma = 1.0 / X.shape[1]

    K = -gamma * manhattan_distances(X, Y)
    np.exp(K, K)  # exponentiate K in-place
    return K


def cosine_similarity(X, Y=None, dense_output=True):
    """Compute cosine similarity between samples in X and Y.

    Cosine similarity, or the cosine kernel, computes similarity as the
    normalized dot product of X and Y:

        K(X, Y) = <X, Y> / (||X||*||Y||)

    On L2-normalized data, this function is equivalent to linear_kernel.

    Read more in the :ref:`User Guide <cosine_similarity>`.

    Parameters
    ----------
    X : ndarray or sparse array, shape: (n_samples_X, n_features)
        Input data.

    Y : ndarray or sparse array, shape: (n_samples_Y, n_features)
        Input data. If ``None``, the output will be the pairwise
        similarities between all samples in ``X``.

    dense_output : boolean (optional), default True
        Whether to return dense output even when the input is sparse. If
        ``False``, the output is sparse if both input arrays are sparse.

        .. versionadded:: 0.17
           parameter ``dense_output`` for dense output.

    Returns
    -------
    kernel matrix : array
        An array with shape (n_samples_X, n_samples_Y).
    """
    # to avoid recursive import

    X, Y = check_pairwise_arrays(X, Y)

    X_normalized = normalize(X, copy=True)
    if X is Y:
        Y_normalized = X_normalized
    else:
        Y_normalized = normalize(Y, copy=True)

    K = safe_sparse_dot(X_normalized, Y_normalized.T,
                        dense_output=dense_output)

    return K


def additive_chi2_kernel(X, Y=None):
    """Computes the additive chi-squared kernel between observations in X and Y

    The chi-squared kernel is computed between each pair of rows in X and Y.  X
    and Y have to be non-negative. This kernel is most commonly applied to
    histograms.

    The chi-squared kernel is given by::

        k(x, y) = -Sum [(x - y)^2 / (x + y)]

    It can be interpreted as a weighted difference per entry.

    Read more in the :ref:`User Guide <chi2_kernel>`.

    Notes
    -----
    As the negative of a distance, this kernel is only conditionally positive
    definite.


    Parameters
    ----------
    X : array-like of shape (n_samples_X, n_features)

    Y : array of shape (n_samples_Y, n_features)

    Returns
    -------
    kernel_matrix : array of shape (n_samples_X, n_samples_Y)

    References
    ----------
    * Zhang, J. and Marszalek, M. and Lazebnik, S. and Schmid, C.
      Local features and kernels for classification of texture and object
      categories: A comprehensive study
      International Journal of Computer Vision 2007
      https://research.microsoft.com/en-us/um/people/manik/projects/trade-off/papers/ZhangIJCV06.pdf


    See also
    --------
    chi2_kernel : The exponentiated version of the kernel, which is usually
        preferable.

    sklearn.kernel_approximation.AdditiveChi2Sampler : A Fourier approximation
        to this kernel.
    """
    if issparse(X) or issparse(Y):
        raise ValueError("additive_chi2 does not support sparse matrices.")
    X, Y = check_pairwise_arrays(X, Y)
    if (X < 0).any():
        raise ValueError("X contains negative values.")
    if Y is not X and (Y < 0).any():
        raise ValueError("Y contains negative values.")

    result = np.zeros((X.shape[0], Y.shape[0]), dtype=X.dtype)
    _chi2_kernel_fast(X, Y, result)
    return result


def chi2_kernel(X, Y=None, gamma=1.):
    """Computes the exponential chi-squared kernel X and Y.

    The chi-squared kernel is computed between each pair of rows in X and Y.  X
    and Y have to be non-negative. This kernel is most commonly applied to
    histograms.

    The chi-squared kernel is given by::

        k(x, y) = exp(-gamma Sum [(x - y)^2 / (x + y)])

    It can be interpreted as a weighted difference per entry.

    Read more in the :ref:`User Guide <chi2_kernel>`.

    Parameters
    ----------
    X : array-like of shape (n_samples_X, n_features)

    Y : array of shape (n_samples_Y, n_features)

    gamma : float, default=1.
        Scaling parameter of the chi2 kernel.

    Returns
    -------
    kernel_matrix : array of shape (n_samples_X, n_samples_Y)

    References
    ----------
    * Zhang, J. and Marszalek, M. and Lazebnik, S. and Schmid, C.
      Local features and kernels for classification of texture and object
      categories: A comprehensive study
      International Journal of Computer Vision 2007
      https://research.microsoft.com/en-us/um/people/manik/projects/trade-off/papers/ZhangIJCV06.pdf

    See also
    --------
    additive_chi2_kernel : The additive version of this kernel

    sklearn.kernel_approximation.AdditiveChi2Sampler : A Fourier approximation
        to the additive version of this kernel.
    """
    K = additive_chi2_kernel(X, Y)
    K *= gamma
    return np.exp(K, K)


# Helper functions - distance
PAIRWISE_DISTANCE_FUNCTIONS = {
    # If updating this dictionary, update the doc in both distance_metrics()
    # and also in pairwise_distances()!
    'cityblock': manhattan_distances,
    'cosine': cosine_distances,
    'euclidean': euclidean_distances,
    'haversine': haversine_distances,
    'l2': euclidean_distances,
    'l1': manhattan_distances,
    'manhattan': manhattan_distances,
    'precomputed': None,  # HACK: precomputed is always allowed, never called
}


def distance_metrics():
    """Valid metrics for pairwise_distances.

    This function simply returns the valid pairwise distance metrics.
    It exists to allow for a description of the mapping for
    each of the valid strings.

    The valid distance metrics, and the function they map to, are:

    ============   ====================================
    metric         Function
    ============   ====================================
    'cityblock'    metrics.pairwise.manhattan_distances
    'cosine'       metrics.pairwise.cosine_distances
    'euclidean'    metrics.pairwise.euclidean_distances
    'haversine'    metrics.pairwise.haversine_distances
    'l1'           metrics.pairwise.manhattan_distances
    'l2'           metrics.pairwise.euclidean_distances
    'manhattan'    metrics.pairwise.manhattan_distances
    ============   ====================================

    Read more in the :ref:`User Guide <metrics>`.

    """
    return PAIRWISE_DISTANCE_FUNCTIONS


def _dist_wrapper(dist_func, dist_matrix, slice_, *args, **kwargs):
    """Write in-place to a slice of a distance matrix"""
    dist_matrix[:, slice_] = dist_func(*args, **kwargs)


def _parallel_pairwise(X, Y, func, n_jobs, **kwds):
    """Break the pairwise matrix in n_jobs even slices
    and compute them in parallel"""

    if Y is None:
        Y = X
    X, Y, dtype = _return_float_dtype(X, Y)

    if effective_n_jobs(n_jobs) == 1:
        return func(X, Y, **kwds)

    # enforce a threading backend to prevent data communication overhead
    fd = delayed(_dist_wrapper)
    ret = np.empty((X.shape[0], Y.shape[0]), dtype=dtype, order='F')
    Parallel(backend="threading", n_jobs=n_jobs)(
        fd(func, ret, s, X, Y[s], **kwds)
        for s in gen_even_slices(_num_samples(Y), effective_n_jobs(n_jobs)))

    if (X is Y or Y is None) and func is euclidean_distances:
        # zeroing diagonal for euclidean norm.
        # TODO: do it also for other norms.
        np.fill_diagonal(ret, 0)

    return ret


def _pairwise_callable(X, Y, metric, **kwds):
    """Handle the callable case for pairwise_{distances,kernels}
    """
    X, Y = check_pairwise_arrays(X, Y)

    if X is Y:
        # Only calculate metric for upper triangle
        out = np.zeros((X.shape[0], Y.shape[0]), dtype='float')
        iterator = itertools.combinations(range(X.shape[0]), 2)
        for i, j in iterator:
            out[i, j] = metric(X[i], Y[j], **kwds)

        # Make symmetric
        # NB: out += out.T will produce incorrect results
        out = out + out.T

        # Calculate diagonal
        # NB: nonzero diagonals are allowed for both metrics and kernels
        for i in range(X.shape[0]):
            x = X[i]
            out[i, i] = metric(x, x, **kwds)

    else:
        # Calculate all cells
        out = np.empty((X.shape[0], Y.shape[0]), dtype='float')
        iterator = itertools.product(range(X.shape[0]), range(Y.shape[0]))
        for i, j in iterator:
            out[i, j] = metric(X[i], Y[j], **kwds)

    return out


_VALID_METRICS = ['euclidean', 'l2', 'l1', 'manhattan', 'cityblock',
                  'braycurtis', 'canberra', 'chebyshev', 'correlation',
                  'cosine', 'dice', 'hamming', 'jaccard', 'kulsinski',
                  'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
                  'russellrao', 'seuclidean', 'sokalmichener',
                  'sokalsneath', 'sqeuclidean', 'yule', 'wminkowski',
                  'haversine']


def _check_chunk_size(reduced, chunk_size):
    """Checks chunk is a sequence of expected size or a tuple of same
    """
    is_tuple = isinstance(reduced, tuple)
    if not is_tuple:
        reduced = (reduced,)
    if any(isinstance(r, tuple) or not hasattr(r, '__iter__')
           for r in reduced):
        raise TypeError('reduce_func returned %r. '
                        'Expected sequence(s) of length %d.' %
                        (reduced if is_tuple else reduced[0], chunk_size))
    if any(_num_samples(r) != chunk_size for r in reduced):
        actual_size = tuple(_num_samples(r) for r in reduced)
        raise ValueError('reduce_func returned object of length %s. '
                         'Expected same length as input: %d.' %
                         (actual_size if is_tuple else actual_size[0],
                          chunk_size))


def _precompute_metric_params(X, Y, metric=None, **kwds):
    """Precompute data-derived metric parameters if not provided
    """
    if metric == "seuclidean" and 'V' not in kwds:
        if X is Y:
            V = np.var(X, axis=0, ddof=1)
        else:
            V = np.var(np.vstack([X, Y]), axis=0, ddof=1)
        return {'V': V}
    if metric == "mahalanobis" and 'VI' not in kwds:
        if X is Y:
            VI = np.linalg.inv(np.cov(X.T)).T
        else:
            VI = np.linalg.inv(np.cov(np.vstack([X, Y]).T)).T
        return {'VI': VI}
    return {}


def pairwise_distances_chunked(X, Y=None, reduce_func=None,
                               metric='euclidean', n_jobs=None,
                               working_memory=None, **kwds):
    """Generate a distance matrix chunk by chunk with optional reduction

    In cases where not all of a pairwise distance matrix needs to be stored at
    once, this is used to calculate pairwise distances in
    ``working_memory``-sized chunks.  If ``reduce_func`` is given, it is run
    on each chunk and its return values are concatenated into lists, arrays
    or sparse matrices.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or,
        [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    Y : array [n_samples_b, n_features], optional
        An optional second feature array. Only allowed if
        metric != "precomputed".

    reduce_func : callable, optional
        The function which is applied on each chunk of the distance matrix,
        reducing it to needed values.  ``reduce_func(D_chunk, start)``
        is called repeatedly, where ``D_chunk`` is a contiguous vertical
        slice of the pairwise distance matrix, starting at row ``start``.
        It should return an array, a list, or a sparse matrix of length
        ``D_chunk.shape[0]``, or a tuple of such objects.

        If None, pairwise_distances_chunked returns a generator of vertical
        chunks of the distance matrix.

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter, or
        a metric listed in pairwise.PAIRWISE_DISTANCE_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a distance matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    working_memory : int, optional
        The sought maximum memory for temporary distance matrix chunks.
        When None (default), the value of
        ``sklearn.get_config()['working_memory']`` is used.

    `**kwds` : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a scipy.spatial.distance metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Yields
    ------
    D_chunk : array or sparse matrix
        A contiguous slice of distance matrix, optionally processed by
        ``reduce_func``.

    Examples
    --------
    Without reduce_func:

    >>> import numpy as np
    >>> from sklearn.metrics import pairwise_distances_chunked
    >>> X = np.random.RandomState(0).rand(5, 3)
    >>> D_chunk = next(pairwise_distances_chunked(X))
    >>> D_chunk  # doctest: +ELLIPSIS
    array([[0.  ..., 0.29..., 0.41..., 0.19..., 0.57...],
           [0.29..., 0.  ..., 0.57..., 0.41..., 0.76...],
           [0.41..., 0.57..., 0.  ..., 0.44..., 0.90...],
           [0.19..., 0.41..., 0.44..., 0.  ..., 0.51...],
           [0.57..., 0.76..., 0.90..., 0.51..., 0.  ...]])

    Retrieve all neighbors and average distance within radius r:

    >>> r = .2
    >>> def reduce_func(D_chunk, start):
    ...     neigh = [np.flatnonzero(d < r) for d in D_chunk]
    ...     avg_dist = (D_chunk * (D_chunk < r)).mean(axis=1)
    ...     return neigh, avg_dist
    >>> gen = pairwise_distances_chunked(X, reduce_func=reduce_func)
    >>> neigh, avg_dist = next(gen)
    >>> neigh
    [array([0, 3]), array([1]), array([2]), array([0, 3]), array([4])]
    >>> avg_dist  # doctest: +ELLIPSIS
    array([0.039..., 0.        , 0.        , 0.039..., 0.        ])

    Where r is defined per sample, we need to make use of ``start``:

    >>> r = [.2, .4, .4, .3, .1]
    >>> def reduce_func(D_chunk, start):
    ...     neigh = [np.flatnonzero(d < r[i])
    ...              for i, d in enumerate(D_chunk, start)]
    ...     return neigh
    >>> neigh = next(pairwise_distances_chunked(X, reduce_func=reduce_func))
    >>> neigh
    [array([0, 3]), array([0, 1]), array([2]), array([0, 3]), array([4])]

    Force row-by-row generation by reducing ``working_memory``:

    >>> gen = pairwise_distances_chunked(X, reduce_func=reduce_func,
    ...                                  working_memory=0)
    >>> next(gen)
    [array([0, 3])]
    >>> next(gen)
    [array([0, 1])]
    """
    n_samples_X = _num_samples(X)
    if metric == 'precomputed':
        slices = (slice(0, n_samples_X),)
    else:
        if Y is None:
            Y = X
        # We get as many rows as possible within our working_memory budget to
        # store len(Y) distances in each row of output.
        #
        # Note:
        #  - this will get at least 1 row, even if 1 row of distances will
        #    exceed working_memory.
        #  - this does not account for any temporary memory usage while
        #    calculating distances (e.g. difference of vectors in manhattan
        #    distance.
        chunk_n_rows = get_chunk_n_rows(row_bytes=8 * _num_samples(Y),
                                        max_n_rows=n_samples_X,
                                        working_memory=working_memory)
        slices = gen_batches(n_samples_X, chunk_n_rows)

    # precompute data-derived metric params
    params = _precompute_metric_params(X, Y, metric=metric, **kwds)
    kwds.update(**params)

    for sl in slices:
        if sl.start == 0 and sl.stop == n_samples_X:
            X_chunk = X  # enable optimised paths for X is Y
        else:
            X_chunk = X[sl]
        D_chunk = pairwise_distances(X_chunk, Y, metric=metric,
                                     n_jobs=n_jobs, **kwds)
        if ((X is Y or Y is None)
                and PAIRWISE_DISTANCE_FUNCTIONS.get(metric, None)
                is euclidean_distances):
            # zeroing diagonal, taking care of aliases of "euclidean",
            # i.e. "l2"
            D_chunk.flat[sl.start::_num_samples(X) + 1] = 0
        if reduce_func is not None:
            chunk_size = D_chunk.shape[0]
            D_chunk = reduce_func(D_chunk, sl.start)
            _check_chunk_size(D_chunk, chunk_size)
        yield D_chunk


def pairwise_distances(X, Y=None, metric="euclidean", n_jobs=None, **kwds):
    """ Compute the distance matrix from a vector array X and optional Y.

    This method takes either a vector array or a distance matrix, and returns
    a distance matrix. If the input is a vector array, the distances are
    computed. If the input is a distances matrix, it is returned instead.

    This method provides a safe way to take a distance matrix as input, while
    preserving compatibility with many other algorithms that take a vector
    array.

    If Y is given (default is None), then the returned matrix is the pairwise
    distance between the arrays from both X and Y.

    Valid values for metric are:

    - From scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
      'manhattan']. These metrics support sparse matrix inputs.

    - From scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
      'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis',
      'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
      'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
      See the documentation for scipy.spatial.distance for details on these
      metrics. These metrics do not support sparse matrix inputs.

    Note that in the case of 'cityblock', 'cosine' and 'euclidean' (which are
    valid scipy.spatial.distance metrics), the scikit-learn implementation
    will be used, which is faster and has support for sparse matrices (except
    for 'cityblock'). For a verbose description of the metrics from
    scikit-learn, see the __doc__ of the sklearn.pairwise.distance_metrics
    function.

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
             [n_samples_a, n_features] otherwise
        Array of pairwise distances between samples, or a feature array.

    Y : array [n_samples_b, n_features], optional
        An optional second feature array. Only allowed if
        metric != "precomputed".

    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter, or
        a metric listed in pairwise.PAIRWISE_DISTANCE_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a distance matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    **kwds : optional keyword parameters
        Any further parameters are passed directly to the distance function.
        If using a scipy.spatial.distance metric, the parameters are still
        metric dependent. See the scipy docs for usage examples.

    Returns
    -------
    D : array [n_samples_a, n_samples_a] or [n_samples_a, n_samples_b]
        A distance matrix D such that D_{i, j} is the distance between the
        ith and jth vectors of the given matrix X, if Y is None.
        If Y is not None, then D_{i, j} is the distance between the ith array
        from X and the jth array from Y.

    See also
    --------
    pairwise_distances_chunked : performs the same calculation as this
        function, but returns a generator of chunks of the distance matrix, in
        order to limit memory usage.
    paired_distances : Computes the distances between corresponding
                       elements of two arrays
    """
    if (metric not in _VALID_METRICS and
            not callable(metric) and metric != "precomputed"):
        raise ValueError("Unknown metric %s. "
                         "Valid metrics are %s, or 'precomputed', or a "
                         "callable" % (metric, _VALID_METRICS))

    if metric == "precomputed":
        X, _ = check_pairwise_arrays(X, Y, precomputed=True)

        whom = ("`pairwise_distances`. Precomputed distance "
                " need to have non-negative values.")
        check_non_negative(X, whom=whom)
        return X
    elif metric in PAIRWISE_DISTANCE_FUNCTIONS:
        func = PAIRWISE_DISTANCE_FUNCTIONS[metric]
    elif callable(metric):
        func = partial(_pairwise_callable, metric=metric, **kwds)
    else:
        if issparse(X) or issparse(Y):
            raise TypeError("scipy distance metrics do not"
                            " support sparse matrices.")

        dtype = bool if metric in PAIRWISE_BOOLEAN_FUNCTIONS else None

        if (dtype == bool and
                (X.dtype != bool or (Y is not None and Y.dtype != bool))):
            msg = "Data was converted to boolean for metric %s" % metric
            warnings.warn(msg, DataConversionWarning)

        X, Y = check_pairwise_arrays(X, Y, dtype=dtype)

        # precompute data-derived metric params
        params = _precompute_metric_params(X, Y, metric=metric, **kwds)
        kwds.update(**params)

        if effective_n_jobs(n_jobs) == 1 and X is Y:
            return distance.squareform(distance.pdist(X, metric=metric,
                                                      **kwds))
        func = partial(distance.cdist, metric=metric, **kwds)

    return _parallel_pairwise(X, Y, func, n_jobs, **kwds)


# These distances recquire boolean arrays, when using scipy.spatial.distance
PAIRWISE_BOOLEAN_FUNCTIONS = [
    'dice',
    'jaccard',
    'kulsinski',
    'matching',
    'rogerstanimoto',
    'russellrao',
    'sokalmichener',
    'sokalsneath',
    'yule',
]

# Helper functions - distance
PAIRWISE_KERNEL_FUNCTIONS = {
    # If updating this dictionary, update the doc in both distance_metrics()
    # and also in pairwise_distances()!
    'additive_chi2': additive_chi2_kernel,
    'chi2': chi2_kernel,
    'linear': linear_kernel,
    'polynomial': polynomial_kernel,
    'poly': polynomial_kernel,
    'rbf': rbf_kernel,
    'laplacian': laplacian_kernel,
    'sigmoid': sigmoid_kernel,
    'cosine': cosine_similarity, }


def kernel_metrics():
    """ Valid metrics for pairwise_kernels

    This function simply returns the valid pairwise distance metrics.
    It exists, however, to allow for a verbose description of the mapping for
    each of the valid strings.

    The valid distance metrics, and the function they map to, are:
      ===============   ========================================
      metric            Function
      ===============   ========================================
      'additive_chi2'   sklearn.pairwise.additive_chi2_kernel
      'chi2'            sklearn.pairwise.chi2_kernel
      'linear'          sklearn.pairwise.linear_kernel
      'poly'            sklearn.pairwise.polynomial_kernel
      'polynomial'      sklearn.pairwise.polynomial_kernel
      'rbf'             sklearn.pairwise.rbf_kernel
      'laplacian'       sklearn.pairwise.laplacian_kernel
      'sigmoid'         sklearn.pairwise.sigmoid_kernel
      'cosine'          sklearn.pairwise.cosine_similarity
      ===============   ========================================

    Read more in the :ref:`User Guide <metrics>`.
    """
    return PAIRWISE_KERNEL_FUNCTIONS


KERNEL_PARAMS = {
    "additive_chi2": (),
    "chi2": frozenset(["gamma"]),
    "cosine": (),
    "linear": (),
    "poly": frozenset(["gamma", "degree", "coef0"]),
    "polynomial": frozenset(["gamma", "degree", "coef0"]),
    "rbf": frozenset(["gamma"]),
    "laplacian": frozenset(["gamma"]),
    "sigmoid": frozenset(["gamma", "coef0"]),
}


def pairwise_kernels(X, Y=None, metric="linear", filter_params=False,
                     n_jobs=None, **kwds):
    """Compute the kernel between arrays X and optional array Y.

    This method takes either a vector array or a kernel matrix, and returns
    a kernel matrix. If the input is a vector array, the kernels are
    computed. If the input is a kernel matrix, it is returned instead.

    This method provides a safe way to take a kernel matrix as input, while
    preserving compatibility with many other algorithms that take a vector
    array.

    If Y is given (default is None), then the returned matrix is the pairwise
    kernel between the arrays from both X and Y.

    Valid values for metric are::
        ['additive_chi2', 'chi2', 'linear', 'poly', 'polynomial', 'rbf',
         'laplacian', 'sigmoid', 'cosine']

    Read more in the :ref:`User Guide <metrics>`.

    Parameters
    ----------
    X : array [n_samples_a, n_samples_a] if metric == "precomputed", or, \
             [n_samples_a, n_features] otherwise
        Array of pairwise kernels between samples, or a feature array.

    Y : array [n_samples_b, n_features]
        A second feature array only if X has shape [n_samples_a, n_features].

    metric : string, or callable
        The metric to use when calculating kernel between instances in a
        feature array. If metric is a string, it must be one of the metrics
        in pairwise.PAIRWISE_KERNEL_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a kernel matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them.

    filter_params : boolean
        Whether to filter invalid parameters or not.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation. This works by breaking
        down the pairwise matrix into n_jobs even slices and computing them in
        parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    **kwds : optional keyword parameters
        Any further parameters are passed directly to the kernel function.

    Returns
    -------
    K : array [n_samples_a, n_samples_a] or [n_samples_a, n_samples_b]
        A kernel matrix K such that K_{i, j} is the kernel between the
        ith and jth vectors of the given matrix X, if Y is None.
        If Y is not None, then K_{i, j} is the kernel between the ith array
        from X and the jth array from Y.

    Notes
    -----
    If metric is 'precomputed', Y is ignored and X is returned.

    """
    # import GPKernel locally to prevent circular imports
    from ..gaussian_process.kernels import Kernel as GPKernel

    if metric == "precomputed":
        X, _ = check_pairwise_arrays(X, Y, precomputed=True)
        return X
    elif isinstance(metric, GPKernel):
        func = metric.__call__
    elif metric in PAIRWISE_KERNEL_FUNCTIONS:
        if filter_params:
            kwds = {k: kwds[k] for k in kwds
                    if k in KERNEL_PARAMS[metric]}
        func = PAIRWISE_KERNEL_FUNCTIONS[metric]
    elif callable(metric):
        func = partial(_pairwise_callable, metric=metric, **kwds)
    else:
        raise ValueError("Unknown kernel %r" % metric)

    return _parallel_pairwise(X, Y, func, n_jobs, **kwds)
