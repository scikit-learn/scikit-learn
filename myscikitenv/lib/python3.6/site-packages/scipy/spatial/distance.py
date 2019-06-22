"""
Distance computations (:mod:`scipy.spatial.distance`)
=====================================================

.. sectionauthor:: Damian Eads

Function Reference
------------------

Distance matrix computation from a collection of raw observation vectors
stored in a rectangular array.

.. autosummary::
   :toctree: generated/

   pdist   -- pairwise distances between observation vectors.
   cdist   -- distances between two collections of observation vectors
   squareform -- convert distance matrix to a condensed one and vice versa
   directed_hausdorff -- directed Hausdorff distance between arrays

Predicates for checking the validity of distance matrices, both
condensed and redundant. Also contained in this module are functions
for computing the number of observations in a distance matrix.

.. autosummary::
   :toctree: generated/

   is_valid_dm -- checks for a valid distance matrix
   is_valid_y  -- checks for a valid condensed distance matrix
   num_obs_dm  -- # of observations in a distance matrix
   num_obs_y   -- # of observations in a condensed distance matrix

Distance functions between two numeric vectors ``u`` and ``v``. Computing
distances over a large collection of vectors is inefficient for these
functions. Use ``pdist`` for this purpose.

.. autosummary::
   :toctree: generated/

   braycurtis       -- the Bray-Curtis distance.
   canberra         -- the Canberra distance.
   chebyshev        -- the Chebyshev distance.
   cityblock        -- the Manhattan distance.
   correlation      -- the Correlation distance.
   cosine           -- the Cosine distance.
   euclidean        -- the Euclidean distance.
   jensenshannon    -- the Jensen-Shannon distance.
   mahalanobis      -- the Mahalanobis distance.
   minkowski        -- the Minkowski distance.
   seuclidean       -- the normalized Euclidean distance.
   sqeuclidean      -- the squared Euclidean distance.
   wminkowski       -- (deprecated) alias of `minkowski`.

Distance functions between two boolean vectors (representing sets) ``u`` and
``v``.  As in the case of numerical vectors, ``pdist`` is more efficient for
computing the distances between all pairs.

.. autosummary::
   :toctree: generated/

   dice             -- the Dice dissimilarity.
   hamming          -- the Hamming distance.
   jaccard          -- the Jaccard distance.
   kulsinski        -- the Kulsinski distance.
   rogerstanimoto   -- the Rogers-Tanimoto dissimilarity.
   russellrao       -- the Russell-Rao dissimilarity.
   sokalmichener    -- the Sokal-Michener dissimilarity.
   sokalsneath      -- the Sokal-Sneath dissimilarity.
   yule             -- the Yule dissimilarity.

:func:`hamming` also operates over discrete numerical vectors.
"""

# Copyright (C) Damian Eads, 2007-2008. New BSD License.

from __future__ import division, print_function, absolute_import

__all__ = [
    'braycurtis',
    'canberra',
    'cdist',
    'chebyshev',
    'cityblock',
    'correlation',
    'cosine',
    'dice',
    'directed_hausdorff',
    'euclidean',
    'hamming',
    'is_valid_dm',
    'is_valid_y',
    'jaccard',
    'jensenshannon',
    'kulsinski',
    'mahalanobis',
    'matching',
    'minkowski',
    'num_obs_dm',
    'num_obs_y',
    'pdist',
    'rogerstanimoto',
    'russellrao',
    'seuclidean',
    'sokalmichener',
    'sokalsneath',
    'sqeuclidean',
    'squareform',
    'wminkowski',
    'yule'
]


import warnings
import numpy as np

from functools import partial
from collections import namedtuple
from scipy._lib.six import callable, string_types
from scipy._lib.six import xrange
from scipy._lib._util import _asarray_validated

from . import _distance_wrap
from . import _hausdorff
from ..linalg import norm
from ..special import rel_entr


def _args_to_kwargs_xdist(args, kwargs, metric, func_name):
    """
    Convert legacy positional arguments to keyword arguments for pdist/cdist.
    """
    if not args:
        return kwargs

    if (callable(metric) and metric not in [
            braycurtis, canberra, chebyshev, cityblock, correlation, cosine,
            dice, euclidean, hamming, jaccard, jensenshannon, kulsinski,
            mahalanobis, matching, minkowski, rogerstanimoto, russellrao,
            seuclidean, sokalmichener, sokalsneath, sqeuclidean, yule,
            wminkowski]):
        raise TypeError('When using a custom metric arguments must be passed'
                        'as keyword (i.e., ARGNAME=ARGVALUE)')

    if func_name == 'pdist':
        old_arg_names = ['p', 'w', 'V', 'VI']
    else:
        old_arg_names = ['p', 'V', 'VI', 'w']

    num_args = len(args)
    warnings.warn('%d metric parameters have been passed as positional.'
                  'This will raise an error in a future version.'
                  'Please pass arguments as keywords(i.e., ARGNAME=ARGVALUE)'
                  % num_args, DeprecationWarning)

    if num_args > 4:
        raise ValueError('Deprecated %s signature accepts only 4'
                         'positional arguments (%s), %d given.'
                         % (func_name, ', '.join(old_arg_names), num_args))

    for old_arg, arg in zip(old_arg_names, args):
        if old_arg in kwargs:
            raise TypeError('%s() got multiple values for argument %s'
                            % (func_name, old_arg))
        kwargs[old_arg] = arg
    return kwargs


def _copy_array_if_base_present(a):
    """Copy the array if its base points to a parent array."""
    if a.base is not None:
        return a.copy()
    return a


def _correlation_cdist_wrap(XA, XB, dm, **kwargs):
    XA = XA - XA.mean(axis=1, keepdims=True)
    XB = XB - XB.mean(axis=1, keepdims=True)
    _distance_wrap.cdist_cosine_double_wrap(XA, XB, dm, **kwargs)


def _correlation_pdist_wrap(X, dm, **kwargs):
    X2 = X - X.mean(axis=1, keepdims=True)
    _distance_wrap.pdist_cosine_double_wrap(X2, dm, **kwargs)


def _convert_to_type(X, out_type):
    return np.ascontiguousarray(X, dtype=out_type)


def _filter_deprecated_kwargs(kwargs, args_blacklist):
    # Filtering out old default keywords
    for k in args_blacklist:
        if k in kwargs:
            del kwargs[k]
            warnings.warn('Got unexpected kwarg %s. This will raise an error'
                          ' in a future version.' % k, DeprecationWarning)


def _nbool_correspond_all(u, v, w=None):
    if u.dtype == v.dtype == bool and w is None:
        not_u = ~u
        not_v = ~v
        nff = (not_u & not_v).sum()
        nft = (not_u & v).sum()
        ntf = (u & not_v).sum()
        ntt = (u & v).sum()
    else:
        dtype = np.find_common_type([int], [u.dtype, v.dtype])
        u = u.astype(dtype)
        v = v.astype(dtype)
        not_u = 1.0 - u
        not_v = 1.0 - v
        if w is not None:
            not_u = w * not_u
            u = w * u
        nff = (not_u * not_v).sum()
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
        ntt = (u * v).sum()
    return (nff, nft, ntf, ntt)


def _nbool_correspond_ft_tf(u, v, w=None):
    if u.dtype == v.dtype == bool and w is None:
        not_u = ~u
        not_v = ~v
        nft = (not_u & v).sum()
        ntf = (u & not_v).sum()
    else:
        dtype = np.find_common_type([int], [u.dtype, v.dtype])
        u = u.astype(dtype)
        v = v.astype(dtype)
        not_u = 1.0 - u
        not_v = 1.0 - v
        if w is not None:
            not_u = w * not_u
            u = w * u
        nft = (not_u * v).sum()
        ntf = (u * not_v).sum()
    return (nft, ntf)


def _validate_cdist_input(XA, XB, mA, mB, n, metric_name, **kwargs):
    if metric_name is not None:
        # get supported types
        types = _METRICS[metric_name].types
        # choose best type
        typ = types[types.index(XA.dtype)] if XA.dtype in types else types[0]
        # validate data
        XA = _convert_to_type(XA, out_type=typ)
        XB = _convert_to_type(XB, out_type=typ)

        # validate kwargs
        _validate_kwargs = _METRICS[metric_name].validator
        if _validate_kwargs:
            kwargs = _validate_kwargs(np.vstack([XA, XB]), mA + mB, n, **kwargs)
    else:
        typ = None
    return XA, XB, typ, kwargs


def _validate_hamming_kwargs(X, m, n, **kwargs):
    w = kwargs.get('w', np.ones((n,), dtype='double'))

    if w.ndim != 1 or w.shape[0] != n:
        raise ValueError("Weights must have same size as input vector. %d vs. %d" % (w.shape[0], n))

    kwargs['w'] = _validate_weights(w)
    return kwargs


def _validate_mahalanobis_kwargs(X, m, n, **kwargs):
    VI = kwargs.pop('VI', None)
    if VI is None:
        if m <= n:
            # There are fewer observations than the dimension of
            # the observations.
            raise ValueError("The number of observations (%d) is too "
                             "small; the covariance matrix is "
                             "singular. For observations with %d "
                             "dimensions, at least %d observations "
                             "are required." % (m, n, n + 1))
        CV = np.atleast_2d(np.cov(X.astype(np.double).T))
        VI = np.linalg.inv(CV).T.copy()
    kwargs["VI"] = _convert_to_double(VI)
    return kwargs


def _validate_minkowski_kwargs(X, m, n, **kwargs):
    if 'p' not in kwargs:
        kwargs['p'] = 2.
    return kwargs


def _validate_pdist_input(X, m, n, metric_name, **kwargs):
    if metric_name is not None:
        # get supported types
        types = _METRICS[metric_name].types
        # choose best type
        typ = types[types.index(X.dtype)] if X.dtype in types else types[0]
        # validate data
        X = _convert_to_type(X, out_type=typ)

        # validate kwargs
        _validate_kwargs = _METRICS[metric_name].validator
        if _validate_kwargs:
            kwargs = _validate_kwargs(X, m, n, **kwargs)
    else:
        typ = None
    return X, typ, kwargs


def _validate_seuclidean_kwargs(X, m, n, **kwargs):
    V = kwargs.pop('V', None)
    if V is None:
        V = np.var(X.astype(np.double), axis=0, ddof=1)
    else:
        V = np.asarray(V, order='c')
        if V.dtype != np.double:
            raise TypeError('Variance vector V must contain doubles.')
        if len(V.shape) != 1:
            raise ValueError('Variance vector V must '
                             'be one-dimensional.')
        if V.shape[0] != n:
            raise ValueError('Variance vector V must be of the same '
                             'dimension as the vectors on which the distances '
                             'are computed.')
    kwargs['V'] = _convert_to_double(V)
    return kwargs


def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u


def _validate_weights(w, dtype=np.double):
    w = _validate_vector(w, dtype=dtype)
    if np.any(w < 0):
        raise ValueError("Input weights should be all non-negative")
    return w


def _validate_wminkowski_kwargs(X, m, n, **kwargs):
    w = kwargs.pop('w', None)
    if w is None:
        raise ValueError('weighted minkowski requires a weight '
                         'vector `w` to be given.')
    kwargs['w'] = _validate_weights(w)
    if 'p' not in kwargs:
        kwargs['p'] = 2.
    return kwargs


def directed_hausdorff(u, v, seed=0):
    """
    Compute the directed Hausdorff distance between two N-D arrays.

    Distances between pairs are calculated using a Euclidean metric.

    Parameters
    ----------
    u : (M,N) ndarray
        Input array.
    v : (O,N) ndarray
        Input array.
    seed : int or None
        Local `numpy.random.RandomState` seed. Default is 0, a random
        shuffling of u and v that guarantees reproducibility.

    Returns
    -------
    d : double
        The directed Hausdorff distance between arrays `u` and `v`,

    index_1 : int
        index of point contributing to Hausdorff pair in `u`

    index_2 : int
        index of point contributing to Hausdorff pair in `v`

    Raises
    ------
    ValueError
        An exception is thrown if `u` and `v` do not have
        the same number of columns.

    Notes
    -----
    Uses the early break technique and the random sampling approach
    described by [1]_. Although worst-case performance is ``O(m * o)``
    (as with the brute force algorithm), this is unlikely in practice
    as the input data would have to require the algorithm to explore
    every single point interaction, and after the algorithm shuffles
    the input points at that. The best case performance is O(m), which
    is satisfied by selecting an inner loop distance that is less than
    cmax and leads to an early break as often as possible. The authors
    have formally shown that the average runtime is closer to O(m).

    .. versionadded:: 0.19.0

    References
    ----------
    .. [1] A. A. Taha and A. Hanbury, "An efficient algorithm for
           calculating the exact Hausdorff distance." IEEE Transactions On
           Pattern Analysis And Machine Intelligence, vol. 37 pp. 2153-63,
           2015.

    See Also
    --------
    scipy.spatial.procrustes : Another similarity test for two data sets

    Examples
    --------
    Find the directed Hausdorff distance between two 2-D arrays of
    coordinates:

    >>> from scipy.spatial.distance import directed_hausdorff
    >>> u = np.array([(1.0, 0.0),
    ...               (0.0, 1.0),
    ...               (-1.0, 0.0),
    ...               (0.0, -1.0)])
    >>> v = np.array([(2.0, 0.0),
    ...               (0.0, 2.0),
    ...               (-2.0, 0.0),
    ...               (0.0, -4.0)])

    >>> directed_hausdorff(u, v)[0]
    2.23606797749979
    >>> directed_hausdorff(v, u)[0]
    3.0

    Find the general (symmetric) Hausdorff distance between two 2-D
    arrays of coordinates:

    >>> max(directed_hausdorff(u, v)[0], directed_hausdorff(v, u)[0])
    3.0

    Find the indices of the points that generate the Hausdorff distance
    (the Hausdorff pair):

    >>> directed_hausdorff(v, u)[1:]
    (3, 3)

    """
    u = np.asarray(u, dtype=np.float64, order='c')
    v = np.asarray(v, dtype=np.float64, order='c')
    if u.shape[1] != v.shape[1]:
        raise ValueError('u and v need to have the same '
                         'number of columns')
    result = _hausdorff.directed_hausdorff(u, v, seed)
    return result


def minkowski(u, v, p=2, w=None):
    """
    Compute the Minkowski distance between two 1-D arrays.

    The Minkowski distance between 1-D arrays `u` and `v`,
    is defined as

    .. math::

       {||u-v||}_p = (\\sum{|u_i - v_i|^p})^{1/p}.


       \\left(\\sum{w_i(|(u_i - v_i)|^p)}\\right)^{1/p}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    p : int
        The order of the norm of the difference :math:`{||u-v||}_p`.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    minkowski : double
        The Minkowski distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.minkowski([1, 0, 0], [0, 1, 0], 1)
    2.0
    >>> distance.minkowski([1, 0, 0], [0, 1, 0], 2)
    1.4142135623730951
    >>> distance.minkowski([1, 0, 0], [0, 1, 0], 3)
    1.2599210498948732
    >>> distance.minkowski([1, 1, 0], [0, 1, 0], 1)
    1.0
    >>> distance.minkowski([1, 1, 0], [0, 1, 0], 2)
    1.0
    >>> distance.minkowski([1, 1, 0], [0, 1, 0], 3)
    1.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if p < 1:
        raise ValueError("p must be at least 1")
    u_v = u - v
    if w is not None:
        w = _validate_weights(w)
        if p == 1:
            root_w = w
        if p == 2:
            # better precision and speed
            root_w = np.sqrt(w)
        else:
            root_w = np.power(w, 1/p)
        u_v = root_w * u_v
    dist = norm(u_v, ord=p)
    return dist


# `minkowski` gained weights in scipy 1.0.  Once we're at say version 1.3,
# deprecated `wminkowski`.  Not done at once because it would be annoying for
# downstream libraries that used `wminkowski` and support multiple scipy
# versions.
def wminkowski(u, v, p, w):
    """
    Compute the weighted Minkowski distance between two 1-D arrays.

    The weighted Minkowski distance between `u` and `v`, defined as

    .. math::

       \\left(\\sum{(|w_i (u_i - v_i)|^p)}\\right)^{1/p}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    p : int
        The order of the norm of the difference :math:`{||u-v||}_p`.
    w : (N,) array_like
        The weight vector.

    Returns
    -------
    wminkowski : double
        The weighted Minkowski distance between vectors `u` and `v`.

    Notes
    -----
    `wminkowski` is DEPRECATED. It implements a definition where weights
    are powered. It is recommended to use the weighted version of `minkowski`
    instead. This function will be removed in a future version of scipy.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.wminkowski([1, 0, 0], [0, 1, 0], 1, np.ones(3))
    2.0
    >>> distance.wminkowski([1, 0, 0], [0, 1, 0], 2, np.ones(3))
    1.4142135623730951
    >>> distance.wminkowski([1, 0, 0], [0, 1, 0], 3, np.ones(3))
    1.2599210498948732
    >>> distance.wminkowski([1, 1, 0], [0, 1, 0], 1, np.ones(3))
    1.0
    >>> distance.wminkowski([1, 1, 0], [0, 1, 0], 2, np.ones(3))
    1.0
    >>> distance.wminkowski([1, 1, 0], [0, 1, 0], 3, np.ones(3))
    1.0

    """
    w = _validate_weights(w)
    return minkowski(u, v, p=p, w=w**p)


def euclidean(u, v, w=None):
    """
    Computes the Euclidean distance between two 1-D arrays.

    The Euclidean distance between 1-D arrays `u` and `v`, is defined as

    .. math::

       {||u-v||}_2

       \\left(\\sum{(w_i |(u_i - v_i)|^2)}\\right)^{1/2}

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    euclidean : double
        The Euclidean distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.euclidean([1, 0, 0], [0, 1, 0])
    1.4142135623730951
    >>> distance.euclidean([1, 1, 0], [0, 1, 0])
    1.0

    """
    return minkowski(u, v, p=2, w=w)


def sqeuclidean(u, v, w=None):
    """
    Compute the squared Euclidean distance between two 1-D arrays.

    The squared Euclidean distance between `u` and `v` is defined as

    .. math::

       {||u-v||}_2^2

       \\left(\\sum{(w_i |(u_i - v_i)|^2)}\\right)

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    sqeuclidean : double
        The squared Euclidean distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.sqeuclidean([1, 0, 0], [0, 1, 0])
    2.0
    >>> distance.sqeuclidean([1, 1, 0], [0, 1, 0])
    1.0

    """
    # Preserve float dtypes, but convert everything else to np.float64
    # for stability.
    utype, vtype = None, None
    if not (hasattr(u, "dtype") and np.issubdtype(u.dtype, np.inexact)):
        utype = np.float64
    if not (hasattr(v, "dtype") and np.issubdtype(v.dtype, np.inexact)):
        vtype = np.float64

    u = _validate_vector(u, dtype=utype)
    v = _validate_vector(v, dtype=vtype)
    u_v = u - v
    u_v_w = u_v  # only want weights applied once
    if w is not None:
        w = _validate_weights(w)
        u_v_w = w * u_v
    return np.dot(u_v, u_v_w)


def correlation(u, v, w=None, centered=True):
    """
    Compute the correlation distance between two 1-D arrays.

    The correlation distance between `u` and `v`, is
    defined as

    .. math::

        1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
                  {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

    where :math:`\\bar{u}` is the mean of the elements of `u`
    and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    correlation : double
        The correlation distance between 1-D array `u` and `v`.

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is not None:
        w = _validate_weights(w)
    if centered:
        umu = np.average(u, weights=w)
        vmu = np.average(v, weights=w)
        u = u - umu
        v = v - vmu
    uv = np.average(u * v, weights=w)
    uu = np.average(np.square(u), weights=w)
    vv = np.average(np.square(v), weights=w)
    dist = 1.0 - uv / np.sqrt(uu * vv)
    return dist


def cosine(u, v, w=None):
    """
    Compute the Cosine distance between 1-D arrays.

    The Cosine distance between `u` and `v`, is defined as

    .. math::

        1 - \\frac{u \\cdot v}
                  {||u||_2 ||v||_2}.

    where :math:`u \\cdot v` is the dot product of :math:`u` and
    :math:`v`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    cosine : double
        The Cosine distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.cosine([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.cosine([100, 0, 0], [0, 1, 0])
    1.0
    >>> distance.cosine([1, 1, 0], [0, 1, 0])
    0.29289321881345254

    """
    # cosine distance is also referred to as 'uncentered correlation',
    #   or 'reflective correlation'
    return correlation(u, v, w=w, centered=False)


def hamming(u, v, w=None):
    """
    Compute the Hamming distance between two 1-D arrays.

    The Hamming distance between 1-D arrays `u` and `v`, is simply the
    proportion of disagreeing components in `u` and `v`. If `u` and `v` are
    boolean vectors, the Hamming distance is

    .. math::

       \\frac{c_{01} + c_{10}}{n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    hamming : double
        The Hamming distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.hamming([1, 0, 0], [0, 1, 0])
    0.66666666666666663
    >>> distance.hamming([1, 0, 0], [1, 1, 0])
    0.33333333333333331
    >>> distance.hamming([1, 0, 0], [2, 0, 0])
    0.33333333333333331
    >>> distance.hamming([1, 0, 0], [3, 0, 0])
    0.33333333333333331

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.shape != v.shape:
        raise ValueError('The 1d arrays must have equal lengths.')
    u_ne_v = u != v
    if w is not None:
        w = _validate_weights(w)
    return np.average(u_ne_v, weights=w)


def jaccard(u, v, w=None):
    """
    Compute the Jaccard-Needham dissimilarity between two boolean 1-D arrays.

    The Jaccard-Needham dissimilarity between 1-D boolean arrays `u` and `v`,
    is defined as

    .. math::

       \\frac{c_{TF} + c_{FT}}
            {c_{TT} + c_{FT} + c_{TF}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    jaccard : double
        The Jaccard distance between vectors `u` and `v`.

    Notes
    -----
    When both `u` and `v` lead to a `0/0` division i.e. there is no overlap
    between the items in the vectors the returned distance is 0. See the
    Wikipedia page on the Jaccard index [1]_, and this paper [2]_.

    .. versionchanged:: 1.2.0
        Previously, when `u` and `v` lead to a `0/0` division, the function
        would return NaN. This was changed to return 0 instead.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Jaccard_index
    .. [2] S. Kosub, "A note on the triangle inequality for the Jaccard
       distance", 2016, Available online: https://arxiv.org/pdf/1612.02696.pdf

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.jaccard([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.jaccard([1, 0, 0], [1, 1, 0])
    0.5
    >>> distance.jaccard([1, 0, 0], [1, 2, 0])
    0.5
    >>> distance.jaccard([1, 0, 0], [1, 1, 1])
    0.66666666666666663

    """
    u = _validate_vector(u)
    v = _validate_vector(v)

    nonzero = np.bitwise_or(u != 0, v != 0)
    unequal_nonzero = np.bitwise_and((u != v), nonzero)
    if w is not None:
        w = _validate_weights(w)
        nonzero = w * nonzero
        unequal_nonzero = w * unequal_nonzero
    a = np.double(unequal_nonzero.sum())
    b = np.double(nonzero.sum())
    return (a / b) if b != 0 else 0


def kulsinski(u, v, w=None):
    """
    Compute the Kulsinski dissimilarity between two boolean 1-D arrays.

    The Kulsinski dissimilarity between two boolean 1-D arrays `u` and `v`,
    is defined as

    .. math::

         \\frac{c_{TF} + c_{FT} - c_{TT} + n}
              {c_{FT} + c_{TF} + n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    kulsinski : double
        The Kulsinski distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.kulsinski([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.kulsinski([1, 0, 0], [1, 1, 0])
    0.75
    >>> distance.kulsinski([1, 0, 0], [2, 1, 0])
    0.33333333333333331
    >>> distance.kulsinski([1, 0, 0], [3, 1, 0])
    -0.5

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is None:
        n = float(len(u))
    else:
        w = _validate_weights(w)
        n = w.sum()
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v, w=w)

    return (ntf + nft - ntt + n) / (ntf + nft + n)


def seuclidean(u, v, V):
    """
    Return the standardized Euclidean distance between two 1-D arrays.

    The standardized Euclidean distance between `u` and `v`.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    V : (N,) array_like
        `V` is an 1-D array of component variances. It is usually computed
        among a larger collection vectors.

    Returns
    -------
    seuclidean : double
        The standardized Euclidean distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.seuclidean([1, 0, 0], [0, 1, 0], [0.1, 0.1, 0.1])
    4.4721359549995796
    >>> distance.seuclidean([1, 0, 0], [0, 1, 0], [1, 0.1, 0.1])
    3.3166247903553998
    >>> distance.seuclidean([1, 0, 0], [0, 1, 0], [10, 0.1, 0.1])
    3.1780497164141406

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    V = _validate_vector(V, dtype=np.float64)
    if V.shape[0] != u.shape[0] or u.shape[0] != v.shape[0]:
        raise TypeError('V must be a 1-D array of the same dimension '
                        'as u and v.')
    return euclidean(u, v, w=1/V)


def cityblock(u, v, w=None):
    """
    Compute the City Block (Manhattan) distance.

    Computes the Manhattan distance between two 1-D arrays `u` and `v`,
    which is defined as

    .. math::

       \\sum_i {\\left| u_i - v_i \\right|}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    cityblock : double
        The City Block (Manhattan) distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.cityblock([1, 0, 0], [0, 1, 0])
    2
    >>> distance.cityblock([1, 0, 0], [0, 2, 0])
    3
    >>> distance.cityblock([1, 0, 0], [1, 1, 0])
    1

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    l1_diff = abs(u - v)
    if w is not None:
        w = _validate_weights(w)
        l1_diff = w * l1_diff
    return l1_diff.sum()


def mahalanobis(u, v, VI):
    """
    Compute the Mahalanobis distance between two 1-D arrays.

    The Mahalanobis distance between 1-D arrays `u` and `v`, is defined as

    .. math::

       \\sqrt{ (u-v) V^{-1} (u-v)^T }

    where ``V`` is the covariance matrix.  Note that the argument `VI`
    is the inverse of ``V``.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    VI : ndarray
        The inverse of the covariance matrix.

    Returns
    -------
    mahalanobis : double
        The Mahalanobis distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> iv = [[1, 0.5, 0.5], [0.5, 1, 0.5], [0.5, 0.5, 1]]
    >>> distance.mahalanobis([1, 0, 0], [0, 1, 0], iv)
    1.0
    >>> distance.mahalanobis([0, 2, 0], [0, 1, 0], iv)
    1.0
    >>> distance.mahalanobis([2, 0, 0], [0, 1, 0], iv)
    1.7320508075688772

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    VI = np.atleast_2d(VI)
    delta = u - v
    m = np.dot(np.dot(delta, VI), delta)
    return np.sqrt(m)


def chebyshev(u, v, w=None):
    """
    Compute the Chebyshev distance.

    Computes the Chebyshev distance between two 1-D arrays `u` and `v`,
    which is defined as

    .. math::

       \\max_i {|u_i-v_i|}.

    Parameters
    ----------
    u : (N,) array_like
        Input vector.
    v : (N,) array_like
        Input vector.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    chebyshev : double
        The Chebyshev distance between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.chebyshev([1, 0, 0], [0, 1, 0])
    1
    >>> distance.chebyshev([1, 1, 0], [0, 1, 0])
    1

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is not None:
        w = _validate_weights(w)
        has_weight = w > 0
        if has_weight.sum() < w.size:
            u = u[has_weight]
            v = v[has_weight]
    return max(abs(u - v))


def braycurtis(u, v, w=None):
    """
    Compute the Bray-Curtis distance between two 1-D arrays.

    Bray-Curtis distance is defined as

    .. math::

       \\sum{|u_i-v_i|} / \\sum{|u_i+v_i|}

    The Bray-Curtis distance is in the range [0, 1] if all coordinates are
    positive, and is undefined if the inputs are of length zero.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    braycurtis : double
        The Bray-Curtis distance between 1-D arrays `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.braycurtis([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.braycurtis([1, 1, 0], [0, 1, 0])
    0.33333333333333331

    """
    u = _validate_vector(u)
    v = _validate_vector(v, dtype=np.float64)
    l1_diff = abs(u - v)
    l1_sum = abs(u + v)
    if w is not None:
        w = _validate_weights(w)
        l1_diff = w * l1_diff
        l1_sum = w * l1_sum
    return l1_diff.sum() / l1_sum.sum()


def canberra(u, v, w=None):
    """
    Compute the Canberra distance between two 1-D arrays.

    The Canberra distance is defined as

    .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    canberra : double
        The Canberra distance between vectors `u` and `v`.

    Notes
    -----
    When `u[i]` and `v[i]` are 0 for given i, then the fraction 0/0 = 0 is
    used in the calculation.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.canberra([1, 0, 0], [0, 1, 0])
    2.0
    >>> distance.canberra([1, 1, 0], [0, 1, 0])
    1.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v, dtype=np.float64)
    if w is not None:
        w = _validate_weights(w)
    olderr = np.seterr(invalid='ignore')
    try:
        abs_uv = abs(u - v)
        abs_u = abs(u)
        abs_v = abs(v)
        d = abs_uv / (abs_u + abs_v)
        if w is not None:
            d = w * d
        d = np.nansum(d)
    finally:
        np.seterr(**olderr)
    return d


def jensenshannon(p, q, base=None):
    """
    Compute the Jensen-Shannon distance (metric) between
    two 1-D probability arrays. This is the square root
    of the Jensen-Shannon divergence.

    The Jensen-Shannon distance between two probability
    vectors `p` and `q` is defined as,

    .. math::

       \\sqrt{\\frac{D(p \\parallel m) + D(q \\parallel m)}{2}}

    where :math:`m` is the pointwise mean of :math:`p` and :math:`q`
    and :math:`D` is the Kullback-Leibler divergence.

    This routine will normalize `p` and `q` if they don't sum to 1.0.

    Parameters
    ----------
    p : (N,) array_like
        left probability vector
    q : (N,) array_like
        right probability vector
    base : double, optional
        the base of the logarithm used to compute the output
        if not given, then the routine uses the default base of
        scipy.stats.entropy.

    Returns
    -------
    js : double
        The Jensen-Shannon distance between `p` and `q`

    .. versionadded:: 1.2.0

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.jensenshannon([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0)
    1.0
    >>> distance.jensenshannon([1.0, 0.0], [0.5, 0.5])
    0.46450140402245893
    >>> distance.jensenshannon([1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
    0.0

    """
    p = np.asarray(p)
    q = np.asarray(q)
    p = p / np.sum(p, axis=0)
    q = q / np.sum(q, axis=0)
    m = (p + q) / 2.0
    left = rel_entr(p, m)
    right = rel_entr(q, m)
    js = np.sum(left, axis=0) + np.sum(right, axis=0)
    if base is not None:
        js /= np.log(base)
    return np.sqrt(js / 2.0)


def yule(u, v, w=None):
    """
    Compute the Yule dissimilarity between two boolean 1-D arrays.

    The Yule dissimilarity is defined as

    .. math::

         \\frac{R}{c_{TT} * c_{FF} + \\frac{R}{2}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2.0 * c_{TF} * c_{FT}`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    yule : double
        The Yule dissimilarity between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.yule([1, 0, 0], [0, 1, 0])
    2.0
    >>> distance.yule([1, 1, 0], [0, 1, 0])
    0.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is not None:
        w = _validate_weights(w)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v, w=w)
    return float(2.0 * ntf * nft / np.array(ntt * nff + ntf * nft))


@np.deprecate(message="spatial.distance.matching is deprecated in scipy 1.0.0; "
                      "use spatial.distance.hamming instead.")
def matching(u, v, w=None):
    """
    Compute the Hamming distance between two boolean 1-D arrays.

    This is a deprecated synonym for :func:`hamming`.
    """
    return hamming(u, v, w=w)


def dice(u, v, w=None):
    """
    Compute the Dice dissimilarity between two boolean 1-D arrays.

    The Dice dissimilarity between `u` and `v`, is

    .. math::

         \\frac{c_{TF} + c_{FT}}
              {2c_{TT} + c_{FT} + c_{TF}}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) ndarray, bool
        Input 1-D array.
    v : (N,) ndarray, bool
        Input 1-D array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    dice : double
        The Dice dissimilarity between 1-D arrays `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.dice([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.dice([1, 0, 0], [1, 1, 0])
    0.3333333333333333
    >>> distance.dice([1, 0, 0], [2, 0, 0])
    -0.3333333333333333

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is not None:
        w = _validate_weights(w)
    if u.dtype == v.dtype == bool and w is None:
        ntt = (u & v).sum()
    else:
        dtype = np.find_common_type([int], [u.dtype, v.dtype])
        u = u.astype(dtype)
        v = v.astype(dtype)
        if w is None:
            ntt = (u * v).sum()
        else:
            ntt = (u * v * w).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v, w=w)
    return float((ntf + nft) / np.array(2.0 * ntt + ntf + nft))


def rogerstanimoto(u, v, w=None):
    """
    Compute the Rogers-Tanimoto dissimilarity between two boolean 1-D arrays.

    The Rogers-Tanimoto dissimilarity between two boolean 1-D arrays
    `u` and `v`, is defined as

    .. math::
       \\frac{R}
            {c_{TT} + c_{FF} + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2(c_{TF} + c_{FT})`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    rogerstanimoto : double
        The Rogers-Tanimoto dissimilarity between vectors
        `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.rogerstanimoto([1, 0, 0], [0, 1, 0])
    0.8
    >>> distance.rogerstanimoto([1, 0, 0], [1, 1, 0])
    0.5
    >>> distance.rogerstanimoto([1, 0, 0], [2, 0, 0])
    -1.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if w is not None:
        w = _validate_weights(w)
    (nff, nft, ntf, ntt) = _nbool_correspond_all(u, v, w=w)
    return float(2.0 * (ntf + nft)) / float(ntt + nff + (2.0 * (ntf + nft)))


def russellrao(u, v, w=None):
    """
    Compute the Russell-Rao dissimilarity between two boolean 1-D arrays.

    The Russell-Rao dissimilarity between two boolean 1-D arrays, `u` and
    `v`, is defined as

    .. math::

      \\frac{n - c_{TT}}
           {n}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    russellrao : double
        The Russell-Rao dissimilarity between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.russellrao([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.russellrao([1, 0, 0], [1, 1, 0])
    0.6666666666666666
    >>> distance.russellrao([1, 0, 0], [2, 0, 0])
    0.3333333333333333

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == v.dtype == bool and w is None:
        ntt = (u & v).sum()
        n = float(len(u))
    elif w is None:
        ntt = (u * v).sum()
        n = float(len(u))
    else:
        w = _validate_weights(w)
        ntt = (u * v * w).sum()
        n = w.sum()
    return float(n - ntt) / n


def sokalmichener(u, v, w=None):
    """
    Compute the Sokal-Michener dissimilarity between two boolean 1-D arrays.

    The Sokal-Michener dissimilarity between boolean 1-D arrays `u` and `v`,
    is defined as

    .. math::

       \\frac{R}
            {S + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n`, :math:`R = 2 * (c_{TF} + c_{FT})` and
    :math:`S = c_{FF} + c_{TT}`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    sokalmichener : double
        The Sokal-Michener dissimilarity between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.sokalmichener([1, 0, 0], [0, 1, 0])
    0.8
    >>> distance.sokalmichener([1, 0, 0], [1, 1, 0])
    0.5
    >>> distance.sokalmichener([1, 0, 0], [2, 0, 0])
    -1.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == v.dtype == bool and w is None:
        ntt = (u & v).sum()
        nff = (~u & ~v).sum()
    elif w is None:
        ntt = (u * v).sum()
        nff = ((1.0 - u) * (1.0 - v)).sum()
    else:
        w = _validate_weights(w)
        ntt = (u * v * w).sum()
        nff = ((1.0 - u) * (1.0 - v) * w).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v)
    return float(2.0 * (ntf + nft)) / float(ntt + nff + 2.0 * (ntf + nft))


def sokalsneath(u, v, w=None):
    """
    Compute the Sokal-Sneath dissimilarity between two boolean 1-D arrays.

    The Sokal-Sneath dissimilarity between `u` and `v`,

    .. math::

       \\frac{R}
            {c_{TT} + R}

    where :math:`c_{ij}` is the number of occurrences of
    :math:`\\mathtt{u[k]} = i` and :math:`\\mathtt{v[k]} = j` for
    :math:`k < n` and :math:`R = 2(c_{TF} + c_{FT})`.

    Parameters
    ----------
    u : (N,) array_like, bool
        Input array.
    v : (N,) array_like, bool
        Input array.
    w : (N,) array_like, optional
        The weights for each value in `u` and `v`. Default is None,
        which gives each value a weight of 1.0

    Returns
    -------
    sokalsneath : double
        The Sokal-Sneath dissimilarity between vectors `u` and `v`.

    Examples
    --------
    >>> from scipy.spatial import distance
    >>> distance.sokalsneath([1, 0, 0], [0, 1, 0])
    1.0
    >>> distance.sokalsneath([1, 0, 0], [1, 1, 0])
    0.66666666666666663
    >>> distance.sokalsneath([1, 0, 0], [2, 1, 0])
    0.0
    >>> distance.sokalsneath([1, 0, 0], [3, 1, 0])
    -2.0

    """
    u = _validate_vector(u)
    v = _validate_vector(v)
    if u.dtype == v.dtype == bool and w is None:
        ntt = (u & v).sum()
    elif w is None:
        ntt = (u * v).sum()
    else:
        w = _validate_weights(w)
        ntt = (u * v * w).sum()
    (nft, ntf) = _nbool_correspond_ft_tf(u, v, w=w)
    denom = np.array(ntt + 2.0 * (ntf + nft))
    if not denom.any():
        raise ValueError('Sokal-Sneath dissimilarity is not defined for '
                         'vectors that are entirely false.')
    return float(2.0 * (ntf + nft)) / denom


_convert_to_double = partial(_convert_to_type, out_type=np.double)
_convert_to_bool = partial(_convert_to_type, out_type=bool)

# adding python-only wrappers to _distance_wrap module
_distance_wrap.pdist_correlation_double_wrap = _correlation_pdist_wrap
_distance_wrap.cdist_correlation_double_wrap = _correlation_cdist_wrap

# Registry of implemented metrics:
# Dictionary with the following structure:
# {
#  metric_name : MetricInfo(aka, types=[double], validator=None)
# }
#
# Where:
# `metric_name` must be equal to python metric name
#
# MetricInfo is a named tuple with fields:
#  'aka' : [list of aliases],
#
#  'validator': f(X, m, n, **kwargs)    # function that check kwargs and
#                                       # computes default values.
#
#  'types': [list of supported types],  # X (pdist) and XA (cdist) are used to
#                                       # choose the type. if there is no match
#                                       # the first type is used. Default double
# }
MetricInfo = namedtuple("MetricInfo", 'aka types validator ')
MetricInfo.__new__.__defaults__ = (['double'], None)

_METRICS = {
    'braycurtis': MetricInfo(aka=['braycurtis']),
    'canberra': MetricInfo(aka=['canberra']),
    'chebyshev': MetricInfo(aka=['chebychev', 'chebyshev', 'cheby', 'cheb', 'ch']),
    'cityblock': MetricInfo(aka=['cityblock', 'cblock', 'cb', 'c']),
    'correlation': MetricInfo(aka=['correlation', 'co']),
    'cosine': MetricInfo(aka=['cosine', 'cos']),
    'dice': MetricInfo(aka=['dice'], types=['bool']),
    'euclidean': MetricInfo(aka=['euclidean', 'euclid', 'eu', 'e']),
    'hamming': MetricInfo(aka=['matching', 'hamming', 'hamm', 'ha', 'h'],
                          types=['double', 'bool'],
                          validator=_validate_hamming_kwargs),
    'jaccard': MetricInfo(aka=['jaccard', 'jacc', 'ja', 'j'],
                          types=['double', 'bool']),
    'jensenshannon': MetricInfo(aka=['jensenshannon', 'js'],
                                types=['double']),
    'kulsinski': MetricInfo(aka=['kulsinski'], types=['bool']),
    'mahalanobis': MetricInfo(aka=['mahalanobis', 'mahal', 'mah'],
                              validator=_validate_mahalanobis_kwargs),
    'minkowski': MetricInfo(aka=['minkowski', 'mi', 'm', 'pnorm'],
                            validator=_validate_minkowski_kwargs),
    'rogerstanimoto': MetricInfo(aka=['rogerstanimoto'], types=['bool']),
    'russellrao': MetricInfo(aka=['russellrao'], types=['bool']),
    'seuclidean': MetricInfo(aka=['seuclidean', 'se', 's'],
                             validator=_validate_seuclidean_kwargs),
    'sokalmichener': MetricInfo(aka=['sokalmichener'], types=['bool']),
    'sokalsneath': MetricInfo(aka=['sokalsneath'], types=['bool']),
    'sqeuclidean': MetricInfo(aka=['sqeuclidean', 'sqe', 'sqeuclid']),
    'wminkowski': MetricInfo(aka=['wminkowski', 'wmi', 'wm', 'wpnorm'],
                             validator=_validate_wminkowski_kwargs),
    'yule': MetricInfo(aka=['yule'], types=['bool']),
}


_METRIC_ALIAS = dict((alias, name)
                     for name, info in _METRICS.items()
                     for alias in info.aka)

_METRICS_NAMES = list(_METRICS.keys())

_TEST_METRICS = {'test_' + name: globals()[name] for name in _METRICS.keys()}


def _select_weighted_metric(mstr, kwargs, out):
    kwargs = dict(kwargs)

    if "w" in kwargs and kwargs["w"] is None:
        # w=None is the same as omitting it
        kwargs.pop("w")

    if mstr.startswith("test_") or mstr in _METRICS['wminkowski'].aka + _METRICS['hamming'].aka:
        # These support weights
        pass
    elif "w" in kwargs:
        if (mstr in _METRICS['seuclidean'].aka or
                mstr in _METRICS['mahalanobis'].aka):
            raise ValueError("metric %s incompatible with weights" % mstr)

        # XXX: C-versions do not support weights
        # need to use python version for weighting
        kwargs['out'] = out
        mstr = "test_%s" % mstr

    return mstr, kwargs


def pdist(X, metric='euclidean', *args, **kwargs):
    """
    Pairwise distances between observations in n-dimensional space.

    See Notes for common calling conventions.

    Parameters
    ----------
    X : ndarray
        An m by n array of m original observations in an
        n-dimensional space.
    metric : str or function, optional
        The distance metric to use. The distance function can
        be 'braycurtis', 'canberra', 'chebyshev', 'cityblock',
        'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
        'jaccard', 'jensenshannon', 'kulsinski', 'mahalanobis', 'matching',
        'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
        'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'.
    *args : tuple. Deprecated.
        Additional arguments should be passed as keyword arguments
    **kwargs : dict, optional
        Extra arguments to `metric`: refer to each metric documentation for a
        list of all possible arguments.

        Some possible arguments:

        p : scalar
        The p-norm to apply for Minkowski, weighted and unweighted.
        Default: 2.

        w : ndarray
        The weight vector for metrics that support weights (e.g., Minkowski).

        V : ndarray
        The variance vector for standardized Euclidean.
        Default: var(X, axis=0, ddof=1)

        VI : ndarray
        The inverse of the covariance matrix for Mahalanobis.
        Default: inv(cov(X.T)).T

        out : ndarray.
        The output array
        If not None, condensed distance matrix Y is stored in this array.
        Note: metric independent, it will become a regular keyword arg in a
        future scipy version

    Returns
    -------
    Y : ndarray
        Returns a condensed distance matrix Y.  For
        each :math:`i` and :math:`j` (where :math:`i<j<m`),where m is the number
        of original observations. The metric ``dist(u=X[i], v=X[j])``
        is computed and stored in entry ``ij``.

    See Also
    --------
    squareform : converts between condensed distance matrices and
                 square distance matrices.

    Notes
    -----
    See ``squareform`` for information on how to calculate the index of
    this entry or to convert the condensed distance matrix to a
    redundant square matrix.

    The following are common calling conventions.

    1. ``Y = pdist(X, 'euclidean')``

       Computes the distance between m points using Euclidean distance
       (2-norm) as the distance metric between the points. The points
       are arranged as m n-dimensional row vectors in the matrix X.

    2. ``Y = pdist(X, 'minkowski', p=2.)``

       Computes the distances using the Minkowski distance
       :math:`||u-v||_p` (p-norm) where :math:`p \\geq 1`.

    3. ``Y = pdist(X, 'cityblock')``

       Computes the city block or Manhattan distance between the
       points.

    4. ``Y = pdist(X, 'seuclidean', V=None)``

       Computes the standardized Euclidean distance. The standardized
       Euclidean distance between two n-vectors ``u`` and ``v`` is

       .. math::

          \\sqrt{\\sum {(u_i-v_i)^2 / V[x_i]}}


       V is the variance vector; V[i] is the variance computed over all
       the i'th components of the points.  If not passed, it is
       automatically computed.

    5. ``Y = pdist(X, 'sqeuclidean')``

       Computes the squared Euclidean distance :math:`||u-v||_2^2` between
       the vectors.

    6. ``Y = pdist(X, 'cosine')``

       Computes the cosine distance between vectors u and v,

       .. math::

          1 - \\frac{u \\cdot v}
                   {{||u||}_2 {||v||}_2}

       where :math:`||*||_2` is the 2-norm of its argument ``*``, and
       :math:`u \\cdot v` is the dot product of ``u`` and ``v``.

    7. ``Y = pdist(X, 'correlation')``

       Computes the correlation distance between vectors u and v. This is

       .. math::

          1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
                   {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

       where :math:`\\bar{v}` is the mean of the elements of vector v,
       and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.

    8. ``Y = pdist(X, 'hamming')``

       Computes the normalized Hamming distance, or the proportion of
       those vector elements between two n-vectors ``u`` and ``v``
       which disagree. To save memory, the matrix ``X`` can be of type
       boolean.

    9. ``Y = pdist(X, 'jaccard')``

       Computes the Jaccard distance between the points. Given two
       vectors, ``u`` and ``v``, the Jaccard distance is the
       proportion of those elements ``u[i]`` and ``v[i]`` that
       disagree.

    10. ``Y = pdist(X, 'chebyshev')``

       Computes the Chebyshev distance between the points. The
       Chebyshev distance between two n-vectors ``u`` and ``v`` is the
       maximum norm-1 distance between their respective elements. More
       precisely, the distance is given by

       .. math::

          d(u,v) = \\max_i {|u_i-v_i|}

    11. ``Y = pdist(X, 'canberra')``

       Computes the Canberra distance between the points. The
       Canberra distance between two points ``u`` and ``v`` is

       .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}


    12. ``Y = pdist(X, 'braycurtis')``

       Computes the Bray-Curtis distance between the points. The
       Bray-Curtis distance between two points ``u`` and ``v`` is


       .. math::

            d(u,v) = \\frac{\\sum_i {|u_i-v_i|}}
                           {\\sum_i {|u_i+v_i|}}

    13. ``Y = pdist(X, 'mahalanobis', VI=None)``

       Computes the Mahalanobis distance between the points. The
       Mahalanobis distance between two points ``u`` and ``v`` is
       :math:`\\sqrt{(u-v)(1/V)(u-v)^T}` where :math:`(1/V)` (the ``VI``
       variable) is the inverse covariance. If ``VI`` is not None,
       ``VI`` will be used as the inverse covariance matrix.

    14. ``Y = pdist(X, 'yule')``

       Computes the Yule distance between each pair of boolean
       vectors. (see yule function documentation)

    15. ``Y = pdist(X, 'matching')``

       Synonym for 'hamming'.

    16. ``Y = pdist(X, 'dice')``

       Computes the Dice distance between each pair of boolean
       vectors. (see dice function documentation)

    17. ``Y = pdist(X, 'kulsinski')``

       Computes the Kulsinski distance between each pair of
       boolean vectors. (see kulsinski function documentation)

    18. ``Y = pdist(X, 'rogerstanimoto')``

       Computes the Rogers-Tanimoto distance between each pair of
       boolean vectors. (see rogerstanimoto function documentation)

    19. ``Y = pdist(X, 'russellrao')``

       Computes the Russell-Rao distance between each pair of
       boolean vectors. (see russellrao function documentation)

    20. ``Y = pdist(X, 'sokalmichener')``

       Computes the Sokal-Michener distance between each pair of
       boolean vectors. (see sokalmichener function documentation)

    21. ``Y = pdist(X, 'sokalsneath')``

       Computes the Sokal-Sneath distance between each pair of
       boolean vectors. (see sokalsneath function documentation)

    22. ``Y = pdist(X, 'wminkowski', p=2, w=w)``

       Computes the weighted Minkowski distance between each pair of
       vectors. (see wminkowski function documentation)

    23. ``Y = pdist(X, f)``

       Computes the distance between all pairs of vectors in X
       using the user supplied 2-arity function f. For example,
       Euclidean distance between the vectors could be computed
       as follows::

         dm = pdist(X, lambda u, v: np.sqrt(((u-v)**2).sum()))

       Note that you should avoid passing a reference to one of
       the distance functions defined in this library. For example,::

         dm = pdist(X, sokalsneath)

       would calculate the pair-wise distances between the vectors in
       X using the Python function sokalsneath. This would result in
       sokalsneath being called :math:`{n \\choose 2}` times, which
       is inefficient. Instead, the optimized C version is more
       efficient, and we call it using the following syntax.::

         dm = pdist(X, 'sokalsneath')

    """
    # You can also call this as:
    #     Y = pdist(X, 'test_abc')
    # where 'abc' is the metric being tested.  This computes the distance
    # between all pairs of vectors in X using the distance metric 'abc' but
    # with a more succinct, verifiable, but less efficient implementation.

    X = _asarray_validated(X, sparse_ok=False, objects_ok=True, mask_ok=True,
                           check_finite=False)
    kwargs = _args_to_kwargs_xdist(args, kwargs, metric, "pdist")

    X = np.asarray(X, order='c')

    s = X.shape
    if len(s) != 2:
        raise ValueError('A 2-dimensional array must be passed.')

    m, n = s
    out = kwargs.pop("out", None)
    if out is None:
        dm = np.empty((m * (m - 1)) // 2, dtype=np.double)
    else:
        if out.shape != (m * (m - 1) // 2,):
            raise ValueError("output array has incorrect shape.")
        if not out.flags.c_contiguous:
            raise ValueError("Output array must be C-contiguous.")
        if out.dtype != np.double:
            raise ValueError("Output array must be double type.")
        dm = out

    # compute blacklist for deprecated kwargs
    if(metric in _METRICS['jensenshannon'].aka
       or metric == 'test_jensenshannon' or metric == jensenshannon):
        kwargs_blacklist = ["p", "w", "V", "VI"]

    elif(metric in _METRICS['minkowski'].aka
         or metric in _METRICS['wminkowski'].aka
         or metric in ['test_minkowski', 'test_wminkowski']
         or metric in [minkowski, wminkowski]):
        kwargs_blacklist = ["V", "VI"]

    elif(metric in _METRICS['seuclidean'].aka or
         metric == 'test_seuclidean' or metric == seuclidean):
        kwargs_blacklist = ["p", "w", "VI"]

    elif(metric in _METRICS['mahalanobis'].aka
         or metric == 'test_mahalanobis' or metric == mahalanobis):
        kwargs_blacklist = ["p", "w", "V"]

    else:
        kwargs_blacklist = ["p", "V", "VI"]

    _filter_deprecated_kwargs(kwargs, kwargs_blacklist)

    if callable(metric):
        mstr = getattr(metric, '__name__', 'UnknownCustomMetric')
        metric_name = _METRIC_ALIAS.get(mstr, None)

        if metric_name is not None:
            X, typ, kwargs = _validate_pdist_input(X, m, n,
                                                   metric_name, **kwargs)

        k = 0
        for i in xrange(0, m - 1):
            for j in xrange(i + 1, m):
                dm[k] = metric(X[i], X[j], **kwargs)
                k = k + 1

    elif isinstance(metric, string_types):
        mstr = metric.lower()

        mstr, kwargs = _select_weighted_metric(mstr, kwargs, out)

        metric_name = _METRIC_ALIAS.get(mstr, None)

        if metric_name is not None:
            X, typ, kwargs = _validate_pdist_input(X, m, n,
                                                   metric_name, **kwargs)

            # get pdist wrapper
            pdist_fn = getattr(_distance_wrap,
                               "pdist_%s_%s_wrap" % (metric_name, typ))
            pdist_fn(X, dm, **kwargs)
            return dm

        elif mstr in ['old_cosine', 'old_cos']:
            warnings.warn('"old_cosine" is deprecated and will be removed in '
                          'a future version. Use "cosine" instead.',
                          DeprecationWarning)
            X = _convert_to_double(X)
            norms = np.einsum('ij,ij->i', X, X, dtype=np.double)
            np.sqrt(norms, out=norms)
            nV = norms.reshape(m, 1)
            # The numerator u * v
            nm = np.dot(X, X.T)
            # The denom. ||u||*||v||
            de = np.dot(nV, nV.T)
            dm = 1.0 - (nm / de)
            dm[xrange(0, m), xrange(0, m)] = 0.0
            dm = squareform(dm)
        elif mstr.startswith("test_"):
            if mstr in _TEST_METRICS:
                dm = pdist(X, _TEST_METRICS[mstr], **kwargs)
            else:
                raise ValueError('Unknown "Test" Distance Metric: %s' % mstr[5:])
        else:
            raise ValueError('Unknown Distance Metric: %s' % mstr)
    else:
        raise TypeError('2nd argument metric must be a string identifier '
                        'or a function.')
    return dm


def squareform(X, force="no", checks=True):
    """
    Convert a vector-form distance vector to a square-form distance
    matrix, and vice-versa.

    Parameters
    ----------
    X : ndarray
        Either a condensed or redundant distance matrix.
    force : str, optional
        As with MATLAB(TM), if force is equal to ``'tovector'`` or
        ``'tomatrix'``, the input will be treated as a distance matrix or
        distance vector respectively.
    checks : bool, optional
        If set to False, no checks will be made for matrix
        symmetry nor zero diagonals. This is useful if it is known that
        ``X - X.T1`` is small and ``diag(X)`` is close to zero.
        These values are ignored any way so they do not disrupt the
        squareform transformation.

    Returns
    -------
    Y : ndarray
        If a condensed distance matrix is passed, a redundant one is
        returned, or if a redundant one is passed, a condensed distance
        matrix is returned.

    Notes
    -----
    1. v = squareform(X)

       Given a square d-by-d symmetric distance matrix X,
       ``v = squareform(X)`` returns a ``d * (d-1) / 2`` (or
       :math:`{n \\choose 2}`) sized vector v.

      :math:`v[{n \\choose 2}-{n-i \\choose 2} + (j-i-1)]` is the distance
      between points i and j. If X is non-square or asymmetric, an error
      is returned.

    2. X = squareform(v)

      Given a ``d*(d-1)/2`` sized v for some integer ``d >= 2`` encoding
      distances as described, ``X = squareform(v)`` returns a d by d distance
      matrix X.  The ``X[i, j]`` and ``X[j, i]`` values are set to
      :math:`v[{n \\choose 2}-{n-i \\choose 2} + (j-i-1)]` and all
      diagonal elements are zero.

    In SciPy 0.19.0, ``squareform`` stopped casting all input types to
    float64, and started returning arrays of the same dtype as the input.

    """

    X = np.ascontiguousarray(X)

    s = X.shape

    if force.lower() == 'tomatrix':
        if len(s) != 1:
            raise ValueError("Forcing 'tomatrix' but input X is not a "
                             "distance vector.")
    elif force.lower() == 'tovector':
        if len(s) != 2:
            raise ValueError("Forcing 'tovector' but input X is not a "
                             "distance matrix.")

    # X = squareform(v)
    if len(s) == 1:
        if s[0] == 0:
            return np.zeros((1, 1), dtype=X.dtype)

        # Grab the closest value to the square root of the number
        # of elements times 2 to see if the number of elements
        # is indeed a binomial coefficient.
        d = int(np.ceil(np.sqrt(s[0] * 2)))

        # Check that v is of valid dimensions.
        if d * (d - 1) != s[0] * 2:
            raise ValueError('Incompatible vector size. It must be a binomial '
                             'coefficient n choose 2 for some integer n >= 2.')

        # Allocate memory for the distance matrix.
        M = np.zeros((d, d), dtype=X.dtype)

        # Since the C code does not support striding using strides.
        # The dimensions are used instead.
        X = _copy_array_if_base_present(X)

        # Fill in the values of the distance matrix.
        _distance_wrap.to_squareform_from_vector_wrap(M, X)

        # Return the distance matrix.
        return M
    elif len(s) == 2:
        if s[0] != s[1]:
            raise ValueError('The matrix argument must be square.')
        if checks:
            is_valid_dm(X, throw=True, name='X')

        # One-side of the dimensions is set here.
        d = s[0]

        if d <= 1:
            return np.array([], dtype=X.dtype)

        # Create a vector.
        v = np.zeros((d * (d - 1)) // 2, dtype=X.dtype)

        # Since the C code does not support striding using strides.
        # The dimensions are used instead.
        X = _copy_array_if_base_present(X)

        # Convert the vector to squareform.
        _distance_wrap.to_vector_from_squareform_wrap(X, v)
        return v
    else:
        raise ValueError(('The first argument must be one or two dimensional '
                          'array. A %d-dimensional array is not '
                          'permitted') % len(s))


def is_valid_dm(D, tol=0.0, throw=False, name="D", warning=False):
    """
    Return True if input array is a valid distance matrix.

    Distance matrices must be 2-dimensional numpy arrays.
    They must have a zero-diagonal, and they must be symmetric.

    Parameters
    ----------
    D : ndarray
        The candidate object to test for validity.
    tol : float, optional
        The distance matrix should be symmetric. `tol` is the maximum
        difference between entries ``ij`` and ``ji`` for the distance
        metric to be considered symmetric.
    throw : bool, optional
        An exception is thrown if the distance matrix passed is not valid.
    name : str, optional
        The name of the variable to checked. This is useful if
        throw is set to True so the offending variable can be identified
        in the exception message when an exception is thrown.
    warning : bool, optional
        Instead of throwing an exception, a warning message is
        raised.

    Returns
    -------
    valid : bool
        True if the variable `D` passed is a valid distance matrix.

    Notes
    -----
    Small numerical differences in `D` and `D.T` and non-zeroness of
    the diagonal are ignored if they are within the tolerance specified
    by `tol`.

    """
    D = np.asarray(D, order='c')
    valid = True
    try:
        s = D.shape
        if len(D.shape) != 2:
            if name:
                raise ValueError(('Distance matrix \'%s\' must have shape=2 '
                                  '(i.e. be two-dimensional).') % name)
            else:
                raise ValueError('Distance matrix must have shape=2 (i.e. '
                                 'be two-dimensional).')
        if tol == 0.0:
            if not (D == D.T).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' must be '
                                     'symmetric.') % name)
                else:
                    raise ValueError('Distance matrix must be symmetric.')
            if not (D[xrange(0, s[0]), xrange(0, s[0])] == 0).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' diagonal must '
                                      'be zero.') % name)
                else:
                    raise ValueError('Distance matrix diagonal must be zero.')
        else:
            if not (D - D.T <= tol).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' must be '
                                      'symmetric within tolerance %5.5f.')
                                     % (name, tol))
                else:
                    raise ValueError('Distance matrix must be symmetric within'
                                     ' tolerance %5.5f.' % tol)
            if not (D[xrange(0, s[0]), xrange(0, s[0])] <= tol).all():
                if name:
                    raise ValueError(('Distance matrix \'%s\' diagonal must be'
                                      ' close to zero within tolerance %5.5f.')
                                     % (name, tol))
                else:
                    raise ValueError(('Distance matrix \'%s\' diagonal must be'
                                      ' close to zero within tolerance %5.5f.')
                                     % tol)
    except Exception as e:
        if throw:
            raise
        if warning:
            warnings.warn(str(e))
        valid = False
    return valid


def is_valid_y(y, warning=False, throw=False, name=None):
    """
    Return True if the input array is a valid condensed distance matrix.

    Condensed distance matrices must be 1-dimensional numpy arrays.
    Their length must be a binomial coefficient :math:`{n \\choose 2}`
    for some positive integer n.

    Parameters
    ----------
    y : ndarray
        The condensed distance matrix.
    warning : bool, optional
        Invokes a warning if the variable passed is not a valid
        condensed distance matrix. The warning message explains why
        the distance matrix is not valid.  `name` is used when
        referencing the offending variable.
    throw : bool, optional
        Throws an exception if the variable passed is not a valid
        condensed distance matrix.
    name : bool, optional
        Used when referencing the offending variable in the
        warning or exception message.

    """
    y = np.asarray(y, order='c')
    valid = True
    try:
        if len(y.shape) != 1:
            if name:
                raise ValueError(('Condensed distance matrix \'%s\' must '
                                  'have shape=1 (i.e. be one-dimensional).')
                                 % name)
            else:
                raise ValueError('Condensed distance matrix must have shape=1 '
                                 '(i.e. be one-dimensional).')
        n = y.shape[0]
        d = int(np.ceil(np.sqrt(n * 2)))
        if (d * (d - 1) / 2) != n:
            if name:
                raise ValueError(('Length n of condensed distance matrix '
                                  '\'%s\' must be a binomial coefficient, i.e.'
                                  'there must be a k such that '
                                  '(k \\choose 2)=n)!') % name)
            else:
                raise ValueError('Length n of condensed distance matrix must '
                                 'be a binomial coefficient, i.e. there must '
                                 'be a k such that (k \\choose 2)=n)!')
    except Exception as e:
        if throw:
            raise
        if warning:
            warnings.warn(str(e))
        valid = False
    return valid


def num_obs_dm(d):
    """
    Return the number of original observations that correspond to a
    square, redundant distance matrix.

    Parameters
    ----------
    d : ndarray
        The target distance matrix.

    Returns
    -------
    num_obs_dm : int
        The number of observations in the redundant distance matrix.

    """
    d = np.asarray(d, order='c')
    is_valid_dm(d, tol=np.inf, throw=True, name='d')
    return d.shape[0]


def num_obs_y(Y):
    """
    Return the number of original observations that correspond to a
    condensed distance matrix.

    Parameters
    ----------
    Y : ndarray
        Condensed distance matrix.

    Returns
    -------
    n : int
        The number of observations in the condensed distance matrix `Y`.

    """
    Y = np.asarray(Y, order='c')
    is_valid_y(Y, throw=True, name='Y')
    k = Y.shape[0]
    if k == 0:
        raise ValueError("The number of observations cannot be determined on "
                         "an empty distance matrix.")
    d = int(np.ceil(np.sqrt(k * 2)))
    if (d * (d - 1) / 2) != k:
        raise ValueError("Invalid condensed distance matrix passed. Must be "
                         "some k where k=(n choose 2) for some n >= 2.")
    return d


def cdist(XA, XB, metric='euclidean', *args, **kwargs):
    """
    Compute distance between each pair of the two collections of inputs.

    See Notes for common calling conventions.

    Parameters
    ----------
    XA : ndarray
        An :math:`m_A` by :math:`n` array of :math:`m_A`
        original observations in an :math:`n`-dimensional space.
        Inputs are converted to float type.
    XB : ndarray
        An :math:`m_B` by :math:`n` array of :math:`m_B`
        original observations in an :math:`n`-dimensional space.
        Inputs are converted to float type.
    metric : str or callable, optional
        The distance metric to use.  If a string, the distance function can be
        'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
        'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon',
        'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
        'wminkowski', 'yule'.
    *args : tuple. Deprecated.
        Additional arguments should be passed as keyword arguments
    **kwargs : dict, optional
        Extra arguments to `metric`: refer to each metric documentation for a
        list of all possible arguments.

        Some possible arguments:

        p : scalar
        The p-norm to apply for Minkowski, weighted and unweighted.
        Default: 2.

        w : ndarray
        The weight vector for metrics that support weights (e.g., Minkowski).

        V : ndarray
        The variance vector for standardized Euclidean.
        Default: var(vstack([XA, XB]), axis=0, ddof=1)

        VI : ndarray
        The inverse of the covariance matrix for Mahalanobis.
        Default: inv(cov(vstack([XA, XB].T))).T

        out : ndarray
        The output array
        If not None, the distance matrix Y is stored in this array.
        Note: metric independent, it will become a regular keyword arg in a
        future scipy version

    Returns
    -------
    Y : ndarray
        A :math:`m_A` by :math:`m_B` distance matrix is returned.
        For each :math:`i` and :math:`j`, the metric
        ``dist(u=XA[i], v=XB[j])`` is computed and stored in the
        :math:`ij` th entry.

    Raises
    ------
    ValueError
        An exception is thrown if `XA` and `XB` do not have
        the same number of columns.

    Notes
    -----
    The following are common calling conventions:

    1. ``Y = cdist(XA, XB, 'euclidean')``

       Computes the distance between :math:`m` points using
       Euclidean distance (2-norm) as the distance metric between the
       points. The points are arranged as :math:`m`
       :math:`n`-dimensional row vectors in the matrix X.

    2. ``Y = cdist(XA, XB, 'minkowski', p=2.)``

       Computes the distances using the Minkowski distance
       :math:`||u-v||_p` (:math:`p`-norm) where :math:`p \\geq 1`.

    3. ``Y = cdist(XA, XB, 'cityblock')``

       Computes the city block or Manhattan distance between the
       points.

    4. ``Y = cdist(XA, XB, 'seuclidean', V=None)``

       Computes the standardized Euclidean distance. The standardized
       Euclidean distance between two n-vectors ``u`` and ``v`` is

       .. math::

          \\sqrt{\\sum {(u_i-v_i)^2 / V[x_i]}}.

       V is the variance vector; V[i] is the variance computed over all
       the i'th components of the points. If not passed, it is
       automatically computed.

    5. ``Y = cdist(XA, XB, 'sqeuclidean')``

       Computes the squared Euclidean distance :math:`||u-v||_2^2` between
       the vectors.

    6. ``Y = cdist(XA, XB, 'cosine')``

       Computes the cosine distance between vectors u and v,

       .. math::

          1 - \\frac{u \\cdot v}
                   {{||u||}_2 {||v||}_2}

       where :math:`||*||_2` is the 2-norm of its argument ``*``, and
       :math:`u \\cdot v` is the dot product of :math:`u` and :math:`v`.

    7. ``Y = cdist(XA, XB, 'correlation')``

       Computes the correlation distance between vectors u and v. This is

       .. math::

          1 - \\frac{(u - \\bar{u}) \\cdot (v - \\bar{v})}
                   {{||(u - \\bar{u})||}_2 {||(v - \\bar{v})||}_2}

       where :math:`\\bar{v}` is the mean of the elements of vector v,
       and :math:`x \\cdot y` is the dot product of :math:`x` and :math:`y`.


    8. ``Y = cdist(XA, XB, 'hamming')``

       Computes the normalized Hamming distance, or the proportion of
       those vector elements between two n-vectors ``u`` and ``v``
       which disagree. To save memory, the matrix ``X`` can be of type
       boolean.

    9. ``Y = cdist(XA, XB, 'jaccard')``

       Computes the Jaccard distance between the points. Given two
       vectors, ``u`` and ``v``, the Jaccard distance is the
       proportion of those elements ``u[i]`` and ``v[i]`` that
       disagree where at least one of them is non-zero.

    10. ``Y = cdist(XA, XB, 'chebyshev')``

       Computes the Chebyshev distance between the points. The
       Chebyshev distance between two n-vectors ``u`` and ``v`` is the
       maximum norm-1 distance between their respective elements. More
       precisely, the distance is given by

       .. math::

          d(u,v) = \\max_i {|u_i-v_i|}.

    11. ``Y = cdist(XA, XB, 'canberra')``

       Computes the Canberra distance between the points. The
       Canberra distance between two points ``u`` and ``v`` is

       .. math::

         d(u,v) = \\sum_i \\frac{|u_i-v_i|}
                              {|u_i|+|v_i|}.

    12. ``Y = cdist(XA, XB, 'braycurtis')``

       Computes the Bray-Curtis distance between the points. The
       Bray-Curtis distance between two points ``u`` and ``v`` is


       .. math::

            d(u,v) = \\frac{\\sum_i (|u_i-v_i|)}
                          {\\sum_i (|u_i+v_i|)}

    13. ``Y = cdist(XA, XB, 'mahalanobis', VI=None)``

       Computes the Mahalanobis distance between the points. The
       Mahalanobis distance between two points ``u`` and ``v`` is
       :math:`\\sqrt{(u-v)(1/V)(u-v)^T}` where :math:`(1/V)` (the ``VI``
       variable) is the inverse covariance. If ``VI`` is not None,
       ``VI`` will be used as the inverse covariance matrix.

    14. ``Y = cdist(XA, XB, 'yule')``

       Computes the Yule distance between the boolean
       vectors. (see `yule` function documentation)

    15. ``Y = cdist(XA, XB, 'matching')``

       Synonym for 'hamming'.

    16. ``Y = cdist(XA, XB, 'dice')``

       Computes the Dice distance between the boolean vectors. (see
       `dice` function documentation)

    17. ``Y = cdist(XA, XB, 'kulsinski')``

       Computes the Kulsinski distance between the boolean
       vectors. (see `kulsinski` function documentation)

    18. ``Y = cdist(XA, XB, 'rogerstanimoto')``

       Computes the Rogers-Tanimoto distance between the boolean
       vectors. (see `rogerstanimoto` function documentation)

    19. ``Y = cdist(XA, XB, 'russellrao')``

       Computes the Russell-Rao distance between the boolean
       vectors. (see `russellrao` function documentation)

    20. ``Y = cdist(XA, XB, 'sokalmichener')``

       Computes the Sokal-Michener distance between the boolean
       vectors. (see `sokalmichener` function documentation)

    21. ``Y = cdist(XA, XB, 'sokalsneath')``

       Computes the Sokal-Sneath distance between the vectors. (see
       `sokalsneath` function documentation)


    22. ``Y = cdist(XA, XB, 'wminkowski', p=2., w=w)``

       Computes the weighted Minkowski distance between the
       vectors. (see `wminkowski` function documentation)

    23. ``Y = cdist(XA, XB, f)``

       Computes the distance between all pairs of vectors in X
       using the user supplied 2-arity function f. For example,
       Euclidean distance between the vectors could be computed
       as follows::

         dm = cdist(XA, XB, lambda u, v: np.sqrt(((u-v)**2).sum()))

       Note that you should avoid passing a reference to one of
       the distance functions defined in this library. For example,::

         dm = cdist(XA, XB, sokalsneath)

       would calculate the pair-wise distances between the vectors in
       X using the Python function `sokalsneath`. This would result in
       sokalsneath being called :math:`{n \\choose 2}` times, which
       is inefficient. Instead, the optimized C version is more
       efficient, and we call it using the following syntax::

         dm = cdist(XA, XB, 'sokalsneath')

    Examples
    --------
    Find the Euclidean distances between four 2-D coordinates:

    >>> from scipy.spatial import distance
    >>> coords = [(35.0456, -85.2672),
    ...           (35.1174, -89.9711),
    ...           (35.9728, -83.9422),
    ...           (36.1667, -86.7833)]
    >>> distance.cdist(coords, coords, 'euclidean')
    array([[ 0.    ,  4.7044,  1.6172,  1.8856],
           [ 4.7044,  0.    ,  6.0893,  3.3561],
           [ 1.6172,  6.0893,  0.    ,  2.8477],
           [ 1.8856,  3.3561,  2.8477,  0.    ]])


    Find the Manhattan distance from a 3-D point to the corners of the unit
    cube:

    >>> a = np.array([[0, 0, 0],
    ...               [0, 0, 1],
    ...               [0, 1, 0],
    ...               [0, 1, 1],
    ...               [1, 0, 0],
    ...               [1, 0, 1],
    ...               [1, 1, 0],
    ...               [1, 1, 1]])
    >>> b = np.array([[ 0.1,  0.2,  0.4]])
    >>> distance.cdist(a, b, 'cityblock')
    array([[ 0.7],
           [ 0.9],
           [ 1.3],
           [ 1.5],
           [ 1.5],
           [ 1.7],
           [ 2.1],
           [ 2.3]])

    """
    # You can also call this as:
    #     Y = cdist(XA, XB, 'test_abc')
    # where 'abc' is the metric being tested.  This computes the distance
    # between all pairs of vectors in XA and XB using the distance metric 'abc'
    # but with a more succinct, verifiable, but less efficient implementation.

    kwargs = _args_to_kwargs_xdist(args, kwargs, metric, "cdist")

    XA = np.asarray(XA, order='c')
    XB = np.asarray(XB, order='c')

    s = XA.shape
    sB = XB.shape

    if len(s) != 2:
        raise ValueError('XA must be a 2-dimensional array.')
    if len(sB) != 2:
        raise ValueError('XB must be a 2-dimensional array.')
    if s[1] != sB[1]:
        raise ValueError('XA and XB must have the same number of columns '
                         '(i.e. feature dimension.)')

    mA = s[0]
    mB = sB[0]
    n = s[1]
    out = kwargs.pop("out", None)
    if out is None:
        dm = np.empty((mA, mB), dtype=np.double)
    else:
        if out.shape != (mA, mB):
            raise ValueError("Output array has incorrect shape.")
        if not out.flags.c_contiguous:
            raise ValueError("Output array must be C-contiguous.")
        if out.dtype != np.double:
            raise ValueError("Output array must be double type.")
        dm = out

    # compute blacklist for deprecated kwargs
    if(metric in _METRICS['minkowski'].aka or
       metric in _METRICS['wminkowski'].aka or
       metric in ['test_minkowski', 'test_wminkowski'] or
       metric in [minkowski, wminkowski]):
        kwargs_blacklist = ["V", "VI"]
    elif(metric in _METRICS['seuclidean'].aka or
         metric == 'test_seuclidean' or metric == seuclidean):
        kwargs_blacklist = ["p", "w", "VI"]
    elif(metric in _METRICS['mahalanobis'].aka or
         metric == 'test_mahalanobis' or metric == mahalanobis):
        kwargs_blacklist = ["p", "w", "V"]
    else:
        kwargs_blacklist = ["p", "V", "VI"]

    _filter_deprecated_kwargs(kwargs, kwargs_blacklist)

    if callable(metric):

        mstr = getattr(metric, '__name__', 'Unknown')
        metric_name = _METRIC_ALIAS.get(mstr, None)

        XA, XB, typ, kwargs = _validate_cdist_input(XA, XB, mA, mB, n,
                                                    metric_name, **kwargs)

        for i in xrange(0, mA):
            for j in xrange(0, mB):
                dm[i, j] = metric(XA[i], XB[j], **kwargs)

    elif isinstance(metric, string_types):
        mstr = metric.lower()

        mstr, kwargs = _select_weighted_metric(mstr, kwargs, out)

        metric_name = _METRIC_ALIAS.get(mstr, None)
        if metric_name is not None:
            XA, XB, typ, kwargs = _validate_cdist_input(XA, XB, mA, mB, n,
                                                        metric_name, **kwargs)
            # get cdist wrapper
            cdist_fn = getattr(_distance_wrap,
                               "cdist_%s_%s_wrap" % (metric_name, typ))
            cdist_fn(XA, XB, dm, **kwargs)
            return dm

        elif mstr.startswith("test_"):
            if mstr in _TEST_METRICS:
                dm = cdist(XA, XB, _TEST_METRICS[mstr], **kwargs)
            else:
                raise ValueError('Unknown "Test" Distance Metric: %s' % mstr[5:])
        else:
            raise ValueError('Unknown Distance Metric: %s' % mstr)
    else:
        raise TypeError('2nd argument metric must be a string identifier '
                        'or a function.')
    return dm
