# cython: language_level=3
# cython: annotate=False
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: profile=False
# cython: linetrace=False
# cython: initializedcheck=False
# cython: binding=False
# distutils: define_macros=CYTHON_TRACE_NOGIL=0


# By Jake Vanderplas (2013) <jakevdp@cs.washington.edu>
# written for the scikit-learn project
# License: BSD

import numpy as np
cimport numpy as np
from cython cimport final

np.import_array()  # required in order to use C-API


# First, define a function to get an ndarray from a memory bufffer
cdef extern from "arrayobject.h":
    object PyArray_SimpleNewFromData(int nd, np.npy_intp* dims,
                                     int typenum, void* data)


cdef inline np.ndarray _buffer_to_ndarray(const DTYPE_t* x, np.npy_intp n):
    # Wrap a memory buffer with an ndarray. Warning: this is not robust.
    # In particular, if x is deallocated before the returned array goes
    # out of scope, this could cause memory errors.  Since there is not
    # a possibility of this for our use-case, this should be safe.

    # Note: this Segfaults unless np.import_array() is called above
    return PyArray_SimpleNewFromData(1, &n, DTYPECODE, <void*>x)


from libc.math cimport fabs, sqrt, exp, pow, cos, sin, asin
cdef DTYPE_t INF = np.inf

from scipy.sparse import csr_matrix, issparse

from ..utils._typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t, DTYPECODE
from ..utils._typedefs import DTYPE, ITYPE
from ..utils import check_array

######################################################################
# newObj function
#  this is a helper function for pickling
def newObj(obj):
    return obj.__new__(obj)


######################################################################
# metric mappings
#  These map from metric id strings to class names
METRIC_MAPPING = {'euclidean': EuclideanDistance,
                  'l2': EuclideanDistance,
                  'minkowski': MinkowskiDistance,
                  'p': MinkowskiDistance,
                  'manhattan': ManhattanDistance,
                  'cityblock': ManhattanDistance,
                  'l1': ManhattanDistance,
                  'chebyshev': ChebyshevDistance,
                  'infinity': ChebyshevDistance,
                  'seuclidean': SEuclideanDistance,
                  'mahalanobis': MahalanobisDistance,
                  'wminkowski': WMinkowskiDistance,
                  'hamming': HammingDistance,
                  'canberra': CanberraDistance,
                  'braycurtis': BrayCurtisDistance,
                  'matching': MatchingDistance,
                  'jaccard': JaccardDistance,
                  'dice': DiceDistance,
                  'kulsinski': KulsinskiDistance,
                  'rogerstanimoto': RogersTanimotoDistance,
                  'russellrao': RussellRaoDistance,
                  'sokalmichener': SokalMichenerDistance,
                  'sokalsneath': SokalSneathDistance,
                  'haversine': HaversineDistance,
                  'pyfunc': PyFuncDistance}

BOOL_METRICS = [
    "matching",
    "jaccard",
    "dice",
    "kulsinski",
    "rogerstanimoto",
    "russellrao",
    "sokalmichener",
    "sokalsneath",
]

def get_valid_metric_ids(L):
    """Given an iterable of metric class names or class identifiers,
    return a list of metric IDs which map to those classes.

    Example:
    >>> L = get_valid_metric_ids([EuclideanDistance, 'ManhattanDistance'])
    >>> sorted(L)
    ['cityblock', 'euclidean', 'l1', 'l2', 'manhattan']
    """
    return [key for (key, val) in METRIC_MAPPING.items()
            if (val.__name__ in L) or (val in L)]


######################################################################
# Distance Metric Classes
cdef class DistanceMetric:
    """DistanceMetric class

    This class provides a uniform interface to fast distance metric
    functions.  The various metrics can be accessed via the :meth:`get_metric`
    class method and the metric string identifier (see below).

    Examples
    --------
    >>> from sklearn.metrics import DistanceMetric
    >>> dist = DistanceMetric.get_metric('euclidean')
    >>> X = [[0, 1, 2],
             [3, 4, 5]]
    >>> dist.pairwise(X)
    array([[ 0.        ,  5.19615242],
           [ 5.19615242,  0.        ]])

    Available Metrics

    The following lists the string metric identifiers and the associated
    distance metric classes:

    **Metrics intended for real-valued vector spaces:**

    ==============  ====================  ========  ===============================
    identifier      class name            args      distance function
    --------------  --------------------  --------  -------------------------------
    "euclidean"     EuclideanDistance     -         ``sqrt(sum((x - y)^2))``
    "manhattan"     ManhattanDistance     -         ``sum(|x - y|)``
    "chebyshev"     ChebyshevDistance     -         ``max(|x - y|)``
    "minkowski"     MinkowskiDistance     p         ``sum(|x - y|^p)^(1/p)``
    "wminkowski"    WMinkowskiDistance    p, w      ``sum(|w * (x - y)|^p)^(1/p)``
    "seuclidean"    SEuclideanDistance    V         ``sqrt(sum((x - y)^2 / V))``
    "mahalanobis"   MahalanobisDistance   V or VI   ``sqrt((x - y)' V^-1 (x - y))``
    ==============  ====================  ========  ===============================

    **Metrics intended for two-dimensional vector spaces:**  Note that the haversine
    distance metric requires data in the form of [latitude, longitude] and both
    inputs and outputs are in units of radians.

    ============  ==================  ===============================================================
    identifier    class name          distance function
    ------------  ------------------  ---------------------------------------------------------------
    "haversine"   HaversineDistance   ``2 arcsin(sqrt(sin^2(0.5*dx) + cos(x1)cos(x2)sin^2(0.5*dy)))``
    ============  ==================  ===============================================================


    **Metrics intended for integer-valued vector spaces:**  Though intended
    for integer-valued vectors, these are also valid metrics in the case of
    real-valued vectors.

    =============  ====================  ========================================
    identifier     class name            distance function
    -------------  --------------------  ----------------------------------------
    "hamming"      HammingDistance       ``N_unequal(x, y) / N_tot``
    "canberra"     CanberraDistance      ``sum(|x - y| / (|x| + |y|))``
    "braycurtis"   BrayCurtisDistance    ``sum(|x - y|) / (sum(|x|) + sum(|y|))``
    =============  ====================  ========================================

    **Metrics intended for boolean-valued vector spaces:**  Any nonzero entry
    is evaluated to "True".  In the listings below, the following
    abbreviations are used:

     - N  : number of dimensions
     - NTT : number of dims in which both values are True
     - NTF : number of dims in which the first value is True, second is False
     - NFT : number of dims in which the first value is False, second is True
     - NFF : number of dims in which both values are False
     - NNEQ : number of non-equal dimensions, NNEQ = NTF + NFT
     - NNZ : number of nonzero dimensions, NNZ = NTF + NFT + NTT

    =================  =======================  ===============================
    identifier         class name               distance function
    -----------------  -----------------------  -------------------------------
    "jaccard"          JaccardDistance          NNEQ / NNZ
    "matching"         MatchingDistance         NNEQ / N
    "dice"             DiceDistance             NNEQ / (NTT + NNZ)
    "kulsinski"        KulsinskiDistance        (NNEQ + N - NTT) / (NNEQ + N)
    "rogerstanimoto"   RogersTanimotoDistance   2 * NNEQ / (N + NNEQ)
    "russellrao"       RussellRaoDistance       (N - NTT) / N
    "sokalmichener"    SokalMichenerDistance    2 * NNEQ / (N + NNEQ)
    "sokalsneath"      SokalSneathDistance      NNEQ / (NNEQ + 0.5 * NTT)
    =================  =======================  ===============================

    **User-defined distance:**

    ===========    ===============    =======
    identifier     class name         args
    -----------    ---------------    -------
    "pyfunc"       PyFuncDistance     func
    ===========    ===============    =======

    Here ``func`` is a function which takes two one-dimensional numpy
    arrays, and returns a distance.  Note that in order to be used within
    the BallTree, the distance must be a true metric:
    i.e. it must satisfy the following properties

    1) Non-negativity: d(x, y) >= 0
    2) Identity: d(x, y) = 0 if and only if x == y
    3) Symmetry: d(x, y) = d(y, x)
    4) Triangle Inequality: d(x, y) + d(y, z) >= d(x, z)

    Because of the Python object overhead involved in calling the python
    function, this will be fairly slow, but it will have the same
    scaling as other distances.
    """
    def __cinit__(self):
        self.p = 2
        self.vec = np.zeros(1, dtype=DTYPE, order='C')
        self.mat = np.zeros((1, 1), dtype=DTYPE, order='C')
        self.size = 1

    def __reduce__(self):
        """
        reduce method used for pickling
        """
        return (newObj, (self.__class__,), self.__getstate__())

    def __getstate__(self):
        """
        get state for pickling
        """
        if self.__class__.__name__ == "PyFuncDistance":
            return (float(self.p), np.asarray(self.vec), np.asarray(self.mat), self.func, self.kwargs)
        return (float(self.p), np.asarray(self.vec), np.asarray(self.mat))

    def __setstate__(self, state):
        """
        set state for pickling
        """
        self.p = state[0]
        self.vec = state[1]
        self.mat = state[2]
        if self.__class__.__name__ == "PyFuncDistance":
            self.func = state[3]
            self.kwargs = state[4]
        self.size = self.vec.shape[0]

    @classmethod
    def get_metric(cls, metric, **kwargs):
        """Get the given distance metric from the string identifier.

        See the docstring of DistanceMetric for a list of available metrics.

        Parameters
        ----------
        metric : string or class name
            The distance metric to use
        **kwargs
            additional arguments will be passed to the requested metric
        """
        if isinstance(metric, DistanceMetric):
            return metric

        if callable(metric):
            return PyFuncDistance(metric, **kwargs)

        # Map the metric string ID to the metric class
        if isinstance(metric, type) and issubclass(metric, DistanceMetric):
            pass
        else:
            try:
                metric = METRIC_MAPPING[metric]
            except:
                raise ValueError("Unrecognized metric '%s'" % metric)

        # In Minkowski special cases, return more efficient methods
        if metric is MinkowskiDistance:
            p = kwargs.pop('p', 2)
            if p == 1:
                return ManhattanDistance(**kwargs)
            elif p == 2:
                return EuclideanDistance(**kwargs)
            elif np.isinf(p):
                return ChebyshevDistance(**kwargs)
            else:
                return MinkowskiDistance(p, **kwargs)
        else:
            return metric(**kwargs)

    def __init__(self):
        if self.__class__ is DistanceMetric:
            raise NotImplementedError("DistanceMetric is an abstract class")

    def _validate_data(self, X):
        """Validate the input data.

        This should be overridden in a base class if a specific input format
        is required.
        """
        return

    cdef DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                      ITYPE_t size) nogil except -1:
        """Compute the distance between vectors x1 and x2

        This should be overridden in a base class.
        """
        return -999

    cdef DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                       ITYPE_t size) nogil except -1:
        """Compute the reduced distance between vectors x1 and x2.

        This can optionally be overridden in a base class.

        The reduced distance is any measure that yields the same rank as the
        distance, but is more efficient to compute.  For example, for the
        Euclidean metric, the reduced distance is the squared-euclidean
        distance.
        """
        return self.dist(x1, x2, size)

    cdef DTYPE_t sparse_dist(
        self,
        const DTYPE_t[:] x1_data,
        const ITYPE_t[:] x1_indices,
        const DTYPE_t[:] x2_data,
        const ITYPE_t[:] x2_indices,
    ) nogil except -1:
        """Compute the reduced distance between vectors x1 and x2
        given non null coordinates and their corresponding indices.

        This should be overridden in a base class.
        """
        return -999

    cdef DTYPE_t sparse_rdist(
        self,
        const DTYPE_t[:] x1_data,
        const ITYPE_t[:] x1_indices,
        const DTYPE_t[:] x2_data,
        const ITYPE_t[:] x2_indices,
    ) nogil except -1:
        """Compute the reduced distance between vectors x1 and x2
        given non null coordinates and their corresponding indices.

        This can optionally be overridden in a base class.

        The reduced distance is any measure that yields the same rank as the
        distance, but is more efficient to compute.  For example, for the
        Euclidean metric, the reduced distance is the squared-euclidean
        distance.
        """
        return self.sparse_dist(x1_data, x1_indices, x2_data, x2_indices)

    cdef int pdist(self, const DTYPE_t[:, ::1] X, DTYPE_t[:, ::1] D) except -1:
        """compute the pairwise distances between points in X"""
        cdef ITYPE_t i1, i2
        for i1 in range(X.shape[0]):
            for i2 in range(i1, X.shape[0]):
                D[i1, i2] = self.dist(&X[i1, 0], &X[i2, 0], X.shape[1])
                D[i2, i1] = D[i1, i2]
        return 0

    cdef int cdist(self, const DTYPE_t[:, ::1] X, const DTYPE_t[:, ::1] Y,
                   DTYPE_t[:, ::1] D) except -1:
        """compute the cross-pairwise distances between arrays X and Y"""
        cdef ITYPE_t i1, i2
        if X.shape[1] != Y.shape[1]:
            raise ValueError('X and Y must have the same second dimension')
        for i1 in range(X.shape[0]):
            for i2 in range(Y.shape[0]):
                D[i1, i2] = self.dist(&X[i1, 0], &Y[i2, 0], X.shape[1])
        return 0

    cdef DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        """Convert the reduced distance to the distance"""
        return rdist

    cdef DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        """Convert the distance to the reduced distance"""
        return dist

    def rdist_to_dist(self, rdist):
        """Convert the Reduced distance to the true distance.

        The reduced distance, defined for some metrics, is a computationally
        more efficient measure which preserves the rank of the true distance.
        For example, in the Euclidean distance metric, the reduced distance
        is the squared-euclidean distance.

        Parameters
        ----------
        rdist : double
            Reduced distance.

        Returns
        -------
        double
            True distance.
        """
        return rdist

    def dist_to_rdist(self, dist):
        """Convert the true distance to the reduced distance.

        The reduced distance, defined for some metrics, is a computationally
        more efficient measure which preserves the rank of the true distance.
        For example, in the Euclidean distance metric, the reduced distance
        is the squared-euclidean distance.

        Parameters
        ----------
        dist : double
            True distance.

        Returns
        -------
        double
            Reduced distance.
        """
        return dist

    def pairwise(self, X, Y=None):
        """Compute the pairwise distances between X and Y

        This is a convenience routine for the sake of testing.  For many
        metrics, the utilities in scipy.spatial.distance.cdist and
        scipy.spatial.distance.pdist will be faster.

        Parameters
        ----------
        X : array-like
            Array of shape (Nx, D), representing Nx points in D dimensions.
        Y : array-like (optional)
            Array of shape (Ny, D), representing Ny points in D dimensions.
            If not specified, then Y=X.

        Returns
        -------
        dist : ndarray
            The shape (Nx, Ny) array of pairwise distances between points in
            X and Y.
        """
        cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] Xarr
        cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] Yarr
        cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] Darr

        Xarr = np.asarray(X, dtype=DTYPE, order='C')
        self._validate_data(Xarr)
        if Y is None:
            Darr = np.zeros((Xarr.shape[0], Xarr.shape[0]),
                         dtype=DTYPE, order='C')
            self.pdist(Xarr, Darr)
        else:
            Yarr = np.asarray(Y, dtype=DTYPE, order='C')
            self._validate_data(Yarr)
            Darr = np.zeros((Xarr.shape[0], Yarr.shape[0]),
                         dtype=DTYPE, order='C')
            self.cdist(Xarr, Yarr, Darr)
        return Darr


#------------------------------------------------------------
# Euclidean Distance
#  d = sqrt(sum(x_i^2 - y_i^2))
cdef class EuclideanDistance(DistanceMetric):
    r"""Euclidean Distance metric

    .. math::
       D(x, y) = \sqrt{ \sum_i (x_i - y_i) ^ 2 }
    """
    def __init__(self):
        self.p = 2

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return euclidean_dist(x1, x2, size)

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        return euclidean_rdist(x1, x2, size)

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return sqrt(rdist)

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        return dist * dist

    def rdist_to_dist(self, rdist):
        return np.sqrt(rdist)

    def dist_to_rdist(self, dist):
        return dist ** 2


#------------------------------------------------------------
# SEuclidean Distance
#  d = sqrt(sum((x_i - y_i2)^2 / v_i))
cdef class SEuclideanDistance(DistanceMetric):
    r"""Standardized Euclidean Distance metric

    .. math::
       D(x, y) = \sqrt{ \sum_i \frac{ (x_i - y_i) ^ 2}{V_i} }
    """
    def __init__(self, V):
        self.vec = np.asarray(V, dtype=DTYPE)
        self.size = self.vec.shape[0]
        self.p = 2

    def _validate_data(self, X):
        if X.shape[1] != self.size:
            raise ValueError('SEuclidean dist: size of V does not match')

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        cdef DTYPE_t tmp, d=0
        cdef np.intp_t j
        for j in range(size):
            tmp = x1[j] - x2[j]
            d += tmp * tmp / self.vec[j]
        return d

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return sqrt(self.rdist(x1, x2, size))

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return sqrt(rdist)

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        return dist * dist

    def rdist_to_dist(self, rdist):
        return np.sqrt(rdist)

    def dist_to_rdist(self, dist):
        return dist ** 2


#------------------------------------------------------------
# Manhattan Distance
#  d = sum(abs(x_i - y_i))
cdef class ManhattanDistance(DistanceMetric):
    r"""Manhattan/City-block Distance metric

    .. math::
       D(x, y) = \sum_i |x_i - y_i|
    """
    def __init__(self):
        self.p = 1

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef DTYPE_t d = 0
        cdef np.intp_t j
        for j in range(size):
            d += fabs(x1[j] - x2[j])
        return d


#------------------------------------------------------------
# Chebyshev Distance
#  d = max_i(abs(x_i - y_i))
cdef class ChebyshevDistance(DistanceMetric):
    """Chebyshev/Infinity Distance

    .. math::
       D(x, y) = max_i (|x_i - y_i|)

    Examples
    --------
    >>> from sklearn.metrics import DistanceMetric
    >>> dist = DistanceMetric.get_metric('chebyshev')
    >>> X = [[0, 1, 2],
    ...      [3, 4, 5]]
    >>> Y = [[-1, 0, 1],
    ...      [3, 4, 5]]
    >>> dist.pairwise(X, Y)
    array([[1.732..., 5.196...],
           [6.928..., 0....   ]])
    """
    def __init__(self):
        self.p = INF

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef DTYPE_t d = 0
        cdef np.intp_t j
        for j in range(size):
            d = fmax(d, fabs(x1[j] - x2[j]))
        return d


#------------------------------------------------------------
# Minkowski Distance
#  d = sum(x_i^p - y_i^p) ^ (1/p)
cdef class MinkowskiDistance(DistanceMetric):
    r"""Minkowski Distance

    .. math::
       D(x, y) = [\sum_i (x_i - y_i)^p] ^ (1/p)

    Minkowski Distance requires p >= 1 and finite. For p = infinity,
    use ChebyshevDistance.
    Note that for p=1, ManhattanDistance is more efficient, and for
    p=2, EuclideanDistance is more efficient.
    """
    def __init__(self, p):
        if p < 1:
            raise ValueError("p must be greater than 1")
        elif np.isinf(p):
            raise ValueError("MinkowskiDistance requires finite p. "
                             "For p=inf, use ChebyshevDistance.")
        self.p = p

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        cdef DTYPE_t d=0
        cdef np.intp_t j
        for j in range(size):
            d += pow(fabs(x1[j] - x2[j]), self.p)
        return d

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return pow(self.rdist(x1, x2, size), 1. / self.p)

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return pow(rdist, 1. / self.p)

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        return pow(dist, self.p)

    def rdist_to_dist(self, rdist):
        return rdist ** (1. / self.p)

    def dist_to_rdist(self, dist):
        return dist ** self.p


#------------------------------------------------------------
# W-Minkowski Distance
#  d = sum(w_i^p * (x_i^p - y_i^p)) ^ (1/p)
cdef class WMinkowskiDistance(DistanceMetric):
    r"""Weighted Minkowski Distance

    .. math::
       D(x, y) = [\sum_i |w_i * (x_i - y_i)|^p] ^ (1/p)

    Weighted Minkowski Distance requires p >= 1 and finite.

    Parameters
    ----------
    p : int
        The order of the norm of the difference :math:`{||u-v||}_p`.
    w : (N,) array-like
        The weight vector.

    """
    def __init__(self, p, w):
        if p < 1:
            raise ValueError("p must be greater than 1")
        elif np.isinf(p):
            raise ValueError("WMinkowskiDistance requires finite p. "
                             "For p=inf, use ChebyshevDistance.")
        self.p = p
        self.vec = np.asarray(w, dtype=DTYPE)
        self.size = self.vec.shape[0]

    def _validate_data(self, X):
        if X.shape[1] != self.size:
            raise ValueError('WMinkowskiDistance dist: '
                             'size of w does not match')

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        cdef DTYPE_t d=0
        cdef np.intp_t j
        for j in range(size):
            d += pow(self.vec[j] * fabs(x1[j] - x2[j]), self.p)
        return d

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return pow(self.rdist(x1, x2, size), 1. / self.p)

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return pow(rdist, 1. / self.p)

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        return pow(dist, self.p)

    def rdist_to_dist(self, rdist):
        return rdist ** (1. / self.p)

    def dist_to_rdist(self, dist):
        return dist ** self.p


#------------------------------------------------------------
# Mahalanobis Distance
#  d = sqrt( (x - y)^T V^-1 (x - y) )
cdef class MahalanobisDistance(DistanceMetric):
    """Mahalanobis Distance

    .. math::
       D(x, y) = \sqrt{ (x - y)^T V^{-1} (x - y) }

    Parameters
    ----------
    V : array-like
        Symmetric positive-definite covariance matrix.
        The inverse of this matrix will be explicitly computed.
    VI : array-like
        optionally specify the inverse directly.  If VI is passed,
        then V is not referenced.
    """
    def __init__(self, V=None, VI=None):
        if VI is None:
            if V is None:
                raise ValueError("Must provide either V or VI "
                                 "for Mahalanobis distance")
            VI = np.linalg.inv(V)
        if VI.ndim != 2 or VI.shape[0] != VI.shape[1]:
            raise ValueError("V/VI must be square")

        self.mat = np.asarray(VI, dtype=float, order='C')

        self.size = self.mat.shape[0]

        # we need vec as a work buffer
        self.vec = np.zeros(self.size, dtype=DTYPE)

    def _validate_data(self, X):
        if X.shape[1] != self.size:
            raise ValueError('Mahalanobis dist: size of V does not match')

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        cdef DTYPE_t tmp, d = 0
        cdef np.intp_t i, j

        # compute (x1 - x2).T * VI * (x1 - x2)
        for i in range(size):
            self.vec[i] = x1[i] - x2[i]

        for i in range(size):
            tmp = 0
            for j in range(size):
                tmp += self.mat[i, j] * self.vec[j]
            d += tmp * self.vec[i]
        return d

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return sqrt(self.rdist(x1, x2, size))

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return sqrt(rdist)

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        return dist * dist

    def rdist_to_dist(self, rdist):
        return np.sqrt(rdist)

    def dist_to_rdist(self, dist):
        return dist ** 2


#------------------------------------------------------------
# Hamming Distance
#  d = N_unequal(x, y) / N_tot
cdef class HammingDistance(DistanceMetric):
    r"""Hamming Distance

    Hamming distance is meant for discrete-valued vectors, though it is
    a valid metric for real-valued vectors.

    .. math::
       D(x, y) = \frac{1}{N} \sum_i \delta_{x_i, y_i}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int n_unequal = 0
        cdef np.intp_t j
        for j in range(size):
            if x1[j] != x2[j]:
                n_unequal += 1
        return float(n_unequal) / size


#------------------------------------------------------------
# Canberra Distance
#  D(x, y) = sum[ abs(x_i - y_i) / (abs(x_i) + abs(y_i)) ]
cdef class CanberraDistance(DistanceMetric):
    r"""Canberra Distance

    Canberra distance is meant for discrete-valued vectors, though it is
    a valid metric for real-valued vectors.

    .. math::
       D(x, y) = \sum_i \frac{|x_i - y_i|}{|x_i| + |y_i|}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef DTYPE_t denom, d = 0
        cdef np.intp_t j
        for j in range(size):
            denom = fabs(x1[j]) + fabs(x2[j])
            if denom > 0:
                d += fabs(x1[j] - x2[j]) / denom
        return d


#------------------------------------------------------------
# Bray-Curtis Distance
#  D(x, y) = sum[abs(x_i - y_i)] / sum[abs(x_i) + abs(y_i)]
cdef class BrayCurtisDistance(DistanceMetric):
    r"""Bray-Curtis Distance

    Bray-Curtis distance is meant for discrete-valued vectors, though it is
    a valid metric for real-valued vectors.

    .. math::
       D(x, y) = \frac{\sum_i |x_i - y_i|}{\sum_i(|x_i| + |y_i|)}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef DTYPE_t num = 0, denom = 0
        cdef np.intp_t j
        for j in range(size):
            num += fabs(x1[j] - x2[j])
            denom += fabs(x1[j]) + fabs(x2[j])
        if denom > 0:
            return num / denom
        else:
            return 0.0


#------------------------------------------------------------
# Jaccard Distance (boolean)
#  D(x, y) = N_unequal(x, y) / N_nonzero(x, y)
cdef class JaccardDistance(DistanceMetric):
    r"""Jaccard Distance

    Jaccard Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{N_{TF} + N_{FT}}{N_{TT} + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, n_eq = 0, nnz = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            nnz += (tf1 or tf2)
            n_eq += (tf1 and tf2)
        # Based on https://github.com/scipy/scipy/pull/7373
        # When comparing two all-zero vectors, scipy>=1.2.0 jaccard metric
        # was changed to return 0, instead of nan.
        if nnz == 0:
            return 0
        return (nnz - n_eq) * 1.0 / nnz


#------------------------------------------------------------
# Matching Distance (boolean)
#  D(x, y) = n_neq / n
cdef class MatchingDistance(DistanceMetric):
    r"""Matching Distance

    Matching Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{N_{TF} + N_{FT}}{N}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, n_neq = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            n_neq += (tf1 != tf2)
        return n_neq * 1. / size


#------------------------------------------------------------
# Dice Distance (boolean)
#  D(x, y) = n_neq / (2 * ntt + n_neq)
cdef class DiceDistance(DistanceMetric):
    r"""Dice Distance

    Dice Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{N_{TF} + N_{FT}}{2 * N_{TT} + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, n_neq = 0, ntt = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            ntt += (tf1 and tf2)
            n_neq += (tf1 != tf2)
        return n_neq / (2.0 * ntt + n_neq)


#------------------------------------------------------------
# Kulsinski Distance (boolean)
#  D(x, y) = (ntf + nft - ntt + n) / (n_neq + n)
cdef class KulsinskiDistance(DistanceMetric):
    r"""Kulsinski Distance

    Kulsinski Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = 1 - \frac{N_{TT}}{N + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, ntt = 0, n_neq = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            n_neq += (tf1 != tf2)
            ntt += (tf1 and tf2)
        return (n_neq - ntt + size) * 1.0 / (n_neq + size)


#------------------------------------------------------------
# Rogers-Tanimoto Distance (boolean)
#  D(x, y) = 2 * n_neq / (n + n_neq)
cdef class RogersTanimotoDistance(DistanceMetric):
    r"""Rogers-Tanimoto Distance

    Rogers-Tanimoto Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{2 (N_{TF} + N_{FT})}{N + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, n_neq = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            n_neq += (tf1 != tf2)
        return (2.0 * n_neq) / (size + n_neq)


#------------------------------------------------------------
# Russell-Rao Distance (boolean)
#  D(x, y) = (n - ntt) / n
cdef class RussellRaoDistance(DistanceMetric):
    r"""Russell-Rao Distance

    Russell-Rao Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{N - N_{TT}}{N}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, ntt = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            ntt += (tf1 and tf2)
        return (size - ntt) * 1. / size


#------------------------------------------------------------
# Sokal-Michener Distance (boolean)
#  D(x, y) = 2 * n_neq / (n + n_neq)
cdef class SokalMichenerDistance(DistanceMetric):
    r"""Sokal-Michener Distance

    Sokal-Michener Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{2 (N_{TF} + N_{FT})}{N + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, n_neq = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            n_neq += (tf1 != tf2)
        return (2.0 * n_neq) / (size + n_neq)


#------------------------------------------------------------
# Sokal-Sneath Distance (boolean)
#  D(x, y) = n_neq / (0.5 * n_tt + n_neq)
cdef class SokalSneathDistance(DistanceMetric):
    r"""Sokal-Sneath Distance

    Sokal-Sneath Distance is a dissimilarity measure for boolean-valued
    vectors. All nonzero entries will be treated as True, zero entries will
    be treated as False.

    .. math::
       D(x, y) = \frac{N_{TF} + N_{FT}}{N_{TT} / 2 + N_{TF} + N_{FT}}
    """
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        cdef int tf1, tf2, ntt = 0, n_neq = 0
        cdef np.intp_t j
        for j in range(size):
            tf1 = x1[j] != 0
            tf2 = x2[j] != 0
            n_neq += (tf1 != tf2)
            ntt += (tf1 and tf2)
        return n_neq / (0.5 * ntt + n_neq)


#------------------------------------------------------------
# Haversine Distance (2 dimensional)
#  D(x, y) = 2 arcsin{sqrt[sin^2 ((x1 - y1) / 2)
#                          + cos(x1) cos(y1) sin^2 ((x2 - y2) / 2)]}
cdef class HaversineDistance(DistanceMetric):
    """Haversine (Spherical) Distance

    The Haversine distance is the angular distance between two points on
    the surface of a sphere.  The first distance of each point is assumed
    to be the latitude, the second is the longitude, given in radians.
    The dimension of the points must be 2:

    .. math::
       D(x, y) = 2\\arcsin[\\sqrt{\\sin^2((x1 - y1) / 2)
                                + \\cos(x1)\\cos(y1)\\sin^2((x2 - y2) / 2)}]
    """

    def _validate_data(self, X):
        if X.shape[1] != 2:
            raise ValueError("Haversine distance only valid "
                             "in 2 dimensions")

    cdef inline DTYPE_t rdist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) nogil except -1:
        cdef DTYPE_t sin_0 = sin(0.5 * (x1[0] - x2[0]))
        cdef DTYPE_t sin_1 = sin(0.5 * (x1[1] - x2[1]))
        return (sin_0 * sin_0 + cos(x1[0]) * cos(x2[0]) * sin_1 * sin_1)

    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return 2 * asin(sqrt(self.rdist(x1, x2, size)))

    cdef inline DTYPE_t _rdist_to_dist(self, DTYPE_t rdist) nogil except -1:
        return 2 * asin(sqrt(rdist))

    cdef inline DTYPE_t _dist_to_rdist(self, DTYPE_t dist) nogil except -1:
        cdef DTYPE_t tmp = sin(0.5 * dist)
        return tmp * tmp

    def rdist_to_dist(self, rdist):
        return 2 * np.arcsin(np.sqrt(rdist))

    def dist_to_rdist(self, dist):
        tmp = np.sin(0.5 * dist)
        return tmp * tmp


#------------------------------------------------------------
# Yule Distance (boolean)
#  D(x, y) = 2 * ntf * nft / (ntt * nff + ntf * nft)
# [This is not a true metric, so we will leave it out.]
#
#cdef class YuleDistance(DistanceMetric):
#    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
#                             ITYPE_t size):
#        cdef int tf1, tf2, ntf = 0, nft = 0, ntt = 0, nff = 0
#        cdef np.intp_t j
#        for j in range(size):
#            tf1 = x1[j] != 0
#            tf2 = x2[j] != 0
#            ntt += tf1 and tf2
#            ntf += tf1 and (tf2 == 0)
#            nft += (tf1 == 0) and tf2
#        nff = size - ntt - ntf - nft
#        return (2.0 * ntf * nft) / (ntt * nff + ntf * nft)


#------------------------------------------------------------
# Cosine Distance
#  D(x, y) = dot(x, y) / (|x| * |y|)
# [This is not a true metric, so we will leave it out.]
#
#cdef class CosineDistance(DistanceMetric):
#    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
#                             ITYPE_t size):
#        cdef DTYPE_t d = 0, norm1 = 0, norm2 = 0
#        cdef np.intp_t j
#        for j in range(size):
#            d += x1[j] * x2[j]
#            norm1 += x1[j] * x1[j]
#            norm2 += x2[j] * x2[j]
#        return 1.0 - d / sqrt(norm1 * norm2)


#------------------------------------------------------------
# Correlation Distance
#  D(x, y) = dot((x - mx), (y - my)) / (|x - mx| * |y - my|)
# [This is not a true metric, so we will leave it out.]
#
#cdef class CorrelationDistance(DistanceMetric):
#    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
#                             ITYPE_t size):
#        cdef DTYPE_t mu1 = 0, mu2 = 0, x1nrm = 0, x2nrm = 0, x1Tx2 = 0
#        cdef DTYPE_t tmp1, tmp2
#
#        cdef np.intp_t i
#        for i in range(size):
#            mu1 += x1[i]
#            mu2 += x2[i]
#        mu1 /= size
#        mu2 /= size
#
#        for i in range(size):
#            tmp1 = x1[i] - mu1
#            tmp2 = x2[i] - mu2
#            x1nrm += tmp1 * tmp1
#            x2nrm += tmp2 * tmp2
#            x1Tx2 += tmp1 * tmp2
#
#        return (1. - x1Tx2) / sqrt(x1nrm * x2nrm)


#------------------------------------------------------------
# User-defined distance
#
cdef class PyFuncDistance(DistanceMetric):
    """PyFunc Distance

    A user-defined distance

    Parameters
    ----------
    func : function
        func should take two numpy arrays as input, and return a distance.
    """
    def __init__(self, func, **kwargs):
        self.func = func
        self.kwargs = kwargs

    # in cython < 0.26, GIL was required to be acquired during definition of
    # the function and inside the body of the function. This behaviour is not
    # allowed in cython >= 0.26 since it is a redundant GIL acquisition. The
    # only way to be back compatible is to inherit `dist` from the base class
    # without GIL and called an inline `_dist` which acquire GIL.
    cdef inline DTYPE_t dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                             ITYPE_t size) nogil except -1:
        return self._dist(x1, x2, size)

    cdef inline DTYPE_t _dist(self, const DTYPE_t* x1, const DTYPE_t* x2,
                              ITYPE_t size) except -1 with gil:
        cdef np.ndarray x1arr
        cdef np.ndarray x2arr
        x1arr = _buffer_to_ndarray(x1, size)
        x2arr = _buffer_to_ndarray(x2, size)
        d = self.func(x1arr, x2arr, **self.kwargs)
        try:
            # Cython generates code here that results in a TypeError
            # if d is the wrong type.
            return d
        except TypeError:
            raise TypeError("Custom distance function must accept two "
                            "vectors and return a float.")


cdef inline double fmax(double a, double b) nogil:
    return max(a, b)


######################################################################
# Datasets Pair Classes
cdef class DatasetsPair:
    """Abstract class which wraps a pair of datasets (X, Y).

    This class allows computing distances between two vectors (X_i, Y_j)
    (rows of X and Y) at a time given the pair of their indices (i, j).

    X and Y can be stored as np.ndarrays or CSR matrices in subclasses.

    This class avoids the overhead of dispatching distance computations
    to :class:`sklearn.metrics.DistanceMetric` based on the physical
    representation of the vectors (sparse vs. dense). It makes use of
    cython.final to remove the overhead of method calls' dispatch.

    Parameters
    ----------
    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """

    # The `distance_metric` attribute is defined in _dist_metrics.pxd

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        str metric="euclidean",
        dict metric_kwargs=dict(),
    ) -> DatasetsPair:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : array-like of shape (n_X, d)
            Input data.

        Y : array-like of shape (n_Y, d)
            Input data.

        metric : str, default='euclidean'
            The distance metric to use for argkmin. The default metric is
            a fast implementation of the standard Euclidean metric.
            For a list of available metrics, see the documentation of
            :class:`~sklearn.metrics.DistanceMetric`.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        Returns
        -------
        datasets_pair: DatasetsPair
            The suited DatasetsPair implementation.
        """
        cdef:
            DistanceMetric distance_metric = DistanceMetric.get_metric(
                metric,
                **metric_kwargs
            )

        # check_array can be expensive, and we prefer to simply coerce from lists
        # to ndarrays eventually to get their dtype itemsize
        X = np.asarray(X) if isinstance(X, (tuple, list)) else X
        Y = np.asarray(Y) if isinstance(Y, (tuple, list)) else Y

        if X.dtype.itemsize != 8 or Y.dtype.itemsize != 8:
            raise ValueError("32bits datasets aren't supported for X and Y yet.")

        X = check_array(X, dtype=DTYPE, accept_sparse='csr')
        Y = check_array(Y, dtype=DTYPE, accept_sparse='csr')

        if X.shape[1] != Y.shape[1]:
            raise ValueError("Vectors of X and Y must have the same "
                              "number of dimensions but currently are "
                              f"respectively {X.shape[1]}-dimensional "
                              f"and {Y.shape[1]}-dimensional.")

        distance_metric._validate_data(X)
        distance_metric._validate_data(Y)

        if not issparse(X) and not issparse(Y):
            return DenseDenseDatasetsPair(X, Y, distance_metric)
        if issparse(X) and not issparse(Y):
            return SparseDenseDatasetsPair(X, Y, distance_metric)
        if not issparse(X) and issparse(Y):
            return DenseSparseDatasetsPair(X, Y, distance_metric)
        return SparseSparseDatasetsPair(X, Y, distance_metric)

    @classmethod
    def unpack_csr_matrix(cls, X: csr_matrix):
        """Ensure getting ITYPE instead of int internally used for CSR matrices."""
        # TODO: this adds another level of checks and conversion, could we remove it?
        X_data = check_array(X.data, dtype=DTYPE, ensure_2d=False)
        X_indices = check_array(X.indices, dtype=ITYPE, ensure_2d=False)
        X_indptr = check_array(X.indptr, dtype=ITYPE, ensure_2d=False)
        return X_data, X_indptr, X_indptr

    def __init__(self, DistanceMetric distance_metric):
        self.distance_metric = distance_metric

    cdef ITYPE_t n_X(self) nogil:
        return -999

    cdef ITYPE_t n_Y(self) nogil:
        return -999

    cdef DTYPE_t proxy_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.dist(i, j)

    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return -1

@final
cdef class DenseDenseDatasetsPair(DatasetsPair):
    """Compute distances between vectors of two arrays.

    Parameters
    ----------
    X: ndarray of shape (n_X, d)
        Rows represent vectors.

    Y: ndarray of shape (n_Y, d)
        Rows represent vectors.

    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """

    # The `X`, `Y` and `d` attributes are defined in _dist_metrics.pxd

    def __cinit__(self):
        # Initializing memory view to prevent memory errors and seg-faults
        # in rare cases where __init__ is not called
        self.X = np.empty((1, 1), dtype=DTYPE, order='C')
        self.Y = np.empty((1, 1), dtype=DTYPE, order='C')

    def __init__(self, X, Y, DistanceMetric distance_metric):
        DatasetsPair.__init__(self, distance_metric)
        self.X = check_array(X, dtype=DTYPE, order='C')
        self.Y = check_array(Y, dtype=DTYPE, order='C')
        self.d = X.shape[1]

    @final
    cdef ITYPE_t n_X(self) nogil:
        return self.X.shape[0]

    @final
    cdef ITYPE_t n_Y(self) nogil:
        return self.Y.shape[0]

    @final
    cdef DTYPE_t proxy_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.distance_metric.rdist(&self.X[i, 0],
                                          &self.Y[j, 0],
                                          self.d)

    @final
    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.distance_metric.dist(&self.X[i, 0],
                                         &self.Y[j, 0],
                                         self.d)

@final
cdef class SparseSparseDatasetsPair(DatasetsPair):
    """Compute distances between vectors of two sparse matrices.

    Parameters
    ----------
    X: sparse matrix of shape (n_X, d)
        Rows represent vectors.

    Y: sparse matrix of shape (n_X, d)
        Rows represent vectors.

    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """
    cdef:
        const DTYPE_t[:] X_data
        const ITYPE_t[:] X_indices,
        const ITYPE_t[:] X_indptr,

        const DTYPE_t[:] Y_data
        const ITYPE_t[:] Y_indices
        const ITYPE_t[:] Y_indptr

    def __cinit__(self):
        # Initializing memory view to prevent memory errors and seg-faults
        # in rare cases where __init__ is not called
        self.X_data = np.empty((1), dtype=DTYPE, order='C')
        self.X_indices = np.empty((1), dtype=ITYPE, order='C')
        self.X_indptr = np.empty((1), dtype=ITYPE, order='C')

        self.Y_data = np.empty((1), dtype=DTYPE, order='C')
        self.Y_indices = np.empty((1), dtype=ITYPE, order='C')
        self.Y_indptr = np.empty((1), dtype=ITYPE, order='C')

    def __init__(self, X, Y, DistanceMetric distance_metric):
        DatasetsPair.__init__(self, distance_metric)

        X = check_array(X, dtype=DTYPE, accept_sparse='csr')
        Y = check_array(Y, dtype=DTYPE, accept_sparse='csr')

        self.X_data, self.X_indices, self.X_indptr = self.unpack_csr_matrix(X)
        self.Y_data, self.Y_indices, self.Y_indptr = self.unpack_csr_matrix(Y)

    @final
    cdef ITYPE_t n_X(self) nogil:
        return self.X_indptr.shape[0] - 1

    @final
    cdef ITYPE_t n_Y(self) nogil:
        return self.Y_indptr.shape[0] -1

    @final
    cdef DTYPE_t proxy_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        cdef:
            ITYPE_t xi_start = self.X_indptr[i]
            ITYPE_t xi_end = self.X_indptr[i + 1]
            ITYPE_t yj_start = self.Y_indptr[j]
            ITYPE_t yj_end = self.Y_indptr[j + 1]

        return self.distance_metric.sparse_rdist(
            self.X_data[xi_start:xi_end],
            self.X_indices[xi_start:xi_end],
            self.Y_data[yj_start:yj_end],
            self.Y_indices[yj_start:yj_end],
        )

    @final
    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        cdef:
            ITYPE_t xi_start = self.X_indptr[i]
            ITYPE_t xi_end = self.X_indptr[i + 1]
            ITYPE_t yj_start = self.Y_indptr[j]
            ITYPE_t yj_end = self.Y_indptr[j + 1]

        return self.distance_metric.sparse_dist(
            self.X_data[xi_start:xi_end],
            self.X_indices[xi_start:xi_end],
            self.Y_data[yj_start:yj_end],
            self.Y_indices[yj_start:yj_end]
        )

@final
cdef class SparseDenseDatasetsPair(DatasetsPair):
    """Compute distances between vectors of a sparse matrix and a dense array.

    Parameters
    ----------
    X: sparse matrix of shape (n_X, d)
        Rows represent vectors.

    Y: ndarray of shape (n_Y, d)
        Rows represent vectors.

    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """
    cdef:
        const DTYPE_t[:] X_data
        const ITYPE_t[:] X_indices,
        const ITYPE_t[:] X_indptr,

        const DTYPE_t[:, ::1] Y  # shape: (n_Y, d)
        const ITYPE_t[:] Y_indices

    def __cinit__(self):
        # Initializing memory view to prevent memory errors and seg-faults
        # in rare cases where __init__ is not called
        self.X_data = np.empty((1), dtype=DTYPE, order='C')
        self.X_indices = np.empty((1), dtype=ITYPE, order='C')
        self.X_indptr = np.empty((1), dtype=ITYPE, order='C')

        self.Y = np.empty((1, 1), dtype=DTYPE, order='C')
        self.Y_indices = np.empty((1), dtype=ITYPE, order='C')


    def __init__(self, X, Y, DistanceMetric distance_metric):
        DatasetsPair.__init__(self, distance_metric)

        X = check_array(X, dtype=DTYPE, accept_sparse='csr')
        self.X_data, self.X_indices, self.X_indptr = self.unpack_csr_matrix(X)

        self.Y = check_array(Y, dtype=DTYPE)
        self.Y_indices = np.arange(self.Y.shape[1], dtype=ITYPE)

    @final
    cdef ITYPE_t n_X(self) nogil:
        return self.X_indptr.shape[0] - 1

    @final
    cdef ITYPE_t n_Y(self) nogil:
        return self.Y.shape[0]

    @final
    cdef DTYPE_t proxy_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        cdef:
            ITYPE_t xi_start = self.X_indptr[i]
            ITYPE_t xi_end = self.X_indptr[i + 1]

        # TODO: the 2D to 1D memory-view conversion might make computation slower, see:
        # https://github.com/scikit-learn/scikit-learn/issues/17299
        # Ideally, we could pass pointers and indices and access elements
        # then in distance_metric.dist
        return self.distance_metric.sparse_rdist(
            self.X_data[xi_start:xi_end],
            self.X_indices[xi_start:xi_end],
            self.Y[j, :],
            self.Y_indices
        )

    @final
    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        cdef:
            ITYPE_t xi_start = self.X_indptr[i]
            ITYPE_t xi_end = self.X_indptr[i + 1]

        # TODO: same as previous comment
        return self.distance_metric.sparse_dist(
            self.X_data[xi_start:xi_end],
            self.X_indices[xi_start:xi_end],
            self.Y[j, :],
            self.Y_indices
        )

@final
cdef class DenseSparseDatasetsPair(DatasetsPair):
    """Compute distances between vectors of a dense array and a sparse matrix.

    Parameters
    ----------
    X: ndarray of shape (n_X, d)
        Rows represent vectors.

    Y: sparse matrix of shape (n_Y, d)
        Rows represent vectors.

    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """
    cdef:
        # As distance metrics are symmetric functions, we can
        # simply rely on the other DatasetsPair and swap arguments.
        DatasetsPair datasets_pair

    def __init__(self, X, Y, DistanceMetric distance_metric):
        DatasetsPair.__init__(self, distance_metric)
        # Swapping arguments on the constructor
        self.datasets_pair = SparseDenseDatasetsPair(Y, X, distance_metric)

    @final
    cdef ITYPE_t n_X(self) nogil:
        # Swapping interface
        return self.datasets_pair.n_Y()

    @final
    cdef ITYPE_t n_Y(self) nogil:
        # Swapping interface
        return self.datasets_pair.n_X()

    @final
    cdef DTYPE_t proxy_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        # Swapping arguments on the same interface
        return self.datasets_pair.proxy_dist(j, i)

    @final
    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        # Swapping arguments on the same interface
        return self.datasets_pair.dist(j, i)
