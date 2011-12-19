"""
Small collection of auxiliary functions that operate on arrays

"""
cimport numpy as np
import  numpy as np

cimport cython

cdef extern from "cblas.h":
   enum CBLAS_ORDER:
       CblasRowMajor=101
       CblasColMajor=102
   enum CBLAS_TRANSPOSE:
       CblasNoTrans=111
       CblasTrans=112
       CblasConjTrans=113
       AtlasConj=114
   enum CBLAS_UPLO:
       CblasUpper=121
       CblasLower=122
   enum CBLAS_DIAG:
       CblasNonUnit=131
       CblasUnit=132

   void cblas_dtrsv(CBLAS_ORDER Order,  CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA,  CBLAS_DIAG Diag,
                 int N, double *A, int lda, double *X,
                 int incX)

   void cblas_strsv(CBLAS_ORDER Order,  CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE TransA,  CBLAS_DIAG Diag,
                 int N, float *A, int lda, float *X,
                 int incX)

cdef extern from "float.h":
   cdef double DBL_MAX
   cdef float FLT_MAX

cdef extern from "src/cholesky_delete.c":
    int double_cholesky_delete (int m, int n, double *L, int go_out)
    int float_cholesky_delete  (int m, int n, float  *L, int go_out)
    
ctypedef np.float64_t DOUBLE


def min_pos(np.ndarray X):
   """
   Find the minimum value of an array over positivie values

   Returns a huge value if none of the values are positive
   """
   if X.dtype.name == 'float32':
      return _float_min_pos(<float *> X.data, X.size)
   elif X.dtype.name == 'float64':
      return _double_min_pos(<double *> X.data, X.size)
   else:
      raise ValueError('Unsupported dtype for array X')


cdef float _float_min_pos (float *X, Py_ssize_t size):
   cdef Py_ssize_t i
   cdef float min_val = DBL_MAX
   for i in range(size):
      if X[i] > 0. and X[i] < min_val:
         min_val = X[i]
   return min_val


cdef double _double_min_pos (double *X, Py_ssize_t size):
   cdef Py_ssize_t i
   cdef np.float64_t min_val = FLT_MAX
   for i in range(size):
      if X[i] > 0. and X[i] < min_val:
         min_val = X[i]
   return min_val

def solve_triangular (np.ndarray X, np.ndarray y):
    """
    Solves a triangular system (overwrites y)

    Note: The lapack function to solve triangular systems was added to
    scipy v0.9. Remove this when we stop supporting earlier versions.
    """
    cdef int lda

    if X.dtype.name == 'float64' and y.dtype.name == 'float64':
       lda = <int> X.strides[0] / sizeof(double)

       cblas_dtrsv (CblasRowMajor, CblasLower, CblasNoTrans,
                    CblasNonUnit, <int> X.shape[0], <double *> X.data,
                    lda, <double *> y.data, 1);

    elif X.dtype.name == 'float32' and y.dtype.name == 'float32':
       lda = <int> X.strides[0] / sizeof(float)

       cblas_strsv (CblasRowMajor, CblasLower, CblasNoTrans,
                    CblasNonUnit, <int> X.shape[0], <float *> X.data,
                    lda, <float *> y.data, 1);
    else:
       raise ValueError ('Unsupported or inconsistent dtype in arrays X, y')


def cholesky_delete (np.ndarray L, int go_out):

    cdef int n = <int> L.shape[0]
    cdef int m

    if L.dtype.name == 'float64':
       m = <int> L.strides[0] / sizeof (double)
       double_cholesky_delete (m, n, <double *> L.data, go_out)
    elif L.dtype.name == 'float32':
       m = <int> L.strides[0] / sizeof (float)
       float_cholesky_delete (m, n, <float *> L.data, go_out)

# TODO: fast_sqdist_foo and fast_pair_sqdist_foo should be templated
# for different dtypes, when the templating infrastructure is finally
# set up.
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fast_sqdist_float32(np.ndarray[np.float32_t, ndim=2] X,
                          np.ndarray[np.float32_t, ndim=2] Y,
                          np.ndarray[np.float32_t, ndim=2] out):
    """
    fast_sqdist_float32(X, Y, out)

    Low-level function for computing squared Euclidean distances between
    two sets of points.

    Parameters
    ----------
    X : ndarray, 2-dimensional, float32
        First set of points, indexed along first axis.
    Y : ndarray, 2-dimensional, float32
        Second set of points, indexed along first axis. Assumed
        that X.shape[1] == Y.shape[1].
    out : ndarray, 2-dimensional, float32
        Array in which to store the output. Assumed to have shape
        (X.shape[0], Y.shape[0]).

    Notes
    -----
    In order to achieve maximal speed this function performs no checks of
    any array metadata. Use the high-level functions defined in the
    `sklearn.metrics.distance` module unless you really know what you're
    doing.
    """
    cdef np.npy_intp i, j, k
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            out[i, j] = 0
            for k in range(X.shape[1]):
                out[i, j] += (X[i, k] - Y[j, k]) ** 2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fast_sqdist_float64(np.ndarray[np.float64_t, ndim=2] X,
                          np.ndarray[np.float64_t, ndim=2] Y,
                          np.ndarray[np.float64_t, ndim=2] out):
    """
    fast_sqdist_float64(X, Y, out)

    Low-level function for computing squared Euclidean distances between
    two sets of points.

    Parameters
    ----------
    X : ndarray, 2-dimensional, float64
        First set of points, indexed along first axis.
    Y : ndarray, 2-dimensional, float64
        Second set of points, indexed along first axis. Assumed
        that X.shape[1] == Y.shape[1].
    out : ndarray, 2-dimensional, float64
        Array in which to store the output. Assumed to have shape
        (X.shape[0], Y.shape[0]).

    Notes
    -----
    In order to achieve maximal speed this function performs no checks of
    any array metadata. Use the high-level functions defined in the
    `sklearn.metrics.distance` module unless you really know what you're
    doing.
    """
    cdef np.npy_intp i, j, k
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            out[i, j] = 0
            for k in range(X.shape[1]):
                out[i, j] += (X[i, k] - Y[j, k]) ** 2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fast_pair_sqdist_float32(np.ndarray[np.float32_t, ndim=2] X,
                               np.ndarray[np.float32_t, ndim=2] out):
    """
    fast_pair_sqdist_float32(X, out)

    Low-level function for computing squared Euclidean distances between
    every pair of points in a set.

    Parameters
    ----------
    X : ndarray, 2-dimensional, float32
        The set of points, indexed along first axis.
    out : ndarray, 2-dimensional, float32
        Array in which to store the output. Assumed to have shape
        (X.shape[0], X.shape[0]).

    Notes
    -----
    In order to achieve maximal speed this function performs no checks of
    any array metadata. Use the high-level functions defined in the
    `sklearn.metrics.distance` module unless you really know what you're
    doing.
    """
    cdef np.npy_intp i, j, k
    for i in range(X.shape[0]):
        out[i, i] = 0
        for j in range(i + 1, X.shape[0]):
            out[i, j] = 0
            for k in range(X.shape[1]):
                out[i, j] += (X[i, k] - X[j, k]) ** 2
            out[j, i] = out[i, j]


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fast_pair_sqdist_float64(np.ndarray[np.float64_t, ndim=2] X,
                               np.ndarray[np.float64_t, ndim=2] out):
    """
    fast_pair_sqdist_float64(X, out)

    Low-level function for computing squared Euclidean distances between
    every pair of points in a set.

    Parameters
    ----------
    X : ndarray, 2-dimensional, float32
        The set of points, indexed along first axis.
    out : ndarray, 2-dimensional, float32
        Array in which to store the output. Assumed to have shape
        (X.shape[0], X.shape[0]).

    Notes
    -----
    In order to achieve maximal speed this function performs no checks of
    any array metadata. Use the high-level functions defined in the
    `sklearn.metrics.distance` module unless you really know what you're
    doing.
    """
    cdef np.npy_intp i, j, k
    for i in range(X.shape[0]):
        out[i, i] = 0
        for j in range(i + 1, X.shape[0]):
            out[i, j] = 0
            for k in range(X.shape[1]):
                out[i, j] += (X[i, k] - X[j, k]) ** 2
            out[j, i] = out[i, j]
