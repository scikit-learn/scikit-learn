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


cdef extern from "float.h":
   cdef double DBL_MAX
   cdef float FLT_MAX

cdef extern from "src/cholesky_delete.c":
    int c_cholesky_delete "cholesky_delete" (int m, int n, double *L, int go_out)
    
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

def solve_triangular (
    np.ndarray[DOUBLE, ndim=2] X,
    np.ndarray[DOUBLE, ndim=1] y
    ):
    """
    Solves a triangular system (overwrites y)

    TODO: option for upper, lower, etc.
    """
    cdef int lda = <int> X.strides[0] / sizeof(double)

    cblas_dtrsv (CblasRowMajor, CblasLower, CblasNoTrans,
                 CblasNonUnit, <int> X.shape[0], <double *> X.data,
                 lda, <double *> y.data, 1);


def cholesky_delete (
    np.ndarray[DOUBLE, ndim=2] L,
    int go_out):

    cdef int n = <int> L.shape[0]
    cdef int m = <int> L.strides[0] / sizeof (double)
    c_cholesky_delete (m, n, <double *> L.data, go_out)
