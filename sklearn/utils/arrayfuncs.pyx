"""
Small collection of auxiliary functions that operate on arrays

"""
cimport numpy as np
import  numpy as np

cimport cython

from libc.float cimport DBL_MAX, FLT_MAX

cdef extern from "src/cholesky_delete.h":
    int cholesky_delete_dbl(int m, int n, double *L, int go_out)
    int cholesky_delete_flt(int m, int n, float  *L, int go_out)

ctypedef np.float64_t DOUBLE


np.import_array()


def min_pos(np.ndarray X):
   """
   Find the minimum value of an array over positive values

   Returns a huge value if none of the values are positive
   """
   if X.dtype.name == 'float32':
      return _float_min_pos(<float *> X.data, X.size)
   elif X.dtype.name == 'float64':
      return _double_min_pos(<double *> X.data, X.size)
   else:
      raise ValueError('Unsupported dtype for array X')


cdef float _float_min_pos(float *X, Py_ssize_t size):
   cdef Py_ssize_t i
   cdef float min_val = DBL_MAX
   for i in range(size):
      if 0. < X[i] < min_val:
         min_val = X[i]
   return min_val


cdef double _double_min_pos(double *X, Py_ssize_t size):
   cdef Py_ssize_t i
   cdef np.float64_t min_val = FLT_MAX
   for i in range(size):
      if 0. < X[i] < min_val:
         min_val = X[i]
   return min_val


# we should be using np.npy_intp or Py_ssize_t for indices, but BLAS wants int
def cholesky_delete(np.ndarray L, int go_out):
    cdef int n = <int> L.shape[0]
    cdef int m = <int> L.strides[0]

    if L.dtype.name == 'float64':
        cholesky_delete_dbl(m / sizeof(double), n, <double *> L.data, go_out)
    elif L.dtype.name == 'float32':
        cholesky_delete_flt(m / sizeof(float),  n, <float *> L.data,  go_out)
    else:
        raise TypeError("unsupported dtype %r." % L.dtype)
