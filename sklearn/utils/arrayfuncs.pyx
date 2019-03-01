# cython: language_level=3
"""
Small collection of auxiliary functions that operate on arrays

"""

cimport numpy as np
import  numpy as np
cimport cython
from cython cimport floating
from libc.math cimport fabs
from libc.float cimport DBL_MAX, FLT_MAX

from ._cython_blas cimport _copy, _rotg, _rot

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


# General Cholesky Delete.
# Remove an element from the cholesky factorization
# m = columns
# n = rows
#
# TODO: put transpose as an option
def cholesky_delete(np.ndarray[floating, ndim=2] L, int go_out):
   cdef:
      int n = L.shape[0]
      int m = L.strides[0]
      floating c, s
      floating *L1
      int i
   
   if floating is float:
      m /= sizeof(float)
   else:
      m /= sizeof(double)

   # delete row go_out
   L1 = &L[0, 0] + (go_out * m)
   for i in range(go_out, n-1):
      _copy(i + 2, L1 + m, 1, L1, 1)
      L1 += m

   L1 = &L[0, 0] + (go_out * m)
   for i in range(go_out, n-1):
      _rotg(L1 + i, L1 + i + 1, &c, &s)
      if L1[i] < 0:
         # Diagonals cannot be negative
         L1[i] = fabs(L1[i])
         c = -c
         s = -s

      L1[i + 1] = 0.  # just for cleanup
      L1 += m

      _rot(n - i - 2, L1 + i, m, L1 + i + 1, m, c, s)
