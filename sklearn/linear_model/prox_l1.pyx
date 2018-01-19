# Synopsis: Soft-thresholding operator
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from libc.math cimport fabs
from cython cimport floating
from utils cimport fmax, fabs
from blas_api cimport fused_scal


cdef inline void prox_l1(int n, floating *w, floating reg,
                         floating ajj) nogil:
    """Computes (in-place)

        argmin .5 * ||z - w / ajj||_2^2 + (reg / ajj) * ||z||_1
          z
    """
    cdef int k
    if ajj == 0.:
        fused_scal(n, ajj, w, 1)
        return
    for k in range(n):
        if w[k] != 0.:
            w[k] = w[k] * fmax(1. - reg / fabs(w[k]), 0.) / ajj
