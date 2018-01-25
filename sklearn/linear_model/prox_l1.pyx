# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from libc.math cimport fabs
from types cimport floating, complexing
from utils cimport fmax, cabs, cabsf
from blas_api cimport fused_scal

cdef inline void prox_l1(int n, complexing *w, floating reg,
                         floating ajj) nogil except *:
    if ajj == 0.:
        fused_scal(n, <complexing>ajj, w, 1)
        return
    if complexing is double or complexing is float:
        abs = fabs
    elif complexing is complex:
        abs = cabs
    else:
        abs = cabsf
    cdef int k
    for k in range(n):
        if w[k] != 0.:
            if ajj != 0.:
                w[k] = w[k] * fmax(1. - reg / abs(w[k]), 0.) / ajj
            else:
                w[k] = 0.
