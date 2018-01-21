# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from types cimport floating, complexing
from blas_api cimport (fused_nrm2, fused_scal)
from utils cimport fmax


cdef inline void prox_l2(int n,
                         complexing *w,
                         floating reg,
                         floating ajj) nogil except *:
    cdef floating scaling

    # N.B.: scaling = ||w||_2
    fused_nrm2(n,
               w,
               1,
               &scaling)
    if scaling != 0.:  # Else w must be the zero vector; do nothing
        scaling = fmax(1. - reg / scaling, 0.)
        if ajj != 0.:
            scaling = scaling / ajj
        else:
            scaling = 0.

        # N.B.: w *= scaling
        if scaling != 1.:
            fused_scal(n,
                       <complexing>scaling,
                       w,
                       1)
