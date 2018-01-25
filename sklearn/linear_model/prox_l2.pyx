# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from types cimport floating, complexing
from blas_api cimport (fused_nrm2, fused_scal)
from utils cimport fmax


cdef inline void prox_l2(int n, complexing *w, floating reg,
                         floating ajj) nogil except *:
    cdef floating scaling

    if ajj == 0.:
        scaling = 0.
    else:
        # N.B.: scaling = ||w||_2
        scaling = <floating>fused_nrm2(n, w, 1)
        if scaling == 0.:
            # w must be the zero vector; do nothing
            return
        scaling = fmax(1. - reg / scaling, 0.) / ajj

    # N.B.: w *= scaling
    fused_scal(n, <complexing>scaling, w, 1)
