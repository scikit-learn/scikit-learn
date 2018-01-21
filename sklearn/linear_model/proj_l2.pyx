# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

cimport numpy as np
from types cimport floating, complexing
from blas_api cimport (fused_nrm2, fused_scal)


cdef inline void proj_l2(int n,
                         complexing *w,
                         floating reg,
                         floating ajj) nogil except *:
    cdef floating scaling

    # N.B.: scaling = ||w||_2
    if ajj == 0.:
        scaling = 0.
    else:
        fused_nrm2(n,
                   w,
                   1,
                   &scaling)
        if scaling > ajj * reg:
            scaling = reg / scaling
        else:
            scaling = 1. / ajj

    # N.B.: w *= scaling
    fused_scal(n,
               <complexing>scaling,
               w,
               1)
