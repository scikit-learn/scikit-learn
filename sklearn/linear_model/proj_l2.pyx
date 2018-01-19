# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

cimport numpy as np
from cython cimport floating
from blas_api cimport fused_nrm2, fused_scal


cdef inline void proj_l2(int n, floating *w, floating reg,
                         floating ajj) nogil:
    """Computes (in-place)

        argmin .5 * ||z - w / ajj||_2^2 subject to ||z||_2 <= reg
    """
    cdef floating scaling

    # N.B.: scaling = ||w||_2
    if ajj == 0.:
        scaling = 0.
    else:
        scaling = fused_nrm2(n, w, 1)
        if scaling > ajj * reg:
            scaling = reg / scaling
        else:
            scaling = 1. / ajj

    # N.B.: w *= scaling
    fused_scal(n, scaling, w, 1)

