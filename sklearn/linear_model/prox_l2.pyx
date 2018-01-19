# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating
from blas_api cimport (fused_nrm2, fused_scal)
from utils cimport fmax


cdef inline void prox_l2(int n, floating *w, floating reg,
                         floating ajj) nogil:
    """Computes (in-place)

        argmin .5 * ||z - w / ajj||_2^2 + (reg / ajj) * ||z||_2
          z
    """
    cdef floating scaling

    if ajj == 0.:
        scaling = 0.
    else:
        # N.B.: scaling = ||w||_2
        scaling = fused_nrm2(n, w, 1)
        if scaling == 0.:
            # w must be the zero vector; do nothing
            return
        scaling = fmax(1. - reg / scaling, 0.) / ajj

    # N.B.: w *= scaling
    fused_scal(n, scaling, w, 1)

