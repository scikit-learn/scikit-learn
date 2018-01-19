# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating

"""In-place L21 Group-Lasso thresholding.
"""
cdef void prox_l2(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil

