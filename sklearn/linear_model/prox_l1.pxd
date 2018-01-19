# Synopsis: Soft-thresholding operator
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating

"""In-place soft-thresholding for real and complex data.
"""
cdef void prox_l1(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil

