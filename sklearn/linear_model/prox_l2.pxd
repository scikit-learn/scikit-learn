# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from types cimport floating, complexing

"""In-place L21 Group-Lasso thresholding for real and complex data.
"""
cdef void prox_l2(int n_tasks,
                  complexing *Wj,
                  floating reg,
                  floating ajj) nogil except *

