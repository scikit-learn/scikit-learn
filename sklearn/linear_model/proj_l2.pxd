# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from types cimport floating, complexing

"""In-place L2 projection for real and complex data.
Projects w / ajj unto the L2 ball defined by ||z||_2 \le reg.
"""
cdef void proj_l2(int n, complexing *w, floating reg,
                  floating ajj) nogil except *

