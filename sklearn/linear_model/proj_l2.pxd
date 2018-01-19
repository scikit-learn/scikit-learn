# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating

"""In-place L2 projection. Projects w / ajj unto the L2 ball defined
by ||z||_2 <= reg.
"""
cdef void proj_l2(int n, floating *w, floating reg,
                  floating ajj) nogil

