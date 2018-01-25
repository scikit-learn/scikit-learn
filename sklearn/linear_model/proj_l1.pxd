# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from types cimport floating, complexing

"""Projection unto l1 ball for real and complex data.
"""
cdef void proj_l1(int n, complexing *w, floating reg,
                  floating ajj) nogil except *

