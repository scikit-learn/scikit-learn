# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating

"""Projection unto l1 ball.
"""
cdef void enet_projection(unsigned int m, floating *v, floating *out, floating radius,
                          floating l1_ratio) nogil
cdef void proj_l1(int n, floating *w, floating reg,
                  floating ajj) nogil

