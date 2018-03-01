# Synopsis: Declaration of proximal operators
# Each solves a problem of the form
#
#      argmin .5 * ||ajj * z - wj||^2 + ajj * g(z)
#         z
# where g is reg times the penalty, or else we are dealing with a constraint
# of the form norm(z) <= reg.
#
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

from cython cimport floating

cdef void prox_l2(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil
cdef void proj_l2(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil
cdef void prox_l1(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil
cdef void proj_l1(int n_tasks, floating *Wj, floating reg,
                  floating ajj) nogil
