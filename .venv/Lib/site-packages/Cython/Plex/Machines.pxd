cimport cython

from .Actions cimport Action
from .Transitions cimport TransitionMap

cdef int maxint


@cython.final
cdef class Machine:
    cdef readonly list states
    cdef readonly dict  initial_states
    cdef readonly Py_ssize_t next_state_number

    cpdef new_state(self)
    cpdef new_initial_state(self, name)
    cpdef make_initial_state(self, name, state)


@cython.final
cdef class Node:
    cdef readonly TransitionMap transitions
    cdef readonly Action action
    cdef public set epsilon_closure
    cdef readonly Py_ssize_t number
    cdef readonly int action_priority


@cython.final
cdef class FastMachine:
    cdef readonly dict initial_states
    cdef readonly dict new_state_template
    cdef readonly list states
    cdef readonly Py_ssize_t next_number

    cpdef make_initial_state(self, name, state)
