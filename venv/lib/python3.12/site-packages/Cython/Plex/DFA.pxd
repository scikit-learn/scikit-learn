cimport cython

from . cimport Machines


@cython.final
cdef class StateMap:
    cdef Machines.FastMachine new_machine
    cdef dict old_to_new_dict
    cdef dict new_to_old_dict

    cdef old_to_new(self, set old_state_set)
    cdef highest_priority_action(self, set state_set)
    cdef make_key(self, set state_set)
