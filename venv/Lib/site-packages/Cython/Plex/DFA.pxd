# cython: auto_pickle=False

cimport cython

from . cimport Machines
from .Transitions cimport TransitionMap


@cython.final
cdef class StateMap:
    cdef Machines.FastMachine new_machine
    cdef dict old_to_new_dict
    cdef dict new_to_old_dict

    cdef old_to_new(self, dict old_state_set)

    @cython.locals(state=Machines.Node)
    cdef highest_priority_action(self, dict state_set)

    cdef make_key(self, dict state_set)


@cython.locals(new_machine=Machines.FastMachine, transitions=TransitionMap)
cpdef nfa_to_dfa(Machines.Machine old_machine, debug=*)

cdef set_epsilon_closure(dict state_set)
cdef dict epsilon_closure(Machines.Node state)

@cython.locals(state_set_2=dict, state2=Machines.Node)
cdef add_to_epsilon_closure(dict state_set, Machines.Node state)
