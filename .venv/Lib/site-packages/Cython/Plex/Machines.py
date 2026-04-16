"""
Python Lexical Analyser

Classes for building NFAs and DFAs
"""

import cython
from .Transitions import TransitionMap

maxint = 2**31-1  # sentinel value

LOWEST_PRIORITY = -maxint


class Machine:
    """A collection of Nodes representing an NFA or DFA."""
    def __init__(self):
        self.states = []  # [Node]
        self.initial_states = {}  # {(name, bol): Node}
        self.next_state_number = 1

    def __del__(self):
        for state in self.states:
            state.destroy()

    def new_state(self):
        """Add a new state to the machine and return it."""
        s = Node()
        n: cython.Py_ssize_t = self.next_state_number
        self.next_state_number = n + 1
        s.number = n
        self.states.append(s)
        return s

    def new_initial_state(self, name):
        state = self.new_state()
        self.make_initial_state(name, state)
        return state

    def make_initial_state(self, name, state):
        self.initial_states[name] = state

    def get_initial_state(self, name):
        return self.initial_states[name]

    def dump(self, file):
        file.write("Plex.Machine:\n")
        if self.initial_states is not None:
            file.write("   Initial states:\n")
            for (name, state) in sorted(self.initial_states.items()):
                file.write("      '%s': %d\n" % (name, state.number))
        for s in self.states:
            s.dump(file)


class Node:
    """A state of an NFA or DFA."""

    def __init__(self):
        # Preinitialise the list of empty transitions, because
        # the nfa-to-dfa algorithm needs it
        self.transitions = TransitionMap()      # TransitionMap
        self.action_priority = LOWEST_PRIORITY  # integer
        self.action = None  # Action
        self.number = 0     # for debug output
        self.epsilon_closure = None  # used by nfa_to_dfa()

    def destroy(self):
        self.transitions = None
        self.action = None
        self.epsilon_closure = None

    def add_transition(self, event, new_state):
        self.transitions.add(event, new_state)

    def link_to(self, state):
        """Add an epsilon-move from this state to another state."""
        self.add_transition('', state)

    def set_action(self, action, priority):
        """Make this an accepting state with the given action. If
        there is already an action, choose the action with highest
        priority."""
        if priority > self.action_priority:
            self.action = action
            self.action_priority = priority

    def get_action(self):
        return self.action

    def get_action_priority(self):
        return self.action_priority

    def is_accepting(self):
        return self.action is not None

    def __str__(self):
        return "State %d" % self.number

    def dump(self, file):
        # Header
        file.write("   State %d:\n" % self.number)
        # Transitions
        #        self.dump_transitions(file)
        self.transitions.dump(file)
        # Action
        action = self.action
        priority = self.action_priority
        if action is not None:
            file.write("      %s [priority %d]\n" % (action, priority))

    def __lt__(self, other):
        return self.number < other.number

    def __hash__(self):
        # Prevent overflowing hash values due to arbitrarily large unsigned addresses.
        return id(self) & maxint


class FastMachine:
    """
    FastMachine is a deterministic machine represented in a way that
    allows fast scanning.
    """
    def __init__(self):
        self.initial_states = {}  # {state_name:state}
        self.states = []          # [state]  where state = {event:state, 'else':state, 'action':Action}
        self.next_number = 1      # for debugging
        self.new_state_template = {
            '': None, 'bol': None, 'eol': None, 'eof': None, 'else': None
        }

    def __del__(self):
        for state in self.states:
            state.clear()

    def new_state(self, action=None):
        number: cython.Py_ssize_t = self.next_number
        self.next_number = number + 1
        result = self.new_state_template.copy()
        result['number'] = number
        result['action'] = action
        self.states.append(result)
        return result

    def make_initial_state(self, name, state):
        self.initial_states[name] = state

    def add_transitions(self, state: dict, event, new_state, maxint: cython.int = maxint):
        code:  cython.int
        code0: cython.int
        code1: cython.int

        if type(event) is tuple:
            code0, code1 = event
            if code0 == -maxint:
                state['else'] = new_state
            elif code1 != maxint:
                for code in range(code0, code1):
                    state[chr(code)] = new_state
        else:
            state[event] = new_state

    def get_initial_state(self, name):
        return self.initial_states[name]

    def dump(self, file):
        file.write("Plex.FastMachine:\n")
        file.write("   Initial states:\n")
        for name, state in sorted(self.initial_states.items()):
            file.write("      %s: %s\n" % (repr(name), state['number']))
        for state in self.states:
            self.dump_state(state, file)

    def dump_state(self, state, file):
        # Header
        file.write("   State %d:\n" % state['number'])
        # Transitions
        self.dump_transitions(state, file)
        # Action
        action = state['action']
        if action is not None:
            file.write("      %s\n" % action)

    def dump_transitions(self, state, file):
        chars_leading_to_state = {}
        special_to_state = {}
        for (c, s) in state.items():
            if len(c) == 1:
                chars = chars_leading_to_state.get(id(s))
                if chars is None:
                    chars = []
                    chars_leading_to_state[id(s)] = chars
                chars.append(c)
            elif len(c) <= 4:
                special_to_state[c] = s
        ranges_to_state = {}
        for state in self.states:
            char_list = chars_leading_to_state.get(id(state))
            if char_list:
                ranges = self.chars_to_ranges(char_list)
                ranges_to_state[ranges] = state
        for ranges in sorted(ranges_to_state):
            key = self.ranges_to_string(ranges)
            state = ranges_to_state[ranges]
            file.write("      %s --> State %d\n" % (key, state['number']))
        for key in ('bol', 'eol', 'eof', 'else'):
            state = special_to_state.get(key)
            if state:
                file.write("      %s --> State %d\n" % (key, state['number']))

    def chars_to_ranges(self, char_list: list) -> tuple:
        char_list.sort()

        c1: cython.Py_UCS4
        c2: cython.Py_UCS4
        i: cython.Py_ssize_t = 0
        n: cython.Py_ssize_t = len(char_list)
        result = []
        while i < n:
            c1 = ord(char_list[i])
            c2 = c1
            i += 1
            while i < n and ord(char_list[i]) == c2 + 1:
                i += 1
                c2 += 1
            result.append((chr(c1), chr(c2)))
        return tuple(result)

    def ranges_to_string(self, range_list) -> str:
        return ','.join(map(self.range_to_string, range_list))

    def range_to_string(self, range_tuple: tuple):
        (c1, c2) = range_tuple
        if c1 == c2:
            return repr(c1)
        else:
            return f"{c1!r}..{c2!r}"
