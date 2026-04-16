from __future__ import absolute_import

import cython

from Cython.Plex.Actions cimport Action

cdef class Scanner:

    cdef public lexicon
    cdef public stream
    cdef public name
    cdef public unicode buffer
    cdef public Py_ssize_t buf_start_pos
    cdef public Py_ssize_t next_pos
    cdef public Py_ssize_t cur_pos
    cdef public Py_ssize_t cur_line
    cdef public Py_ssize_t cur_line_start
    cdef public Py_ssize_t start_pos
    cdef tuple current_scanner_position_tuple
    cdef public tuple last_token_position_tuple
    cdef public str text
    cdef public initial_state # int?
    cdef public state_name
    cdef public list queue
    cdef public bint trace
    cdef public cur_char
    cdef public long input_state

    cdef public level

    @cython.locals(action=Action)
    cpdef tuple read(self)
    cdef inline unread(self, token, value, position)
    cdef inline get_current_scan_pos(self)
    cdef inline tuple scan_a_token(self)

    # This could be @final and @cfunc, but that would hide it from Parsing.py
    # unless that's compiled as well (which it isn't with "minimal" compilation).
    #@cython.final
    cpdef tuple position(self)  # used frequently by Parsing.py

    cdef run_machine_inlined(self)

    cdef inline begin(self, state)
    cdef inline produce(self, value, text = *)
