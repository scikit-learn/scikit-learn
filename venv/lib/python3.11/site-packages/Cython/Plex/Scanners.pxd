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
    cdef public text
    cdef public initial_state # int?
    cdef public state_name
    cdef public list queue
    cdef public bint trace
    cdef public cur_char
    cdef public long input_state

    cdef public level

    @cython.locals(input_state=long)
    cdef inline next_char(self)
    @cython.locals(action=Action)
    cpdef tuple read(self)
    cdef inline unread(self, token, value, position)
    cdef inline get_current_scan_pos(self)
    cdef inline tuple scan_a_token(self)
    ##cdef tuple position(self)  # used frequently by Parsing.py

    @cython.final
    @cython.locals(cur_pos=Py_ssize_t, cur_line=Py_ssize_t, cur_line_start=Py_ssize_t,
                   input_state=long, next_pos=Py_ssize_t, state=dict,
                   buf_start_pos=Py_ssize_t, buf_len=Py_ssize_t, buf_index=Py_ssize_t,
                   trace=bint, discard=Py_ssize_t, data=unicode, buffer=unicode)
    cdef run_machine_inlined(self)

    cdef inline begin(self, state)
    cdef inline produce(self, value, text = *)
