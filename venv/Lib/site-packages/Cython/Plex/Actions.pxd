# cython: language_level=3

cdef class Action:
    cdef perform(self, token_stream, text)

cdef class Return(Action):
    cdef object value
    cdef perform(self, token_stream, text)

cdef class Call(Action):
    cdef object function
    cdef perform(self, token_stream, text)

cdef class Method(Action):
    cdef str name
    cdef dict kwargs

cdef class Begin(Action):
    cdef object state_name
    cdef perform(self, token_stream, text)

cdef class Ignore(Action):
    cdef perform(self, token_stream, text)

cdef class Text(Action):
    cdef perform(self, token_stream, text)
