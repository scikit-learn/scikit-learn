cimport cython

cdef class Shrubbery:
    cdef int width, height

    def __init__(self, int w, int h):
        self.width = w
        self.height = h