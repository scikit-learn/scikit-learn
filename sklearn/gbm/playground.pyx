cimport cython

cdef class MyClass:
    cdef int width, height

    def __init__(self, int w, int h):
        self.width = w
        self.height = h

def hello():
    o = MyClass(9, 5)
    return zob(o)

cdef int zob (MyClass o) nogil:
    return o.width