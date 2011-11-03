cdef extern from "_fast_math.h":
    cdef extern float fast_log(float x)

def py_fast_log(x):
    return fast_log(x)
