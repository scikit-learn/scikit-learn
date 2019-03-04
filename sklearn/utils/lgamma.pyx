# cython: language_level=3

cdef extern from "src/gamma.h":
    cdef double sklearn_lgamma(double x)


cdef double lgamma(double x):
    if x <= 0:
        raise ValueError("x must be strictly positive, got %f" % x)
    return sklearn_lgamma(x)
