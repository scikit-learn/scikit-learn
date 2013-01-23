cimport cython

cdef extern from "src/Gamma.h":
    cdef double LogGamma(double x) except +

cdef double lgamma(double x):
    return LogGamma(x)
