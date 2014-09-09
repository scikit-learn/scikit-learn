# cdef extern from "stdlib.h":
#     int strcmp(const char*, const char*)
cdef extern from "math.h":
    double pow(double, double)
cimport cython
cdef double numerator = 1.0

cdef LearningRate get_learning_rate(int learning_rate):
    if learning_rate == 1:
        return Constant.__new__(Constant)
    elif learning_rate == 2:
        return Optimal.__new__(Optimal)
    elif learning_rate == 3:
        return InvScaling.__new__(InvScaling)

cdef class LearningRate:
    cdef double step(self, double eta0, double alpha, int t, int power_t):
        pass

cdef class Constant(LearningRate):
    cdef double step(self, double eta0, double alpha, int t, int power_t):
        return eta0

cdef class Optimal(LearningRate):
    @cython.cdivision(True)
    cdef double step(self, double eta0, double alpha, int t, int power_t):
        return 1.0 / alpha * t

cdef class InvScaling(LearningRate):
    @cython.cdivision(True)
    cdef double step(self, double eta0, double alpha, int t, int power_t):
        return eta0 / pow(t, power_t)
