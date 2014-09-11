# cdef extern from "stdlib.h":
#     int strcmp(const char*, const char*)
cdef extern from "math.h":
    double pow(double, double) nogil
cdef extern from "math.h":
    double fmin(double, double) nogil
cimport cython

cdef LearningRate get_learning_rate(int learning_rate):
    if learning_rate == 1:
        return Constant.__new__(Constant)
    elif learning_rate == 2:
        return Optimal.__new__(Optimal)
    elif learning_rate == 3:
        return InvScaling.__new__(InvScaling)
    elif learning_rate == 4:
        return PA1.__new__(PA1)
    elif learning_rate == 5:
        return PA2.__new__(PA2)

cdef class LearningRate:
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        pass
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double y, double p, int is_hinge) nogil:
        return -eta * gradient

cdef class Constant(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0

cdef class Optimal(LearningRate):
    @cython.cdivision(True)
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return 1.0 / (alpha * t)

cdef class InvScaling(LearningRate):
    @cython.cdivision(True)
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0 / pow(t, power_t)

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p) nogil:
        if is_hinge:
            # classification
            update *= y
        elif y - p < 0:
            # regression
            update *= -1

cdef class PA1(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double y, double p, int is_hinge) nogil:
        update = loss / norm
        update = fmin(C, update)
        update *= self._get_multiplier(is_hinge, y, p)
        return update

cdef class PA2(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double y, double p, int is_hinge) nogil:
        update = loss / (norm + 0.5 / C)
        update *= self._get_multiplier(is_hinge, y, p)
        return update
