cdef extern from "string.h":
    int strcmp(char*, char*)
cdef extern from "math.h":
    double pow(double, double) nogil
cdef extern from "math.h":
    double fmin(double, double) nogil
cimport cython

cdef LearningRate get_learning_rate(char* learning_rate):
    if strcmp(learning_rate, "constant") == 0:
        return Constant.__new__(Constant)
    elif strcmp(learning_rate, "optimal") == 0:
        return Optimal.__new__(Optimal)
    elif strcmp(learning_rate, "invscaling") == 0:
        return InvScaling.__new__(InvScaling)
    elif strcmp(learning_rate, "pa1") == 0:
        return PA1.__new__(PA1)
    elif strcmp(learning_rate, "pa2") == 0:
        return PA2.__new__(PA2)

cdef class LearningRate:
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        pass
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil:
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
    cdef double _get_multiplier(self, int is_hinge, double p, double y) nogil:
        if is_hinge:
            # classification
            return y
        elif y - p < 0:
            # regression
            return -1
        else:
            return 1

cdef class PA1(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil:
        update = loss / norm
        update = fmin(C, update)
        update *= self._get_multiplier(is_hinge, p, y)
        return update

cdef class PA2(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil:
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil:
        update = loss / (norm + 0.5 / C)
        update *= self._get_multiplier(is_hinge, p, y)
        return update
