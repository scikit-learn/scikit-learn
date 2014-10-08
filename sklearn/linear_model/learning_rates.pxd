cdef LearningRate get_learning_rate(char* learning_rate)

cdef class LearningRate(object):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil

cdef class Constant(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil

cdef class Optimal(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil

cdef class InvScaling(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p) nogil

cdef class PA1(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil

cdef class PA2(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t) nogil
    cdef double update(self, double gradient, double loss, double eta, double norm, double C, double p, double y, int is_hinge) nogil
