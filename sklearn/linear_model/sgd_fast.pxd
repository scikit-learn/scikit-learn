cdef class LossFunction:
    cdef double loss(self, double p, double y) nogil
    cdef double _dloss(self, double p, double y) nogil

cdef class Classification(LossFunction):
    pass

cdef class Regression(LossFunction):
    pass
