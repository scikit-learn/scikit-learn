# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.



################################################################################
# Extension types for various classification and regression loss functions


cdef class LossFunction:
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Regression(LossFunction):
    cpdef double loss(self,double p, double y)
    cpdef double dloss(self,double p, double y)

cdef class Classification(LossFunction):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class ModifiedHuber(Classification):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Hinge(Classification):
    cdef double threshold
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Log(Classification):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class SquaredLoss(Regression):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Huber(Regression):
    cdef double c
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)