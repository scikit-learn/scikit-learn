# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.

cimport numpy as np

cdef extern from "math.h":
    cdef extern double exp(double x)
    cdef extern double log(double x)
    cdef extern double sqrt(double x)
    cdef extern double pow(double x, double y)

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

# -----------------------------------------
# Headers for Loss Function Extension Types
# -----------------------------------------

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
