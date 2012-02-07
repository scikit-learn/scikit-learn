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

# -----------------------------------------
# Headers for Weight Vector Extension Types
# -----------------------------------------

cdef class WeightVector:
    cdef w
    cdef double *w_data_ptr
    cdef double wscale
    cdef unsigned int n_features
    cdef double add(self, double *X_data_ptr, unsigned int offset,
                    unsigned int n_features, double c)
    cdef double dot(self, double *X_data_ptr, unsigned int offset,
                    unsigned int n_features)
    cdef void scale(self, double c)
    cdef void reset_wscale(self)
    cdef double norm(self)