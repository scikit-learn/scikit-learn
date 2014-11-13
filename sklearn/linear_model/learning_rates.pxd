cimport numpy as np
from sklearn.utils.weight_vector cimport WeightVector

cdef class LearningRate(object):
    cdef double eta(self, double *eta_ptr, double eta0, double alpha, double t,
                    double power_t, double gradient, int n_features,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_intercept) nogil
    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge) nogil

cdef class Constant(LearningRate):
    pass

cdef class Optimal(LearningRate):
    pass

cdef class InvScaling(LearningRate):
    pass

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p) nogil

cdef class Adaptive(LearningRate):
    cdef double _compute_eta(self, double full_gradient,
                             double val, int idx, double eta0, int n_features) nogil
    cdef double _compute_intercept_eta(self, double gradient, double eta0) nogil

cdef class AdaGrad(Adaptive):
    cdef double* accugrad
    cdef double intercept_accugrad
    cdef double eps0

cdef class AdaDelta(Adaptive):
    cdef double rho0
    cdef double eps0
    cdef double intercept_accugrad
    cdef double intercept_accudelta
    cdef double* accugrad
    cdef double* accudelta

cdef class PA1(PA):
    pass

cdef class PA2(PA):
    pass


