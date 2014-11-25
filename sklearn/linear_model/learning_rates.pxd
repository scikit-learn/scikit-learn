cimport numpy as np
from sklearn.utils.weight_vector cimport WeightVector

cdef class LearningRate(object):
    cdef double eta0
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_intercept) nogil
    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge) nogil

cdef class Constant(LearningRate):
    pass

cdef class Optimal(LearningRate):
    cdef double alpha

cdef class InvScaling(LearningRate):
    cdef double power_t

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p) nogil

cdef class Adaptive(LearningRate):
    cdef double alpha
    cdef int n_features
    cdef double eps0
    cdef double intercept_accugrad
    cdef double* accugrad
    cdef np.ndarray accugrad_array
    cdef double _compute_eta(self, double full_gradient,
                             double val, int idx) nogil
    cdef double _compute_intercept_eta(self, double gradient) nogil

cdef class AdaGrad(Adaptive):
    pass

cdef class AdaDelta(Adaptive):
    cdef double rho0
    cdef double intercept_accudelta
    cdef double* accudelta
    cdef np.ndarray accudelta_array

cdef class PA1(PA):
    pass

cdef class PA2(PA):
    pass


