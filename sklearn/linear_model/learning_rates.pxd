cimport numpy as np

cdef class LearningRate(object):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr, int* x_ind_ptr,
                       int xnnz, int n_features)

cdef class Constant(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)

cdef class Optimal(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)

cdef class InvScaling(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p)

cdef class AdaGrad(LearningRate):
    cdef double sum_squared_grad
    cdef np.ndarray sum_squared_grad_vector
    cdef double eps0
    cdef double eta(self, double eta0, double alpha, double t, double power_t)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr, int* x_ind_ptr,
                       int xnnz, int n_features)

cdef class AdaDelta(LearningRate):
    cdef double rho0
    cdef double eps0
    cdef double accugrad
    cdef double accudelta
    cdef double eta(self, double eta0, double alpha, double t, double power_t)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr, int* x_ind_ptr,
                       int xnnz, int n_features)

cdef class PA1(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr, int* x_ind_ptr,
                       int xnnz, int n_features)

cdef class PA2(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr, int* x_ind_ptr,
                       int xnnz, int n_features)


