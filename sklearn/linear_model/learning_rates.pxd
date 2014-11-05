cimport numpy as np

cdef class LearningRate(object):
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge)

cdef class Constant(LearningRate):
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz)

cdef class Optimal(LearningRate):
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz)

cdef class InvScaling(LearningRate):
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz)

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double y, double p)

cdef class AdaGrad(LearningRate):
    cdef double sum_squared_grad
    cdef np.ndarray sum_squared_grad_vector
    cdef double eps0

cdef class AdaDelta(LearningRate):
    cdef double rho0
    cdef double eps0
    cdef double accugrad
    cdef double accudelta

cdef class PA1(PA):
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge)

cdef class PA2(PA):
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge)


