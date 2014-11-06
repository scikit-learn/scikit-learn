cimport numpy as np
from sklearn.utils.weight_vector cimport WeightVector

cdef class LearningRate(object):
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, int n_features, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz, WeightVector w)
    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge)

cdef class Constant(LearningRate):
    pass

cdef class Optimal(LearningRate):
    pass

cdef class InvScaling(LearningRate):
    pass

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
    pass

cdef class PA2(PA):
    pass


