cdef extern from "math.h":
    double pow(double, double) nogil
    double fmin(double, double) nogil
    double sqrt(double) nogil
import numpy as np
cimport cython
cimport numpy as np

cdef class LearningRate:
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz):
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            eta_ptr[idx] = eta0
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge):
        return gradient

cdef class Constant(LearningRate):

    def __reduce__(self):
        return Constant, ()

cdef class Optimal(LearningRate):
    @cython.cdivision(True)
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz):
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            eta_ptr[idx] = 1.0 / (alpha * (t - 1))

    def __reduce__(self):
        return Optimal, ()

cdef class InvScaling(LearningRate):
    @cython.cdivision(True)
    cdef void eta(self, double *eta_ptr, double eta0, double alpha, double t,
                  double power_t, double gradient, double* x_data_ptr,
                  int* x_ind_ptr, int xnnz):
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            eta_ptr[idx] = eta0 / pow(t, power_t)

    def __reduce__(self):
        return InvScaling, ()

cdef class AdaGrad(LearningRate):

    def __cinit__(self, double sum_squared_grad, double eps0, double rho0):
        self.sum_squared_grad = sum_squared_grad
        self.sum_squared_grad_vector = None
        self.eps0 = eps0

    def __reduce__(self):
        return AdaGrad, (self.sum_squared_grad, self.eps0, self.rho0)


cdef class AdaDelta(LearningRate):
    def __cinit__(self, double sum_squared_grad, double eps0, double rho0):
        self.rho0 = rho0
        self.eps0 = eps0
        self.accugrad = 0
        self.accudelta = 0

    def __reduce__(self):
        return AdaDelta, (self.sum_squared_grad, self.eps0, self.rho0)

cdef class PA(LearningRate):
    cdef double _get_multiplier(self, int is_hinge, double p, double y):
        if is_hinge:
            # classification
            return y
        elif y - p < 0.0:
            # regression
            return -1.0
        else:
            return 1.0

cdef class PA1(PA):

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge):
        update = loss / norm
        update = fmin(C, update)
        update *= self._get_multiplier(is_hinge, p, y)
        return update

    def __reduce__(self):
        return PA1, ()

cdef class PA2(PA):

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge):
        update = loss / (norm + 0.5 / C)
        update *= self._get_multiplier(is_hinge, p, y)
        return update

    def __reduce__(self):
        return PA2, ()
