cdef extern from "math.h":
    double pow(double, double) nogil
    double fmin(double, double) nogil
    double sqrt(double) nogil
import numpy as np
cimport cython
cimport numpy as np

cdef class LearningRate:
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        pass
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr,
                       int* x_ind_ptr, int xnnz, int n_features):
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            gradient_ptr[idx] = -eta * gradient
        return -eta * gradient

cdef class Constant(LearningRate):
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0

    def __reduce__(self):
        return Constant, ()

cdef class Optimal(LearningRate):
    @cython.cdivision(True)
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return 1.0 / (alpha * (t - 1))

    def __reduce__(self):
        return Optimal, ()

cdef class InvScaling(LearningRate):
    @cython.cdivision(True)
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0 / pow(t, power_t)

    def __reduce__(self):
        return InvScaling, ()

cdef class AdaGrad(LearningRate):

    def __cinit__(self, double sum_squared_grad, double eps0, double rho0):
        self.sum_squared_grad = sum_squared_grad
        self.sum_squared_grad_vector = None
        self.eps0 = eps0

    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr,
                       int* x_ind_ptr, int xnnz, int n_features):
        if self.sum_squared_grad_vector is None:
            self.sum_squared_grad_vector = np.zeros((n_features,), dtype=np.float64, order="c")
        self.sum_squared_grad += gradient ** 2.0 + self.eps0
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            self.sum_squared_grad_vector[idx] += gradient * val * gradient * val + self.eps0
            gradient_ptr[idx] = -(eta / sqrt(self.sum_squared_grad_vector[idx])) * gradient
        return -(eta / sqrt(self.sum_squared_grad)) * gradient

    def __reduce__(self):
        return AdaGrad, (self.sum_squared_grad, self.eps0, self.rho0)


cdef class AdaDelta(LearningRate):
    def __cinit__(self, double sum_squared_grad, double eps0, double rho0):
        self.rho0 = rho0
        self.eps0 = eps0
        self.accugrad = 0
        self.accudelta = 0

    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr,
                       int* x_ind_ptr, int xnnz, int n_features):
        agrad = self.rho0 * self.accugrad + \
            (1. - self.rho0) * gradient * gradient
        dx = - sqrt((self.accudelta + self.eps0) /
                    (agrad + self.eps0)) * gradient
        self.accudelta = self.rho0 * self.accudelta +\
            (1. - self.rho0) * dx * dx
        self.accugrad = agrad
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            gradient_ptr[idx] = dx
        return dx

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
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr,
                       int* x_ind_ptr, int xnnz, int n_features):
        update = loss / norm
        update = fmin(C, update)
        update *= self._get_multiplier(is_hinge, p, y)
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            gradient_ptr[idx] = update
        return update

    def __reduce__(self):
        return PA1, ()

cdef class PA2(PA):
    cdef double eta(self, double eta0, double alpha, double t, double power_t):
        return eta0

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss, double eta,
                       double norm, double C, double p, double y,
                       int is_hinge, double* gradient_ptr, double* x_data_ptr,
                       int* x_ind_ptr, int xnnz, int n_features):
        update = loss / (norm + 0.5 / C)
        update *= self._get_multiplier(is_hinge, p, y)
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            gradient_ptr[idx] = update
        return update

    def __reduce__(self):
        return PA2, ()
