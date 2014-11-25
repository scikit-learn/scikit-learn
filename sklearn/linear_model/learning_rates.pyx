cdef extern from "math.h":
    double pow(double, double) nogil
    double fmin(double, double) nogil
    double sqrt(double) nogil
import numpy as np
from libc.stdlib cimport malloc
cimport cython
cimport numpy as np
from sklearn.utils.weight_vector cimport WeightVector

cdef class LearningRate:
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_itercept) nogil:
        eta_ptr[0] = self.eta0
        return eta_ptr[0]

    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge) nogil:
        return gradient

cdef class Constant(LearningRate):
    def __cinit__(self, double eta0, double alpha, double power_t,
                  int n_features, double eps0, double rho0):
        self.eta0 = eta0

    def __reduce__(self):
        return Constant, (self.eta0, 0, 0, 0, 0, 0)

cdef class Optimal(LearningRate):
    def __cinit__(self, double eta0, double alpha, double power_t,
                  int n_features, double eps0, double rho0):
        self.eta0 = eta0
        self.alpha = alpha

    @cython.cdivision(True)
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_itercept) nogil:
        cdef double alpha = self.alpha
        eta_ptr[0] = 1.0 / (alpha * (t - 1))
        return eta_ptr[0]

    def __reduce__(self):
        return Optimal, (self.eta0, self.alpha, 0, 0, 0, 0)

cdef class InvScaling(LearningRate):
    def __cinit__(self, double eta0, double alpha, double power_t,
                  int n_features, double eps0, double rho0):
        self.eta0 = eta0
        self.power_t = power_t

    @cython.cdivision(True)
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_itercept) nogil:
        cdef double eta0 = self.eta0
        cdef double power_t = self.power_t

        eta_ptr[0] = eta0 / pow(t, power_t)
        return eta_ptr[0]

    def __reduce__(self):
        return InvScaling, (self.eta0, 0, self.power_t, 0, 0, 0)

cdef class Adaptive(LearningRate):
    cdef double _compute_eta(self, double full_gradient,
                             double val, int idx) nogil:
        pass

    cdef double _compute_intercept_eta(self, double full_gradient) nogil:
        pass

    @cython.cdivision(True)
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_itercept) nogil:

        cdef int counter = 0
        cdef int j = 0
        cdef double val
        cdef double intercept_eta = 0.0
        cdef double alpha = self.alpha

        # we can use a sparse trick here
        for j in range(xnnz):
            idx = x_ind_ptr[j]
            val = x_data_ptr[j]
            full_gradient = gradient * val
            eta_ptr[idx] = self._compute_eta(full_gradient, val, idx)

        if fit_itercept == 1:
            intercept_eta = self._compute_intercept_eta(gradient)

        # TODO: remove this
        with gil:
            for j in range(2):
                print(eta_ptr[j])

        return intercept_eta

cdef class AdaGrad(Adaptive):

    def __cinit__(self, double eta0, double alpha, double power_t,
                  int n_features, double eps0, double rho0):
        self.eta0 = eta0
        self.alpha = alpha
        self.n_features = n_features
        self.eps0 = eps0
        self.intercept_accugrad = 0.0
        self.accugrad_array = np.zeros((n_features,),
                                       dtype=np.float64,
                                       order="c")
        self.accugrad = <double *> self.accugrad_array.data

        # self.accugrad = <double*> malloc(n_features * sizeof(double))
        # for i in range(n_features):
        #     self.accugrad[i] = 0.0

    @cython.cdivision(True)
    cdef double _compute_eta(self, double full_gradient,
                             double val, int idx) nogil:
        cdef double eps0 = self.eps0
        cdef double eta0 = self.eta0
        cdef double* accugrad = self.accugrad
        accugrad[idx] += full_gradient * full_gradient

        return eta0 / sqrt(accugrad[idx] + eps0)

    @cython.cdivision(True)
    cdef double _compute_intercept_eta(self, double gradient) nogil:
        cdef double eps0 = self.eps0
        cdef double eta0 = self.eta0
        self.intercept_accugrad += gradient * gradient
        return eta0 / sqrt(self.intercept_accugrad + eps0)

    def __reduce__(self):
        return AdaGrad, (self.eta0, self.alpha, 0,
                         self.n_features, self.eps0, 0)


cdef class AdaDelta(Adaptive):
    def __cinit__(self, double eta0, double alpha, double power_t,
                  int n_features, double eps0, double rho0):
        self.eta0 = eta0
        self.alpha = alpha
        self.n_features = n_features
        self.rho0 = rho0
        self.eps0 = eps0
        self.intercept_accugrad = 0.0
        self.intercept_accudelta = 0.0
        self.accugrad_array = np.zeros((n_features,),
                                       dtype=np.float64,
                                       order="c")
        self.accudelta_array = np.zeros((n_features,),
                                        dtype=np.float64,
                                        order="c")
        self.accugrad = <double *> self.accugrad_array.data
        self.accudelta = <double *> self.accudelta_array.data

    @cython.cdivision(True)
    cdef double _compute_eta(self, double full_gradient,
                             double val, int idx) nogil:
        cdef double eta0 = self.eta0
        cdef double eps0 = self.eps0
        cdef double rho0 = self.rho0
        cdef double* accugrad = self.accugrad
        cdef double* accudelta = self.accudelta

        accugrad[idx] *= rho0
        accugrad[idx] += (1. - rho0) * full_gradient * full_gradient
        dx = sqrt((accudelta[idx] + eps0) / (accugrad[idx] + eps0))
        accudelta[idx] *= rho0
        accudelta[idx] += (1. - rho0) * dx * dx

        return eta0 * dx

    @cython.cdivision(True)
    cdef double _compute_intercept_eta(self, double gradient) nogil:
        cdef double eta0 = self.eta0
        cdef double eps0 = self.eps0
        cdef double rho0 = self.rho0

        self.intercept_accugrad *= rho0
        self.intercept_accugrad += (1. - rho0) * gradient * gradient
        dx = sqrt((self.intercept_accudelta + eps0) /
                  (self.intercept_accugrad + eps0))
        self.intercept_accudelta *= rho0
        self.intercept_accudelta += (1. - rho0) * dx * dx

        return eta0 * dx

    def __reduce__(self):
        return AdaDelta, (self.eta0, self.alpha, 0, self.n_features,
                          self.eps0, self.rho0)

cdef class PA(LearningRate):
    cdef double eta(self, double *eta_ptr, double t, double gradient,
                    double* x_data_ptr, int* x_ind_ptr, int xnnz,
                    WeightVector w, double intercept, int fit_itercept) nogil:
        eta_ptr[0] = -1.0

    cdef double _get_multiplier(self, int is_hinge, double p, double y) nogil:
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
    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge) nogil:
        update = loss / norm
        update = fmin(C, update)
        update *= self._get_multiplier(is_hinge, p, y)
        return update

    def __reduce__(self):
        return PA1, (0, 0, 0, 0, 0, 0)

cdef class PA2(PA):

    @cython.cdivision(True)
    cdef double update(self, double gradient, double loss,
                       double norm, double C, double p, double y,
                       int is_hinge) nogil:
        update = loss / (norm + 0.5 / C)
        update *= self._get_multiplier(is_hinge, p, y)
        return update

    def __reduce__(self):
        return PA2, (0, 0, 0, 0, 0, 0)
