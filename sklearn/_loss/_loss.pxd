# cython: language_level=3

import numpy as np
cimport numpy as np

np.import_array()


# Fused types for y_true, y_pred, raw_prediction
ctypedef fused Y_DTYPE_C:
    np.npy_float64
    np.npy_float32


# Fused types for gradient and hessian
ctypedef fused G_DTYPE_C:
    np.npy_float64
    np.npy_float32


# Struct to return 2 doubles
ctypedef struct double_pair:
   double val1
   double val2


# C base class for loss functions
cdef class CyLossFunction:
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfSquaredError(CyLossFunction):
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyAbsoluteError(CyLossFunction):
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyPinballLoss(CyLossFunction):
    cdef readonly double quantile  # readonly makes it accessible from Python
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfPoissonLoss(CyLossFunction):
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfGammaLoss(CyLossFunction):
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfTweedieLoss(CyLossFunction):
    cdef readonly double power  # readonly makes it accessible from Python
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfTweedieLossIdentity(CyLossFunction):
    cdef readonly double power  # readonly makes it accessible from Python
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil


cdef class CyHalfBinomialLoss(CyLossFunction):
    cdef double cy_loss(self, double y_true, double raw_prediction) nogil
    cdef double cy_gradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cy_grad_hess(self, double y_true, double raw_prediction) nogil
