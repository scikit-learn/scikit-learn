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
cdef class cLossFunction:
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cHalfSquaredError(cLossFunction):
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cAbsoluteError(cLossFunction):
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cPinballLoss(cLossFunction):
    cdef readonly double quantile  # readonly makes it inherited by children
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cHalfPoissonLoss(cLossFunction):
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cHalfGammaLoss(cLossFunction):
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cHalfTweedieLoss(cLossFunction):
    cdef readonly double power  # readonly makes it inherited by children
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil


cdef class cBinaryCrossEntropy(cLossFunction):
    cdef double closs(self, double y_true, double raw_prediction) nogil
    cdef double cgradient(self, double y_true, double raw_prediction) nogil
    cdef double_pair cgrad_hess(self, double y_true, double raw_prediction) nogil
