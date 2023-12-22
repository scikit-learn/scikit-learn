# License: BSD 3 clause
"""Helper to load LossFunction from sgd_fast.pyx to sag_fast.pyx"""

from sklearn._loss._loss cimport CyLossFunction

cdef class Regression(CyLossFunction):
    cdef double cy_loss(self, double y, double p) noexcept nogil
    cdef double cy_gradient(self, double y, double p) noexcept nogil


cdef class Classification(CyLossFunction):
    cdef double cy_loss(self, double y, double p) noexcept nogil
    cdef double cy_gradient(self, double y, double p) noexcept nogil


cdef class Log(Classification):
    cdef double cy_loss(self, double y, double p) noexcept nogil
    cdef double cy_gradient(self, double y, double p) noexcept nogil


cdef class SquaredLoss(Regression):
    cdef double cy_loss(self, double y, double p) noexcept nogil
    cdef double cy_gradient(self, double y, double p) noexcept nogil
