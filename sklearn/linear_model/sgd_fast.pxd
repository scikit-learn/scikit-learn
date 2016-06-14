"""Helper to load LossFunction from sgd_fast.pyx to sag_fast.pyx"""
# Licence: BSD 3 clause

from cython cimport floating

cdef class LossFunction:
    cdef floating loss(self, floating p, floating y) nogil
    cdef floating _dloss(self, floating p, floating y) nogil


cdef class Regression(LossFunction):
    cdef floating loss(self, floating p, floating y) nogil
    cdef floating _dloss(self, floating p, floating y) nogil


cdef class Classification(LossFunction):
    cdef floating loss(self, floating p, floating y) nogil
    cdef floating _dloss(self, floating p, floating y) nogil


cdef class Log(Classification):
    cdef floating loss(self, floating p, floating y) nogil
    cdef floating _dloss(self, floating p, floating y) nogil


cdef class SquaredLoss(Regression):
    cdef floating loss(self, floating p, floating y) nogil
    cdef floating _dloss(self, floating p, floating y) nogil
