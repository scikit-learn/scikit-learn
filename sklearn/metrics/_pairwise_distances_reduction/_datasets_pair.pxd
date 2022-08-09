cimport numpy as cnp
from ...metrics._dist_metrics cimport DistanceMetric


cdef class DatasetsPair:
    cdef DistanceMetric distance_metric

    cdef cnp.intp_t n_samples_X(self) nogil

    cdef cnp.intp_t n_samples_Y(self) nogil

    cdef cnp.float64_t dist(self, cnp.intp_t i, cnp.intp_t j) nogil

    cdef cnp.float64_t surrogate_dist(self, cnp.intp_t i, cnp.intp_t j) nogil


cdef class DenseDenseDatasetsPair(DatasetsPair):
    cdef:
        const cnp.float64_t[:, ::1] X
        const cnp.float64_t[:, ::1] Y
        cnp.intp_t d
