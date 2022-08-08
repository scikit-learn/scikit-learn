from ...utils._typedefs cimport DTYPE_t, ITYPE_t
from ...metrics._dist_metrics cimport DistanceMetric


cdef class DatasetsPair:
    cdef DistanceMetric distance_metric

    cdef ITYPE_t n_samples_X(self) nogil

    cdef ITYPE_t n_samples_Y(self) nogil

    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil

    cdef DTYPE_t surrogate_dist(self, ITYPE_t i, ITYPE_t j) nogil


cdef class DenseDenseDatasetsPair(DatasetsPair):
    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y
        ITYPE_t d
