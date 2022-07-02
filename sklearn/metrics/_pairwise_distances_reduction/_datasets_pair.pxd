from ...utils._typedefs cimport DTYPE_t, ITYPE_t, SPARSE_INDEX_TYPE_t
from ...metrics._dist_metrics cimport DistanceMetric


cdef class DatasetsPair:
    cdef:
        DistanceMetric distance_metric
        ITYPE_t n_features

    cdef ITYPE_t n_samples_X(self) nogil

    cdef ITYPE_t n_samples_Y(self) nogil

    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil

    cdef DTYPE_t surrogate_dist(self, ITYPE_t i, ITYPE_t j) nogil


cdef class DenseDenseDatasetsPair(DatasetsPair):
    cdef:
        const DTYPE_t[:, ::1] X
        const DTYPE_t[:, ::1] Y


cdef class SparseSparseDatasetsPair(DatasetsPair):
    cdef:
        const DTYPE_t[:] X_data
        const SPARSE_INDEX_TYPE_t[:] X_indices
        const SPARSE_INDEX_TYPE_t[:] X_indptr

        const DTYPE_t[:] Y_data
        const SPARSE_INDEX_TYPE_t[:] Y_indices
        const SPARSE_INDEX_TYPE_t[:] Y_indptr


cdef class SparseDenseDatasetsPair(DatasetsPair):
    cdef:
        const DTYPE_t[:] X_data
        const SPARSE_INDEX_TYPE_t[:] X_indices
        const SPARSE_INDEX_TYPE_t[:] X_indptr

        const DTYPE_t[:] Y_data
        const SPARSE_INDEX_TYPE_t[:] Y_indices
        ITYPE_t n_Y


cdef class DenseSparseDatasetsPair(DatasetsPair):
    cdef:
        # As distance metrics are commutative, we can simply rely
        # on the implementation of SparseDenseDatasetsPair and
        # swap arguments.
        DatasetsPair datasets_pair
