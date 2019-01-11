import numpy as np
cimport numpy as np

from .types import HISTOGRAM_DTYPE

ctypedef np.npy_uint8 NPY_X_BINNED_DTYPE
ctypedef np.npy_float32 NPY_Y_DTYPE

cdef packed struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count

cpdef void _subtract_histograms(unsigned int n_bins,
                                hist_struct [:] hist_a,
                                hist_struct [:] hist_b,
                                hist_struct [:] out) nogil

cpdef void _build_histogram(unsigned int n_bins,
                            unsigned int [:] sample_indices,
                            NPY_X_BINNED_DTYPE [:] binned_feature,
                            NPY_Y_DTYPE [:] ordered_gradients,
                            NPY_Y_DTYPE [:] ordered_hessians,
                            hist_struct [:] out) nogil

cpdef void _build_histogram_no_hessian(
    unsigned int n_bins,
    unsigned int [:] sample_indices,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] ordered_gradients,
    hist_struct [:] out) nogil

cpdef void _build_histogram_root_no_hessian(
    unsigned int n_bins,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] all_gradients,
    hist_struct [:] out) nogil

cpdef void _build_histogram_root(
    unsigned int n_bins,
    NPY_X_BINNED_DTYPE [:] binned_feature,
    NPY_Y_DTYPE [:] all_gradients,
    NPY_Y_DTYPE [:] all_hessians,
    hist_struct [:] out) nogil