import numpy as np
cimport numpy as np

from .types import HISTOGRAM_DTYPE
from .types cimport X_BINNED_DTYPE_C
from .types cimport Y_DTYPE_C
from .types cimport hist_struct

# See histogram.pyx for docstrings and details

cpdef void _subtract_histograms(
    unsigned int n_bins,
    hist_struct [:] hist_a,
    hist_struct [:] hist_b,
    hist_struct [:] out) nogil

cpdef void _build_histogram(
    unsigned int n_bins,
    unsigned int [:] sample_indices,
    X_BINNED_DTYPE_C [:] binned_feature,
    Y_DTYPE_C [:] ordered_gradients,
    Y_DTYPE_C [:] ordered_hessians,
    hist_struct [:] out) nogil

cpdef void _build_histogram_no_hessian(
    unsigned int n_bins,
    unsigned int [:] sample_indices,
    X_BINNED_DTYPE_C [:] binned_feature,
    Y_DTYPE_C [:] ordered_gradients,
    hist_struct [:] out) nogil

cpdef void _build_histogram_root_no_hessian(
    unsigned int n_bins,
    X_BINNED_DTYPE_C [:] binned_feature,
    Y_DTYPE_C [:] all_gradients,
    hist_struct [:] out) nogil

cpdef void _build_histogram_root(
    unsigned int n_bins,
    X_BINNED_DTYPE_C [:] binned_feature,
    Y_DTYPE_C [:] all_gradients,
    Y_DTYPE_C [:] all_hessians,
    hist_struct [:] out) nogil
