# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Authors: Danny Sullivan <dbsullivan23@gmail.com>
#          Tom Dupre la Tour <tom.dupre-la-tour@m4x.org>
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
import scipy.sparse as sp
from libc.math cimport fabs
from libc.time cimport time, time_t

from ..utils.seq_dataset cimport SequentialDataset
from .sgd_fast cimport LossFunction

cdef extern from "sgd_fast_helpers.h":
    bint skl_isfinite(double) nogil

# This sparse implementation is taken from section 4.3 of "Minimizing Finite
# Sums with the Stochastic Average Gradient" by
# Mark Schmidt, Nicolas Le Roux, Francis Bach. 2013. <hal-00860051>
#
# https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf


cdef inline double fmax(double x, double y) nogil:
    if x > y:
        return x
    return y


def sag(SequentialDataset dataset,
        np.ndarray[double, ndim=1, mode='c'] weights_array,
        double intercept,
        int n_samples,
        int n_features,
        double tol,
        int max_iter,
        LossFunction loss,
        double step_size,
        double alpha,
        np.ndarray[double, ndim=1, mode='c'] sum_gradient_init,
        np.ndarray[double, ndim=1, mode='c'] gradient_memory_init,
        np.ndarray[bint, ndim=1, mode='c'] seen_init,
        int num_seen,
        bint fit_intercept,
        double intercept_sum_gradient,
        double intercept_decay,
        bint verbose):

    # true if the weights or intercept are NaN or infinity
    cdef bint infinity = False
    # the pointer to the coef_ or weights
    cdef double* weights = <double * >weights_array.data
    # the data pointer for X, the training set
    cdef double *x_data_ptr
    # the index pointer for the column of the data
    cdef int *x_ind_ptr
    # the label for the sample
    cdef double y
    # the prediction for y
    cdef double p
    # the sample weight
    cdef double sample_weight
    # the number of non-zero features for this sample
    cdef int xnnz
    # helper variable for indexes
    cdef int idx
    # the number of pass through all samples
    cdef int n_iter = 0
    # helper to track iterations through samples
    cdef int itr
    # the index (row number) of the current sample
    cdef int current_index
    # the maximum change in weights, used to compute stopping criterea
    cdef double max_change
    # a holder variable for the max weight, used to compute stopping criterea
    cdef double max_weight
    # the start time of the fit
    cdef time_t start_time
    # the end time of the fit
    cdef time_t end_time
    # precomputation since the step size does not change in this implementation
    cdef double wscale_update = 1.0 - step_size * alpha

    # vector of booleans indicating whether this sample has been seen
    cdef bint* seen = <bint*> seen_init.data
    # the sum of gradients for each feature
    cdef double* sum_gradient = <double*> sum_gradient_init.data
    # the previously seen gradient for each sample
    cdef double* gradient_memory = <double*> gradient_memory_init.data

    # the cumulative sums needed for JIT params
    cdef np.ndarray[double, ndim=1] cumulative_sums_array = \
        np.empty(n_samples, dtype=np.double, order="c")
    cdef double* cumulative_sums = <double*> cumulative_sums_array.data
    # the index for the last time this feature was updated
    cdef np.ndarray[int, ndim=1] feature_hist_array = \
        np.zeros(n_features, dtype=np.int32, order="c")
    cdef int* feature_hist = <int*> feature_hist_array.data
    # the previous weights to use to compute stopping criteria
    cdef np.ndarray[double, ndim=1] previous_weights_array = \
        np.zeros(n_features, dtype=np.double, order="c")
    cdef double* previous_weights = <double*> previous_weights_array.data

    # the scalar used for multiplying z
    cdef double wscale = 1.0

    # helpers
    cdef double val, cum_sum, gradient_change
    cdef int j

    # the cumulative sums for each iteration for the sparse implementation
    cumulative_sums[0] = 0.0

    with nogil:
        start_time = time(NULL)
        while True:
            for itr in range(n_samples):

                # extract a random sample
                current_index = dataset.random(&x_data_ptr,
                                               &x_ind_ptr,
                                               &xnnz,
                                               &y,
                                               &sample_weight)

                # update the number of samples seen and the seen array
                if seen[current_index] == 0:
                    num_seen += 1
                    seen[current_index] = 1

                # make the weight updates
                if itr > 0:
                    for j in range(xnnz):
                        idx = x_ind_ptr[j]

                        cum_sum = cumulative_sums[itr - 1]
                        if feature_hist[idx] != 0:
                            cum_sum -= cumulative_sums[feature_hist[idx] - 1]

                        weights[idx] -= cum_sum * sum_gradient[idx]

                        feature_hist[idx] = itr

                        # check to see that the weight is not inf or NaN
                        if not skl_isfinite(weights[idx]):
                            infinity = True
                            break

                # check to see that the intercept is not inf or NaN
                if infinity or not skl_isfinite(intercept):
                    infinity = True
                    break

                # find the current prediction, gradient
                p = (wscale * sparse_dot(x_data_ptr, x_ind_ptr, xnnz,
                                         weights)) + intercept
                gradient = loss._dloss(p, y) * sample_weight

                # make the updates to the sum of gradients
                gradient_change = gradient - gradient_memory[current_index]
                for j in range(xnnz):
                    idx = x_ind_ptr[j]
                    val = x_data_ptr[j]
                    sum_gradient[idx] += gradient_change * val

                if fit_intercept:
                    intercept_sum_gradient += gradient_change
                    intercept -= (step_size * intercept_sum_gradient
                                  / num_seen * intercept_decay)

                # update the gradient memory for this sample
                gradient_memory[current_index] = gradient

                # L2 regularization by simply rescaling the weights
                wscale *= wscale_update

                if itr == 0:
                    cumulative_sums[0] = step_size / (wscale * num_seen)
                else:
                    cumulative_sums[itr] = (cumulative_sums[itr - 1]
                                            + step_size / (wscale * num_seen))

                # if wscale gets too small, we need to reset the scale
                if wscale < 1e-9:
                    if verbose:
                        with gil:
                            print("rescaling...")
                    scale_weights(weights, wscale, n_features, n_samples,
                                  itr, cumulative_sums, feature_hist,
                                  sum_gradient)
                    wscale = 1.0

            if infinity:
                break

            # check if the stopping criteria is reached
            scale_weights(weights, wscale, n_features, n_samples, n_samples - 1,
                          cumulative_sums, feature_hist, sum_gradient)
            wscale = 1.0
            
            max_change = 0.0
            max_weight = 0.0
            for j in range(n_features):
                max_weight = fmax(max_weight, fabs(weights[j]))
            for j in range(n_features):
                max_change = fmax(max_change,
                                  fabs(weights[j] - previous_weights[j]))
                previous_weights[j] = weights[j]

            n_iter += 1
            if max_change / max_weight <= tol:
                if verbose:
                    end_time = time(NULL)
                    with gil:
                        print("convergence after %d epochs took %d seconds" %
                              (n_iter, end_time - start_time))
                break

            if n_iter >= max_iter:
                if verbose:
                    end_time = time(NULL)
                    with gil:
                        print(("max_iter reached after %d seconds") %
                              (end_time - start_time))
                break

    if infinity:
        raise ValueError(("Floating-point under-/overflow occurred at epoch"
                          " #%d. Lowering the step_size or scaling the input"
                          " data with StandardScaler or"
                          " MinMaxScaler might help.") % (n_iter + 1))


    return intercept, num_seen, n_iter, intercept_sum_gradient


cdef void scale_weights(double* weights, double wscale, int n_features,
                        int n_samples, int itr, double* cumulative_sums,
                        int* feature_hist, double* sum_gradient) nogil:
    cdef int j
    cdef double cum_sum
    for j in range(n_features):
        cum_sum = cumulative_sums[itr]
        if feature_hist[j] != 0:
            cum_sum -= cumulative_sums[feature_hist[j] - 1]

        weights[j] -= cum_sum * sum_gradient[j]

        weights[j] *= wscale
        feature_hist[j] = (itr + 1) % n_samples

    cumulative_sums[itr % n_samples] = 0.0


cdef double sparse_dot(double* x_data_ptr, int* x_ind_ptr, int xnnz,
                       double* w_data_ptr) nogil:
    """Compute dot product between sparse vector x and dense vector w."""
    cdef int j, idx
    cdef double innerprod = 0.0

    # only consider nonzero values of x
    for j in range(xnnz):
        idx = x_ind_ptr[j]
        innerprod += w_data_ptr[idx] * x_data_ptr[j]
    return innerprod


def get_max_squared_sum(X):
    """Maximum squared sum of X over samples. 

    Used in ``get_auto_step_size()``, for SAG solver.

    Parameter
    ---------
    X : {numpy array, scipy CSR sparse matrix}, shape (n_samples, n_features)
        Training vector. X must be in C order.

    Returns
    -------
    max_squared_sum : double
        Maximum squared sum of X over samples.
    """

    # CSR sparse matrix X
    cdef np.ndarray[double] X_data
    cdef np.ndarray[int] X_indptr
    cdef double *X_data_ptr
    cdef int *X_indptr_ptr
    cdef int offset

    # Dense numpy array X
    cdef np.ndarray[double, ndim=2] X_ndarray
    cdef int stride

    # Both cases
    cdef bint sparse = sp.issparse(X)
    cdef double *x_data_ptr
    cdef double max_squared_sum = 0.0
    cdef double current_squared_sum = 0.0
    cdef int n_samples = X.shape[0]
    cdef int nnz = X.shape[1]
    cdef int i, j
    cdef double val

    if sparse:
        X_data = X.data
        X_indptr = X.indptr
        X_data_ptr = <double *>X_data.data
        X_indptr_ptr = <int *>X_indptr.data
    else:
        X_ndarray = X
        stride = X_ndarray.strides[0] / X_ndarray.itemsize
        x_data_ptr = <double *>X_ndarray.data - stride

    with nogil:
        for i in range(n_samples):
            # find next sample data
            if sparse:
                offset = X_indptr_ptr[i]
                nnz = X_indptr_ptr[i + 1] - offset
                x_data_ptr = X_data_ptr + offset
            else:
                x_data_ptr += stride

            # sum of squared non-zero features
            for j in range(nnz):
                val = x_data_ptr[j]
                current_squared_sum += val * val

            if current_squared_sum > max_squared_sum:
                max_squared_sum = current_squared_sum
            current_squared_sum = 0.0


    if not skl_isfinite(max_squared_sum):
        raise ValueError("Floating-point under-/overflow occurred")

    return max_squared_sum
