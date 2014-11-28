# cython: cdivision=True

import numpy as np
cimport numpy as np
import warnings
from libc.math cimport fmax, fabs
from libc.time cimport time, time_t

cdef extern from "sgd_fast_helpers.h":
    bint skl_isfinite(double) nogil

from sklearn.utils.seq_dataset cimport SequentialDataset
from .sgd_fast cimport LossFunction, Classification

# This sparse implementation is taken from section 4.3 of "Minimizing Finite
# Sums with the Stochastic Average Gradient" by
# Mark Schmidt, Nicolas Le Roux, Francis Bach. 2013. <hal-00860051>
#
# https://hal.inria.fr/hal-00860051/PDF/sag_journal.pdf

def sag_sparse(SequentialDataset dataset,
               np.ndarray[double, ndim=1, mode='c'] weights_array,
               double intercept_init,
               int n_samples,
               int n_features,
               double tol,
               int max_iter,
               LossFunction loss,
               double eta,
               double alpha,
               np.ndarray[double, ndim=1, mode='c'] sum_gradient_init,
               np.ndarray[double, ndim=1, mode='c'] gradient_memory_init,
               np.ndarray[bint, ndim=1, mode='c'] seen_init,
               int num_seen_init,
               double weight_pos,
               double weight_neg,
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
    # the sample weight
    cdef double sample_weight
    # the number of non-zero features for this sample
    cdef int xnnz
    # helper variable for indexes
    cdef int idx
    # the total number of interations through the data
    cdef int total_iter = 0
    # the index (row number) of the current sample
    cdef int current_index
    # helper variable for the weight of a pos/neg class
    cdef double class_weight
    # the maximum change in weights, used to compute stopping criterea
    cdef double max_change
    # a holder variable for the max weight, used to compute stopping criterea
    cdef double max_weight
    # whether or not the max iter has been reached
    cdef bint max_iter_reached = False
    # the start time of the fit
    cdef time_t start_time
    # the end time of the fit
    cdef time_t end_time

    # the total number of samples seen
    cdef double num_seen = num_seen_init

    # vector of booleans indicating whether this sample has been seen
    cdef bint* seen = <bint*> seen_init.data

    # the sum of gradients for each feature
    cdef double* sum_gradient = <double*> sum_gradient_init.data

    # the previously seen gradient for each sample
    cdef double* gradient_memory = <double*> gradient_memory_init.data

    # the cumulative sums needed for JIT params
    cdef np.ndarray[double, ndim=1] cumulative_sums_array = \
        np.empty(max_iter * n_samples + 1,
                 dtype=np.double,
                 order="c")
    cdef double* cumulative_sums = <double*> cumulative_sums_array.data

    # the index for the last time this feature was updated
    cdef np.ndarray[int, ndim=1] feature_hist_array = \
        np.zeros(n_features,
                 dtype=np.int32,
                 order="c")
    cdef int* feature_hist = <int*> feature_hist_array.data

    # the previous weights to use to compute stopping criteria
    cdef np.ndarray[double, ndim=1] previous_weights_array = \
        np.zeros(n_features,
                 dtype=np.double,
                 order="c")
    cdef double* previous_weights = <double*> previous_weights_array.data

    # the scalar used for multiplying z
    cdef double wscale = 1.0

    # the cumulative sums for each iteration for the sparse implementation
    cumulative_sums[0] = 0.0

    cdef double intercept = intercept_init

    with nogil:
        start_time = time(NULL)
        while True:
            for k in range(n_samples):

                # extract a random sample
                current_index = dataset.random(&x_data_ptr,
                                               &x_ind_ptr,
                                               &xnnz,
                                               &y,
                                               &sample_weight)

                # update the number of samples seen and the seen array
                if seen[current_index] == 0:
                    num_seen += 1.0
                    seen[current_index] = 1

                # make the weight updates
                for j in range(xnnz):
                    idx = x_ind_ptr[j]
                    weights[idx] -= ((cumulative_sums[total_iter] -
                                      cumulative_sums[feature_hist[idx]]) *
                                     sum_gradient[idx])
                    feature_hist[idx] = total_iter

                    # check to see that the weight is not inf or NaN
                    if not skl_isfinite(weights[idx]):
                        infinity = True
                        break

                # check to see if we have already encountered a bad weight or
                # that the intercept is not inf or NaN
                if infinity or not skl_isfinite(intercept):
                    infinity = True
                    break

                # find the current prediction, gradient
                p = (wscale * dot(x_data_ptr, x_ind_ptr,
                     weights, xnnz)) + intercept
                gradient = loss._dloss(p, y)

                # find the class_weight
                if y > 0.0:
                    class_weight = weight_pos
                else:
                    class_weight = weight_neg

                gradient *= sample_weight * class_weight

                # make the updates to the sum of gradients
                for j in range(xnnz):
                    idx = x_ind_ptr[j]
                    val = x_data_ptr[j]
                    update = val * gradient
                    sum_gradient[idx] += (update -
                                          gradient_memory[current_index] * val)

                gradient_memory[current_index] = gradient

                intercept -= eta * gradient * intercept_decay

                wscale *= 1.0 - eta * alpha

                # if wscale gets too small, we need to reset the scale
                if wscale < 1e-9:
                    scale_weights(weights, wscale, n_features, total_iter,
                                  cumulative_sums, feature_hist, sum_gradient)
                    wscale = 1.0

                cumulative_sums[total_iter + 1] = (cumulative_sums[total_iter]
                                                   + eta / (wscale * num_seen))

                total_iter += 1

            # check if the stopping criteria is reached
            scale_weights(weights, wscale, n_features, total_iter,
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

            if max_change / max_weight <= tol:
                if verbose:
                    end_time = time(NULL)
                    with gil:
                        print("convergence after %d epochs took %d seconds" %
                              ((total_iter / n_samples),
                              (end_time - start_time)))
                break
            if total_iter / n_samples >= max_iter:
                if verbose:
                    end_time = time(NULL)
                    with gil:
                        print(("max_iter reached after %d seconds") %
                              (end_time - start_time))
                max_iter_reached = True
                break


    if infinity:
        raise ValueError(("Floating-point under-/overflow occurred at epoch"
                          " #%d. Lowering the eta0 or scaling the input data"
                          " with StandardScaler or"
                          " MinMaxScaler might help.") % (total_iter + 1))

    if max_iter_reached:
        warnings.warn("The max_iter was reached which means "
                      "the coef_ did not converge")

    return intercept, num_seen


cdef void scale_weights(double* weights, double wscale, int n_features,
                        int total_iter, double* cumulative_sums,
                        int* feature_hist, double* sum_gradient) nogil:
    for j in range(n_features):
        weights[j] -= ((cumulative_sums[total_iter] -
                            cumulative_sums[feature_hist[j]]) *
                           sum_gradient[j])
        weights[j] *= wscale

def get_auto_eta(SequentialDataset dataset, double alpha,
                 int n_samples, LossFunction loss):
    cdef double *x_data_ptr
    cdef int *x_ind_ptr
    cdef double y
    cdef double sample_weight
    cdef int xnnz
    cdef double max_squared_sum = 0.0
    cdef double current_squared_sum = 0.0

    with nogil:
        for i in range(n_samples):
            dataset.next(&x_data_ptr,
                         &x_ind_ptr,
                         &xnnz,
                         &y,
                         &sample_weight)
            for j in range(xnnz):
                val = x_data_ptr[j]
                current_squared_sum += val * val

            if current_squared_sum > max_squared_sum:
                max_squared_sum = current_squared_sum
            current_squared_sum = 0.0

    if isinstance(loss, Classification):
        # Lipschitz for log loss
        return 4.0 / (max_squared_sum + 4.0 * alpha)
    else:
        # Lipschitz for squared loss
        return 1.0 / (max_squared_sum + alpha)

cdef double dot(double* x_data_ptr, int* x_ind_ptr, double* w_data_ptr,
                int xnnz) nogil:
        cdef int j
        cdef int idx
        cdef double innerprod = 0.0

        for j in range(xnnz):
            idx = x_ind_ptr[j]
            innerprod += w_data_ptr[idx] * x_data_ptr[j]
        return innerprod

