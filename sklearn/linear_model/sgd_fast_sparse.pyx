# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.

import numpy as np
import sys
from time import time

cimport numpy as np
cimport cython
cimport sgd_fast

from sgd_fast cimport LossFunction, log, sqrt, pow

# Penalty constants
DEF L1 = 1
DEF L2 = 2
DEF ELASTICNET = 3

# Learning rate constants
DEF CONSTANT = 1
DEF OPTIMAL = 2
DEF INVSCALING = 3


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def plain_sgd(np.ndarray[double, ndim=1] w,
              double intercept,
              LossFunction loss,
              int penalty_type,
              double alpha, double rho,
              np.ndarray[double, ndim=1] X_data,
              np.ndarray[int, ndim=1] X_indices,
              np.ndarray[int, ndim=1] X_indptr,
              np.ndarray[double, ndim=1] Y,
              int n_iter, int fit_intercept,
              int verbose, int shuffle, int seed,
              double weight_pos, double weight_neg,
              np.ndarray[double, ndim=1] sample_weight,
              int learning_rate, double eta0,
              double power_t):
    """Cython impl. of SGD with different loss functions and penalties

    This representation assumes X represented using the Compressed Sparse Row
    representation of scipy.sparse.

    Parameters
    ----------
    w : ndarray[double, ndim=1]
        The allocated coef_ vector.
    intercept : double
        The initial intercept
    loss : LossFunction
        A concrete LossFunction object.
    penalty_type : int
        The penalty 2 for L2, 1 for L1, and 3 for Elastic-Net.
    alpha : float
        The regularization parameter.
    rho : float
        The elastic net hyperparameter.
    X : csr_matrix[double, ndim=2]
        The dataset as a Compressed Sparse Row matrix
        (see scipy.sparse.csr_matrix).
    Y : ndarray[double, ndim=1]
        The labels.
    n_iter : int
        The number of iterations (epochs).
    fit_intercept : int
        Whether or not to fit the intercept (1 or 0).
    verbose : int
        Print verbose output; 0 for quite.
    shuffle : int
        Whether to shuffle the training data before each epoch.
    weight_pos : float
        The weight of the positive class.
    weight_neg : float
        The weight of the negative class.
    seed : int
        The seed of the pseudo random number generator to use when
        shuffling the data
    sample_weight : array, shape = [n_samples]
        The importance weight of each sample.
    learning_rate : int
        The learning rate:
        (1) constant, eta = eta0
        (2) optimal, eta = 1.0/(t+t0)
        (3) inverse scaling, eta = eta0 / pow(t, power_t)
    eta0 : double
        The initial learning rate.
    power_t : double
        The exponent for inverse scaling learning rate.

    Returns
    -------
    w : array, shape [n_features]
        The fitted weight vector.
    intercept : float
        The fitted intercept term.
    """
    # get the data information into easy vars
    cdef unsigned int n_samples = Y.shape[0]
    cdef unsigned int n_features = w.shape[0]

    cdef double *w_data_ptr = <double *>w.data
    cdef double *X_data_ptr = <double *>X_data.data
    cdef int *X_indptr_ptr = <int *>X_indptr.data
    cdef int *X_indices_ptr = <int *>X_indices.data
    cdef double *Y_data_ptr = <double *>Y.data

    cdef double *sample_weight_data = <double *>sample_weight.data

    cdef np.ndarray[int, ndim=1, mode="c"] index = np.arange(n_samples,
                                                             dtype=np.int32)
    cdef int *index_data_ptr = <int *>index.data
    cdef int offset = 0
    cdef int xnnz = 0
    cdef double wscale = 1.0
    cdef double eta = 0.0
    cdef double p = 0.0
    cdef double update = 0.0
    cdef double sumloss = 0.0
    cdef double wnorm = 0.0
    cdef double t = 0.0
    cdef double y = 0.0
    cdef double class_weight = 1.0
    cdef unsigned int count = 0
    cdef unsigned int epoch = 0
    cdef unsigned int i = 0
    cdef int sample_idx = 0

    # q vector is only used for L1 regularization
    cdef np.ndarray[double, ndim=1, mode="c"] q = None
    cdef double *q_data_ptr = NULL
    if penalty_type != L2:
        q = np.zeros((n_features,), dtype=np.float64, order="c")
        q_data_ptr = <double *> q.data
    cdef double u = 0.0

    if penalty_type == L2:
        rho = 1.0
    elif penalty_type == L1:
        rho = 0.0
    
    cdef double typw = sqrt(1.0 / sqrt(alpha))

    if learning_rate == OPTIMAL:
        # computing eta0, the initial learning rate
        eta0 = typw / max(1.0, loss.dloss(-typw, 1.0))
    else:
        eta = eta0

    if learning_rate == OPTIMAL:
        # initialize t such that eta at first example equals eta0
        t = 1.0 / (eta0 * alpha)
    else:
        t = 1.0

    t_start = time()
    for epoch in xrange(n_iter):
        if verbose > 0:
            print("-- Epoch %d" % (epoch + 1))
        if shuffle:
            np.random.RandomState(seed).shuffle(index)
        for i in xrange(n_samples):
            sample_idx = index_data_ptr[i]
            offset = X_indptr_ptr[sample_idx]
            xnnz = X_indptr_ptr[sample_idx + 1] - offset
            y = Y_data_ptr[sample_idx]
            if learning_rate == OPTIMAL:
                eta = 1.0 / (alpha * t)
            elif learning_rate == INVSCALING:
                eta = eta0 / pow(t, power_t)
            p = (dot(w_data_ptr, X_data_ptr, X_indices_ptr,
                     offset, xnnz) * wscale) + intercept
            sumloss += loss.loss(p, y)
            if y > 0:
                class_weight = weight_pos
            else:
                class_weight = weight_neg
            update = eta * loss.dloss(p, y) * class_weight * \
                sample_weight_data[sample_idx]
            if update != 0.0:
                add(w_data_ptr, wscale, X_data_ptr, X_indices_ptr,
                    offset, xnnz, -update)
                if fit_intercept == 1:
                    intercept -= update * 0.01
            if penalty_type != L1:
                wscale *= (1.0 - (rho * eta * alpha))
                if wscale < 1e-9:
                    w *= wscale
                    wscale = 1.0
            if penalty_type == L1 or penalty_type == ELASTICNET:
                u += ((1.0 - rho) * eta * alpha)
                l1penalty(w_data_ptr, wscale, q_data_ptr,
                          X_indices_ptr, offset, xnnz, u)
            t += 1
            count += 1

        # report epoche information
        if verbose > 0:
            wnorm = sqrt(np.dot(w, w) * wscale * wscale)
            print("Norm: %.2f, NNZs: %d, "\
            "Bias: %.6f, T: %d, Avg. loss: %.6f" % (wnorm,
                                                    w.nonzero()[0].shape[0],
                                                    intercept, count,
                                                    sumloss / count))
            print("Total training time: %.2f seconds." % (time() - t_start))

        # floating-point under-/overflow check.
        if np.any(np.isinf(w)) or np.any(np.isnan(w)) \
           or np.isnan(intercept) or np.isinf(intercept):
            raise ValueError("floating-point under-/overflow occured.")

    w *= wscale
    return w, intercept


cdef inline double max(double a, double b):
    return a if a >= b else b


cdef inline double min(double a, double b):
    return a if a <= b else b


cdef double dot(double *w_data_ptr, double *X_data_ptr, int *X_indices_ptr,
                int offset, int xnnz):
    cdef double sum = 0.0
    cdef int j
    for j in xrange(xnnz):
        sum += w_data_ptr[X_indices_ptr[offset + j]] * X_data_ptr[offset + j]
    return sum


cdef double add(double *w_data_ptr, double wscale, double *X_data_ptr,
                int *X_indices_ptr, int offset, int xnnz, double c):
    """Scales example x by constant c and adds it to the weight vector w"""
    cdef int j
    cdef int idx
    cdef double val
    cdef double innerprod = 0.0
    cdef double xsqnorm = 0.0
    for j in xrange(xnnz):
        idx = X_indices_ptr[offset + j]
        val = X_data_ptr[offset + j]
        innerprod += (w_data_ptr[idx] * val)
        xsqnorm += (val * val)
        w_data_ptr[idx] += val * (c / wscale)
    return (xsqnorm * c * c) + (2.0 * innerprod * wscale * c)


cdef void l1penalty(double *w_data_ptr, double wscale, double *q_data_ptr,
                    int *X_indices_ptr, int offset, int xnnz, double u):
    """Apply the L1 penalty to each updated feature

    This implements the truncated gradient approach by
    [Tsuruoka, Y., Tsujii, J., and Ananiadou, S., 2009].
    """
    cdef double z = 0.0
    cdef int j = 0
    cdef int idx = 0
    for j in xrange(xnnz):
        idx = X_indices_ptr[offset + j]
        z = w_data_ptr[idx]
        if (wscale * w_data_ptr[idx]) > 0.0:
            w_data_ptr[idx] = max(
                0.0, w_data_ptr[idx] - ((u + q_data_ptr[idx]) / wscale))

        elif (wscale * w_data_ptr[idx]) < 0.0:
            w_data_ptr[idx] = min(
                0.0, w_data_ptr[idx] + ((u - q_data_ptr[idx]) / wscale))

        q_data_ptr[idx] += (wscale * (w_data_ptr[idx] - z))
