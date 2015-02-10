# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>,
#          Artem Sobolev
# License: BSD 3 clause

from __future__ import print_function

import cython
import numpy as np
cimport numpy as np


np.import_array()

from libc.math cimport sqrt

def factorize_matrix_sgd(np.ndarray[double, ndim=1] data,
                         np.ndarray[int, ndim=1] row,
                         np.ndarray[int, ndim=1] col,
                         np.ndarray[double, ndim=2] L,
                         np.ndarray[double, ndim=2] R,
                         np.ndarray[double, ndim=1] bias_samples,
                         np.ndarray[double, ndim=1] bias_features,
                         np.npy_intp n_samples,
                         np.npy_intp n_features,
                         np.npy_intp n_components,
                         int n_iter,
                         int fit_intercept,
                         double regularization,
                         double learning_rate,
                         object random_state,
                         int verbose,
                         int adagrad):

    cdef np.npy_intp i, j, k
    cdef double value, estimation, error

    if adagrad:
        feature_grad_norms = np.zeros(n_features)
        sample_grad_norms = np.zeros(n_samples)

        feature_bias_grad_norms = np.zeros(n_features)
        sample_bias_grad_norms = np.zeros(n_samples)

    indices = range(data.shape[0])

    for current_iter in range(n_iter):
        random_state.shuffle(indices)

        for value_index in indices:
            value = data[value_index]
            i = row[value_index]
            j = col[value_index]
            
            estimation = bias_samples[i] + bias_features[j] + np.dot(L[i, :], R[:, j])
            error = estimation - value

            L_grad = error * R[:, j] + regularization * L[i, :]
            R_grad = error *  L[i, :] + regularization * R[:, j]

            r_discount = 1
            l_discount = 1

            if adagrad:
                feature_grad_norms[j] += np.dot(R_grad, R_grad)
                sample_grad_norms[i] += np.dot(L_grad, L_grad)

                r_discount = sqrt(feature_grad_norms[j]) + 1e-6
                l_discount = sqrt(sample_grad_norms[i]) + 1e-6

            R[:, j] -= learning_rate / r_discount * R_grad
            L[i, :] -= learning_rate / l_discount * L_grad

            if fit_intercept:
                s_discount = 1
                f_discount = 1

                if adagrad:
                    feature_bias_grad_norms[j] += error ** 2
                    sample_bias_grad_norms[i] += error ** 2

                    s_discount = sqrt(sample_bias_grad_norms[i]) + 1e-6
                    f_discount = sqrt(feature_bias_grad_norms[j]) + 1e-6

                bias_samples[i] -= learning_rate / s_discount * error
                bias_features[j] -= learning_rate / f_discount * error


        if verbose > 0:
            _log_step(L, R, bias_features, bias_samples,
                      row, col, data, n_iter, current_iter,
                      regularization)


cdef rr1(double* bias,
         np.ndarray[double, ndim=1] w,
         np.ndarray[double, ndim=1] b,
         np.ndarray[double, ndim=2] X,
         np.ndarray[double, ndim=1] y,
         double lam, int n_iter,
         int fit_intercept):

    cdef int n_samples = X.shape[0]
    cdef int n_components = X.shape[1]

    errors = y - np.dot(X, w) - b - bias[0]

    for current_iter in range(n_iter):
        for k in range(n_components):
            data = X[:, k]
            errors += w[k] * data
            w[k] = 0

            a = np.dot(data, data)
            d = np.dot(data, errors)

            w[k] = d / (lam + a)
            errors -= w[k] * data

        if fit_intercept:
            bias[0] = (y - np.dot(X, w) - b).mean()


cdef rr(double* bias,
         np.ndarray[double, ndim=1] w,
         np.ndarray[double, ndim=1] b,
         np.ndarray[double, ndim=2] X,
         np.ndarray[double, ndim=1] y,
         double lam,
         int n_iter,
         int fit_intercept):

    cdef int n_samples = X.shape[0]
    cdef int n_components = X.shape[1]

    A = np.linalg.inv(np.dot(X.T, X) + n_samples * lam * np.eye(n_components))
    D = np.dot(A, X.T)

    w *= 0
    w += np.dot(D, y - bias[0] - b)

    if fit_intercept:
        bias[0] = (y - np.dot(X, w) - b).mean()

cdef _rr(int rr_algo,
        double* bias,
        np.ndarray[double, ndim=1] w,
        np.ndarray[double, ndim=1] b,
        np.ndarray[double, ndim=2] X,
        np.ndarray[double, ndim=1] y,
        double lam, int n_iter,
        int fit_intercept):

    if rr_algo == 1:
        rr1(bias, w, b, X, y, lam, n_iter, fit_intercept)
    else:
        rr(bias, w, b, X, y, lam, n_iter, fit_intercept)


cdef _log_step(L, R, bias_features, bias_samples, row, col, data, n_iter, current_iter, regularization):
    output_step = max(1, n_iter / 100)

    if current_iter % output_step == 0 or current_iter == n_iter - 1:
        reg = regularization * ((L * L).sum() + (R * R).sum())
        loss = _rmse_loss(data, row, col, L, R, bias_samples, bias_features) + reg
        print("Epoch %d / %d. Loss: %.5f" % (current_iter + 1, n_iter, loss))

def factorize_matrix_als(np.ndarray[double, ndim=1] data,
                           np.ndarray[int, ndim=1] row,
                           np.ndarray[int, ndim=1] col,
                           np.ndarray[double, ndim=2] L,
                           np.ndarray[double, ndim=2] R,
                           np.ndarray[double, ndim=1] bias_samples,
                           np.ndarray[double, ndim=1] bias_features,
                           np.npy_intp n_samples,
                           np.npy_intp n_features,
                           np.npy_intp n_components,
                           int n_iter,
                           int fit_intercept,
                           double regularization,
                           int verbose,
                           int rr_algo):

    cdef np.npy_intp i, j, k
    cdef double value, estimation, error, newL

    cdef double loss = _rmse(data, row, col, L, R)

    # List of samples (features) that share a value with
    # the given feature (sample)
    # TODO: std::map?
    cdef dict feature_to_samples = {}
    cdef dict sample_to_features = {}

    for idx in range(data.shape[0]):
        i = row[idx]
        j = col[idx]

        if i not in sample_to_features:
            sample_to_features[i] = ([], [])

        if j not in feature_to_samples:
            feature_to_samples[j] = ([], [])

        samples, indices = feature_to_samples[j]
        samples.append(i)
        indices.append(idx)

        features, indices = sample_to_features[i]
        features.append(j)
        indices.append(idx)


    cdef int output_step = max(1, n_iter / 100)

    for current_iter in range(n_iter):
        for feature in range(n_features):
            rel_samples, indices = feature_to_samples[feature]
            _rr(rr_algo, &(bias_features[feature]), R[:, feature],
                bias_samples[rel_samples], L[rel_samples, :],
                data[indices], regularization, 1, fit_intercept)

        for sample in range(n_samples):
            rel_features, indices = sample_to_features[sample]
            _rr(rr_algo, &(bias_samples[sample]), L[sample, :],
                bias_features[rel_features], R[:, rel_features].T,
                data[indices], regularization, 1, fit_intercept)

        if verbose > 0:
            _log_step(L, R, bias_features, bias_samples,
                      row, col, data, n_iter, current_iter,
                      regularization)


cdef _rmse_loss(np.ndarray[double, ndim=1] data,
          np.ndarray[int, ndim=1] row,
          np.ndarray[int, ndim=1] col,
          np.ndarray[double, ndim=2] L,
          np.ndarray[double, ndim=2] R,
          np.ndarray[double, ndim=1] bias_samples,
          np.ndarray[double, ndim=1] bias_features):

    cdef np.npy_intp i, j, k
    cdef double value, estimation, error
    cdef total_error = 0
    for value_index in range(data.shape[0]):
        value = data[value_index]
        i = row[value_index]
        j = col[value_index]

        estimation = bias_samples[i] + bias_features[j]
        for k in range(L.shape[1]):
            estimation += L[i, k] * R[k, j]
        error = value - estimation

        total_error += error * error

    return sqrt(total_error / data.shape[0])


def _rmse(np.ndarray[double, ndim=1] data,
          np.ndarray[int, ndim=1] row,
          np.ndarray[int, ndim=1] col,
          np.ndarray[double, ndim=2] L,
          np.ndarray[double, ndim=2] R):

    # FIXME: reuse _rmse_loss
    cdef np.npy_intp i, j, k
    cdef double value, estimation, error
    cdef total_error = 0
    for value_index in range(data.shape[0]):
        value = data[value_index]
        i = row[value_index]
        j = col[value_index]
        
        estimation = 0
        for k in range(L.shape[1]):
            estimation += L[i, k] * R[k, j]
        error = value - estimation
        
        total_error += error * error
        
    return sqrt(total_error / data.shape[0])
