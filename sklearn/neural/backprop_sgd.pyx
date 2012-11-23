cimport cython
from libc.math cimport exp, fabs, log, sqrt
from libc.stdint cimport uint64_t

cimport numpy as np
import numpy as np

from ..utils.extmath import logsumexp, safe_sparse_dot

np.import_array()


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
cdef void logistic(double[:] z) nogil:
    """Logistic sigmoid: 1. / (1. + np.exp(-x)), computed in-place."""
    #cdef np.ndarray[np.float64_t, ndim=1] x = z.ravel()
    for i in xrange(z.shape[0]):
        z[i] = 1. / (1. + exp(-z[i]))


@cython.boundscheck(False)
cdef inline double _logsumexp(double[:] x) nogil:
    """Fast implementation of 1-d logsumexp for double precision floats.

    Computes log(sum(exp(x))) and stores the result in out.
    """
    # XXX this might be interesting for sklearn.utils.extmath,
    # and maybe even scipy.misc

    cdef double total, vmax

    vmax = x[0]
    for i in range(1, x.shape[0]):
        vmax = max(x[i], vmax)
    total = 0.
    for i in range(0, x.shape[0]):
        total += exp(x[i] - vmax)
    return log(total) + vmax


@cython.boundscheck(False)
cdef void log_softmax(double[:, :] X, double[:, :] log_y_output) nogil:
    """Logistic K-way softmax (exp(X).T / exp(X).sum(axis=1)).T, in log domain.

    In-place implementation.
    """
    cdef double logsum

    for i in range(X.shape[0]):
        logsum = _logsumexp(X[i, :])
        for j in range(X.shape[1]):
            log_y_output[i, j] = X[i, j] - logsum


@cython.boundscheck(False)
@cython.cdivision(True)
cdef double log_loss_multiclass(double[:, :] output, uint64_t[:, :] T) nogil:
    """Multiclass cross-entropy loss function."""
    #loss = -np.mean(np.sum(log_y_output * T, axis=1))
    cdef double total

    log_softmax(output, output)

    total = 0.
    for i in range(output.shape[0]):
        for j in range(output.shape[1]):
            total += output[i, j] * T[i, j]

    return -(total / output.shape[0])


@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def backprop_sgd(self, X, np.ndarray[np.int64_t, ndim=2] Y):
    cdef double alpha, alpha2
    cdef double loss, penalty, prev_loss, tol
    cdef double lr, momentum
    cdef double sqrt_n_features, sqrt_n_hidden
    cdef bint multiclass, verbose
    cdef np.int32_t i, n_features, n_hidden, n_samples, n_targets

    cdef np.ndarray[np.float64_t, ndim=2] \
        w_hidden, w_output, output_grad, \
        v_hidden, v_output, log_y_output, y_output, z_output
    # XXX fix mode on these
    cdef np.ndarray[np.float64_t, ndim=2] hidden_grad
    cdef np.ndarray[np.float64_t, ndim=1] \
        bias_hidden, bias_output, bias_hidden_grad, bias_output_grad, \
        v_bias_hidden, v_bias_output

    n_samples, n_features = X.shape
    n_hidden = self.n_hidden
    n_targets = Y.shape[1]

    multiclass = n_targets > 1 and not self._lbin.multilabel

    w_hidden = self.coef_hidden_
    w_output = self.coef_output_
    bias_hidden = self.intercept_hidden_
    bias_output = self.intercept_output_

    alpha = self.alpha
    alpha2 = .5 * alpha
    lr = self.learning_rate

    momentum = self.momentum
    prev_loss = np.inf
    tol = self.tol
    verbose = self.verbose

    # scale the learning rate by the number of features for faster convergence
    sqrt_n_features = sqrt(n_features)
    sqrt_n_hidden = sqrt(n_hidden)
    lr_hidden = lr * sqrt_n_features
    lr_output = lr * sqrt_n_hidden

    # set up output arrays to prevent allocation inside the loop
    z_output = np.empty((n_samples, n_targets), dtype=np.float64)
    y_output = np.empty((n_samples, n_targets), dtype=np.float64)
    output_grad = np.empty((n_targets, n_hidden), dtype=np.float64)

    # velocities for momentum method
    v_hidden = np.zeros(self.coef_hidden_.shape)
    v_output = np.zeros(self.coef_output_.shape)
    v_bias_hidden = np.zeros(self.intercept_hidden_.shape)
    v_bias_output = np.zeros(self.intercept_output_.shape)

    for epoch in xrange(self.max_iter):
        # make predictions
        y_hidden = safe_sparse_dot(X, w_hidden.T) + bias_hidden
        logistic(y_hidden)
        #z_output = np.dot(y_hidden, w_output.T) + bias_output
        np.dot(y_hidden, w_output.T, out=z_output)

        for i in xrange(n_samples):
            for j in xrange(n_targets):
                z_output[i, j] += bias_output[j]
        if multiclass:
            #log_y_output = log_softmax(z_output)
            log_softmax(z_output, log_y_output)
        else:
            logistic(z_output)
            log_y_output = z_output
        for i in xrange(n_samples):
            for j in xrange(n_targets):
                y_output[i, j] = exp(log_y_output[i, j])

        delta_output = y_output - Y
        delta_hidden = (np.dot(delta_output, w_output)
                       * (y_hidden * (1 - y_hidden)))

        #output_grad = np.dot(delta_output.T, y_hidden)
        np.dot(delta_output.T, y_hidden, out=output_grad)
        for i in xrange(n_targets):
            for j in xrange(n_hidden):
                output_grad[i, j] /= n_samples
        hidden_grad = safe_sparse_dot(delta_hidden.T, X)
        for i in xrange(n_hidden):
            for j in xrange(n_features):
                hidden_grad[i, j] /= n_samples

        if alpha > 0:
            hidden_grad += alpha * w_hidden
            output_grad += alpha * w_output

        bias_hidden_grad = delta_hidden.sum(axis=0) + alpha * bias_hidden
        bias_hidden_grad /= n_samples
        bias_output_grad = delta_output.sum(axis=0) + alpha * bias_output
        bias_output_grad /= n_samples

        for i in xrange(n_hidden):
            v_bias_hidden[i] *= momentum
            v_bias_hidden[i] -= lr * bias_hidden_grad[i]
            bias_hidden[i] += v_bias_hidden[i]
            for j in xrange(n_features):
                v_hidden[i, j] *= momentum
                v_hidden[i, j] -= lr_hidden * hidden_grad[i, j]
                w_hidden[i, j] += v_hidden[i, j]
        for i in xrange(n_targets):
            v_bias_output[i] *= momentum
            v_bias_output[i] -= lr * bias_output_grad[i]
            bias_output[i] += v_bias_output[i]
            for j in xrange(n_hidden):
                v_output[i, j] *= momentum
                v_output[i, j] -= lr_output * output_grad[i, j]
                w_output[i, j] += v_output[i, j]

        if n_targets == 1:
            loss = -np.mean(Y * log_y_output +
                            (1 - Y) * np.log(1 - y_output))
        else:
            loss = -np.mean(np.sum(log_y_output * Y, axis=1))

        # penalize loss
        penalty = 0
        if alpha2 > 0:
            for i in xrange(n_hidden):
                penalty += bias_hidden[i]
                for j in xrange(n_features):
                    penalty += w_hidden[i, j] ** 2
            for i in xrange(n_targets):
                penalty += bias_output[i]
                for j in xrange(n_hidden):
                    penalty += w_output[i, j] ** 2
        loss += alpha2 * penalty / n_samples

        if verbose:
            print("iteration %d, avg. loss = %.5f" % (epoch + 1, loss))
        if fabs(loss - prev_loss) < tol:
            break
        prev_loss = loss

    return
