# encoding: utf-8
# filename: mlp_fast.pyx
# cython: profile=False
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
import numpy as np
cimport numpy as np
cimport cython

import warnings

from ..utils import shuffle

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
ctypedef np.uint64_t uint64_t

cdef extern from "math.h":
    DTYPE_t log "log"(DTYPE_t) nogil
    DTYPE_t exp "exp"(DTYPE_t) nogil

cdef extern from "cblas.h":
    enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
        AtlasConj=114

    void daxpy "cblas_daxpy"(int N, double alpha, double *X, int incX,
                             double *Y, int incY)
    double ddot "cblas_ddot"(int N, double *X, int incX, double *Y, int incY)
    void dger "cblas_dger"(CBLAS_ORDER Order, int M, int N, double alpha,
                double *X, int incX, double *Y, int incY, double *A, int lda)
    void dgemm "cblas_dgemm"(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                    CBLAS_TRANSPOSE TransB, int M, int N,
                    int K, double alpha, double *A,
                    int lda, double *B, int ldb,
                    double beta, double *C, int ldc)
    void dgemv "cblas_dgemv"(CBLAS_ORDER Order,
                      CBLAS_TRANSPOSE TransA, int M, int N,
                      double alpha, double *A, int lda,
                      double *X, int incX, double beta,
                      double *Y, int incY)
    double dnrm2 "cblas_dnrm2"(int N, double *X, int incX)
    void dcopy "cblas_dcopy"(int N, double *X, int incX, double *Y, int incY)
    void dscal "cblas_dscal"(int N, double alpha, double *X, int incX)



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


cdef void log_softmax(double[:, :] X, double[:, :] log_y_output) nogil:
    """Logistic K-way softmax (exp(X).T / exp(X).sum(axis=1)).T, in log domain.

    In-place implementation.
    """
    cdef double logsum

    for i in range(X.shape[0]):
        logsum = _logsumexp(X[i, :])
        for j in range(X.shape[1]):
            log_y_output[i, j] = X[i, j] - logsum


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


cdef class LossFunction:
    """Base class for loss functions"""

    def loss(self,
             np.ndarray[DTYPE_t] y,
             int batch_start, int batch_end,
             np.ndarray[DTYPE_t] p,
             np.ndarray[DTYPE_t, ndim=2] out):
        """Evaluate the loss function.

        Parameters
        ----------
        y : np.ndarray, shape = [n_samples]
            The true value (aka target).
        p : np.ndarray, shape = [n_samples]
            The prediction.

        Returns
        -------

        np.ndarray, shape = [n_samples]
            The loss evaluated at `p` and `y`.
        """
        raise NotImplementedError()

    def dloss(self,
              np.ndarray[DTYPE_t, ndim=2] y,
              int batch_start, int batch_end,
              np.ndarray[DTYPE_t, ndim=2] p,
              np.ndarray[DTYPE_t, ndim=2] out):
        """Evaluate the derivative of the loss function with respect to
        the prediction `p`.

        Parameters
        ----------
        y : np.ndarray, shape = [n_samples]
            The true value (aka target).
        p : np.ndarray, shape = [n_samples]
            The prediction.

        Returns
        -------

        np.ndarray, shape = [n_samples]
            The derivative of the loss function at `p` and `y`.
        """
        raise NotImplementedError()


cdef class SquaredLoss(LossFunction):
    """Squared loss function."""

    def loss(self,
             np.ndarray[DTYPE_t, ndim=2] y,
             int batch_start, int batch_end,
             np.ndarray[DTYPE_t, ndim=2] p,
             np.ndarray[DTYPE_t, ndim=2] out):
        #out[:] = 0.5 * (y - p) * (y - p)
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = 0.5 * (yval - p_ptr[j]) * (yval - p_ptr[j])
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]

    def dloss(self,
              np.ndarray[DTYPE_t, ndim=2] y,
              int batch_start, int batch_end,
              np.ndarray[DTYPE_t, ndim=2] p,
              np.ndarray[DTYPE_t, ndim=2] out):
        #np.subtract(y, p, out)
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = (yval - p_ptr[j])
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]

    def __reduce__(self):
        return SquaredLoss, ()


cdef class CrossEntropyLoss(LossFunction):
    """Cross entropy loss function"""
    def loss(self,
             np.ndarray[DTYPE_t, ndim=2] y,
             int batch_start, int batch_end,
             np.ndarray[DTYPE_t, ndim=2] p,
             np.ndarray[DTYPE_t, ndim=2] out):
        #out[:] = -y * np.log(p) - (1 - y) * np.log(1 - p)
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = - yval * log(p_ptr[j]) - (1 - yval) * log(1 - p_ptr[j])
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]

    def dloss(self,
              np.ndarray[DTYPE_t, ndim=2] y,
             int batch_start, int batch_end,
              np.ndarray[DTYPE_t, ndim=2] p,
              np.ndarray[DTYPE_t, ndim=2] out):
        #out[:] = y - p
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = yval - p_ptr[j]
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]


cdef class MultiCrossEntropyLoss(LossFunction):
    """Multinomial cross entropy loss function."""

    def loss(self,
             np.ndarray[DTYPE_t, ndim=2] y,
             int batch_start, int batch_end,
             np.ndarray[DTYPE_t, ndim=2] p,
             np.ndarray[DTYPE_t, ndim=2] out):
        #out[:] = -y * np.log(p)
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = -yval * log(p_ptr[j])
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]

    def dloss(self,
              np.ndarray[DTYPE_t, ndim=2] y,
              int batch_start, int batch_end,
              np.ndarray[DTYPE_t, ndim=2] p,
              np.ndarray[DTYPE_t, ndim=2] out):
        #out[:] = y - p
        cdef DTYPE_t *y_ptr, *p_ptr, *o_ptr
        cdef DTYPE_t  yval
        cdef int i, j
        y_ptr = <DTYPE_t*>y.data + batch_start * y.shape[1]
        p_ptr = <DTYPE_t*>p.data
        o_ptr = <DTYPE_t*>out.data
        for i in xrange(0, p.shape[0]):
            for j in xrange(0, p.shape[1]):
                yval = y_ptr[j]
                o_ptr[j] = yval - p_ptr[j]
            y_ptr += y.shape[1]
            p_ptr += y.shape[1]
            o_ptr += y.shape[1]


cdef class OutputFunction:
    """Base class for ouput functions"""

    def output(self,
               np.ndarray[DTYPE_t, ndim=2] x,
               np.ndarray[DTYPE_t, ndim=2] out):
        raise NotImplementedError()

    def doutput(self,
                np.ndarray[DTYPE_t, ndim=2] x,
                np.ndarray[DTYPE_t, ndim=2] out):
        raise NotImplementedError()


cdef class Tanh(OutputFunction):

    def output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        np.tanh(x, out)

    def doutput(self,
                np.ndarray[DTYPE_t, ndim=2] x,
                np.ndarray[DTYPE_t, ndim=2] out):
        np.multiply(-x, x, out)
        out += 1

cdef class SoftMax(OutputFunction):

    def output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        # TODO: get rid of this allocation
        log_loss_multiclass(x, out)
        #r = np.logaddexp.reduce(x, axis=1)[:, np.newaxis]
        #np.subtract(x, r, out)
        #np.exp(out, out)


cdef class LogSig(OutputFunction):

    def output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        np.exp(-x, out)
        out += 1
        np.reciprocal(out, out)


cpdef forward(np.ndarray[DTYPE_t, ndim=2] X,
             int batch_start, int batch_end,
             np.ndarray[DTYPE_t, ndim=2] weights_hidden,
             np.ndarray[DTYPE_t, ndim=1] bias_hidden,
             np.ndarray[DTYPE_t, ndim=2] weights_output,
             np.ndarray[DTYPE_t, ndim=1] bias_output,
             np.ndarray[DTYPE_t, ndim=2] x_hidden,
             np.ndarray[DTYPE_t, ndim=2] x_output,
             OutputFunction output,
             OutputFunction hidden):

    cdef int batch_size = batch_end - batch_start
    # Hidden layer
    if batch_size == 1:
        dgemv(
            CblasRowMajor,
            CblasTrans,
            weights_hidden.shape[0], weights_hidden.shape[1],  # M, N
            1.0,  # alpha
            <DTYPE_t*>weights_hidden.data, weights_hidden.shape[1],
            <DTYPE_t*>X.data + batch_start * X.shape[1],
            1, 0.0,  # incX, beta
            <DTYPE_t*>x_hidden.data, 1)  #incY
    else:
        dgemm(
            CblasRowMajor,
            CblasNoTrans, CblasNoTrans,
            batch_size, weights_hidden.shape[1], X.shape[1],   # M, N, K
            1.0,
            <DTYPE_t*>X.data + batch_start*X.shape[1], X.shape[1],
            <DTYPE_t*>weights_hidden.data, weights_hidden.shape[1],
            0.0,
            <DTYPE_t*>x_hidden.data, x_hidden.shape[1])
    x_hidden += bias_hidden
    hidden.output(x_hidden, x_hidden)

    # Output layer
    #np.dot(x_hidden, weights_output, x_output)
    if batch_size == 1:
        dgemv(
            CblasRowMajor,
            CblasTrans,
            weights_output.shape[0], weights_output.shape[1],  # M, N
            1.0,  # alpha
            <DTYPE_t*>weights_output.data, weights_output.shape[1],
            <DTYPE_t*>(x_hidden.data),
            1, 0.0,  # incX, beta
            <DTYPE_t*>x_output.data, 1)  #incY
    else:
        dgemm(
            CblasRowMajor,
            CblasNoTrans, CblasNoTrans,
            batch_size, weights_output.shape[1], x_hidden.shape[1],  # M, N, K
            1.0,
            <DTYPE_t*>x_hidden.data, x_hidden.shape[1],
            <DTYPE_t*>weights_output.data, weights_output.shape[1],
            0.0, # beta
            <DTYPE_t*>x_output.data, x_output.shape[1])
    x_output += bias_output
    output.output(x_output, x_output)


cpdef backward(np.ndarray[DTYPE_t, ndim=2] X,
              int batch_start, int batch_end,
              np.ndarray[DTYPE_t, ndim=2] x_output,
              np.ndarray[DTYPE_t, ndim=2] x_hidden,
              np.ndarray[DTYPE_t, ndim=2] y,
              np.ndarray[DTYPE_t, ndim=2] weights_output,
              np.ndarray[DTYPE_t, ndim=1] bias_output,
              np.ndarray[DTYPE_t, ndim=2] weights_hidden,
              np.ndarray[DTYPE_t, ndim=1] bias_hidden,
              np.ndarray[DTYPE_t, ndim=2] delta_o,
              np.ndarray[DTYPE_t, ndim=2] delta_h,
              np.ndarray[DTYPE_t, ndim=2] dx_output,
              np.ndarray[DTYPE_t, ndim=2] dx_hidden,
              np.ndarray[DTYPE_t, ndim=2] weights_moment_o,
              np.ndarray[DTYPE_t, ndim=2] weights_moment_h,
              LossFunction loss,
              OutputFunction output,
              OutputFunction hidden,
              np.float64_t lr,
              np.float64_t lr_moment):

    cdef int batch_size = batch_end - batch_start

    # Output layer
    if isinstance(loss, (CrossEntropyLoss, MultiCrossEntropyLoss)):
        loss.dloss(y, batch_start, batch_end, x_output, delta_o)
    else:
        loss.dloss(y, batch_start, batch_end, x_output, delta_o)
        output.doutput(x_output, dx_output)
        delta_o *= dx_output

    # Hidden layer
    #np.dot(delta_o, weights_output.T, delta_h)
    if batch_size == 1:
       dgemv(
           CblasRowMajor,
           CblasNoTrans,
           weights_output.shape[0], weights_output.shape[1],  # M, N
           1.0,  # alpha
           <DTYPE_t*>weights_output.data, weights_output.shape[1],
           <DTYPE_t*>delta_o.data, 1,  # x, incX
           0.0,  # beta
           <DTYPE_t*>delta_h.data, 1)  # y, incY
    else:
       dgemm(
           CblasRowMajor,
           CblasNoTrans, CblasTrans,
           batch_size, delta_h.shape[1], delta_o.shape[1],  # M, N, K
           1.0,
           <DTYPE_t*>delta_o.data, delta_o.shape[1],
           <DTYPE_t*>weights_output.data, weights_output.shape[1],
           0.0, # beta
           <DTYPE_t*>delta_h.data, delta_h.shape[1])
    hidden.doutput(x_hidden, dx_hidden)
    delta_h *= dx_hidden

    # Update weights
    #weights_output += lr / batch_size * np.dot(x_hidden.T, delta_o)
    dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
          weights_output.shape[0], weights_output.shape[1], delta_o.shape[0],  # M, N, K
          lr / batch_size,  # alpha
          <DTYPE_t*>x_hidden.data, x_hidden.shape[1],
          <DTYPE_t*>delta_o.data, delta_o.shape[1],
          1.0,  # beta: add to previous value
          <DTYPE_t*>weights_output.data, weights_output.shape[1])

    weights_output += lr_moment * weights_moment_o

    bias_output += lr * np.mean(delta_o, axis=0)

    #weights_hidden += lr / batch_size * np.dot(X[batch_start:batch_end,:].T, delta_h)
    dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
          weights_hidden.shape[0], weights_hidden.shape[1], delta_h.shape[0],
          lr / batch_size,  # alpha
          <DTYPE_t*>X.data + batch_start * X.shape[1], X.shape[1],
          <DTYPE_t*>delta_h.data, delta_h.shape[1],
          1.0,  # beta: add to previous value
          <DTYPE_t*>weights_hidden.data, weights_hidden.shape[1])
    weights_hidden += lr_moment * weights_moment_h

    bias_hidden += lr * np.mean(delta_h, axis=0)


def sgd(np.ndarray[DTYPE_t, ndim=2] X not None,
        np.ndarray[DTYPE_t, ndim=2] y not None,
        LossFunction loss not None,
        OutputFunction output not None,
        OutputFunction hidden not None,
        np.ndarray[DTYPE_t, ndim=2] weights_hidden not None,
        np.ndarray[DTYPE_t, ndim=2] weights_output not None,
        np.ndarray[DTYPE_t, ndim=1] bias_hidden not None,
        np.ndarray[DTYPE_t, ndim=1] bias_output not None,
        np.float64_t lr,
        np.float64_t lr_moment,
        int n_hidden,
        int max_epochs,
        int batch_size,
        int shuffle_data):
    """Stochastic gradient descent for multilayer perceptron

    Parameters
    ----------

    Returns
    -------

    """
    cdef int n_samples = X.shape[0]
    cdef int n_features = X.shape[1]
    cdef int n_outs = y.shape[1]
    cdef int n_batches = n_samples / batch_size
    cdef int i
    cdef int j

    cdef np.ndarray[DTYPE_t, ndim=2] x_hidden = np.empty((batch_size, n_hidden))
    cdef np.ndarray[DTYPE_t, ndim=2] delta_h = np.empty((batch_size, n_hidden))
    cdef np.ndarray[DTYPE_t, ndim=2] dx_hidden = np.empty((batch_size, n_hidden))
    cdef np.ndarray[DTYPE_t, ndim=2] x_output = np.empty((batch_size, n_outs))
    cdef np.ndarray[DTYPE_t, ndim=2] delta_o = np.empty((batch_size, n_outs))
    cdef np.ndarray[DTYPE_t, ndim=2] dx_output = np.empty((batch_size, n_outs))

    cdef np.ndarray[DTYPE_t, ndim=2] loss_output = np.empty((batch_size, n_outs))

    cdef np.ndarray[DTYPE_t, ndim=2] weights_prev_o
    cdef np.ndarray[DTYPE_t, ndim=2] weights_prev_prev_o
    cdef np.ndarray[DTYPE_t, ndim=2] weights_moment_o

    cdef np.ndarray[DTYPE_t, ndim=2] weights_prev_h
    cdef np.ndarray[DTYPE_t, ndim=2] weights_prev_prev_h
    cdef np.ndarray[DTYPE_t, ndim=2] weights_moment_h

    if y.shape[0] != n_samples:
        raise ValueError("Shapes of X and y don't fit.")

    if n_samples % batch_size != 0:
        warnings.warn("Discarding some samples: \
            sample size not divisible by chunk size.")

    # Generate weights
    # TODO: Use Nguyen-Widrow initialization
    weights_hidden[:] = np.random.uniform(size=(n_features, n_hidden)) \
        / np.sqrt(n_features)
    bias_hidden[:] = np.zeros(n_hidden)
    weights_output[:] = np.random.uniform(size=(n_hidden, n_outs)) \
        / np.sqrt(n_hidden)
    bias_output[:] = np.zeros(n_outs)

    weights_prev_o = np.array(weights_output)
    weights_prev_prev_o = np.array(weights_output)
    weights_moment_o = np.zeros((n_hidden, n_outs))

    weights_prev_h = np.array(weights_hidden)
    weights_prev_prev_h = np.array(weights_hidden)
    weights_moment_h = np.zeros((n_features, n_hidden))


    for i in range(max_epochs):
        if shuffle_data:
            X, y = shuffle(X, y)
        for j in range(batch_size, n_samples + 1, batch_size):

            forward(X,
                    j - batch_size, j,
                    weights_hidden,
                    bias_hidden,
                    weights_output,
                    bias_output,
                    x_hidden,
                    x_output,
                    output,
                    hidden)

            #loss.loss(y, j - batch_size, j, x_output, loss_output)
            #print loss_output.mean(axis=0).sum()


            weights_moment_o[:] = weights_prev_o
            weights_moment_o -= weights_prev_prev_o
            weights_moment_h[:] = weights_prev_h
            weights_moment_h -= weights_prev_prev_h

            backward(X,
                     j - batch_size, j,
                     x_output,
                     x_hidden,
                     y,
                     weights_output,
                     bias_output,
                     weights_hidden,
                     bias_hidden,
                     delta_o,
                     delta_h,
                     dx_output,
                     dx_hidden,
                     weights_moment_o,
                     weights_moment_h,
                     loss,
                     output,
                     hidden,
                     lr,
                     lr_moment)

            weights_prev_prev_o[:] = weights_prev_o
            weights_prev_o[:] = weights_output
            weights_prev_prev_h[:] = weights_prev_h
            weights_prev_h[:] = weights_hidden

    return weights_hidden, bias_hidden, weights_output, bias_output


def predict(np.ndarray[DTYPE_t, ndim=2] X not None,
            np.ndarray[DTYPE_t, ndim=2] weights_hidden not None,
            np.ndarray[DTYPE_t, ndim=1] bias_hidden not None,
            np.ndarray[DTYPE_t, ndim=2] weights_output not None,
            np.ndarray[DTYPE_t, ndim=1] bias_output not None,
            OutputFunction output not None,
            OutputFunction hidden not None):

    cdef int n_samples = X.shape[0]
    cdef int n_hidden = weights_hidden.shape[1]
    cdef int n_outs = weights_output.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] x_hidden = np.empty((n_samples, n_hidden))
    cdef np.ndarray[DTYPE_t, ndim=2] x_output = np.empty((n_samples, n_outs))

    forward(X,
            0, X.shape[0],
            weights_hidden,
            bias_hidden,
            weights_output,
            bias_output,
            x_hidden,
            x_output,
            output,
            hidden)

    return x_output
