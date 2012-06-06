# encoding: utf-8
# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
import numpy as np
cimport numpy as np

import warnings

from ..utils import shuffle

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef class LossFunction:
    """Base class for loss functions"""

    cdef inline loss(self,
                     np.ndarray[DTYPE_t] y,
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

    cdef inline dloss(self,
                      np.ndarray[DTYPE_t] y,
                      np.ndarray[DTYPE_t] p,
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

    cdef inline loss(self,
                     np.ndarray[DTYPE_t, ndim=2] y,
                     np.ndarray[DTYPE_t, ndim=2] p,
                     np.ndarray[DTYPE_t, ndim=2] out):
        out[:] = 0.5 * (y - p) * (y - p)

    cdef inline dloss(self,
                      np.ndarray[DTYPE_t, ndim=2] y,
                      np.ndarray[DTYPE_t, ndim=2] p,
                      np.ndarray[DTYPE_t, ndim=2] out):
        out[:] = y - p

    def __reduce__(self):
        return SquaredLoss, ()


cdef class CrossEntropyLoss(LossFunction):
    """Cross entropy loss function"""
    cdef inline loss(self,
                     np.ndarray[DTYPE_t, ndim=2] y,
                     np.ndarray[DTYPE_t, ndim=2] p,
                     np.ndarray[DTYPE_t, ndim=2] out):
        out[:] -y * np.log(p) - (1 - y) * np.log(1 - p)

    cdef inline dloss(self,
                      np.ndarray[DTYPE_t, ndim=2] y,
                      np.ndarray[DTYPE_t, ndim=2] p,
                      np.ndarray[DTYPE_t, ndim=2] out):
        out[:] = y - p


cdef class MultiCrossEntropyLoss(LossFunction):
    """Multinomial cross entropy loss function."""

    cdef inline loss(self,
                     np.ndarray[DTYPE_t, ndim=2] y,
                     np.ndarray[DTYPE_t, ndim=2] p,
                     np.ndarray[DTYPE_t, ndim=2] out):
        out[:] = -y * np.log(p)

    cdef inline dloss(self,
                      np.ndarray[DTYPE_t, ndim=2] y,
                      np.ndarray[DTYPE_t, ndim=2] p,
                      np.ndarray[DTYPE_t, ndim=2] out):
        out[:] = y - p


cdef class OutputFunction:
    """Base class for ouput functions"""

    cdef inline output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        raise NotImplementedError()

    cdef inline np.ndarray[DTYPE_t, ndim=2] doutput(self, np.ndarray[DTYPE_t, ndim=2] x):
        raise NotImplementedError()


cdef class Tanh(OutputFunction):

    cdef inline output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        np.tanh(x, out)

    cdef inline np.ndarray[DTYPE_t, ndim=2] doutput(self, np.ndarray[DTYPE_t, ndim=2] x):
        return -x * x + 1


cdef class SoftMax(OutputFunction):

    cdef inline output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        np.exp(x, out)
        out /= np.sum(out, axis=1)[:, np.newaxis]


cdef class LogSig(OutputFunction):

    cdef inline output(self, np.ndarray[DTYPE_t, ndim=2] x, np.ndarray[DTYPE_t, ndim=2] out):
        np.exp(-x, out)
        out[:] = 1 / (1 + out)


cdef forward(np.ndarray[DTYPE_t, ndim=2] X,
             np.ndarray[DTYPE_t, ndim=2] weights_hidden,
             np.ndarray[DTYPE_t, ndim=1] bias_hidden,
             np.ndarray[DTYPE_t, ndim=2] weights_output,
             np.ndarray[DTYPE_t, ndim=1] bias_output,
             np.ndarray[DTYPE_t, ndim=2] x_hidden,
             np.ndarray[DTYPE_t, ndim=2] x_output,
             OutputFunction output,
             OutputFunction hidden):

    # Hidden layer
    np.dot(X, weights_hidden, x_hidden)
    x_hidden += bias_hidden
    hidden.output(x_hidden, x_hidden)

    # Output layer
    np.dot(x_hidden, weights_output, x_output)
    x_output += bias_output
    output.output(x_output, x_output)


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
    cdef np.ndarray[DTYPE_t, ndim=2] x_output = np.empty((batch_size, n_outs))
    cdef np.ndarray[DTYPE_t, ndim=2] delta_o = np.empty((batch_size, n_outs))

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

    if shuffle_data:
        X, y = shuffle(X, y)

    for i in range(max_epochs):
        for j in range(batch_size, n_samples + 1, batch_size):

            forward(X[j - batch_size:j],
                    weights_hidden,
                    bias_hidden,
                    weights_output,
                    bias_output,
                    x_hidden,
                    x_output,
                    output,
                    hidden)

            # BACKWARD

            # Output layer
            if isinstance(loss, (CrossEntropyLoss, MultiCrossEntropyLoss)):
                loss.dloss(y[j - batch_size:j], x_output, delta_o)
            else:
                loss.dloss(y[j - batch_size:j], x_output, delta_o)
                delta_o *= output.doutput(x_output)

            # Hidden layer
            np.dot(delta_o, weights_output.T, delta_h)
            delta_h *= hidden.doutput(x_hidden)

            # Update weights
            weights_output += lr / batch_size * np.dot(x_hidden.T, delta_o)
            bias_output += lr * np.mean(delta_o, axis=0)
            weights_hidden += lr / batch_size * np.dot(X[j - batch_size:j].T, delta_h)
            bias_hidden += lr * np.mean(delta_h, axis=0)

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
            weights_hidden,
            bias_hidden,
            weights_output,
            bias_output,
            x_hidden,
            x_output,
            output,
            hidden)

    return x_output
