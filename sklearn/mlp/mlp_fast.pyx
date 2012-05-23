# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
import numpy as np
cimport numpy as np

import warnings
import pdb

from ..utils import shuffle

cdef class LossFunction:
    """Base class for loss functions"""

    cpdef np.ndarray loss(self, np.ndarray y, np.ndarray p):
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

    cpdef np.ndarray dloss(self, np.ndarray y, np.ndarray p):
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

    cpdef np.ndarray loss(self, np.ndarray y, np.ndarray p):
        return 0.5 * (y - p) * (y - p)

    cpdef np.ndarray dloss(self, np.ndarray y, np.ndarray p):
        return y - p

    def __reduce__(self):
        return SquaredLoss, ()


cdef class OutputFunction:
    """Base class for ouput functions"""

    cpdef np.ndarray output(self, np.ndarray x, np.ndarray out=None):
        raise NotImplementedError()

    cpdef np.ndarray doutput(self, np.ndarray x, np.ndarray out=None):
        raise NotImplementedError()


cdef class Tanh(OutputFunction):

    cpdef np.ndarray output(self, np.ndarray x, np.ndarray out=None):
        if out is not None:
            return np.tanh(x, out)
        else:
            return np.tanh(x)

    cpdef np.ndarray doutput(self, np.ndarray x, np.ndarray out=None):
        if out is not None:
            out[:] = -x * x + 1
        else:
            return -x * x + 1


cpdef forward(np.ndarray X,
              np.ndarray weights_hidden,
              np.ndarray bias_hidden,
              np.ndarray weights_output,
              np.ndarray bias_output,
              np.ndarray x_hidden,
              np.ndarray x_output,
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


cpdef backward(np.ndarray y,
               np.ndarray weights_hidden,
               np.ndarray weights_output,
               np.ndarray x_hidden,
               np.ndarray x_output,
               OutputFunction output,
               OutputFunction hidden,
               LossFunction loss,
               np.ndarray delta_h,
               np.ndarray delta_o):

    # Output layer
    delta_o[:] = loss.dloss(y, x_output) * output.doutput(x_output)

    # Hidden layer
    delta_h[:] = np.dot(delta_o, weights_output.T) * hidden.doutput(x_hidden)


cpdef sgd(np.ndarray X,
          np.ndarray y,
          LossFunction loss,
          OutputFunction output,
          OutputFunction hidden,
          np.ndarray weights_hidden,
          np.ndarray weights_output,
          np.ndarray bias_hidden,
          np.ndarray bias_output,
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

    cdef np.ndarray x_hidden = np.empty((batch_size, n_hidden))
    cdef np.ndarray delta_h = np.empty((batch_size, n_hidden))
    cdef np.ndarray x_output = np.empty((batch_size, n_outs))
    cdef np.ndarray delta_o = np.empty((batch_size, n_outs))

    if y.shape[0] != n_samples:
        raise ValueError("Shapes of X and y don't fit.")

    if n_samples % batch_size != 0:
        warnings.warn("Discarding some samples: \
            sample size not divisible by chunk size.")

    # Generate weights
    # TODO: Use Nguyen-Widrow initialization
    weights_hidden[:] = np.random.standard_normal(size=(n_features, n_hidden)) \
        / np.sqrt(n_features)
    bias_hidden[:] = np.zeros(n_hidden)
    weights_output[:] = np.random.standard_normal(size=(n_hidden, n_outs)) \
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

            backward(y[j - batch_size:j],
                     weights_hidden,
                     weights_output,
                     x_hidden,
                     x_output,
                     output,
                     hidden,
                     loss,
                     delta_h,
                     delta_o)

            # Update weights
            weights_output += lr / batch_size * np.dot(x_hidden.T, delta_o)
            bias_output += lr * np.mean(delta_o, axis=0)
            weights_hidden += lr / batch_size * np.dot(X[j - batch_size:j].T, delta_h)
            bias_hidden += lr * np.mean(delta_h, axis=0)

    return weights_hidden, bias_hidden, weights_output, bias_output

cpdef predict(np.ndarray X,
              np.ndarray weights_hidden,
              np.ndarray bias_hidden,
              np.ndarray weights_output,
              np.ndarray bias_output,
              OutputFunction output,
              OutputFunction hidden):

    cdef int n_samples = X.shape[0]
    cdef int n_hidden = weights_hidden.shape[1]
    cdef int n_outs = weights_output.shape[1]
    cdef np.ndarray x_hidden = np.empty((n_samples, n_hidden))
    cdef np.ndarray x_output = np.empty((n_samples, n_outs))

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
