"""Utilities for the neural network modules
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause

import numpy as np

from ..utils.fixes import expit as logistic_sigmoid
from ..utils import check_random_state


def identity(X):
    """Return the input array."""
    return X


def logistic(X):
    """Compute the logistic function."""
    return logistic_sigmoid(X, out=X)


def tanh(X):
    """Compute the hyperbolic tan function."""
    return np.tanh(X, out=X)


def relu(X):
    """Compute the rectified linear unit function."""
    np.clip(X, 0, np.finfo(X.dtype).max, out=X)
    return X


def softmax(X):
    """Compute the K-way softmax function. """
    tmp = X - X.max(axis=1)[:, np.newaxis]
    np.exp(tmp, out=X)
    X /= X.sum(axis=1)[:, np.newaxis]

    return X


ACTIVATIONS = {'identity': identity, 'tanh': tanh, 'logistic': logistic,
               'relu': relu, 'softmax': softmax}


def _d_logistic(Z):
    """Compute the derivative of the logistic function."""
    return Z * (1 - Z)


def _d_tanh(Z):
    """Compute the derivative of the hyperbolic tan function."""
    return 1 - (Z ** 2)


def _d_relu(Z):
    """Compute the derivative of the rectified linear unit function."""
    return (Z > 0).astype('b')


DERIVATIVE_FUNCTIONS = {'tanh': _d_tanh, 'logistic': _d_logistic,
                        'relu': _d_relu}


def _squared_loss(y_true, y_pred):
    """Compute the square loss for regression."""
    return np.sum((y_true - y_pred) ** 2) / (2 * len(y_true))


def _log_loss(y_true, y_prob):
    """Compute Logistic loss for binary class.

    Max/Min clipping is enabled to prevent
    invalid zero value in log computation.
    """
    y_prob = np.clip(y_prob, 1e-10, 1 - 1e-10)

    return -np.sum(y_true * np.log(y_prob) +
                  (1 - y_true) * np.log(1 - y_prob)) / y_prob.shape[0]


LOSS_FUNCTIONS = {'squared_loss': _squared_loss, 'log_loss': _log_loss}


def init_weights(weight_scale, n_features, n_outputs, random_state):
    """Initialize the parameter weights."""
    rng = check_random_state(random_state)

    coef = rng.uniform(-1, 1, (n_features, n_outputs))
    intercept = rng.uniform(-1, 1, n_outputs)

    if weight_scale != 1:
        coef *= weight_scale
        intercept *= weight_scale

    return coef, intercept


def clear_layer_lists(*args):
    """Clear layer list from memory."""
    for layer_list in args:
        del layer_list[:]
