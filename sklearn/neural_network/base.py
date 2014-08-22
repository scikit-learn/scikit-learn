"""Utilities for the neural network modules
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause

import numpy as np

from ..utils.fixes import expit as logistic_sigmoid


def identity(X):
    """Simply return the input array.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Data, where n_samples is the number of samples
        and n_features is the number of features.

    Returns
    -------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Same as the input data.
    """
    return X


def logistic(X):
    """Compute the logistic function inplace.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    return logistic_sigmoid(X, out=X)


def tanh(X):
    """Compute the hyperbolic tan function inplace.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    return np.tanh(X, out=X)


def relu(X):
    """Compute the rectified linear unit function inplace.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    np.clip(X, 0, np.finfo(X.dtype).max, out=X)
    return X


def softmax(X):
    """Compute the  K-way softmax function inplace.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    X_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    tmp = X - X.max(axis=1)[:, np.newaxis]
    np.exp(tmp, out=X)
    X /= X.sum(axis=1)[:, np.newaxis]

    return X


ACTIVATIONS = {'identity': identity, 'tanh': tanh, 'logistic': logistic,
               'relu': relu, 'softmax': softmax}


def logistic_derivative(Z):
    """Compute the derivative of the logistic function.

    Parameters
    ----------
    Z : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    Z_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    return Z * (1 - Z)


def tanh_derivative(Z):
    """Compute the derivative of the hyperbolic tan function.

    Parameters
    ----------
    Z : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    Z_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    return 1 - (Z ** 2)


def relu_derivative(Z):
    """Compute the derivative of the rectified linear unit function.

    Parameters
    ----------
    Z : {array-like, sparse matrix}, shape (n_samples, n_features)
        The input data.

    Returns
    -------
    Z_new : {array-like, sparse matrix}, shape (n_samples, n_features)
        The transformed data.
    """
    return (Z > 0).astype(Z.dtype)


DERIVATIVES = {'tanh': tanh_derivative, 'logistic': logistic_derivative,
               'relu': relu_derivative}


def squared_loss(y_true, y_pred):
    """Compute the squared loss for regression.

    Parameters
    ----------
    y_true : array-like or label indicator matrix
        Ground truth (correct) labels.

    y_pred : array-like or label indicator matrix
        Predicted labels, as returned by a regression estimator.

    Returns
    -------
    score : float
        The degree to which the samples are correctly predicted.
    """
    return ((y_true - y_pred) ** 2).mean() / 2


def log_loss(y_true, y_prob):
    """Compute Logistic loss for classification.

    Parameters
    ----------
    y_true : array-like or label indicator matrix
        Ground truth (correct) labels.

    y_pred : array-like of float, shape = (n_samples, n_classes)
        Predicted probabilities, as returned by a classifier's
        predict_proba method.

    Returns
    -------
    score : float
        The degree to which the samples are correctly predicted.
    """
    y_prob = np.clip(y_prob, 1e-10, 1 - 1e-10)

    return -np.sum(y_true * np.log(y_prob) +
                  (1 - y_true) * np.log(1 - y_prob)) / y_prob.shape[0]


LOSS_FUNCTIONS = {'squared_loss': squared_loss, 'log_loss': log_loss}