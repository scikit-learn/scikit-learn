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

class learning_rate_class():
    def __init__(self, learning_rate='rprop', 
                       learning_rate_init=0.1,
                       momentum=0.1,
                       nesterovs_momentum=False):
        self.learning_rate = learning_rate
        self.learning_rate_init = learning_rate_init
        self.momentum = momentum
        self.nesterovs_momentum = nesterovs_momentum

        self.learning_rate_list = None
        self.update_value_list = None

        self.velocity = None
        self.r = None
        self.historical_grad = 0
        self.count = 0

        

    def get_update(self, grad_list, theta, nerf=0.5, buff=1.2):
        n_parameters = grad_list.shape
        self.count += 1

        if self.learning_rate == 'constant':
            # momentum
            if self.velocity is None:
                self.velocity = np.zeros(n_parameters)

            self.velocity = (self.momentum * self.velocity
                             - self.learning_rate_init * grad_list)


            if self.nesterovs_momentum:
                return - (self.momentum * self.velocity 
                          - self.learning_rate_init * grad_list)

            return - self.velocity

        elif self.learning_rate == 'constant_avg':
            self.historical_grad += grad_list

            if self.count > 1:
                self.historical_grad *= (self.count - 1)
                self.historical_grad /= self.count

            return self.learning_rate_init * self.historical_grad

        elif self.learning_rate == 'rprop':
            grad_sign_list = np.sign(grad_list)

            if self.learning_rate_list is None:
                # Initialize learning rates
                self.learning_rate_list = np.ones(n_parameters) * self.learning_rate_init

            else:
                # Update learning rates
                result_grad_sign_list = grad_sign_list * self.old_grad_sign_list

                self.learning_rate_list[result_grad_sign_list > 0] *= buff
                self.learning_rate_list[result_grad_sign_list < 0] *= nerf

            self.old_grad_sign_list = grad_sign_list.copy()

            return self.learning_rate_list * grad_sign_list 

        elif self.learning_rate == 'bordes':
            if self.r is None:
                self.r = 1
                self.t0 = self.learning_rate_init
                self.t = 0.

                self.hess_approx = np.ones(n_parameters)

            else:
                self.r += 1
                self.t += 1

                self.hess_approx += (2./ self.r) * ((theta - self.old_theta) / 
                                                   (grad_list -  self.old_grad_list) -  self.hess_approx)

            self.old_theta = theta.copy()
            self.old_grad_list = grad_list.copy()

            return self.hess_approx * (1. / (self.t0 + self.t)) * grad_list

        elif self.learning_rate == 'adagrad':
            fudge_factor = 1e-6 
            
            self.historical_grad += grad_list ** 2

            adjusted_grad = grad_list / (fudge_factor + np.sqrt(self.historical_grad))
         
            return self.learning_rate_init * adjusted_grad

        elif self.learning_rate == 'rmsprop':
            if self.r is None:
                self.decay = 0.8
                self.r = 0

            self.r = (1 - self.decay) * grad_list ** 2 + self.decay * self.r

            update = (self.learning_rate_init / np.sqrt(self.r)) * grad_list
         
            return update