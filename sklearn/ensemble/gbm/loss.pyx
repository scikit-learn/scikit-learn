# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
"""
This module contains the loss classes.

Specific losses are used for regression, binary classification or multiclass
classification.
"""
from abc import ABC, abstractmethod

cimport cython

import numpy as np
cimport numpy as np

from scipy.special import expit, logsumexp

from .types import Y_DTYPE

ctypedef np.npy_float32 NPY_Y_DTYPE


cdef get_threads_chunks(unsigned int total_size):
    """Get start and end indices of threads in an array of size total_size.

    The interval [0, total_size - 1] is divided into n_threads contiguous
    regions, and the starts and ends of each region are returned. Used to
    simulate a 'static' scheduling.
    """
    cdef:
        np.ndarray[np.uint32_t] sizes
        np.ndarray[np.uint32_t] starts
        np.ndarray[np.uint32_t] ends
        unsigned int n_threads

    n_threads = 1  # TODO: change this
    sizes = np.full(n_threads, total_size // n_threads, dtype=np.uint32)
    sizes[:total_size % n_threads] += 1
    starts = np.zeros(n_threads, dtype=np.uint32)
    starts[1:] = np.cumsum(sizes[:-1])
    ends = starts + sizes

    return starts, ends, n_threads


class BaseLoss(ABC):
    """Base class for a loss."""

    def init_gradients_and_hessians(self, n_samples, prediction_dim):
        """Return initial gradients and hessians.

        Unless hessians are constant, arrays are initialized with undefined
        values.

        Parameters
        ----------
        n_samples : int
            The number of samples passed to `fit()`
        prediction_dim : int
            The dimension of a raw prediction, i.e. the number of trees
            built at each iteration. Equals 1 for regression and binary
            classification, or K where K is the number of classes for
            multiclass classification.

        Returns
        -------
        gradients : array-like, shape=(n_samples * prediction_dim)
        hessians : array-like, shape=(n_samples * prediction_dim).
            If hessians are constant (e.g. for ``LeastSquares`` loss, shape
            is (1,) and the array is initialized to ``1``.
        """
        shape = n_samples * prediction_dim
        gradients = np.empty(shape=shape, dtype=Y_DTYPE)
        if self.hessian_is_constant:
            hessians = np.ones(shape=1, dtype=Y_DTYPE)
        else:
            hessians = np.empty(shape=shape, dtype=Y_DTYPE)

        return gradients, hessians

    @abstractmethod
    def get_baseline_prediction(self, y_train, prediction_dim):
        """Return initial predictions (before the first iteration).

        Parameters
        ----------
        y_train : array-like, shape=(n_samples,)
            The target training values.
        prediction_dim : int
            The dimension of one prediction: 1 for binary classification and
            regression, n_classes for multiclass classification.

        Returns
        -------
        baseline_prediction: float or array of shape (1, prediction_dim)
            The baseline prediction.
        """
        pass

    @abstractmethod
    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        """Update gradients and hessians arrays, inplace.

        The gradients (resp. hessians) are the first (resp. second) order
        derivatives of the loss for each sample with respect to the
        predictions of model, evaluated at iteration ``i - 1``.

        Parameters
        ----------
        gradients : array-like, shape=(n_samples * prediction_dim)
            The gradients (treated as OUT array).
        hessians : array-like, shape=(n_samples * prediction_dim) or \
            (1,)
            The hessians (treated as OUT array).
        y_true : array-like, shape=(n_samples,)
            The true target values or each training sample.
        raw_predictions : array-like, shape=(n_samples, prediction_dim)
            The raw_predictions (i.e. values from the trees) of the tree
            ensemble at iteration ``i - 1``.
        """
        pass


class LeastSquares(BaseLoss):
    """Least squares loss, for regression.

    For a given sample x_i, least squares loss is defined as::

        loss(x_i) = (y_true_i - raw_pred_i)**2
    """

    hessian_is_constant = True

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        loss = np.power(y_true - raw_predictions, 2)
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        return np.mean(y_train)

    def inverse_link_function(self, raw_predictions):
        return raw_predictions

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        raw_predictions = raw_predictions.reshape(-1)
        return _update_gradients_least_squares(gradients, y_true,
                                               raw_predictions)


def _update_gradients_least_squares(NPY_Y_DTYPE[:] gradients, NPY_Y_DTYPE[:] y_true, NPY_Y_DTYPE[:] raw_predictions):
    cdef:
        unsigned int n_samples
        unsigned int i

    n_samples = raw_predictions.shape[0]
    for i in range(n_samples):
        # Note: a more correct exp is 2 * (raw_predictions - y_true) but
        # since we use 1 for the constant hessian value (and not 2) this
        # is strictly equivalent for the leaves values.
        gradients[i] = raw_predictions[i] - y_true[i]


## class BinaryCrossEntropy(BaseLoss):
##     """Binary cross-entropy loss, for binary classification.
## 
##     For a given sample x_i, the binary cross-entropy loss is defined as the
##     negative log-likelihood of the model which can be expressed as::
## 
##         loss(x_i) = log(1 + exp(raw_pred_i)) - y_true_i * raw_pred_i
## 
##     See The Elements of Statistical Learning, by Hastie, Tibshirani, Friedman.
##     """
## 
##     hessian_is_constant = False
##     inverse_link_function = staticmethod(expit)
## 
##     def __call__(self, y_true, raw_predictions, average=True):
##         # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
##         # return a view.
##         raw_predictions = raw_predictions.reshape(-1)
##         # logaddexp(0, x) = log(1 + exp(x))
##         loss = np.logaddexp(0, raw_predictions) - y_true * raw_predictions
##         return loss.mean() if average else loss
## 
##     def get_baseline_prediction(self, y_train, prediction_dim):
##         proba_positive_class = np.mean(y_train)
##         eps = np.finfo(y_train.dtype).eps
##         proba_positive_class = np.clip(proba_positive_class, eps, 1 - eps)
##         # log(x / 1 - x) is the anti function of sigmoid, or the link function
##         # of the Binomial model.
##         return np.log(proba_positive_class / (1 - proba_positive_class))
## 
##     def update_gradients_and_hessians(self, gradients, hessians, y_true,
##                                       raw_predictions):
##         raw_predictions = raw_predictions.reshape(-1)
##         return _update_gradients_hessians_binary_crossentropy(
##             gradients, hessians, y_true, raw_predictions)
## 
##     def predict_proba(self, raw_predictions):
##         # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
##         # return a view.
##         raw_predictions = raw_predictions.reshape(-1)
##         proba = np.empty((raw_predictions.shape[0], 2), dtype=np.float32)
##         proba[:, 1] = expit(raw_predictions)
##         proba[:, 0] = 1 - proba[:, 1]
##         return proba
## 
## 
## def _update_gradients_hessians_binary_crossentropy(float [:] gradients,
## float [:] hessians, float_or_double [:] y_true, double [:] raw_predictions):
##     cdef:
##         unsigned int n_samples
##         unsigned int i
##         unsigned int thread_idx
##         unsigned int n_threads
##         unsigned int [:] starts
##         unsigned int [:] ends
##     n_samples = raw_predictions.shape[0]
##     starts, ends, n_threads = get_threads_chunks(total_size=n_samples)
##     for thread_idx in range(n_threads):
##         for i in range(starts[thread_idx], ends[thread_idx]):
##             gradients[i] = <float>expit(raw_predictions[i]) - y_true[i]
##             gradient_abs = np.abs(gradients[i])
##             hessians[i] = gradient_abs * (1. - gradient_abs)
## 
## 
## class CategoricalCrossEntropy(BaseLoss):
##     """Categorical cross-entropy loss, for multiclass classification.
## 
##     For a given sample x_i, the categorical cross-entropy loss is defined as
##     the negative log-likelihood of the model and generalizes the binary
##     cross-entropy to more than 2 classes.
##     """
## 
##     hessian_is_constant = False
## 
##     def __call__(self, y_true, raw_predictions, average=True):
##         one_hot_true = np.zeros_like(raw_predictions)
##         prediction_dim = raw_predictions.shape[1]
##         for k in range(prediction_dim):
##             one_hot_true[:, k] = (y_true == k)
## 
##         loss = (logsumexp(raw_predictions, axis=1) -
##                 (one_hot_true * raw_predictions).sum(axis=1))
##         return loss.mean() if average else loss
## 
##     def get_baseline_prediction(self, y_train, prediction_dim):
##         init_value = np.zeros(
##             shape=(1, prediction_dim),
##             dtype=np.float32
##         )
##         eps = np.finfo(y_train.dtype).eps
##         for k in range(prediction_dim):
##             proba_kth_class = np.mean(y_train == k)
##             proba_kth_class = np.clip(proba_kth_class, eps, 1 - eps)
##             init_value[:, k] += np.log(proba_kth_class)
## 
##         return init_value
## 
##     def update_gradients_and_hessians(self, gradients, hessians, y_true,
##                                       raw_predictions):
##         return _update_gradients_hessians_categorical_crossentropy(
##             gradients, hessians, y_true, raw_predictions)
## 
##     def predict_proba(self, raw_predictions):
##         # TODO: This could be done in parallel
##         # compute softmax (using exp(log(softmax)))
##         return np.exp(raw_predictions -
##                       logsumexp(raw_predictions, axis=1)[:, np.newaxis])
## 
## 
## def _update_gradients_hessians_categorical_crossentropy(
##         float [:] gradients, float [:] hessians, float_or_double [:] y_true,
##         float_or_double [:, :] raw_predictions):
##     # Here gradients and hessians are of shape
##     # (n_samples * prediction_dim,).
##     # y_true is of shape (n_samples,).
##     # raw_predictions is of shape (n_samples, raw_predictions)
##     cdef:
##         unsigned int n_samples
##         unsigned int prediction_dim
##         unsigned int i
##         unsigned int k
##         unsigned int thread_idx
##         unsigned int n_threads
##         unsigned int [:] starts
##         unsigned int [:] ends
##         float p_k
## 
##     n_samples = raw_predictions.shape[0]
##     prediction_dim = raw_predictions.shape[1]
##     starts, ends, n_threads = get_threads_chunks(total_size=n_samples)
##     for k in range(prediction_dim):
##         gradients_at_k = gradients[n_samples * k:n_samples * (k + 1)]
##         hessians_at_k = hessians[n_samples * k:n_samples * (k + 1)]
##         for thread_idx in range(n_threads):
##             for i in range(starts[thread_idx], ends[thread_idx]):
##                 # p_k is the probability that class(ith sample) == k.
##                 # This is a regular softmax.
##                 p_k = np.exp(raw_predictions[i, k] -
##                              logsumexp(raw_predictions[i, :]))
##                 gradients_at_k[i] = p_k - (y_true[i] == k)
##                 hessians_at_k[i] = p_k * (1. - p_k)


_LOSSES = {'least_squares': LeastSquares}
