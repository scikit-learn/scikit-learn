# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""
This module contains the loss classes.

Specific losses are used for regression, binary classification or multiclass
classification.
"""
from abc import ABC, abstractmethod

cimport cython
from cython.parallel import prange
from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np
from scipy.special import expit
try:
    from scipy.special import logsumexp
except ImportError:
    from scipy.misc import logsumexp

from libc.math cimport fabs, exp, log

from .types import Y_DTYPE
from .types cimport Y_DTYPE_C
from .types import G_H_DTYPE
from .types cimport G_H_DTYPE_C


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
        gradients = np.empty(shape=shape, dtype=G_H_DTYPE)
        if self.hessians_are_constant:
            # if the hessians are constant, we consider they are equal to 1.
            # this is correct as long as we adjust the gradients. See e.g. LS
            # loss
            hessians = np.ones(shape=1, dtype=G_H_DTYPE)
        else:
            hessians = np.empty(shape=shape, dtype=G_H_DTYPE)

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
        baseline_prediction: float or array of shape (prediction_dim, 1)
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

    hessians_are_constant = True

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        loss = np.power(y_true - raw_predictions, 2)
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        return np.mean(y_train).astype(Y_DTYPE)

    @staticmethod
    def inverse_link_function(raw_predictions):
        return raw_predictions

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        raw_predictions = raw_predictions.reshape(-1)
        return _update_gradients_least_squares(gradients, y_true,
                                               raw_predictions)


cdef void _update_gradients_least_squares(
        G_H_DTYPE_C [::1] gradients,
        const Y_DTYPE_C [::1] y_true,
        const Y_DTYPE_C [::1] raw_predictions) nogil:
    cdef:
        int n_samples
        int i

    n_samples = raw_predictions.shape[0]
    for i in prange(n_samples, schedule='static'):
        # Note: a more correct exp is 2 * (raw_predictions - y_true) but
        # since we use 1 for the constant hessian value (and not 2) this
        # is strictly equivalent for the leaves values.
        gradients[i] = raw_predictions[i] - y_true[i]


class BinaryCrossEntropy(BaseLoss):
    """Binary cross-entropy loss, for binary classification.

    For a given sample x_i, the binary cross-entropy loss is defined as the
    negative log-likelihood of the model which can be expressed as::

        loss(x_i) = log(1 + exp(raw_pred_i)) - y_true_i * raw_pred_i

    See The Elements of Statistical Learning, by Hastie, Tibshirani, Friedman.
    """

    hessians_are_constant = False
    inverse_link_function = staticmethod(expit)

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        # logaddexp(0, x) = log(1 + exp(x))
        loss = np.logaddexp(0, raw_predictions) - y_true * raw_predictions
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        proba_positive_class = np.mean(y_train)
        eps = np.finfo(y_train.dtype).eps
        proba_positive_class = np.clip(proba_positive_class, eps, 1 - eps)
        # log(x / 1 - x) is the anti function of sigmoid, or the link function
        # of the Binomial model.
        return np.log(proba_positive_class / (1 - proba_positive_class))

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        return _update_gradients_hessians_binary_crossentropy(
            gradients, hessians, y_true, raw_predictions)

    def predict_proba(self, raw_predictions):
        # shape (n_samples, 1) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        proba = np.empty((raw_predictions.shape[0], 2), dtype=Y_DTYPE)
        proba[:, 1] = expit(raw_predictions)
        proba[:, 0] = 1 - proba[:, 1]
        return proba


cdef void _update_gradients_hessians_binary_crossentropy(
        G_H_DTYPE_C [::1] gradients,
        G_H_DTYPE_C [::1] hessians,
        const Y_DTYPE_C [::1] y_true,
        const Y_DTYPE_C [::1] raw_predictions) nogil:
    cdef:
        int n_samples
        G_H_DTYPE_C gradient_abs
        int i

    n_samples = raw_predictions.shape[0]
    for i in prange(n_samples, schedule='static'):
        gradients[i] = cexpit(raw_predictions[i]) - y_true[i]
        gradient_abs = fabs(gradients[i])
        hessians[i] = gradient_abs * (1. - gradient_abs)


class CategoricalCrossEntropy(BaseLoss):
    """Categorical cross-entropy loss, for multiclass classification.

    For a given sample x_i, the categorical cross-entropy loss is defined as
    the negative log-likelihood of the model and generalizes the binary
    cross-entropy to more than 2 classes.
    """

    hessians_are_constant = False

    def __call__(self, y_true, raw_predictions, average=True):
        one_hot_true = np.zeros_like(raw_predictions)
        prediction_dim = raw_predictions.shape[0]
        for k in range(prediction_dim):
            one_hot_true[k, :] = (y_true == k)

        loss = (logsumexp(raw_predictions, axis=0) -
                (one_hot_true * raw_predictions).sum(axis=0))
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        init_value = np.zeros(shape=(prediction_dim, 1), dtype=Y_DTYPE)
        eps = np.finfo(y_train.dtype).eps
        for k in range(prediction_dim):
            proba_kth_class = np.mean(y_train == k)
            proba_kth_class = np.clip(proba_kth_class, eps, 1 - eps)
            init_value[k, :] += np.log(proba_kth_class)

        return init_value

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        return _update_gradients_hessians_categorical_crossentropy(
            gradients, hessians, y_true, raw_predictions)

    def predict_proba(self, raw_predictions):
        # TODO: This could be done in parallel
        # compute softmax (using exp(log(softmax)))
        proba = np.exp(raw_predictions -
                       logsumexp(raw_predictions, axis=0)[np.newaxis, :])
        return proba.T


cdef void _update_gradients_hessians_categorical_crossentropy(
        G_H_DTYPE_C [::1] gradients,  # shape (n_samples * pred_dim,), OUT
        G_H_DTYPE_C [::1] hessians,  # shape (n_samples * pred_dim,), OUT
        const Y_DTYPE_C [::1] y_true,  # shape (n_samples,), IN
        # shape (n_samples, n_tree_per_iter), IN
        const Y_DTYPE_C [:, ::1] raw_predictions) nogil:
    cdef:
        unsigned int prediction_dim = raw_predictions.shape[0]
        int n_samples = raw_predictions.shape[1]
        unsigned int k
        int i
        Y_DTYPE_C * p = <Y_DTYPE_C *> malloc(sizeof(Y_DTYPE_C) *
                                             (prediction_dim * n_samples))
        Y_DTYPE_C p_k
        G_H_DTYPE_C [::1] gradients_at_k,
        G_H_DTYPE_C [::1] hessians_at_k,

    for i in prange(n_samples, schedule='static'):
        # first compute softmaxes of sample i for each class
        for k in range(prediction_dim):
            p[i * prediction_dim + k] = raw_predictions[k, i]
        compute_softmax(p + (i * prediction_dim), prediction_dim)
        # then update gradients and hessians
        for k in range(prediction_dim):
            # p_k is the probability that class(ith sample) == k.
            p_k = p[i * prediction_dim + k]
            gradients[n_samples * k + i] = p_k - (y_true[i] == k)
            hessians[n_samples * k + i] = p_k * (1. - p_k)
    free(p)


cdef inline void compute_softmax(
        Y_DTYPE_C * p,  # IN OUT, treated as array with <pred_dim> entries
        const unsigned int prediction_dim) nogil:
    """Compute softmaxes of values in p."""

    cdef:
        Y_DTYPE_C max_value = p[0]
        Y_DTYPE_C sum_exps = 0.
        unsigned int k

    # Compute max value of array for numerical stability
    for k in range(1, prediction_dim):
        if max_value < p[k]:
            max_value = p[k]

    for k in range(prediction_dim):
        p[k] = exp(p[k] - max_value)
        sum_exps += p[k]

    for k in range(prediction_dim):
        p[k] /= sum_exps


cdef inline Y_DTYPE_C cexpit(const Y_DTYPE_C x) nogil:
    """Custom expit (logistic sigmoid function)"""
    return 1. / (1. + exp(-x))


_LOSSES = {
    'least_squares': LeastSquares,
    'binary_crossentropy': BinaryCrossEntropy,
    'categorical_crossentropy': CategoricalCrossEntropy
}
