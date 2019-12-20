"""
This module contains the loss classes.

Specific losses are used for regression, binary classification or multiclass
classification.
"""
# Author: Nicolas Hug

from abc import ABC, abstractmethod

import numpy as np
from scipy.special import expit
try:  # logsumexp was moved from mist to special in 0.19
    from scipy.special import logsumexp
except ImportError:
    from scipy.misc import logsumexp

from .common import Y_DTYPE
from .common import G_H_DTYPE
from ._loss import _update_gradients_least_squares
from ._loss import _update_gradients_least_absolute_deviation
from ._loss import _update_gradients_hessians_binary_crossentropy
from ._loss import _update_gradients_hessians_categorical_crossentropy


class BaseLoss(ABC):
    """Base class for a loss."""

    # This variable indicates whether the loss requires the leaves values to
    # be updated once the tree has been trained. The trees are trained to
    # predict a Newton-Raphson step (see grower._finalize_leaf()). But for
    # some losses (e.g. least absolute deviation) we need to adjust the tree
    # values to account for the "line search" of the gradient descent
    # procedure. See the original paper Greedy Function Approximation: A
    # Gradient Boosting Machine by Friedman
    # (https://statweb.stanford.edu/~jhf/ftp/trebst.pdf) for the theory.
    need_update_leaves_values = False

    def init_gradients_and_hessians(self, n_samples, prediction_dim):
        """Return initial gradients and hessians.

        Unless hessians are constant, arrays are initialized with undefined
        values.

        Parameters
        ----------
        n_samples : int
            The number of samples passed to `fit()`.
        prediction_dim : int
            The dimension of a raw prediction, i.e. the number of trees
            built at each iteration. Equals 1 for regression and binary
            classification, or K where K is the number of classes for
            multiclass classification.

        Returns
        -------
        gradients : ndarray, shape (prediction_dim, n_samples)
            The initial gradients. The array is not initialized.
        hessians : ndarray, shape (prediction_dim, n_samples)
            If hessians are constant (e.g. for `LeastSquares` loss, the
            array is initialized to ``1``. Otherwise, the array is allocated
            without being initialized.
        """
        shape = (prediction_dim, n_samples)
        gradients = np.empty(shape=shape, dtype=G_H_DTYPE)
        if self.hessians_are_constant:
            # If the hessians are constant, we consider they are equal to 1.
            # - This is correct for the half LS loss
            # - For LAD loss, hessians are actually 0, but they are always
            #   ignored anyway.
            hessians = np.ones(shape=(1, 1), dtype=G_H_DTYPE)
        else:
            hessians = np.empty(shape=shape, dtype=G_H_DTYPE)

        return gradients, hessians

    @abstractmethod
    def get_baseline_prediction(self, y_train, prediction_dim):
        """Return initial predictions (before the first iteration).

        Parameters
        ----------
        y_train : ndarray, shape (n_samples,)
            The target training values.
        prediction_dim : int
            The dimension of one prediction: 1 for binary classification and
            regression, n_classes for multiclass classification.

        Returns
        -------
        baseline_prediction : float or ndarray, shape (1, prediction_dim)
            The baseline prediction.
        """

    @abstractmethod
    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        """Update gradients and hessians arrays, inplace.

        The gradients (resp. hessians) are the first (resp. second) order
        derivatives of the loss for each sample with respect to the
        predictions of model, evaluated at iteration ``i - 1``.

        Parameters
        ----------
        gradients : ndarray, shape (prediction_dim, n_samples)
            The gradients (treated as OUT array).
        hessians : ndarray, shape (prediction_dim, n_samples) or \
            (1,)
            The hessians (treated as OUT array).
        y_true : ndarray, shape (n_samples,)
            The true target values or each training sample.
        raw_predictions : ndarray, shape (prediction_dim, n_samples)
            The raw_predictions (i.e. values from the trees) of the tree
            ensemble at iteration ``i - 1``.
        """


class LeastSquares(BaseLoss):
    """Least squares loss, for regression.

    For a given sample x_i, least squares loss is defined as::

        loss(x_i) = 0.5 * (y_true_i - raw_pred_i)**2

    This actually computes the half least squares loss to optimize simplify
    the computation of the gradients and get a unit hessian (and be consistent
    with what is done in LightGBM).
    """

    hessians_are_constant = True

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        loss = 0.5 * np.power(y_true - raw_predictions, 2)
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        return np.mean(y_train)

    @staticmethod
    def inverse_link_function(raw_predictions):
        return raw_predictions

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        gradients = gradients.reshape(-1)
        _update_gradients_least_squares(gradients, y_true, raw_predictions)


class LeastAbsoluteDeviation(BaseLoss):
    """Least absolute deviation, for regression.

    For a given sample x_i, the loss is defined as::

        loss(x_i) = |y_true_i - raw_pred_i|
    """

    hessians_are_constant = True
    # This variable indicates whether the loss requires the leaves values to
    # be updated once the tree has been trained. The trees are trained to
    # predict a Newton-Raphson step (see grower._finalize_leaf()). But for
    # some losses (e.g. least absolute deviation) we need to adjust the tree
    # values to account for the "line search" of the gradient descent
    # procedure. See the original paper Greedy Function Approximation: A
    # Gradient Boosting Machine by Friedman
    # (https://statweb.stanford.edu/~jhf/ftp/trebst.pdf) for the theory.
    need_update_leaves_values = True

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        loss = np.abs(y_true - raw_predictions)
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        return np.median(y_train)

    @staticmethod
    def inverse_link_function(raw_predictions):
        return raw_predictions

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        gradients = gradients.reshape(-1)
        _update_gradients_least_absolute_deviation(gradients, y_true,
                                                   raw_predictions)

    def update_leaves_values(self, grower, y_true, raw_predictions):
        # Update the values predicted by the tree with
        # median(y_true - raw_predictions).
        # See note about need_update_leaves_values in BaseLoss.

        # TODO: ideally this should be computed in parallel over the leaves
        # using something similar to _update_raw_predictions(), but this
        # requires a cython version of median()
        for leaf in grower.finalized_leaves:
            indices = leaf.sample_indices
            median_res = np.median(y_true[indices] - raw_predictions[indices])
            leaf.value = grower.shrinkage * median_res
            # Note that the regularization is ignored here


class BinaryCrossEntropy(BaseLoss):
    """Binary cross-entropy loss, for binary classification.

    For a given sample x_i, the binary cross-entropy loss is defined as the
    negative log-likelihood of the model which can be expressed as::

        loss(x_i) = log(1 + exp(raw_pred_i)) - y_true_i * raw_pred_i

    See The Elements of Statistical Learning, by Hastie, Tibshirani, Friedman,
    section 4.4.1 (about logistic regression).
    """

    hessians_are_constant = False
    inverse_link_function = staticmethod(expit)

    def __call__(self, y_true, raw_predictions, average=True):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        # logaddexp(0, x) = log(1 + exp(x))
        loss = np.logaddexp(0, raw_predictions) - y_true * raw_predictions
        return loss.mean() if average else loss

    def get_baseline_prediction(self, y_train, prediction_dim):
        if prediction_dim > 2:
            raise ValueError(
                "loss='binary_crossentropy' is not defined for multiclass"
                " classification with n_classes=%d, use"
                " loss='categorical_crossentropy' instead" % prediction_dim)
        proba_positive_class = np.mean(y_train)
        eps = np.finfo(y_train.dtype).eps
        proba_positive_class = np.clip(proba_positive_class, eps, 1 - eps)
        # log(x / 1 - x) is the anti function of sigmoid, or the link function
        # of the Binomial model.
        return np.log(proba_positive_class / (1 - proba_positive_class))

    def update_gradients_and_hessians(self, gradients, hessians, y_true,
                                      raw_predictions):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        gradients = gradients.reshape(-1)
        hessians = hessians.reshape(-1)
        _update_gradients_hessians_binary_crossentropy(
            gradients, hessians, y_true, raw_predictions)

    def predict_proba(self, raw_predictions):
        # shape (1, n_samples) --> (n_samples,). reshape(-1) is more likely to
        # return a view.
        raw_predictions = raw_predictions.reshape(-1)
        proba = np.empty((raw_predictions.shape[0], 2), dtype=Y_DTYPE)
        proba[:, 1] = expit(raw_predictions)
        proba[:, 0] = 1 - proba[:, 1]
        return proba


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
        _update_gradients_hessians_categorical_crossentropy(
            gradients, hessians, y_true, raw_predictions)

    def predict_proba(self, raw_predictions):
        # TODO: This could be done in parallel
        # compute softmax (using exp(log(softmax)))
        proba = np.exp(raw_predictions -
                       logsumexp(raw_predictions, axis=0)[np.newaxis, :])
        return proba.T


_LOSSES = {
    'least_squares': LeastSquares,
    'least_absolute_deviation': LeastAbsoluteDeviation,
    'binary_crossentropy': BinaryCrossEntropy,
    'categorical_crossentropy': CategoricalCrossEntropy
}
