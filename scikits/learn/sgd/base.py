# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Stochastic Gradient Descent (SGD) abstract base class"""

import numpy as np

from ..base import BaseEstimator, ClassifierMixin
from .sgd_fast import Hinge, Log, ModifiedHuber

class BaseSGD(BaseEstimator, ClassifierMixin):
    """Base class for dense and sparse SGD"""

    def __init__(self, loss="hinge", penalty='l2', alpha=0.0001,
                 rho=0.85, coef_=None, intercept_=None,
                 fit_intercept=True, n_iter=5, shuffle=False,
                 verbose=0, n_jobs=-1):
        self.loss = loss
        self.penalty = penalty
        self.alpha = alpha
        self.rho = rho
        self.coef_ = np.asarray(coef_) if coef_ is not None else None
        self.intercept_ = intercept_
        if self.intercept_ is not None:
            self.intercept_ = np.asarray(intercept_)
        self.fit_intercept = fit_intercept
        self.n_iter = int(n_iter)
        if self.n_iter <= 0:
            raise ValueError("n_iter must be greater than zero.")
        if not isinstance(shuffle, bool):
            raise ValueError("shuffle must be either True or False")
        self.shuffle = shuffle
        self.verbose = verbose
        self.n_jobs = n_jobs
        self._get_loss_function()
        self._get_penalty_type()

    def _get_loss_function(self):
        """Get concrete LossFunction"""
        loss_functions = {
            "hinge" : Hinge(),
            "log" : Log(),
            "modifiedhuber" : ModifiedHuber(),
        }
        try:
            self.loss_function = loss_functions[self.loss]
        except KeyError:
            raise ValueError("The loss %s is not supported. " % self.loss)

    def _get_penalty_type(self):
        penalty_types = {"l2": 2, "l1": 1, "elasticnet": 3}
        try:
            self.penalty_type = penalty_types[self.penalty]
            if self.penalty_type == 2:
                self.rho = 1.0
            elif self.penalty_type == 1:
                self.rho = 0.0
        except KeyError:
            raise ValueError("Penalty %s is not supported. " % self.penalty)

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : array or scipy.sparse matrix of shape [n_samples, n_features]
           Whether the numpy.array or scipy.sparse matrix is accepted dependes
           on the actual implementation

        Returns
        -------
        array, shape = [n_samples]
           Array containing the predicted class labels.
        """
        scores = self.predict_margin(X)
        if self.classes.shape[0] == 2:
            indices = np.array(scores > 0, dtype=np.int)
        else:
            indices = scores.argmax(axis=1)
        return self.classes[np.ravel(indices)]

    def predict_proba(self, X):
        """Predict class membership probability

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] if n_classes == 2 else [n_samples, n_classes]
            Contains the membership probabilities of the positive class.

        """
        if (isinstance(self.loss_function, Log) and
            self.classes.shape[0] == 2):
            return 1.0 / (1.0 + np.exp(-self.predict_margin(X)))
        else:
            raise NotImplementedError("%s loss does not provide "
                                      "this functionality" % self.loss)

