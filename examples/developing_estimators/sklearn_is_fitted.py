"""
========================================
`__sklearn_is_fitted__` as Developer API
========================================

The `__sklearn_is_fitted__` method is a convention used in scikit-learn for
checking whether an estimator object has been fitted or not. This method is
typically implemented in custom estimator classes that are built on top of
scikit-learn's base classes like `BaseEstimator` or its subclasses.

Developers should use :func:`~sklearn.utils.validation.check_is_fitted`
at the beginning of all methods except `fit`. If they need to customize or
speed-up the check, they can implement the `__sklearn_is_fitted__` method as
shown below.

In this example the custom estimator showcases the usage of the
`__sklearn_is_fitted__` method and the `check_is_fitted` utility function
as developer APIs. The `__sklearn_is_fitted__` method checks fitted status
by verifying the presence of the `_is_fitted` attribute.
"""

# %%
# An example custom estimator implementing a simple classifier
# ------------------------------------------------------------
# This code snippet defines a custom estimator class called `CustomEstimator`
# that extends both the `BaseEstimator` and `ClassifierMixin` classes from
# scikit-learn and showcases the usage of the `__sklearn_is_fitted__` method
# and the `check_is_fitted` utility function.

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted


class CustomEstimator(BaseEstimator, ClassifierMixin):
    def __init__(self, parameter=1):
        self.parameter = parameter

    def fit(self, X, y):
        """
        Fit the estimator to the training data.
        """
        self.classes_ = sorted(set(y))
        # Custom attribute to track if the estimator is fitted
        self._is_fitted = True
        return self

    def predict(self, X):
        """
        Perform Predictions

        If the estimator is not fitted, then raise NotFittedError
        """
        check_is_fitted(self)
        # Perform prediction logic
        predictions = [self.classes_[0]] * len(X)
        return predictions

    def score(self, X, y):
        """
        Calculate Score

        If the estimator is not fitted, then raise NotFittedError
        """
        check_is_fitted(self)
        # Perform scoring logic
        return 0.5

    def __sklearn_is_fitted__(self):
        """
        Check fitted status and return a Boolean value.
        """
        return hasattr(self, "_is_fitted") and self._is_fitted
