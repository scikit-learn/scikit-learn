"""
===================================================
``__sklearn_is_fitted__`` as Internal Developer API
===================================================

Developers should use :func:`~sklearn.sklearn.utils.validation.check_is_fitted`
at the beginning of all methods except ``fit``. If they need to customize or
speed-up the check, the can implement the ``__sklearn_is_fitted__`` method as
shown bellow.

In this example the custom estimator showcases the usage of the
``__sklearn_is_fitted__`` method and ``check_is_fitted`` utility
function as developer APIs. The ``__sklearn_is_fitted__``
method checks whether the estimator has been fitted by verifying
the presence of the `_is_fitted` attribute.
"""

# %%
# An example custom estimator implementing a simple classifier
# ------------------------------------------------------------
# This code snippet defines a custom estimator class called `CustomEstimator`
# that extends the `BaseEstimator` and `ClassifierMixin` classes from
# scikit-learn. It showcases the usage of the ``__sklearn_is_fitted__`` method
# and the `check_is_fitted` utility function as developer APIs.

# Author: Kushan <kushansharma1@gmail.com>
#
# License: BSD 3 clause

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted


class CustomEstimator(BaseEstimator, ClassifierMixin):
    def __init__(self):
        self._is_fitted = False

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
        
        If estimator is not fitted, then raise NotFittedError
        """
        check_is_fitted(self)
        # Perform prediction logic
        predictions = [self.classes_[0]] * len(X)
        return predictions

    def score(self, X, y):
        """
        Calculate Score
        
        If estimator is not fitted, then raise NotFittedError
        """
        check_is_fitted(self)
        # Perform scoring logic
        return 0.5

    def __sklearn_is_fitted__(self):
        """Returns True if the estimator is fitted, False otherwise."""
        return hasattr(self, "_is_fitted") and self._is_fitted
