"""
===============================================
__sklearn_is_fitted__ as Internal Developer API
===============================================

The `__sklearn_is_fitted__` method is like a internal developer API which
used internally by scikit-learn to check if an estimator object has been
fitted or not. `check_is_fitted` utilizes the `__sklearn_is_fitted__`
method of estimator, if it exists, to check the fitted status.

Users are encouraged to use the
:func:`~sklearn.sklearn.utils.validation.check_is_fitted` function
provided by scikit-learn. But by relying on `__sklearn_is_fitted__` when
available, `check_is_fitted` allows the custom estimator to define their
own convention for indicating fitted status.

In this example the custom estimator showcases the usage of the
`__sklearn_is_fitted__` method and `check_is_fitted` utility
function as internal developer APIs. The `__sklearn_is_fitted__`
method checks whether the estimator has been fitted by verifying
the presence of the `_is_fitted` attribute.
"""

# %%
# An example custom estimator implementing a simple classifier
# ------------------------------------------------------------
# This code snippet defines a custom estimator class called `CustomEstimator`
# that extends the `BaseEstimator` and `ClassifierMixin` classes from
# scikit-learn. It showcases the usage of the `__sklearn_is_fitted__` method
# and the `check_is_fitted` utility function as internal developer APIs.

# Author: Kushan <kushansharma1@gmail.com>
#
# License: BSD 3 clause

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted


class CustomEstimator(BaseEstimator, ClassifierMixin):
    def __init__(self, parameter=1):
        self._is_fitted = False
        self.parameter = parameter

    def fit(self, X, y):
        """
        Fit the estimator to the training data.
        Parameters:
        -----------
        - X : array-like of shape (n_samples, n_features)
            The input samples.
        - y : array-like of shape (n_samples,)
            The target values.
        Returns:
        --------
        - self : object
            Returns self.
        """
        self.classes_ = sorted(set(y))
        # Custom attribute to track if the estimator is fitted
        self._is_fitted = True
        return self

    def predict(self, X):
        """
        Perform predictions using the fitted estimator.
        Parameters:
        -----------
        - X : array-like of shape (n_samples, n_features)
            The input samples.
        Returns:
        --------
        - predictions : list
            The predicted labels for the input samples.
        Raises:
        -------
        - NotFittedError: If the estimator is not fitted.
        """
        check_is_fitted(self)
        # Perform prediction logic
        predictions = [self.classes_[0]] * len(X)
        return predictions

    def score(self, X, y):
        """
        Calculate the score of the estimator.
        Parameters:
        -----------
        - X : array-like of shape (n_samples, n_features)
            The input samples.
        - y : array-like of shape (n_samples,)
            The target values.
        Returns:
        --------
        - score : float
            The calculated score.
        Raises:
        -------
        - NotFittedError: If the estimator is not fitted.
        """
        check_is_fitted(self)
        # Perform scoring logic
        return 0.5

    def __sklearn_is_fitted__(self):
        """Internal developer API. Returns True if the estimator is fitted,
        False otherwise."""
        return hasattr(self, '_is_fitted') and self._is_fitted


# %%
# Example usage
# -------------
# The provided code demonstrates the usage of the `CustomEstimator`
# class by fitting the estimator to a sample dataset, making predictions,
# and calculating the score.

X = [[1, 2, 3], [4, 5, 6]]
y = [0, 1]

# %%
# Create an instance of the CustomEstimator

estimator = CustomEstimator()

# %%
# try to predict before data was fitted
#
# `estimator.predict(X)` Output: NotFittedError
# check if the estimater is fitted or not

estimator.__sklearn_is_fitted__()

# %%
# Now fit the estimator to the training data

estimator.fit(X, y)

# %%
# Perform predictions on the same dataset
# and calculate the score of the estimator on the training data


predictions = estimator.predict(X)
score = estimator.score(X, y)

# this time it will return True
estimator.__sklearn_is_fitted__()

# %%
# Get the sample results

print("Predictions:", predictions)
print("Score:", score)
