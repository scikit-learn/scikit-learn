"""
=================================================================
Internal Developer API
=================================================================

The `__sklearn_is_fitted__` method is like a internal developer API used
internally by scikit-learn to check if an estimator object has been
fitted or not. Users are encouraged to use the `check_is_fitted` function 
provided by scikit-learn, which internally relies on `__sklearn_is_fitted__`,
to check the fitted status of an estimator in their code.

In this example the custom estimator showcases the usage of the 
`__sklearn_is_fitted__` method and `check_is_fitted` utility 
function as internal developer APIs. The `__sklearn_is_fitted__` 
method checks whether the estimator has been fitted by verifying 
the presence of the `_is_fitted` attribute. It can be used internally
 within the estimator's methods to ensure that the estimator is fitted
 before performing predictions or scoring.

@author: Kushan
"""

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted

class CustomEstimator(BaseEstimator, ClassifierMixin):
    """Custom estimator implementing a simple classifier."""
    
    def __init__(self, parameter=1):
        self.parameter = parameter

    def fit(self, X, y):
        """
        Fit the estimator to the training data.
        
        Parameters:
        - X : array-like of shape (n_samples, n_features)
            The input samples.
        - y : array-like of shape (n_samples,)
            The target values.
            
        Returns:
        - self : object
            Returns self.
        """
        self.classes_ = sorted(set(y))
        self._is_fitted = True  # Custom attribute to track if the estimator is fitted
        return self

    def predict(self, X):
        """
        Perform predictions using the fitted estimator.
        
        Parameters:
        - X : array-like of shape (n_samples, n_features)
            The input samples.
            
        Returns:
        - predictions : list
            The predicted labels for the input samples.
        
        Raises:
        - NotFittedError: If the estimator is not fitted.
        """
        check_is_fitted(self, '_is_fitted')
        predictions = [self.classes_[0]] * len(X)  # Perform prediction logic
        return predictions

    def score(self, X, y):
        """
        Calculate the score of the estimator.
        
        Parameters:
        - X : array-like of shape (n_samples, n_features)
            The input samples.
        - y : array-like of shape (n_samples,)
            The target values.
            
        Returns:
        - score : float
            The calculated score.
        
        Raises:
        - NotFittedError: If the estimator is not fitted.
        """
        check_is_fitted(self, '_is_fitted')
        # Perform scoring logic
        return 0.5

    def __sklearn_is_fitted__(self):
        """
        Internal developer API. Returns True if the estimator is fitted, False otherwise.
        """
        return hasattr(self, '_is_fitted')
