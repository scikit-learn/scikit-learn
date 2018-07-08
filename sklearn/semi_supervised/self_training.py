import numpy as np

from ..base import BaseEstimator
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils import safe_mask


def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if not hasattr(estimator, "predict_proba"):
        raise ValueError("The base estimator should implement predict_proba!")


class SelfTraining(BaseEstimator):

    """Self-Training classifier

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and `predict_proba`.

    threshold : float
        Threshold above which predictions are added to the labeled dataset

    max_iter : integer
        Change maximum number of iterations allowed

    """
    def __init__(self, estimator, threshold=0.7, max_iter=500):
        self.estimator = estimator
        self.threshold = threshold
        self.max_iter = max_iter

    def fit(self, X, y):
        """
        Fits SelfTraining Estimator to dataset

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data
        y : array-like, shape = (n_samples, 1)
            array representing the labels

        Returns
        -------
        self: returns an instance of self.
        """
        X, y = check_X_y(X, y)
        _check_estimator(self.estimator)

        # Data usable for supervised training
        X_labeled = X[safe_mask(X, np.where(y != -1))][0]
        y_labeled = y[safe_mask(y, np.where(y != -1))][0]

        # Unlabeled data
        X_unlabeled = X[safe_mask(X, np.where(y == -1))][0]
        y_unlabeled = y[safe_mask(y, np.where(y == -1))][0]

        iter = 0
        while (len(X_labeled) < len(X) and iter < self.max_iter):
            iter += 1
            self.estimator.fit(X_labeled, y_labeled)

            # Select prediction where confidence is above the threshold
            pred = self.predict(X_unlabeled)
            max_proba = np.max(self.predict_proba(X_unlabeled), axis=1)
            confident = np.where(max_proba > self.threshold)[0]

            # Add newly labeled confident predictions to the dataset
            X_labeled = np.append(X_labeled, X_unlabeled[confident], axis=0)
            y_labeled = np.append(y_labeled, pred[confident], axis=0)

            # Remove already labeled data from unlabeled dataset
            X_unlabeled = np.delete(X_unlabeled, confident, axis=0)
            y_unlabeled = np.delete(y_unlabeled, confident, axis=0)

        self.estimator.fit(X_labeled, y_labeled)
        return self.estimator

    def predict(self, X):
        """Predict on a dataset.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, 1)
            array with predicted labels
        """
        check_is_fitted(self, 'estimator')
        X = check_array(X)
        return self.estimator.predict(X)

    def predict_proba(self, X):
        """Predict probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, n_features)
            array with prediction probabilities
        """
        _check_estimator(self.estimator)
        check_is_fitted(self, 'estimator')
        return self.estimator.predict_proba(X)
