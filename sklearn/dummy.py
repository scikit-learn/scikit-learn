
# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .utils import check_random_state
from .utils.fixes import unique
from .utils.validation import safe_asarray


class DummyClassifier(BaseEstimator, ClassifierMixin):
    """
    DummyClassifier is a classifier that makes predictions using simple rules.

    This classifier is useful as a simple baseline to compare with other
    (real) classifiers. Do not use it for real problems.

    Parameters
    ----------
    strategy: str
        Strategy to use to generate predictions.
            * "stratified": generates predictions by respecting the training
              set's class distribution.
            * "most_frequent": always predicts the most frequent label in the
              training set.
            * "uniform": generates predictions uniformly at random.

    random_state: int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use.

    Attributes
    ----------
    `classes_` : array, shape = [n_classes]
        Class labels.

    `class_prior_` : array, shape = [n_classes]
        Probability of each class.
    """

    def __init__(self, strategy="stratified", random_state=None):
        self.strategy = strategy
        self.random_state = random_state

    def fit(self, X, y):
        """Fit the random classifier.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """
        if self.strategy not in ("most_frequent", "stratified", "uniform"):
            raise ValueError("Unknown strategy type.")

        self.classes_, y = unique(y, return_inverse=True)
        self.class_prior_ = np.bincount(y) / float(y.shape[0])
        return self

    def predict(self, X):
        """
        Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array, shape = [n_samples]
            Predicted target values for X.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]
        rs = check_random_state(self.random_state)

        if self.strategy == "most_frequent":
            ret = np.ones(n_samples, dtype=int) * self.class_prior_.argmax()
        elif self.strategy == "stratified":
            ret = self.predict_proba(X).argmax(axis=1)
        elif self.strategy == "uniform":
            ret = rs.randint(len(self.classes_), size=n_samples)

        return self.classes_[ret]

    def predict_proba(self, X):
        """
        Return probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered arithmetically.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]
        n_classes = len(self.classes_)
        rs = check_random_state(self.random_state)

        if self.strategy == "most_frequent":
            ind = np.ones(n_samples, dtype=int) * self.class_prior_.argmax()
            out = np.zeros((n_samples, n_classes), dtype=np.float64)
            out[:, ind] = 1.0
        elif self.strategy == "stratified":
            out = rs.multinomial(1, self.class_prior_, size=n_samples)
        elif self.strategy == "uniform":
            out = np.ones((n_samples, n_classes), dtype=np.float64)
            out /= n_classes

        return out

    def predict_log_proba(self, X):
        """
        Return log probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the log probability of the sample for each class in
            the model, where classes are ordered arithmetically.
        """
        return np.log(self.predict_proba(X))


class DummyRegressor(BaseEstimator, RegressorMixin):
    """
    DummyRegressor is a regressor that always predicts the mean of the training
    targets.

    This regressor is useful as a simple baseline to compare with other
    (real) regressors. Do not use it for real problems.

    Attributes
    ----------
    `y_mean_` : float
        Mean of the training targets.
    """

    def fit(self, X, y):
        """Fit the random regressor.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """
        self.y_mean_ = np.mean(y)
        return self

    def predict(self, X):
        """
        Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array, shape = [n_samples]
            Predicted target values for X.
        """
        if not hasattr(self, "y_mean_"):
            raise ValueError("DummyRegressor not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]

        return np.ones(n_samples) * self.y_mean_
