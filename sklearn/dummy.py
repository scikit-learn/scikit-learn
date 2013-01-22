
# Author: Mathieu Blondel <mathieu@mblondel.org>
#         Arnaud Joly <a.joly@ulg.ac.be>
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
    `classes_` : array or list of array of shape = [n_classes]
        Class labels for each output.

    `n_classes_` : array or list of array of shape = [n_classes]
        Number of label for each output.

    `class_prior_` : array or list of array of shape = [n_classes]
        Probability of each class for each output.

    `n_outputs_` : int,
        Number of outputs.

    `outputs_2d_` : bool,
        True if the output at fit is 2d, else false.
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

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """
        if self.strategy not in ("most_frequent", "stratified", "uniform"):
            raise ValueError("Unknown strategy type.")

        y = np.atleast_1d(y)
        self.output_2d_ = y.ndim == 2

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]
        self.classes_ = []
        self.n_classes_ = []
        self.class_prior_ = []

        for k in xrange(self.n_outputs_):
            classes, y_k = unique(y[:, k], return_inverse=True)
            self.classes_.append(classes)
            self.n_classes_.append(classes.shape[0])
            self.class_prior_.append(np.bincount(y_k) / float(y_k.shape[0]))

        if self.n_outputs_ == 1 and not self.output_2d_:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]
            self.class_prior_ = self.class_prior_[0]

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
        y : array, shape = [n_samples] or [n_samples, n_outputs]
            Predicted target values for X.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]
        rs = check_random_state(self.random_state)

        n_classes_ = self.n_classes_
        classes_ = self.classes_
        class_prior_ = self.class_prior_
        if self.n_outputs_ == 1:
            # Get same type even for self.n_outputs_ == 1
            n_classes_ = [n_classes_]
            classes_ = [classes_]
            class_prior_ = [class_prior_]

        # Compute probability only once
        if self.strategy == "stratified":
            proba = self.predict_proba(X)
            if self.n_outputs_ == 1:
                proba = [proba]

        y = []
        for k in xrange(self.n_outputs_):
            if self.strategy == "most_frequent":
                ret = np.ones(n_samples, dtype=int) * class_prior_[k].argmax()

            elif self.strategy == "stratified":
                ret = proba[k].argmax(axis=1)

            elif self.strategy == "uniform":
                ret = rs.randint(n_classes_[k], size=n_samples)

            y.append(classes_[k][ret])

        y = np.vstack(y).T
        if self.n_outputs_ == 1 and not self.output_2d_:
            y = np.ravel(y)

        return y

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
        P : array-like or list of array-lke of shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered arithmetically, for each
            output.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]
        rs = check_random_state(self.random_state)

        n_classes_ = self.n_classes_
        classes_ = self.classes_
        class_prior_ = self.class_prior_
        if self.n_outputs_ == 1 and not self.output_2d_:
            # Get same type even for self.n_outputs_ == 1
            n_classes_ = [n_classes_]
            classes_ = [classes_]
            class_prior_ = [class_prior_]

        P = []
        for k in xrange(self.n_outputs_):
            if self.strategy == "most_frequent":
                ind = np.ones(n_samples, dtype=int) * class_prior_[k].argmax()
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0

            elif self.strategy == "stratified":
                out = rs.multinomial(1, class_prior_[k], size=n_samples)

            elif self.strategy == "uniform":
                out = np.ones((n_samples, n_classes_[k]), dtype=np.float64)
                out /= n_classes_[k]

            P.append(out)

        if self.n_outputs_ == 1 and not self.output_2d_:
            P = P[0]

        return P

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
        P : array-like or list of array-like of shape = [n_samples, n_classes]
            Returns the log probability of the sample for each class in
            the model, where classes are ordered arithmetically for each
            output.
        """
        proba = self.predict_proba(X)
        if self.n_outputs_ == 1:
            return np.log(proba)
        else:
            return [np.log(p) for p in proba]


class DummyRegressor(BaseEstimator, RegressorMixin):
    """
    DummyRegressor is a regressor that always predicts the mean of the training
    targets.

    This regressor is useful as a simple baseline to compare with other
    (real) regressors. Do not use it for real problems.

    Attributes
    ----------
    `y_mean_` : float or array of shape [n_outputs]
        Mean of the training targets.

    `n_outputs_` : int,
        Number of outputs.

    `outputs_2d_` : bool,
        True if the output at fit is 2d, else false.
    """

    def fit(self, X, y):
        """Fit the random regressor.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """
        y = safe_asarray(y)
        self.y_mean_ = np.reshape(np.mean(y, axis=0), (1, -1))
        self.n_outputs_ = np.size(self.y_mean_)  # y.shape[1] is not safe
        self.output_2d_ = (y.ndim == 2)
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
        y : array, shape = [n_samples]  or [n_samples, n_outputs]
            Predicted target values for X.
        """
        if not hasattr(self, "y_mean_"):
            raise ValueError("DummyRegressor not fitted.")

        X = safe_asarray(X)
        n_samples = X.shape[0]
        y = np.ones((n_samples, 1)) * self.y_mean_

        if self.n_outputs_ == 1 and not self.output_2d_:
            y = np.ravel(y)

        return y
