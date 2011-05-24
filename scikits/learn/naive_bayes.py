""" Naives Bayes classifiers.
"""

# Author: Vincent Michel <vincent.michel@inria.fr>
#         Minor fixes by Fabian Pedregosa
#
# License: BSD Style.
import numpy as np

from .base import BaseEstimator, ClassifierMixin


class GNB(BaseEstimator, ClassifierMixin):
    """
    Gaussian Naive Bayes (GNB)

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    y : array, shape = [n_samples]
        Target vector relative to X

    Attributes
    ----------
    proba_y : array, shape = [n_classes]
        probability of each class.

    theta : array, shape [n_classes * n_features]
        mean of each feature for the different class

    sigma : array, shape [n_classes * n_features]
        variance of each feature for the different class


    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    predict_proba(X) : array
        Predict the probability of each class using the model.

    predict_log_proba(X) : array
        Predict the log-probability of each class using the model.


    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> Y = np.array([1, 1, 1, 2, 2, 2])
    >>> from scikits.learn.naive_bayes import GNB
    >>> clf = GNB()
    >>> clf.fit(X, Y)
    GNB()
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------

    """

    def fit(self, X, y):
        """Fit Gaussian Naive Bayes according to X, y

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.


        Returns
        -------
        self : object
            Returns self.
        """

        X = np.asanyarray(X)
        y = np.asanyarray(y)

        theta = []
        sigma = []
        proba_y = []
        unique_y = np.unique(y)
        for yi in unique_y:
            theta.append(np.mean(X[y == yi, :], 0))
            sigma.append(np.var(X[y == yi, :], 0))
            proba_y.append(np.float(np.sum(y == yi)) / np.size(y))
        self.theta = np.array(theta)
        self.sigma = np.array(sigma)
        self.proba_y = np.array(proba_y)
        self.unique_y = unique_y
        return self

    def predict(self, X):
        """
        Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        X = np.asanyarray(X)
        y_pred = self.unique_y[np.argmax(self.predict_proba(X), 1)]
        return y_pred

    def _joint_log_likelihood(self, X):
        joint_log_likelihood = []
        for i in range(np.size(self.unique_y)):
            jointi = np.log(self.proba_y[i])
            n_ij = - 0.5 * np.sum(np.log(np.pi * self.sigma[i, :]))
            n_ij -= 0.5 * np.sum(((X - self.theta[i, :]) ** 2) / \
                                    (self.sigma[i, :]), 1)
            joint_log_likelihood.append(jointi + n_ij)
        joint_log_likelihood = np.array(joint_log_likelihood).T
        return joint_log_likelihood

    def predict_proba(self, X):
        """
        Return probability estimates for the test vector X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.
        """
        X = np.asanyarray(X)
        joint_log_likelihood = self._joint_log_likelihood(X)
        proba = np.exp(joint_log_likelihood)
        proba = proba / np.sum(proba, 1)[:, np.newaxis]
        return proba

    def predict_log_proba(self, X):
        """
        Return log-probability estimates for the test vector X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array-like, shape = [n_samples, n_classes]
            Returns the log-probability of the sample for each class
            in the model, where classes are ordered by arithmetical
            order.
        """
        log_proba = self._joint_log_likelihood(X)
        # Compute a sum of logs without underflow. Equivalent to:
        # log_proba -= np.log(np.sum(np.exp(log_proba), axis=1))[:, np.newaxis]
        B = np.max(log_proba, axis=1)[:, np.newaxis]
        logaB = log_proba - B
        sup = logaB > -np.inf
        aB = np.zeros_like(logaB)
        aB[sup] = np.exp(logaB[sup])
        log_proba -= np.log(np.sum(aB, axis=1))[:, np.newaxis] + B
        return log_proba
