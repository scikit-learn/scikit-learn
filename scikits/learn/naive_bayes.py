""" Naives Bayes classifiers.
"""

# Author: Vincent Michel <vincent.michel@inria.fr>
#         Amit Aides <amitibo@tx.technion.ac.il>
#
# License: BSD Style.
import numpy as np

from .base import BaseEstimator, ClassifierMixin

eps = np.finfo(np.float).eps

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
    proba_y : array, shape = nb of classes
              probability of each class.
    theta : array of shape nb_class*nb_features
            mean of each feature for the different class
    sigma : array of shape nb_class*nb_features
            variance of each feature for the different class


    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    predict_proba(X) : array
        Predict the probability of each class using the model.

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

    def __init__(self):
        pass

    def fit(self, X, y):
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
        joint_log_likelihood = self._joint_log_likelihood(X)
        proba = np.exp(joint_log_likelihood)
        proba = proba / np.sum(proba, 1)[:, np.newaxis]
        return proba

    def predict_log_proba(self, X):
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


class MNNB(BaseEstimator, ClassifierMixin):
    """
    Multinomial Naive Bayes (MNNB)

    The Multinomial Naive Bayes classifier is suitable for text classification.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.
    y : array, shape = [n_samples]
        Target vector relative to X

    Parameters
    ----------
    alpha_i: float, optional (default=1.0)
        smoothing constant.

    alpha_ratio: float, optional (default=1.0)
        smoothing ratio.

    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    predict_proba(X) : array
        Predict the probability of each class using the model.

    predict_log_proba(X) : array
        Predict the log probability of each class using the model.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randint( 5, size=(6, 100) )
    >>> Y = np.array([1, 2, 3, 4, 5, 6])
    >>> from scikits.learn.naive_bayes import MNNB
    >>> clf = MNNB()
    >>> clf.fit(X, Y)
    MNNB(alpha_ratio=1.0, alpha_i=1.0)
    >>> print clf.predict(X[2])
    3

    See also
    --------

    """

    def __init__(self, alpha_i=1.0, alpha_ratio=1.0):

        self.alpha_i = alpha_i
        self.alpha_ratio = alpha_ratio

    def fit(self, X, y):
        """Fit the Multinomial distribution"""

        #
        # N_c is the count of all words in all documents of class c.
        # N_c_i is the a count of word i in all documents of class c.
        # theta_c is the prior empirical probability of a document of class c.
        # theta_c_i is the (smoothened) empirical likelihood of word i
        # given a document of class c.
        #
        N_c_i_temp = []
        theta_c = []
        self.unique_y = np.unique(y)

        for yi in self.unique_y:
            N_c_i_temp.append(np.sum(X[y == yi, :], 0))
            theta_c.append(np.float(np.sum(y == yi)) / y.size)

        N_c_i = np.array(N_c_i_temp)
        N_c = np.sum(N_c_i, axis=1)

        #
        # Smoothing coefficients
        #
        alpha_i = self.alpha_i
        alpha = self.alpha_ratio * alpha_i * X.shape[1]

        #
        # Estimate the parameters of the distribution
        #
        self.theta_c_i = (N_c_i + alpha_i) / (N_c.reshape(-1, 1) + alpha)
        self.theta_c = np.array(theta_c)

        return self

    def predict(self, X):
        """Predict the classification of samples X"""

        joint_log_likelihood = self._joint_log_likelihood(X)
        y_pred = self.unique_y[np.argmax(joint_log_likelihood, axis=0)]

        return y_pred

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""

        joint_log_likelihood = []
        for i in range(self.unique_y.size):
            jointi = np.log(self.theta_c[i])
            n_ij = np.dot(np.log(self.theta_c_i[i]), X.T)
            joint_log_likelihood.append(jointi + n_ij)

        joint_log_likelihood = np.array(joint_log_likelihood)

        return joint_log_likelihood

    def _mininf(self, X, axis=None):
        """Calculate the minimum of a matrix ignoring -inf values"""
        
        A = X.copy()
        A[np.isinf(X)] = np.inf
        return np.min(X, axis=axis)
        
    def predict_proba(self, X):
        """Predict the posterior probability of samples X"""

        joint_log_likelihood = self._joint_log_likelihood(X)
        
        #
        # The _joint_log_likelihood has very low values that create underflow
        # in the computation of the exponent. Therefore I 'fix' it by adding
        # a minimal value.
        #
        fix = self._mininf(joint_log_likelihood, axis=1)[:, np.newaxis]
        loga_fix = joint_log_likelihood - fix
        proba_fix = np.exp(loga_fix)
        proba = proba_fix / np.sum(proba_fix, 1)[:, np.newaxis]

        return proba

    def predict_log_proba(self, X):
        """Predict the posterior log probability of samples X"""

        joint_log_likelihood = self._joint_log_likelihood(X)
        
        #
        # The _joint_log_likelihood has very low values that create underflow
        # in the computation of the exponent. Therefore I 'fix' it by adding
        # a minimal value.
        #
        fix = self._mininf(joint_log_likelihood, axis=1)[:, np.newaxis]
        loga_fix = joint_log_likelihood - fix
        proba_fix = np.exp(loga_fix)
        log_proba = loga_fix - np.log(np.sum(proba_fix, axis=1))[:, np.newaxis]
        
        return log_proba
