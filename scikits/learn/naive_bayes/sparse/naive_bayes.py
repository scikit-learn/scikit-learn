""" Naive Bayes classifiers for sparse data.
"""

# Author: Amit Aides <amitibo@tx.technion.ac.il>
#
# License: BSD Style.
import numpy as np

from ...base import BaseEstimator, ClassifierMixin


class MultinomialNB(BaseEstimator, ClassifierMixin):
    """
    Multinomial Naive Bayes for sparse matrices

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
    alpha: float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).

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
    >>> from scikits.learn.naive_bayes import MultinomialNB
    >>> clf = MultinomialNB()
    >>> clf.fit(X, Y)
    MultinomialNB(alpha=1.0)
    >>> print clf.predict(X[2])
    3

    """

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    def fit(self, X, y):
        """Fit Multinomial Naive Bayes according to X, y

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
            row_ind = np.nonzero(y == yi)[0]
            N_c_i_temp.append(np.array(X[row_ind, :].sum(axis=0)).ravel())
            theta_c.append(np.float(np.sum(y == yi)) / y.size)

        N_c_i = np.array(N_c_i_temp)
        N_c = np.sum(N_c_i, axis=1)

        #
        # Smoothing coefficients
        #
        alpha_i = self.alpha
        alpha = alpha_i * X.shape[1]

        #
        # Estimate the parameters of the distribution
        #
        self.theta_c_i = (N_c_i + alpha_i) / (N_c.reshape(-1, 1) + alpha)
        self.theta_c = np.array(theta_c)

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

        joint_log_likelihood = self._joint_log_likelihood(X)
        y_pred = self.unique_y[np.argmax(joint_log_likelihood, axis=0)]

        return y_pred

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""

        joint_log_likelihood = []
        for i in range(self.unique_y.size):
            jointi = np.log(self.theta_c[i])
            n_ij = np.log(self.theta_c_i[i]) * X.T
            joint_log_likelihood.append(jointi + n_ij)

        joint_log_likelihood = np.array(joint_log_likelihood)

        return joint_log_likelihood

    def _mininf(self, X, axis=None):
        """Calculate the minimum of a matrix ignoring -inf values"""

        A = X.copy()
        A[np.isinf(X)] = np.inf
        return np.min(X, axis=axis)

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
