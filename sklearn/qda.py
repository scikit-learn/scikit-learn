"""
Quadratic Discriminant Analysis
"""

# Author: Matthieu Perrot <matthieu.perrot@gmail.com>
#
# License: BSD Style.

import warnings

import numpy as np

from .base import BaseEstimator, ClassifierMixin
from .utils.fixes import unique
from .utils import check_arrays

__all__ = ['QDA']


class QDA(BaseEstimator, ClassifierMixin):
    """
    Quadratic Discriminant Analysis (QDA)

    A classifier with a quadratic decision boundary, generated
    by fitting class conditional densities to the data
    and using Bayes' rule.

    The model fits a Gaussian density to each class.

    Parameters
    ----------
    priors : array, optional, shape = [n_classes]
        Priors on classes

    Attributes
    ----------
    `means_` : array-like, shape = [n_classes, n_features]
        Class means
    `priors_` : array-like, shape = [n_classes]
        Class priors (sum to 1)
    `covariances_` : list of array-like, shape = [n_features, n_features]
        Covariance matrices of each class

    Examples
    --------
    >>> from sklearn.qda import QDA
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = QDA()
    >>> clf.fit(X, y)
    QDA(priors=None)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    sklearn.lda.LDA: Linear discriminant analysis
    """

    def __init__(self, priors=None):
        self.priors = np.asarray(priors) if priors is not None else None

    def fit(self, X, y, store_covariances=False, tol=1.0e-4):
        """
        Fit the QDA model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers)
        store_covariances : boolean
            If True the covariance matrices are computed and stored in the
            `self.covariances_` attribute.
        """
        X, y = check_arrays(X, y)
        self.classes_, y = unique(y, return_inverse=True)
        n_samples, n_features = X.shape
        n_classes = len(self.classes_)
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')
        if self.priors is None:
            self.priors_ = np.bincount(y) / float(n_samples)
        else:
            self.priors_ = self.priors

        cov = None
        if store_covariances:
            cov = []
        means = []
        scalings = []
        rotations = []
        for ind in xrange(n_classes):
            Xg = X[y == ind, :]
            meang = Xg.mean(0)
            means.append(meang)
            Xgc = Xg - meang
            # Xgc = U * S * V.T
            U, S, Vt = np.linalg.svd(Xgc, full_matrices=False)
            rank = np.sum(S > tol)
            if rank < n_features:
                warnings.warn("Variables are collinear")
            S2 = (S ** 2) / (len(Xg) - 1)
            if store_covariances:
                # cov = V * (S^2 / (n-1)) * V.T
                cov.append(np.dot(S2 * Vt.T, Vt))
            scalings.append(S2)
            rotations.append(Vt.T)
        if store_covariances:
            self.covariances_ = cov
        self.means_ = np.asarray(means)
        self.scalings = np.asarray(scalings)
        self.rotations = rotations
        return self

    @property
    def classes(self):
        warnings.warn("QDA.classes is deprecated and will be removed in 0.14. "
                      "Use QDA.classes_ instead.", DeprecationWarning,
                      stacklevel=2)
        return self.classes_

    def _decision_function(self, X):
        X = np.asarray(X)
        norm2 = []
        for i in range(len(self.classes_)):
            R = self.rotations[i]
            S = self.scalings[i]
            Xm = X - self.means_[i]
            X2 = np.dot(Xm, R * (S ** (-0.5)))
            norm2.append(np.sum(X2 ** 2, 1))
        norm2 = np.array(norm2).T   # shape = [len(X), n_classes]
        return (-0.5 * (norm2 + np.sum(np.log(self.scalings), 1))
                + np.log(self.priors_))

    def decision_function(self, X):
        """Apply decision function to an array of samples.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Array of samples (test vectors).

        Returns
        -------
        C : array, shape = [n_samples, n_classes] or [n_samples,]
            Decision function values related to each class, per sample.
            In the two-class case, the shape is [n_samples,], giving the
            log likelihood ratio of the positive class.
        """
        dec_func = self._decision_function(X)
        # handle special case of two classes
        if len(self.classes_) == 2:
            return dec_func[:, 1] - dec_func[:, 0]
        return dec_func

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        The predicted class C for each sample in X is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """
        d = self._decision_function(X)
        y_pred = self.classes_.take(d.argmax(1))
        return y_pred

    def predict_proba(self, X):
        """Return posterior probabilities of classification.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Array of samples/test vectors.

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
            Posterior probabilities of classification per class.
        """
        values = self._decision_function(X)
        # compute the likelihood of the underlying gaussian models
        # up to a multiplicative constant.
        likelihood = np.exp(values - values.min(axis=1)[:, np.newaxis])
        # compute posterior probabilities
        return likelihood / likelihood.sum(axis=1)[:, np.newaxis]

    def predict_log_proba(self, X):
        """Return posterior probabilities of classification.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Array of samples/test vectors.

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
            Posterior log-probabilities of classification per class.
        """
        # XXX : can do better to avoid precision overflows
        probas_ = self.predict_proba(X)
        return np.log(probas_)
