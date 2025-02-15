"""
Quadratic Discriminant Analysis
"""
import warnings
import numpy as np
import scipy.ndimage as ndimage
from .base import BaseEstimator, ClassifierMixin

class QDA(BaseEstimator, ClassifierMixin):
    """
    Quadratic Discriminant Analysis (QDA)

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.
    y : array, shape = [n_samples]
        Target vector relative to X

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
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------
    LDA
    """

    def __init__(self, priors=None):
        """Auto-generated docstring for function __init__."""
        self.priors = np.asarray(priors) if priors is not None else None

    def fit(self, X, y, store_covariances=False, tol=0.0001):
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
            self.covariances_ attribute.
        """
        X = np.asarray(X)
        y = np.asarray(y)
        if X.ndim != 2:
            raise ValueError('X must be a 2D array')
        if X.shape[0] != y.shape[0]:
            raise ValueError('Incompatible shapes: X has %s samples, while y has %s' % (X.shape[0], y.shape[0]))
        if y.dtype.char.lower() not in ('b', 'h', 'i'):
            y = y.astype(np.int32)
        n_samples, n_features = X.shape
        classes = np.unique(y)
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')
        classes_indices = [(y == c).ravel() for c in classes]
        if self.priors is None:
            counts = np.array(ndimage.measurements.sum(np.ones(n_samples, dtype=y.dtype), y, index=classes))
            self.priors_ = counts / float(n_samples)
        else:
            self.priors_ = self.priors
        cov = None
        if store_covariances:
            cov = []
        means = []
        scalings = []
        rotations = []
        for group_indices in classes_indices:
            Xg = X[group_indices, :]
            meang = Xg.mean(0)
            means.append(meang)
            Xgc = Xg - meang
            U, S, Vt = np.linalg.svd(Xgc, full_matrices=False)
            rank = np.sum(S > tol)
            if rank < n_features:
                warnings.warn('Variables are collinear')
            S2 = S ** 2 / (len(Xg) - 1)
            if store_covariances:
                cov.append(np.dot(S2 * Vt.T, Vt))
            scalings.append(S2)
            rotations.append(Vt.T)
        if store_covariances:
            self.covariances_ = cov
        self.means_ = np.asarray(means)
        self.scalings = np.asarray(scalings)
        self.rotations = rotations
        self.classes = classes
        return self

    def decision_function(self, X):
        """Apply decision function to an array of samples.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Array of samples (test vectors).

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
            Decision function values related to each class, per sample.
        """
        X = np.asarray(X)
        norm2 = []
        for i in range(len(self.classes)):
            R = self.rotations[i]
            S = self.scalings[i]
            Xm = X - self.means_[i]
            X2 = np.dot(Xm, R * S ** (-0.5))
            norm2.append(np.sum(X2 ** 2, 1))
        norm2 = np.array(norm2).T
        return -0.5 * (norm2 + np.sum(np.log(self.scalings), 1)) + np.log(self.priors_)

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
        d = self.decision_function(X)
        y_pred = self.classes[d.argmax(1)]
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
        values = self.decision_function(X)
        likelihood = np.exp(values - values.min(axis=1)[:, np.newaxis])
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
        probas_ = self.predict_proba(X)
        return np.log(probas_)