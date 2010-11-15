"""
LDA: Linear Discriminant Analysis
"""
# Author: Matthieu Perrot

import warnings

import numpy as np
from scipy import linalg, ndimage

from .base import BaseEstimator, ClassifierMixin


class LDA(BaseEstimator, ClassifierMixin):
    """
    Linear Discriminant Analysis (LDA)

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
    `xbar_` : float, shape = [n_features]
        Over all mean
    `priors_` : array-like, shape = [n_classes]
        Class priors (sum to 1)
    `covariance_` : array-like, shape = [n_features, n_features]
        Covariance matrix (shared by all classes)

    Methods
    -------
    fit(X, y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.lda import LDA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = LDA()
    >>> clf.fit(X, y)
    LDA(priors=None)
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------
    QDA

    """

    def __init__(self, priors=None):
        self.priors = np.asarray(priors) if priors is not None else None

    def fit(self, X, y, store_covariance=False, tol=1.0e-4, **params):
        """
        Fit the LDA model according to the given training data and parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers)
        store_covariance : boolean
            If True the covariance matrix (shared by all classes) is computed
            and stored in self.covariance_ attribute.
        """
        self._set_params(**params)
        X = np.asanyarray(X)
        y = np.asanyarray(y)
        if X.ndim != 2:
            raise ValueError('X must be a 2D array')
        n_samples = X.shape[0]
        n_features = X.shape[1]
        # We need int32 to be able to use ndimage.measurements
        classes = np.unique(y).astype(np.int32)
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')
        classes_indices = [(y == c).ravel() for c in classes]
        if self.priors is None:
            counts = np.array(ndimage.measurements.sum(np.ones(len(y)),
                                                    y, index=classes))
            self.priors_ = counts / float(n_samples)
        else:
            self.priors_ = self.priors

        # Group means n_classes*n_features matrix
        means = []
        Xc = []
        cov = None
        if store_covariance:
            cov = np.zeros((n_features, n_features))
        for group_indices in classes_indices:
            Xg = X[group_indices, :]
            meang = Xg.mean(0)
            means.append(meang)
            # centered group data
            Xgc = Xg - meang
            Xc.append(Xgc)
            if store_covariance:
                cov += np.dot(Xgc.T, Xgc)
        if store_covariance:
            cov /= (n_samples - n_classes)
            self.covariance_ = cov

        means = np.asarray(means)
        Xc = np.concatenate(Xc, 0)

        # ----------------------------
        # 1) within (univariate) scaling by with classes std-dev
        scaling = 1. / Xc.std(0)
        fac = float(1) / (n_samples - n_classes)
        # ----------------------------
        # 2) Within variance scaling
        X = np.sqrt(fac) * (Xc * scaling)
        # SVD of centered (within)scaled data
        U, S, V = linalg.svd(X, full_matrices=0)

        rank = np.sum(S > tol)
        if rank < n_features:
            warnings.warn("Variables are collinear")
        # Scaling of within covariance is: V' 1/S
        scaling = (scaling * V.T[:, :rank].T).T / S[:rank]
        ## ----------------------------
        ## 3) Between variance scaling
        # Overall mean
        xbar = np.dot(self.priors_, means)
        # Scale weighted centers
        X = np.dot(((np.sqrt((n_samples * self.priors_)*fac)) *
                          (means - xbar).T).T, scaling)
        # Centers are living in a space with n_classes-1 dim (maximum)
        # Use svd to find projection in the space spamed by the
        # (n_classes) centers
        _, S, V = linalg.svd(X, full_matrices=0)

        rank = np.sum(S > tol*S[0])
        # compose the scalings
        scaling = np.dot(scaling, V.T[:, :rank])
        self.scaling = scaling
        self.means_ = means
        self.xbar_ = xbar
        self.classes = classes
        return self

    def decision_function(self, X):
        """
        This function return the decision function values related to each
        class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
        """
        X = np.asanyarray(X)
        scaling = self.scaling
        # Remove overall mean (center) and scale
        # a) data
        X = np.dot(X - self.xbar_, scaling)
        # b) centers
        dm = np.dot(self.means_ - self.xbar_, scaling)
        # for each class k, compute the linear discrinant function
        # (p. 87 Hastie) of sphered (scaled data)
        return -0.5 * np.sum(dm ** 2, 1) + \
                np.log(self.priors_) + np.dot(X, dm.T)

    def predict(self, X):
        """
        This function does classification on an array of test vectors X.

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
        """
        This function return posterior probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
        """
        values = self.decision_function(X)
        # compute the likelihood of the underlying gaussian models
        # up to a multiplicative constant.
        likelihood = np.exp(values - values.min(axis=1)[:, np.newaxis])
        # compute posterior probabilities
        return likelihood / likelihood.sum(axis=1)[:, np.newaxis]

