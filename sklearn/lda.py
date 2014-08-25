"""
The :mod:`sklearn.lda` module implements Linear Discriminant Analysis (LDA).
"""
from __future__ import print_function
# Authors: Matthieu Perrot
#          Mathieu Blondel

import warnings

import numpy as np
from scipy import linalg

from .base import BaseEstimator, ClassifierMixin, TransformerMixin
from .utils.extmath import logsumexp
from .utils import check_array, check_X_y

__all__ = ['LDA']


class LDA(BaseEstimator, ClassifierMixin, TransformerMixin):
    """
    Linear Discriminant Analysis (LDA)

    A classifier with a linear decision boundary, generated
    by fitting class conditional densities to the data
    and using Bayes' rule.

    The model fits a Gaussian density to each class, assuming that
    all classes share the same covariance matrix.

    The fitted model can also be used to reduce the dimensionality
    of the input, by projecting it to the most discriminative
    directions.

    Parameters
    ----------
    n_components: int
        Number of components (< n_classes - 1) for dimensionality reduction

    priors : array, optional, shape = [n_classes]
        Priors on classes

    Attributes
    ----------
    coef_ : array-like, shape = [rank, n_classes - 1]
        Coefficients of the features in the linear decision
        function. rank is min(rank_features, n_classes) where
        rank_features is the dimensionality of the spaces spanned
        by the features (i.e. n_features excluding redundant features).

    covariance_ : array-like, shape = [n_features, n_features]
        Covariance matrix (shared by all classes).

    means_ : array-like, shape = [n_classes, n_features]
        Class means.

    priors_ : array-like, shape = [n_classes]
        Class priors (sum to 1).

    scalings_ : array-like, shape = [rank, n_classes - 1]
        Scaling of the features in the space spanned by the class
        centroids.

    xbar_ : float, shape = [n_features]
        Overall mean.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.lda import LDA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = LDA()
    >>> clf.fit(X, y)
    LDA(n_components=None, priors=None)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    sklearn.qda.QDA: Quadratic discriminant analysis

    """

    def __init__(self, n_components=None, priors=None):
        self.n_components = n_components
        self.priors = np.asarray(priors) if priors is not None else None

        if self.priors is not None:
            if (self.priors < 0).any():
                raise ValueError('priors must be non-negative')
            if self.priors.sum() != 1:
                print('warning: the priors do not sum to 1. Renormalizing')
                self.priors = self.priors / self.priors.sum()

    def fit(self, X, y, store_covariance=False, tol=1.0e-4):
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
            and stored in `self.covariance_` attribute.
        """
        X, y = check_X_y(X, y)
        self.classes_, y = np.unique(y, return_inverse=True)
        n_samples, n_features = X.shape
        n_classes = len(self.classes_)
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')
        if self.priors is None:
            self.priors_ = np.bincount(y) / float(n_samples)
        else:
            self.priors_ = self.priors

        # Group means n_classes*n_features matrix
        means = []
        Xc = []
        cov = None
        if store_covariance:
            cov = np.zeros((n_features, n_features))
        for ind in range(n_classes):
            Xg = X[y == ind, :]
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

        self.means_ = np.asarray(means)
        Xc = np.concatenate(Xc, axis=0)

        # ----------------------------
        # 1) within (univariate) scaling by with classes std-dev
        std = Xc.std(axis=0)
        # avoid division by zero in normalization
        std[std == 0] = 1.
        fac = 1. / (n_samples - n_classes)
        # ----------------------------
        # 2) Within variance scaling
        X = np.sqrt(fac) * (Xc / std)
        # SVD of centered (within)scaled data
        U, S, V = linalg.svd(X, full_matrices=False)

        rank = np.sum(S > tol)
        if rank < n_features:
            warnings.warn("Variables are collinear")
        # Scaling of within covariance is: V' 1/S
        scalings = (V[:rank] / std).T / S[:rank]

        ## ----------------------------
        ## 3) Between variance scaling
        # Overall mean
        xbar = np.dot(self.priors_, self.means_)
        # Scale weighted centers
        X = np.dot(((np.sqrt((n_samples * self.priors_) * fac)) *
                    (means - xbar).T).T, scalings)
        # Centers are living in a space with n_classes-1 dim (maximum)
        # Use svd to find projection in the space spanned by the
        # (n_classes) centers
        _, S, V = linalg.svd(X, full_matrices=0)

        rank = np.sum(S > tol * S[0])
        # compose the scalings
        self.scalings_ = np.dot(scalings, V.T[:, :rank])
        self.xbar_ = xbar
        # weight vectors / centroids
        self.coef_ = np.dot(self.means_ - self.xbar_, self.scalings_)
        self.intercept_ = (-0.5 * np.sum(self.coef_ ** 2, axis=1) +
                           np.log(self.priors_))
        return self

    def _decision_function(self, X):
        X = check_array(X)
        # center and scale data
        X = np.dot(X - self.xbar_, self.scalings_)
        return np.dot(X, self.coef_.T) + self.intercept_

    def decision_function(self, X):
        """
        This function returns the decision function values related to each
        class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, n_classes] or [n_samples,]
            Decision function values related to each class, per sample.
            In the two-class case, the shape is [n_samples,], giving the
            log likelihood ratio of the positive class.
        """
        dec_func = self._decision_function(X)
        if len(self.classes_) == 2:
            return dec_func[:, 1] - dec_func[:, 0]
        return dec_func

    def transform(self, X):
        """
        Project the data so as to maximize class separation (large separation
        between projected class means and small variance within each class).

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        X_new : array, shape = [n_samples, n_components]
        """
        X = check_array(X)
        # center and scale data
        X_new = np.dot(X - self.xbar_, self.scalings_)
        n_components = X.shape[1] if self.n_components is None \
            else self.n_components
        return X_new[:, :n_components]

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
        d = self._decision_function(X)
        y_pred = self.classes_.take(d.argmax(1))
        return y_pred

    def predict_proba(self, X):
        """
        This function returns posterior probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
        """
        values = self._decision_function(X)
        # compute the likelihood of the underlying gaussian models
        # up to a multiplicative constant.
        likelihood = np.exp(values - values.max(axis=1)[:, np.newaxis])
        # compute posterior probabilities
        return likelihood / likelihood.sum(axis=1)[:, np.newaxis]

    def predict_log_proba(self, X):
        """
        This function returns posterior log-probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, n_classes]
        """
        values = self._decision_function(X)
        loglikelihood = (values - values.max(axis=1)[:, np.newaxis])
        normalization = logsumexp(loglikelihood, axis=1)
        return loglikelihood - normalization[:, np.newaxis]
