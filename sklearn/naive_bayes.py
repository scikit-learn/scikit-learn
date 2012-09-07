# -*- coding: utf-8 -*-

"""
The :mod:`sklearn.naive_bayes` module implements Naive Bayes algorithms. These
are supervised learning methods based on applying Bayes' theorem with strong
(naive) feature independence assumptions.
"""

# Author: Vincent Michel <vincent.michel@inria.fr>
#         Minor fixes by Fabian Pedregosa
#         Amit Aides <amitibo@tx.technion.ac.il>
#         Yehuda Finkelstein <yehudaf@tx.technion.ac.il>
#         Lars Buitinck <L.J.Buitinck@uva.nl>
#         (parts based on earlier work by Mathieu Blondel)
#
# License: BSD Style.

from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import issparse

from .base import BaseEstimator, ClassifierMixin
from .preprocessing import binarize, LabelBinarizer
from .utils import array2d, atleast2d_or_csr
from .utils.extmath import safe_sparse_dot, logsumexp
from .utils import check_arrays

__all__ = ['BernoulliNB', 'GaussianNB', 'MultinomialNB']


class BaseNB(BaseEstimator, ClassifierMixin):
    """Abstract base class for naive Bayes estimators"""

    __metaclass__ = ABCMeta

    @abstractmethod
    def _joint_log_likelihood(self, X):
        """Compute the unnormalized posterior log probability of X

        I.e. ``log P(c) + log P(x|c)`` for all rows x of X, as an array-like of
        shape [n_classes, n_samples].

        Input is passed to _joint_log_likelihood as-is by predict,
        predict_proba and predict_log_proba.
        """

    def predict(self, X):
        """
        Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Predicted target values for X
        """
        jll = self._joint_log_likelihood(X)
        return self.classes_[np.argmax(jll, axis=1)]

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
            in the model, where classes are ordered arithmetically.
        """
        jll = self._joint_log_likelihood(X)
        # normalize by P(x) = P(f_1, ..., f_n)
        log_prob_x = logsumexp(jll, axis=1)
        return jll - np.atleast_2d(log_prob_x).T

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
            the model, where classes are ordered arithmetically.
        """
        return np.exp(self.predict_log_proba(X))


class GaussianNB(BaseNB):
    """
    Gaussian Naive Bayes (GaussianNB)

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    y : array, shape = [n_samples]
        Target vector relative to X

    Attributes
    ----------
    `class_prior_` : array, shape = [n_classes]
        probability of each class.

    `theta_` : array, shape = [n_classes, n_features]
        mean of each feature per class

    `sigma_` : array, shape = [n_classes, n_features]
        variance of each feature per class

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> Y = np.array([1, 1, 1, 2, 2, 2])
    >>> from sklearn.naive_bayes import GaussianNB
    >>> clf = GaussianNB()
    >>> clf.fit(X, Y)
    GaussianNB()
    >>> print(clf.predict([[-0.8, -1]]))
    [1]
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

        X, y = check_arrays(X, y, sparse_format='dense')

        n_samples, n_features = X.shape

        if n_samples != y.shape[0]:
            raise ValueError("X and y have incompatible shapes")

        self.classes_ = unique_y = np.unique(y)
        n_classes = unique_y.shape[0]

        self.theta_ = np.zeros((n_classes, n_features))
        self.sigma_ = np.zeros((n_classes, n_features))
        self.class_prior_ = np.zeros(n_classes)
        epsilon = 1e-9
        for i, y_i in enumerate(unique_y):
            self.theta_[i, :] = np.mean(X[y == y_i, :], axis=0)
            self.sigma_[i, :] = np.var(X[y == y_i, :], axis=0) + epsilon
            self.class_prior_[i] = np.float(np.sum(y == y_i)) / n_samples
        return self

    def _joint_log_likelihood(self, X):
        X = array2d(X)
        joint_log_likelihood = []
        for i in xrange(np.size(self.classes_)):
            jointi = np.log(self.class_prior_[i])
            n_ij = - 0.5 * np.sum(np.log(np.pi * self.sigma_[i, :]))
            n_ij -= 0.5 * np.sum(((X - self.theta_[i, :]) ** 2) / \
                                 (self.sigma_[i, :]), 1)
            joint_log_likelihood.append(jointi + n_ij)

        joint_log_likelihood = np.array(joint_log_likelihood).T
        return joint_log_likelihood


class BaseDiscreteNB(BaseNB):
    """Abstract base class for naive Bayes on discrete/categorical data

    Any estimator based on this class should provide:

    __init__
    _joint_log_likelihood(X) as per BaseNB
    """

    def fit(self, X, y, sample_weight=None, class_prior=None):
        """Fit Naive Bayes classifier according to X, y

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        class_prior : array, shape [n_classes]
            Custom prior probability per class.
            Overrides the fit_prior parameter.

        Returns
        -------
        self : object
            Returns self.
        """
        X = atleast2d_or_csr(X)

        labelbin = LabelBinarizer()
        Y = labelbin.fit_transform(y)
        self.classes_ = labelbin.classes_
        n_classes = len(self.classes_)
        if Y.shape[1] == 1:
            Y = np.concatenate((1 - Y, Y), axis=1)

        if X.shape[0] != Y.shape[0]:
            msg = "X and y have incompatible shapes."
            if issparse(X):
                msg += "\nNote: Sparse matrices cannot be indexed w/ boolean \
                masks (use `indices=True` in CV)."
            raise ValueError(msg)

        if sample_weight is not None:
            Y *= array2d(sample_weight).T

        if class_prior:
            if len(class_prior) != n_classes:
                raise ValueError(
                        "Number of priors must match number of classes")
            self.class_log_prior_ = np.log(class_prior)
        elif self.fit_prior:
            # empirical prior, with sample_weight taken into account
            y_freq = Y.sum(axis=0)
            self.class_log_prior_ = np.log(y_freq) - np.log(y_freq.sum())
        else:
            self.class_log_prior_ = np.zeros(n_classes) - np.log(n_classes)

        N_c, N_c_i = self._count(X, Y)

        self.feature_log_prob_ = (np.log(N_c_i + self.alpha)
                                - np.log(N_c.reshape(-1, 1)
                                       + self.alpha * X.shape[1]))

        return self

    @staticmethod
    def _count(X, Y):
        """Count feature occurrences.

        Returns (N_c, N_c_i), where
            N_c is the count of all features in all samples of class c;
            N_c_i is the count of feature i in all samples of class c.
        """
        if np.any((X.data if issparse(X) else X) < 0):
            raise ValueError("Input X must be non-negative.")
        N_c_i = safe_sparse_dot(Y.T, X)
        N_c = np.sum(N_c_i, axis=1)

        return N_c, N_c_i

    # XXX The following is a stopgap measure; we need to set the dimensions
    # of class_log_prior_ and feature_log_prob_ correctly.
    def _get_coef(self):
        return self.feature_log_prob_[1] if len(self.classes_) == 2 \
                                         else self.feature_log_prob_

    def _get_intercept(self):
        return self.class_log_prior_[1] if len(self.classes_) == 2 \
                                        else self.class_log_prior_

    coef_ = property(_get_coef)
    intercept_ = property(_get_intercept)


class MultinomialNB(BaseDiscreteNB):
    """
    Naive Bayes classifier for multinomial models

    The multinomial Naive Bayes classifier is suitable for classification with
    discrete features (e.g., word counts for text classification). The
    multinomial distribution normally requires integer feature counts. However,
    in practice, fractional counts such as tf-idf may also work.

    Parameters
    ----------
    alpha: float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).
    fit_prior: boolean
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    Attributes
    ----------
    `intercept_`, `class_log_prior_` : array, shape = [n_classes]
        Smoothed empirical log probability for each class.

    `feature_log_prob_`, `coef_` : array, shape = [n_classes, n_features]
        Empirical log probability of features
        given a class, P(x_i|y).

        (`intercept_` and `coef_` are properties
        referring to `class_log_prior_` and
        `feature_log_prob_`, respectively.)

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randint(5, size=(6, 100))
    >>> Y = np.array([1, 2, 3, 4, 5, 6])
    >>> from sklearn.naive_bayes import MultinomialNB
    >>> clf = MultinomialNB()
    >>> clf.fit(X, Y)
    MultinomialNB(alpha=1.0, fit_prior=True)
    >>> print(clf.predict(X[2]))
    [3]

    Notes
    -----
    For the rationale behind the names `coef_` and `intercept_`, i.e.
    naive Bayes as a linear classifier, see J. Rennie et al. (2003),
    Tackling the poor assumptions of naive Bayes text classifiers, ICML.
    """

    def __init__(self, alpha=1.0, fit_prior=True):
        self.alpha = alpha
        self.fit_prior = fit_prior

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""
        X = atleast2d_or_csr(X)
        return (safe_sparse_dot(X, self.feature_log_prob_.T)
               + self.class_log_prior_)


class BernoulliNB(BaseDiscreteNB):
    """Naive Bayes classifier for multivariate Bernoulli models.

    Like MultinomialNB, this classifier is suitable for discrete data. The
    difference is that while MultinomialNB works with occurrence counts,
    BernoulliNB is designed for binary/boolean features.

    Parameters
    ----------
    alpha: float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).
    binarize: float or None, optional
        Threshold for binarizing (mapping to booleans) of sample features.
        If None, input is presumed to already consist of binary vectors.
    fit_prior: boolean
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    Attributes
    ----------
    `class_log_prior_` : array, shape = [n_classes]
        Log probability of each class (smoothed).

    `feature_log_prob_` : array, shape = [n_classes, n_features]
        Empirical log probability of features given a class, P(x_i|y).

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randint(2, size=(6, 100))
    >>> Y = np.array([1, 2, 3, 4, 4, 5])
    >>> from sklearn.naive_bayes import BernoulliNB
    >>> clf = BernoulliNB()
    >>> clf.fit(X, Y)
    BernoulliNB(alpha=1.0, binarize=0.0, fit_prior=True)
    >>> print(clf.predict(X[2]))
    [3]

    References
    ----------

    C.D. Manning, P. Raghavan and H. Schütze (2008). Introduction to
    Information Retrieval. Cambridge University Press, pp. 234–265.

    A. McCallum and K. Nigam (1998). A comparison of event models for naive
    Bayes text classification. Proc. AAAI/ICML-98 Workshop on Learning for
    Text Categorization, pp. 41–48.

    V. Metsis, I. Androutsopoulos and G. Paliouras (2006). Spam filtering with
    naive Bayes -- Which naive Bayes? 3rd Conf. on Email and Anti-Spam (CEAS).
    """

    def __init__(self, alpha=1.0, binarize=.0, fit_prior=True):
        self.alpha = alpha
        self.binarize = binarize
        self.fit_prior = fit_prior

    def _count(self, X, Y):
        if self.binarize is not None:
            X = binarize(X, threshold=self.binarize)
        return super(BernoulliNB, self)._count(X, Y)

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""

        X = atleast2d_or_csr(X)

        if self.binarize is not None:
            X = binarize(X, threshold=self.binarize)

        n_classes, n_features = self.feature_log_prob_.shape
        n_samples, n_features_X = X.shape

        if n_features_X != n_features:
            raise ValueError("Expected input with %d features, got %d instead"
                             % (n_features, n_features_X))

        neg_prob = np.log(1 - np.exp(self.feature_log_prob_))
        # Compute  neg_prob · (1 - X).T  as  ∑neg_prob - X · neg_prob
        X_neg_prob = (neg_prob.sum(axis=1)
                    - safe_sparse_dot(X, neg_prob.T))
        jll = safe_sparse_dot(X, self.feature_log_prob_.T) + X_neg_prob

        return jll + self.class_log_prior_
