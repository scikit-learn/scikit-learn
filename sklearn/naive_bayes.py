# -*- coding: utf-8 -*-

"""
The :mod:`sklearn.naive_bayes` module implements Naive Bayes algorithms. These
are supervised learning methods based on applying Bayes' theorem with strong
(naive) feature independence assumptions.

It also implements routines for generative classification with fewer
assumptions, i.e. non-naive Bayesian classification.
"""

# Author: Vincent Michel <vincent.michel@inria.fr>
#         Minor fixes by Fabian Pedregosa
#         Amit Aides <amitibo@tx.technion.ac.il>
#         Yehuda Finkelstein <yehudaf@tx.technion.ac.il>
#         Lars Buitinck <L.J.Buitinck@uva.nl>
#         (parts based on earlier work by Mathieu Blondel)
#         Generative Classification by
#         Jake Vanderplas <jakevdp@cs.washington.edu>
#
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.sparse import issparse
import warnings

from .base import BaseEstimator, ClassifierMixin, clone
from .preprocessing import binarize
from .preprocessing import LabelBinarizer
from .preprocessing import label_binarize
from .utils import (array2d, atleast2d_or_csr, column_or_1d, check_arrays,
                    check_random_state, check_arrays)
from .utils.extmath import safe_sparse_dot, logsumexp
from .utils.multiclass import _check_partial_fit_first_call
from .externals import six
from .neighbors import KernelDensity
from .mixture import GMM

__all__ = ['BernoulliNB', 'GaussianNB', 'MultinomialNB', 'GenerativeBayes']


class BaseNB(six.with_metaclass(ABCMeta, BaseEstimator, ClassifierMixin)):
    """Abstract base class for naive Bayes estimators"""

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
        y = column_or_1d(y, warn=True)

        n_samples, n_features = X.shape

        self.classes_ = unique_y = np.unique(y)
        n_classes = unique_y.shape[0]

        self.theta_ = np.zeros((n_classes, n_features))
        self.sigma_ = np.zeros((n_classes, n_features))
        self.class_prior_ = np.zeros(n_classes)
        epsilon = 1e-9
        for i, y_i in enumerate(unique_y):
            Xi = X[y == y_i, :]
            self.theta_[i, :] = np.mean(Xi, axis=0)
            self.sigma_[i, :] = np.var(Xi, axis=0) + epsilon
            self.class_prior_[i] = np.float(Xi.shape[0]) / n_samples
        return self

    def _joint_log_likelihood(self, X):
        X = array2d(X)
        joint_log_likelihood = []
        for i in range(np.size(self.classes_)):
            jointi = np.log(self.class_prior_[i])
            n_ij = - 0.5 * np.sum(np.log(np.pi * self.sigma_[i, :]))
            n_ij -= 0.5 * np.sum(((X - self.theta_[i, :]) ** 2) /
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

    def _update_class_log_prior(self, class_prior=None):
        n_classes = len(self.classes_)
        if class_prior is not None:
            if len(class_prior) != n_classes:
                raise ValueError("Number of priors must match number of"
                                 " classes.")
            self.class_log_prior_ = np.log(class_prior)
        elif self.fit_prior:
            # empirical prior, with sample_weight taken into account
            self.class_log_prior_ = (np.log(self.class_count_)
                                     - np.log(self.class_count_.sum()))
        else:
            self.class_log_prior_ = np.zeros(n_classes) - np.log(n_classes)

    def partial_fit(self, X, y, classes=None, sample_weight=None):
        """Incremental fit on a batch of samples.

        This method is expected to be called several times consecutively
        on different chunks of a dataset so as to implement out-of-core
        or online learning.

        This is especially useful when the whole dataset is too big to fit in
        memory at once.

        This method has some performance overhead hence it is better to call
        partial_fit on chunks of data that are as large as possible
        (as long as fitting in the memory budget) to hide the overhead.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        classes : array-like, shape = [n_classes]
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like, shape = [n_samples], optional
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns self.
        """
        X = atleast2d_or_csr(X, dtype=np.float64)
        _, n_features = X.shape

        if _check_partial_fit_first_call(self, classes):
            # This is the first call to partial_fit:
            # initialize various cumulative counters
            n_effective_classes = len(classes) if len(classes) > 1 else 2
            self.class_count_ = np.zeros(n_effective_classes, dtype=np.float64)
            self.feature_count_ = np.zeros((n_effective_classes, n_features),
                                           dtype=np.float64)

        Y = label_binarize(y, classes=self.classes_)
        if Y.shape[1] == 1:
            Y = np.concatenate((1 - Y, Y), axis=1)

        n_samples, n_classes = Y.shape

        if X.shape[0] != Y.shape[0]:
            msg = "X.shape[0]=%d and y.shape[0]=%d are incompatible."
            raise ValueError(msg % (X.shape[0], y.shape[0]))

        # convert to float to support sample weight consistently
        Y = Y.astype(np.float64)
        if sample_weight is not None:
            Y *= array2d(sample_weight).T

        # Count raw events from data before updating the class log prior
        # and feature log probas
        self._count(X, Y)

        # XXX: OPTIM: we could introduce a public finalization method to
        # be called by the user explicitly just once after several consecutive
        # calls to partial_fit and prior any call to predict[_[log_]proba]
        # to avoid computing the smooth log probas at each call to partial fit
        self._update_feature_log_prob()
        self._update_class_log_prior()
        return self

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

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_arrays(X, y, sparse_format='csr')
        X = X.astype(np.float)
        y = column_or_1d(y, warn=True)
        _, n_features = X.shape

        labelbin = LabelBinarizer()
        Y = labelbin.fit_transform(y)
        self.classes_ = labelbin.classes_
        if Y.shape[1] == 1:
            Y = np.concatenate((1 - Y, Y), axis=1)

        # convert to float to support sample weight consistently
        Y = Y.astype(np.float64)
        if sample_weight is not None:
            Y *= array2d(sample_weight).T

        if class_prior is not None:
            warnings.warn('class_prior has been made an ``__init__`` parameter'
                          ' and will be removed from fit in version 0.15.',
                          DeprecationWarning)
        else:
            class_prior = self.class_prior

        # Count raw events from data before updating the class log prior
        # and feature log probas
        n_effective_classes = Y.shape[1]
        self.class_count_ = np.zeros(n_effective_classes, dtype=np.float64)
        self.feature_count_ = np.zeros((n_effective_classes, n_features),
                                       dtype=np.float64)
        self._count(X, Y)
        self._update_feature_log_prob()
        self._update_class_log_prior(class_prior=class_prior)
        return self

    # XXX The following is a stopgap measure; we need to set the dimensions
    # of class_log_prior_ and feature_log_prob_ correctly.
    def _get_coef(self):
        return (self.feature_log_prob_[1:]
                if len(self.classes_) == 2 else self.feature_log_prob_)

    def _get_intercept(self):
        return (self.class_log_prior_[1:]
                if len(self.classes_) == 2 else self.class_log_prior_)

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
    alpha : float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).

    fit_prior : boolean
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    class_prior : array-like, size (n_classes,)
        Prior probabilities of the classes. If specified the priors are not
        adjusted according to the data.

    Attributes
    ----------
    `class_log_prior_` : array, shape (n_classes, )
        Smoothed empirical log probability for each class.

    `intercept_` : property
        Mirrors ``class_log_prior_`` for interpreting MultinomialNB
        as a linear model.

    `feature_log_prob_`: array, shape (n_classes, n_features)
        Empirical log probability of features
        given a class, ``P(x_i|y)``.

    `coef_` : property
        Mirrors ``feature_log_prob_`` for interpreting MultinomialNB
        as a linear model.

    `class_count_` : array, shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    `feature_count_` : array, shape (n_classes, n_features)
        Number of samples encountered for each (class, feature)
        during fitting. This value is weighted by the sample weight when
        provided.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randint(5, size=(6, 100))
    >>> y = np.array([1, 2, 3, 4, 5, 6])
    >>> from sklearn.naive_bayes import MultinomialNB
    >>> clf = MultinomialNB()
    >>> clf.fit(X, y)
    MultinomialNB(alpha=1.0, class_prior=None, fit_prior=True)
    >>> print(clf.predict(X[2]))
    [3]

    Notes
    -----
    For the rationale behind the names `coef_` and `intercept_`, i.e.
    naive Bayes as a linear classifier, see J. Rennie et al. (2003),
    Tackling the poor assumptions of naive Bayes text classifiers, ICML.

    References
    ----------
    C.D. Manning, P. Raghavan and H. Schuetze (2008). Introduction to
    Information Retrieval. Cambridge University Press, pp. 234-265.
    http://nlp.stanford.edu/IR-book/html/htmledition/
        naive-bayes-text-classification-1.html
    """

    def __init__(self, alpha=1.0, fit_prior=True, class_prior=None):
        self.alpha = alpha
        self.fit_prior = fit_prior
        self.class_prior = class_prior

    def _count(self, X, Y):
        """Count and smooth feature occurrences."""
        if np.any((X.data if issparse(X) else X) < 0):
            raise ValueError("Input X must be non-negative")
        self.feature_count_ += safe_sparse_dot(Y.T, X)
        self.class_count_ += Y.sum(axis=0)

    def _update_feature_log_prob(self):
        """Apply smoothing to raw counts and recompute log probabilities"""
        smoothed_fc = self.feature_count_ + self.alpha
        smoothed_cc = smoothed_fc.sum(axis=1)

        self.feature_log_prob_ = (np.log(smoothed_fc)
                                  - np.log(smoothed_cc.reshape(-1, 1)))

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
    alpha : float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).

    binarize : float or None, optional
        Threshold for binarizing (mapping to booleans) of sample features.
        If None, input is presumed to already consist of binary vectors.

    fit_prior : boolean
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    class_prior : array-like, size=[n_classes,]
        Prior probabilities of the classes. If specified the priors are not
        adjusted according to the data.

    Attributes
    ----------
    `class_log_prior_` : array, shape = [n_classes]
        Log probability of each class (smoothed).

    `feature_log_prob_` : array, shape = [n_classes, n_features]
        Empirical log probability of features given a class, P(x_i|y).

    `class_count_` : array, shape = [n_classes]
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    `feature_count_` : array, shape = [n_classes, n_features]
        Number of samples encountered for each (class, feature)
        during fitting. This value is weighted by the sample weight when
        provided.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randint(2, size=(6, 100))
    >>> Y = np.array([1, 2, 3, 4, 4, 5])
    >>> from sklearn.naive_bayes import BernoulliNB
    >>> clf = BernoulliNB()
    >>> clf.fit(X, Y)
    BernoulliNB(alpha=1.0, binarize=0.0, class_prior=None, fit_prior=True)
    >>> print(clf.predict(X[2]))
    [3]

    References
    ----------

    C.D. Manning, P. Raghavan and H. Schuetze (2008). Introduction to
    Information Retrieval. Cambridge University Press, pp. 234-265.
    http://nlp.stanford.edu/IR-book/html/htmledition/the-bernoulli-model-1.html

    A. McCallum and K. Nigam (1998). A comparison of event models for naive
    Bayes text classification. Proc. AAAI/ICML-98 Workshop on Learning for
    Text Categorization, pp. 41-48.

    V. Metsis, I. Androutsopoulos and G. Paliouras (2006). Spam filtering with
    naive Bayes -- Which naive Bayes? 3rd Conf. on Email and Anti-Spam (CEAS).
    """

    def __init__(self, alpha=1.0, binarize=.0, fit_prior=True,
                 class_prior=None):
        self.alpha = alpha
        self.binarize = binarize
        self.fit_prior = fit_prior
        self.class_prior = class_prior

    def _count(self, X, Y):
        """Count and smooth feature occurrences."""
        if self.binarize is not None:
            X = binarize(X, threshold=self.binarize)
        self.feature_count_ += safe_sparse_dot(Y.T, X)
        self.class_count_ += Y.sum(axis=0)

    def _update_feature_log_prob(self):
        """Apply smoothing to raw counts and recompute log probabilities"""
        n_classes = len(self.classes_)
        smoothed_fc = self.feature_count_ + self.alpha
        smoothed_cc = self.class_count_ + self.alpha * n_classes

        self.feature_log_prob_ = (np.log(smoothed_fc)
                                  - np.log(smoothed_cc.reshape(-1, 1)))

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
        jll = safe_sparse_dot(X, (self.feature_log_prob_ - neg_prob).T)
        jll += self.class_log_prior_ + neg_prob.sum(axis=1)

        return jll


class _NormalApproximation(BaseEstimator):
    """Fit an axis-aligned normal approximation to the data.

    Parameters
    ----------
    epsilon : float
        The minimum variance along any dimension.  Default = 1E-9.
    """
    def __init__(self, epsilon=1E-9):
        self.epsilon = epsilon

    def fit(self, X):
        """Fit the Normal Approximation to data

        Parameters
        ----------
        X: array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        """
        X = array2d(X)
        self.mean = X.mean(0)
        self.var = X.var(0) + self.epsilon
        return self

    def score_samples(self, X):
        """Evaluate the model on the data

        Parameters
        ----------
        X : array_like
            An array of points to query.  Last dimension should match dimension
            of training data (n_features)

        Returns
        -------
        sample scores : ndarray
            The array of log density.  This has shape X.shape[:-1]
        """
        X = array2d(X)
        if X.shape[-1] != self.mean.shape[0]:
            raise ValueError("dimension of X must match that of training data")
        log_norm = -0.5 * (X.shape[-1] * np.log(2 * np.pi)
                           + np.log(self.var).sum())
        return log_norm - 0.5 * ((X - self.mean) ** 2 / self.var).sum(1)

    def score(self, X):
        """Compute the log probability under the model.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : float
            log-likelihood of the data X under the model.
        """
        return np.sum(self.score_samples(X))

    def sample(self, n_samples=1, random_state=None):
        """Generate random samples from the model.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.

        random_state: RandomState or an int seed (0 by default)
            A random number generator instance

        Returns
        -------
        X : array_like, shape (n_samples, n_features)
            List of samples
        """
        rng = check_random_state(random_state)
        return rng.normal(self.mean, np.sqrt(self.var),
                          size=(n_samples, len(self.mean)))


DENSITY_MODELS = {'normal_approximation': _NormalApproximation,
                  'gmm': GMM,
                  'kde': KernelDensity}


class GenerativeBayes(BaseNB):
    """
    Generative Bayes Classifier

    This is a meta-estimator which performs generative Bayesian classification
    using flexible underlying density models.

    Parameters
    ----------
    density_estimator : str or instance
        The density estimator to use for each class.  Options are
            - 'normal_approximation' : Axis-aligned Normal Approximation
                                       (i.e. Gaussian Naive Bayes)
            - 'gmm' : Gaussian Mixture Model
            - 'kde' : Kernel Density Estimate
        The default is 'normal_approximation'.
        Alternatively, a scikit-learn estimator instance can be specified.
        The estimator should contain a ``score_samples`` method with semantics
        similar to those in :class:`sklearn.neighbors.KDE`.
    model_kwds : dict or None
        Additional keyword arguments to be passed to the constructor
        specified by density_estimator.  Ignored if density_estimator is
        a class instance.  Default=None.

    Attributes
    ----------
    `classes_` : array, shape = [n_classes]
        the sorted list of classes

    `class_prior_` : array, shape = [n_classes]
        probability of each class.

    `estimators_` : list, length = [n_classes]
        the density estimator associated with each class
    """
    def __init__(self, density_estimator='normal_approximation',
                 model_kwds=None):
        self.density_estimator = density_estimator
        self.model_kwds = model_kwds

        # run this here to check for any exceptions; we avoid assigning
        # the result here so that the estimator can be cloned.
        self._choose_estimator(density_estimator, self.model_kwds)

    def _choose_estimator(self, density_estimator, kwargs=None):
        """Choose the estimator based on the input"""
        dclass = DENSITY_MODELS.get(density_estimator)

        if dclass is not None:
            if kwargs is None:
                kwargs = {}
            density_estimator = dclass(**kwargs)

        if not hasattr(dclass, 'score_samples'):
            raise ValueError('invalid density_estimator')

        return density_estimator

    def fit(self, X, y):
        """Fit the model using X as training data and y as target values

        Parameters
        ----------
        X : array-like
            Training data. shape = [n_samples, n_features]

        y : array-like
            Target values, array of float values, shape = [n_samples]
        """
        X, y = check_arrays(X, y, sparse_format='dense')
        y = column_or_1d(y, warn=True)

        estimator = self._choose_estimator(self.density_estimator,
                                           self.model_kwds)

        self.classes_ = np.sort(np.unique(y))
        n_classes = len(self.classes_)
        n_samples, self.n_features_ = X.shape

        class_membership_masks = [(y == c) for c in self.classes_]

        self.class_prior_ = np.array([np.float(mask.sum()) / n_samples
                                      for mask in class_membership_masks])
        self.estimators_ = [clone(estimator).fit(X[mask])
                            for mask in class_membership_masks]
        return self

    def _joint_log_likelihood(self, X):
        """Compute the per-class log likelihood of each sample

        Parameters
        ----------
        X : array_like
            Array of samples on which to compute likelihoods.  Shape is
            (n_samples, n_features)

        Returns
        -------
        logL : array_like
            The log likelihood under each class.
            Shape is (n_samples, n_classes).  logL[i, j] gives the log
            likelihood of X[i] within the model representing the class
            self.classes_[j].
        """
        X = array2d(X)

        # GMM API, in particular score() and score_samples(), is
        # not consistent with the rest of the package.  This needs
        # to be addressed eventually...
        if isinstance(self.estimators_[0], GMM):
            return np.array([np.log(prior) + dens.score(X)
                             for (prior, dens)
                             in zip(self.class_prior_,
                                    self.estimators_)]).T
        else:
            return np.array([np.log(prior) + dens.score_samples(X)
                             for (prior, dens)
                             in zip(self.class_prior_,
                                    self.estimators_)]).T

    def sample(self, n_samples=1, random_state=None):
        """Generate random samples from the model.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.

        Returns
        -------
        X : array_like, shape (n_samples, n_features)
            List of samples
        y : array_like, shape (n_samples,)
            List of class labels for the generated samples
        """
        random_state = check_random_state(random_state)
        X = np.empty((n_samples, self.n_features_))
        rand = random_state.rand(n_samples)

        # split samples by class
        prior_cdf = np.cumsum(self.class_prior_)
        labels = prior_cdf.searchsorted(rand)

        # for each class, generate all needed samples
        for i, model in enumerate(self.estimators_):
            model_mask = (labels == i)
            N_model = model_mask.sum()
            if N_model > 0:
                X[model_mask] = model.sample(N_model,
                                             random_state=random_state)

        return X, self.classes_[labels]
