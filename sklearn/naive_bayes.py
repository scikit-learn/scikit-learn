"""
The :mod:`sklearn.naive_bayes` module implements Naive Bayes algorithms. These
are supervised learning methods based on applying Bayes' theorem with strong
(naive) feature independence assumptions.
"""

# Author: Vincent Michel <vincent.michel@inria.fr>
#         Minor fixes by Fabian Pedregosa
#         Amit Aides <amitibo@tx.technion.ac.il>
#         Yehuda Finkelstein <yehudaf@tx.technion.ac.il>
#         Lars Buitinck
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         (parts based on earlier work by Mathieu Blondel)
#         Andrey V. Melnik <andrey.melnik.maths@gmail.com>
#
# License: BSD 3 clause
import warnings

from abc import ABCMeta, abstractmethod
from numbers import Real, Integral

import numpy as np
from scipy.special import logsumexp
from joblib import Parallel

from .base import BaseEstimator, ClassifierMixin
from .base import clone
from .preprocessing import binarize
from .preprocessing import LabelBinarizer
from .preprocessing import label_binarize
from .utils import deprecated
from .utils.extmath import safe_sparse_dot
from .utils.multiclass import _check_partial_fit_first_call
from .utils.validation import check_is_fitted, check_non_negative
from .utils.validation import _check_sample_weight
from .utils.validation import column_or_1d
from .utils.metaestimators import _BaseComposition
from .utils import _safe_indexing, _get_column_indices
from .utils import _print_elapsed_time
from .utils import Bunch
from .utils.fixes import delayed
from .utils._estimator_html_repr import _VisualBlock
from .compose._column_transformer import _is_empty_column_selection

from .utils._param_validation import Interval, Hidden, StrOptions

__all__ = [
    "BernoulliNB",
    "GaussianNB",
    "MultinomialNB",
    "ComplementNB",
    "CategoricalNB",
    "ColumnwiseNB",
]


class _BaseNB(ClassifierMixin, BaseEstimator, metaclass=ABCMeta):
    """Abstract base class for naive Bayes estimators"""

    @abstractmethod
    def _joint_log_likelihood(self, X):
        """Compute the unnormalized posterior log probability of X

        I.e. ``log P(c) + log P(x|c)`` for all rows x of X, as an array-like of
        shape (n_samples, n_classes).

        Public methods predict, predict_proba, predict_log_proba, and
        predict_joint_log_proba pass the input through _check_X before handing it
        over to _joint_log_likelihood. The term "joint log likelihood" is used
        interchangibly with "joint log probability".
        """

    @abstractmethod
    def _check_X(self, X):
        """To be overridden in subclasses with the actual checks.

        Only used in predict* methods.
        """

    def predict_joint_log_proba(self, X):
        """Return joint log probability estimates for the test vector X.

        For each row x of X and class y, the joint log probability is given by
        ``log P(x, y) = log P(y) + log P(x|y),``
        where ``log P(y)`` is the class prior probability and ``log P(x|y)`` is
        the class-conditional probability.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        C : ndarray of shape (n_samples, n_classes)
            Returns the joint log-probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        X = self._check_X(X)
        return self._joint_log_likelihood(X)

    def predict(self, X):
        """
        Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            Predicted target values for X.
        """
        check_is_fitted(self)
        X = self._check_X(X)
        jll = self._joint_log_likelihood(X)
        return self.classes_[np.argmax(jll, axis=1)]

    def predict_log_proba(self, X):
        """
        Return log-probability estimates for the test vector X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        C : array-like of shape (n_samples, n_classes)
            Returns the log-probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        X = self._check_X(X)
        jll = self._joint_log_likelihood(X)
        # normalize by P(x) = P(f_1, ..., f_n)
        log_prob_x = logsumexp(jll, axis=1)
        return jll - np.atleast_2d(log_prob_x).T

    def predict_proba(self, X):
        """
        Return probability estimates for the test vector X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        C : array-like of shape (n_samples, n_classes)
            Returns the probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute :term:`classes_`.
        """
        return np.exp(self.predict_log_proba(X))


class GaussianNB(_BaseNB):
    """
    Gaussian Naive Bayes (GaussianNB).

    Can perform online updates to model parameters via :meth:`partial_fit`.
    For details on algorithm used to update feature means and variance online,
    see Stanford CS tech report STAN-CS-79-773 by Chan, Golub, and LeVeque:

        http://i.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf

    Read more in the :ref:`User Guide <gaussian_naive_bayes>`.

    Parameters
    ----------
    priors : array-like of shape (n_classes,), default=None
        Prior probabilities of the classes. If specified, the priors are not
        adjusted according to the data.

    var_smoothing : float, default=1e-9
        Portion of the largest variance of all features that is added to
        variances for calculation stability.

        .. versionadded:: 0.20

    Attributes
    ----------
    class_count_ : ndarray of shape (n_classes,)
        number of training samples observed in each class.

    class_prior_ : ndarray of shape (n_classes,)
        probability of each class.

    classes_ : ndarray of shape (n_classes,)
        class labels known to the classifier.

    epsilon_ : float
        absolute additive value to variances.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    sigma_ : ndarray of shape (n_classes, n_features)
        Variance of each feature per class.

        .. deprecated:: 1.0
           `sigma_` is deprecated in 1.0 and will be removed in 1.2.
           Use `var_` instead.

    var_ : ndarray of shape (n_classes, n_features)
        Variance of each feature per class.

        .. versionadded:: 1.0

    theta_ : ndarray of shape (n_classes, n_features)
        mean of each feature per class.

    See Also
    --------
    BernoulliNB : Naive Bayes classifier for multivariate Bernoulli models.
    CategoricalNB : Naive Bayes classifier for categorical features.
    ComplementNB : Complement Naive Bayes classifier.
    MultinomialNB : Naive Bayes classifier for multinomial models.

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
    >>> clf_pf = GaussianNB()
    >>> clf_pf.partial_fit(X, Y, np.unique(Y))
    GaussianNB()
    >>> print(clf_pf.predict([[-0.8, -1]]))
    [1]
    """

    _parameter_constraints: dict = {
        "priors": ["array-like", None],
        "var_smoothing": [Interval(Real, 0, None, closed="left")],
    }

    def __init__(self, *, priors=None, var_smoothing=1e-9):
        self.priors = priors
        self.var_smoothing = var_smoothing

    def fit(self, X, y, sample_weight=None):
        """Fit Gaussian Naive Bayes according to X, y.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

            .. versionadded:: 0.17
               Gaussian Naive Bayes supports fitting with *sample_weight*.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._validate_params()
        y = self._validate_data(y=y)
        return self._partial_fit(
            X, y, np.unique(y), _refit=True, sample_weight=sample_weight
        )

    def _check_X(self, X):
        """Validate X, used only in predict* methods."""
        return self._validate_data(X, reset=False)

    @staticmethod
    def _update_mean_variance(n_past, mu, var, X, sample_weight=None):
        """Compute online update of Gaussian mean and variance.

        Given starting sample count, mean, and variance, a new set of
        points X, and optionally sample weights, return the updated mean and
        variance. (NB - each dimension (column) in X is treated as independent
        -- you get variance, not covariance).

        Can take scalar mean and variance, or vector mean and variance to
        simultaneously update a number of independent Gaussians.

        See Stanford CS tech report STAN-CS-79-773 by Chan, Golub, and LeVeque:

        http://i.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf

        Parameters
        ----------
        n_past : int
            Number of samples represented in old mean and variance. If sample
            weights were given, this should contain the sum of sample
            weights represented in old mean and variance.

        mu : array-like of shape (number of Gaussians,)
            Means for Gaussians in original set.

        var : array-like of shape (number of Gaussians,)
            Variances for Gaussians in original set.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        total_mu : array-like of shape (number of Gaussians,)
            Updated mean for each Gaussian over the combined set.

        total_var : array-like of shape (number of Gaussians,)
            Updated variance for each Gaussian over the combined set.
        """
        if X.shape[0] == 0:
            return mu, var

        # Compute (potentially weighted) mean and variance of new datapoints
        if sample_weight is not None:
            n_new = float(sample_weight.sum())
            new_mu = np.average(X, axis=0, weights=sample_weight)
            new_var = np.average((X - new_mu) ** 2, axis=0, weights=sample_weight)
        else:
            n_new = X.shape[0]
            new_var = np.var(X, axis=0)
            new_mu = np.mean(X, axis=0)

        if n_past == 0:
            return new_mu, new_var

        n_total = float(n_past + n_new)

        # Combine mean of old and new data, taking into consideration
        # (weighted) number of observations
        total_mu = (n_new * new_mu + n_past * mu) / n_total

        # Combine variance of old and new data, taking into consideration
        # (weighted) number of observations. This is achieved by combining
        # the sum-of-squared-differences (ssd)
        old_ssd = n_past * var
        new_ssd = n_new * new_var
        total_ssd = old_ssd + new_ssd + (n_new * n_past / n_total) * (mu - new_mu) ** 2
        total_var = total_ssd / n_total

        return total_mu, total_var

    def partial_fit(self, X, y, classes=None, sample_weight=None):
        """Incremental fit on a batch of samples.

        This method is expected to be called several times consecutively
        on different chunks of a dataset so as to implement out-of-core
        or online learning.

        This is especially useful when the whole dataset is too big to fit in
        memory at once.

        This method has some performance and numerical stability overhead,
        hence it is better to call partial_fit on chunks of data that are
        as large as possible (as long as fitting in the memory budget) to
        hide the overhead.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

            .. versionadded:: 0.17

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._validate_params()

        return self._partial_fit(
            X, y, classes, _refit=False, sample_weight=sample_weight
        )

    def _partial_fit(self, X, y, classes=None, _refit=False, sample_weight=None):
        """Actual implementation of Gaussian NB fitting.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        _refit : bool, default=False
            If true, act as though this were the first time we called
            _partial_fit (ie, throw away any past fitting and start over).

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
        """
        if _refit:
            self.classes_ = None

        first_call = _check_partial_fit_first_call(self, classes)
        X, y = self._validate_data(X, y, reset=first_call)
        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        # If the ratio of data variance between dimensions is too small, it
        # will cause numerical errors. To address this, we artificially
        # boost the variance by epsilon, a small fraction of the standard
        # deviation of the largest dimension.
        self.epsilon_ = self.var_smoothing * np.var(X, axis=0).max()

        if first_call:
            # This is the first call to partial_fit:
            # initialize various cumulative counters
            n_features = X.shape[1]
            n_classes = len(self.classes_)
            self.theta_ = np.zeros((n_classes, n_features))
            self.var_ = np.zeros((n_classes, n_features))

            self.class_count_ = np.zeros(n_classes, dtype=np.float64)

            # Initialise the class prior
            # Take into account the priors
            if self.priors is not None:
                priors = np.asarray(self.priors)
                # Check that the provided prior matches the number of classes
                if len(priors) != n_classes:
                    raise ValueError("Number of priors must match number of classes.")
                # Check that the sum is 1
                if not np.isclose(priors.sum(), 1.0):
                    raise ValueError("The sum of the priors should be 1.")
                # Check that the priors are non-negative
                if (priors < 0).any():
                    raise ValueError("Priors must be non-negative.")
                self.class_prior_ = priors
            else:
                # Initialize the priors to zeros for each class
                self.class_prior_ = np.zeros(len(self.classes_), dtype=np.float64)
        else:
            if X.shape[1] != self.theta_.shape[1]:
                msg = "Number of features %d does not match previous data %d."
                raise ValueError(msg % (X.shape[1], self.theta_.shape[1]))
            # Put epsilon back in each time
            self.var_[:, :] -= self.epsilon_

        classes = self.classes_

        unique_y = np.unique(y)
        unique_y_in_classes = np.in1d(unique_y, classes)

        if not np.all(unique_y_in_classes):
            raise ValueError(
                "The target label(s) %s in y do not exist in the initial classes %s"
                % (unique_y[~unique_y_in_classes], classes)
            )

        for y_i in unique_y:
            i = classes.searchsorted(y_i)
            X_i = X[y == y_i, :]

            if sample_weight is not None:
                sw_i = sample_weight[y == y_i]
                N_i = sw_i.sum()
            else:
                sw_i = None
                N_i = X_i.shape[0]

            new_theta, new_sigma = self._update_mean_variance(
                self.class_count_[i], self.theta_[i, :], self.var_[i, :], X_i, sw_i
            )

            self.theta_[i, :] = new_theta
            self.var_[i, :] = new_sigma
            self.class_count_[i] += N_i

        self.var_[:, :] += self.epsilon_

        # Update if only no priors is provided
        if self.priors is None:
            # Empirical prior, with sample_weight taken into account
            self.class_prior_ = self.class_count_ / self.class_count_.sum()

        return self

    def _joint_log_likelihood(self, X):
        joint_log_likelihood = []
        for i in range(np.size(self.classes_)):
            jointi = np.log(self.class_prior_[i])
            n_ij = -0.5 * np.sum(np.log(2.0 * np.pi * self.var_[i, :]))
            n_ij -= 0.5 * np.sum(((X - self.theta_[i, :]) ** 2) / (self.var_[i, :]), 1)
            joint_log_likelihood.append(jointi + n_ij)

        joint_log_likelihood = np.array(joint_log_likelihood).T
        return joint_log_likelihood

    @deprecated(  # type: ignore
        "Attribute `sigma_` was deprecated in 1.0 and will be removed in"
        "1.2. Use `var_` instead."
    )
    @property
    def sigma_(self):
        return self.var_


class _BaseDiscreteNB(_BaseNB):
    """Abstract base class for naive Bayes on discrete/categorical data

    Any estimator based on this class should provide:

    __init__
    _joint_log_likelihood(X) as per _BaseNB
    _update_feature_log_prob(alpha)
    _count(X, Y)
    """

    _parameter_constraints: dict = {
        "alpha": [Interval(Real, 0, None, closed="left"), "array-like"],
        "fit_prior": ["boolean"],
        "class_prior": ["array-like", None],
        "force_alpha": ["boolean", Hidden(StrOptions({"warn"}))],
    }

    def __init__(self, alpha=1.0, fit_prior=True, class_prior=None, force_alpha="warn"):
        self.alpha = alpha
        self.fit_prior = fit_prior
        self.class_prior = class_prior
        self.force_alpha = force_alpha

    @abstractmethod
    def _count(self, X, Y):
        """Update counts that are used to calculate probabilities.

        The counts make up a sufficient statistic extracted from the data.
        Accordingly, this method is called each time `fit` or `partial_fit`
        update the model. `class_count_` and `feature_count_` must be updated
        here along with any model specific counts.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.
        Y : ndarray of shape (n_samples, n_classes)
            Binarized class labels.
        """

    @abstractmethod
    def _update_feature_log_prob(self, alpha):
        """Update feature log probabilities based on counts.

        This method is called each time `fit` or `partial_fit` update the
        model.

        Parameters
        ----------
        alpha : float
            smoothing parameter. See :meth:`_check_alpha`.
        """

    def _check_X(self, X):
        """Validate X, used only in predict* methods."""
        return self._validate_data(X, accept_sparse="csr", reset=False)

    def _check_X_y(self, X, y, reset=True):
        """Validate X and y in fit methods."""
        return self._validate_data(X, y, accept_sparse="csr", reset=reset)

    def _update_class_log_prior(self, class_prior=None):
        """Update class log priors.

        The class log priors are based on `class_prior`, class count or the
        number of classes. This method is called each time `fit` or
        `partial_fit` update the model.
        """
        n_classes = len(self.classes_)
        if class_prior is not None:
            if len(class_prior) != n_classes:
                raise ValueError("Number of priors must match number of classes.")
            self.class_log_prior_ = np.log(class_prior)
        elif self.fit_prior:
            with warnings.catch_warnings():
                # silence the warning when count is 0 because class was not yet
                # observed
                warnings.simplefilter("ignore", RuntimeWarning)
                log_class_count = np.log(self.class_count_)

            # empirical prior, with sample_weight taken into account
            self.class_log_prior_ = log_class_count - np.log(self.class_count_.sum())
        else:
            self.class_log_prior_ = np.full(n_classes, -np.log(n_classes))

    def _check_alpha(self):
        alpha = (
            np.asarray(self.alpha) if not isinstance(self.alpha, Real) else self.alpha
        )
        alpha_min = np.min(alpha)
        if isinstance(alpha, np.ndarray):
            if not alpha.shape[0] == self.n_features_in_:
                raise ValueError(
                    "When alpha is an array, it should contains `n_features`. "
                    f"Got {alpha.shape[0]} elements instead of {self.n_features_in_}."
                )
            # check that all alpha are positive
            if alpha_min < 0:
                raise ValueError("All values in alpha must be greater than 0.")
        alpha_lower_bound = 1e-10
        # TODO(1.4): Replace w/ deprecation of self.force_alpha
        # See gh #22269
        _force_alpha = self.force_alpha
        if _force_alpha == "warn" and alpha_min < alpha_lower_bound:
            _force_alpha = False
            warnings.warn(
                "The default value for `force_alpha` will change to `True` in 1.4. To"
                " suppress this warning, manually set the value of `force_alpha`.",
                FutureWarning,
            )
        if alpha_min < alpha_lower_bound and not _force_alpha:
            warnings.warn(
                "alpha too small will result in numeric errors, setting alpha ="
                f" {alpha_lower_bound:.1e}. Use `force_alpha=True` to keep alpha"
                " unchanged."
            )
            return np.maximum(alpha, alpha_lower_bound)
        return alpha

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
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        first_call = not hasattr(self, "classes_")

        if first_call:
            self._validate_params()

        X, y = self._check_X_y(X, y, reset=first_call)
        _, n_features = X.shape

        if _check_partial_fit_first_call(self, classes):
            # This is the first call to partial_fit:
            # initialize various cumulative counters
            n_classes = len(classes)
            self._init_counters(n_classes, n_features)

        Y = label_binarize(y, classes=self.classes_)
        if Y.shape[1] == 1:
            if len(self.classes_) == 2:
                Y = np.concatenate((1 - Y, Y), axis=1)
            else:  # degenerate case: just one class
                Y = np.ones_like(Y)

        if X.shape[0] != Y.shape[0]:
            msg = "X.shape[0]=%d and y.shape[0]=%d are incompatible."
            raise ValueError(msg % (X.shape[0], y.shape[0]))

        # label_binarize() returns arrays with dtype=np.int64.
        # We convert it to np.float64 to support sample_weight consistently
        Y = Y.astype(np.float64, copy=False)
        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)
            sample_weight = np.atleast_2d(sample_weight)
            Y *= sample_weight.T

        class_prior = self.class_prior

        # Count raw events from data before updating the class log prior
        # and feature log probas
        self._count(X, Y)

        # XXX: OPTIM: we could introduce a public finalization method to
        # be called by the user explicitly just once after several consecutive
        # calls to partial_fit and prior any call to predict[_[log_]proba]
        # to avoid computing the smooth log probas at each call to partial fit
        alpha = self._check_alpha()
        self._update_feature_log_prob(alpha)
        self._update_class_log_prior(class_prior=class_prior)
        return self

    def fit(self, X, y, sample_weight=None):
        """Fit Naive Bayes classifier according to X, y.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._validate_params()
        X, y = self._check_X_y(X, y)
        _, n_features = X.shape

        labelbin = LabelBinarizer()
        Y = labelbin.fit_transform(y)
        self.classes_ = labelbin.classes_
        if Y.shape[1] == 1:
            if len(self.classes_) == 2:
                Y = np.concatenate((1 - Y, Y), axis=1)
            else:  # degenerate case: just one class
                Y = np.ones_like(Y)

        # LabelBinarizer().fit_transform() returns arrays with dtype=np.int64.
        # We convert it to np.float64 to support sample_weight consistently;
        # this means we also don't have to cast X to floating point
        if sample_weight is not None:
            Y = Y.astype(np.float64, copy=False)
            sample_weight = _check_sample_weight(sample_weight, X)
            sample_weight = np.atleast_2d(sample_weight)
            Y *= sample_weight.T

        class_prior = self.class_prior

        # Count raw events from data before updating the class log prior
        # and feature log probas
        n_classes = Y.shape[1]
        self._init_counters(n_classes, n_features)
        self._count(X, Y)
        alpha = self._check_alpha()
        self._update_feature_log_prob(alpha)
        self._update_class_log_prior(class_prior=class_prior)
        return self

    def _init_counters(self, n_classes, n_features):
        self.class_count_ = np.zeros(n_classes, dtype=np.float64)
        self.feature_count_ = np.zeros((n_classes, n_features), dtype=np.float64)

    def _more_tags(self):
        return {"poor_score": True}

    # TODO: Remove in 1.2
    # mypy error: Decorated property not supported
    @deprecated(  # type: ignore
        "Attribute `n_features_` was deprecated in version 1.0 and will be "
        "removed in 1.2. Use `n_features_in_` instead."
    )
    @property
    def n_features_(self):
        return self.n_features_in_


class MultinomialNB(_BaseDiscreteNB):
    """
    Naive Bayes classifier for multinomial models.

    The multinomial Naive Bayes classifier is suitable for classification with
    discrete features (e.g., word counts for text classification). The
    multinomial distribution normally requires integer feature counts. However,
    in practice, fractional counts such as tf-idf may also work.

    Read more in the :ref:`User Guide <multinomial_naive_bayes>`.

    Parameters
    ----------
    alpha : float or array-like of shape (n_features,), default=1.0
        Additive (Laplace/Lidstone) smoothing parameter
        (set alpha=0 and force_alpha=True, for no smoothing).

    force_alpha : bool, default=False
        If False and alpha is less than 1e-10, it will set alpha to
        1e-10. If True, alpha will remain unchanged. This may cause
        numerical errors if alpha is too close to 0.

        .. versionadded:: 1.2
        .. deprecated:: 1.2
           The default value of `force_alpha` will change to `True` in v1.4.

    fit_prior : bool, default=True
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    class_prior : array-like of shape (n_classes,), default=None
        Prior probabilities of the classes. If specified, the priors are not
        adjusted according to the data.

    Attributes
    ----------
    class_count_ : ndarray of shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    class_log_prior_ : ndarray of shape (n_classes,)
        Smoothed empirical log probability for each class.

    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier

    feature_count_ : ndarray of shape (n_classes, n_features)
        Number of samples encountered for each (class, feature)
        during fitting. This value is weighted by the sample weight when
        provided.

    feature_log_prob_ : ndarray of shape (n_classes, n_features)
        Empirical log probability of features
        given a class, ``P(x_i|y)``.

    n_features_ : int
        Number of features of each sample.

        .. deprecated:: 1.0
            Attribute `n_features_` was deprecated in version 1.0 and will be
            removed in 1.2. Use `n_features_in_` instead.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    BernoulliNB : Naive Bayes classifier for multivariate Bernoulli models.
    CategoricalNB : Naive Bayes classifier for categorical features.
    ComplementNB : Complement Naive Bayes classifier.
    GaussianNB : Gaussian Naive Bayes.

    References
    ----------
    C.D. Manning, P. Raghavan and H. Schuetze (2008). Introduction to
    Information Retrieval. Cambridge University Press, pp. 234-265.
    https://nlp.stanford.edu/IR-book/html/htmledition/naive-bayes-text-classification-1.html

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(1)
    >>> X = rng.randint(5, size=(6, 100))
    >>> y = np.array([1, 2, 3, 4, 5, 6])
    >>> from sklearn.naive_bayes import MultinomialNB
    >>> clf = MultinomialNB(force_alpha=True)
    >>> clf.fit(X, y)
    MultinomialNB(force_alpha=True)
    >>> print(clf.predict(X[2:3]))
    [3]
    """

    def __init__(
        self, *, alpha=1.0, force_alpha="warn", fit_prior=True, class_prior=None
    ):
        super().__init__(
            alpha=alpha,
            fit_prior=fit_prior,
            class_prior=class_prior,
            force_alpha=force_alpha,
        )

    def _more_tags(self):
        return {"requires_positive_X": True}

    def _count(self, X, Y):
        """Count and smooth feature occurrences."""
        check_non_negative(X, "MultinomialNB (input X)")
        self.feature_count_ += safe_sparse_dot(Y.T, X)
        self.class_count_ += Y.sum(axis=0)

    def _update_feature_log_prob(self, alpha):
        """Apply smoothing to raw counts and recompute log probabilities"""
        smoothed_fc = self.feature_count_ + alpha
        smoothed_cc = smoothed_fc.sum(axis=1)

        self.feature_log_prob_ = np.log(smoothed_fc) - np.log(
            smoothed_cc.reshape(-1, 1)
        )

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""
        return safe_sparse_dot(X, self.feature_log_prob_.T) + self.class_log_prior_


class ComplementNB(_BaseDiscreteNB):
    """The Complement Naive Bayes classifier described in Rennie et al. (2003).

    The Complement Naive Bayes classifier was designed to correct the "severe
    assumptions" made by the standard Multinomial Naive Bayes classifier. It is
    particularly suited for imbalanced data sets.

    Read more in the :ref:`User Guide <complement_naive_bayes>`.

    .. versionadded:: 0.20

    Parameters
    ----------
    alpha : float or array-like of shape (n_features,), default=1.0
        Additive (Laplace/Lidstone) smoothing parameter
        (set alpha=0 and force_alpha=True, for no smoothing).

    force_alpha : bool, default=False
        If False and alpha is less than 1e-10, it will set alpha to
        1e-10. If True, alpha will remain unchanged. This may cause
        numerical errors if alpha is too close to 0.

        .. versionadded:: 1.2
        .. deprecated:: 1.2
           The default value of `force_alpha` will change to `True` in v1.4.

    fit_prior : bool, default=True
        Only used in edge case with a single class in the training set.

    class_prior : array-like of shape (n_classes,), default=None
        Prior probabilities of the classes. Not used.

    norm : bool, default=False
        Whether or not a second normalization of the weights is performed. The
        default behavior mirrors the implementations found in Mahout and Weka,
        which do not follow the full algorithm described in Table 9 of the
        paper.

    Attributes
    ----------
    class_count_ : ndarray of shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    class_log_prior_ : ndarray of shape (n_classes,)
        Smoothed empirical log probability for each class. Only used in edge
        case with a single class in the training set.

    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier

    feature_all_ : ndarray of shape (n_features,)
        Number of samples encountered for each feature during fitting. This
        value is weighted by the sample weight when provided.

    feature_count_ : ndarray of shape (n_classes, n_features)
        Number of samples encountered for each (class, feature) during fitting.
        This value is weighted by the sample weight when provided.

    feature_log_prob_ : ndarray of shape (n_classes, n_features)
        Empirical weights for class complements.

    n_features_ : int
        Number of features of each sample.

        .. deprecated:: 1.0
            Attribute `n_features_` was deprecated in version 1.0 and will be
            removed in 1.2. Use `n_features_in_` instead.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    BernoulliNB : Naive Bayes classifier for multivariate Bernoulli models.
    CategoricalNB : Naive Bayes classifier for categorical features.
    GaussianNB : Gaussian Naive Bayes.
    MultinomialNB : Naive Bayes classifier for multinomial models.

    References
    ----------
    Rennie, J. D., Shih, L., Teevan, J., & Karger, D. R. (2003).
    Tackling the poor assumptions of naive bayes text classifiers. In ICML
    (Vol. 3, pp. 616-623).
    https://people.csail.mit.edu/jrennie/papers/icml03-nb.pdf

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(1)
    >>> X = rng.randint(5, size=(6, 100))
    >>> y = np.array([1, 2, 3, 4, 5, 6])
    >>> from sklearn.naive_bayes import ComplementNB
    >>> clf = ComplementNB(force_alpha=True)
    >>> clf.fit(X, y)
    ComplementNB(force_alpha=True)
    >>> print(clf.predict(X[2:3]))
    [3]
    """

    _parameter_constraints: dict = {
        **_BaseDiscreteNB._parameter_constraints,
        "norm": ["boolean"],
    }

    def __init__(
        self,
        *,
        alpha=1.0,
        force_alpha="warn",
        fit_prior=True,
        class_prior=None,
        norm=False,
    ):
        super().__init__(
            alpha=alpha,
            force_alpha=force_alpha,
            fit_prior=fit_prior,
            class_prior=class_prior,
        )
        self.norm = norm

    def _more_tags(self):
        return {"requires_positive_X": True}

    def _count(self, X, Y):
        """Count feature occurrences."""
        check_non_negative(X, "ComplementNB (input X)")
        self.feature_count_ += safe_sparse_dot(Y.T, X)
        self.class_count_ += Y.sum(axis=0)
        self.feature_all_ = self.feature_count_.sum(axis=0)

    def _update_feature_log_prob(self, alpha):
        """Apply smoothing to raw counts and compute the weights."""
        comp_count = self.feature_all_ + alpha - self.feature_count_
        logged = np.log(comp_count / comp_count.sum(axis=1, keepdims=True))
        # _BaseNB.predict uses argmax, but ComplementNB operates with argmin.
        if self.norm:
            summed = logged.sum(axis=1, keepdims=True)
            feature_log_prob = logged / summed
        else:
            feature_log_prob = -logged
        self.feature_log_prob_ = feature_log_prob

    def _joint_log_likelihood(self, X):
        """Calculate the class scores for the samples in X."""
        jll = safe_sparse_dot(X, self.feature_log_prob_.T)
        if len(self.classes_) == 1:
            jll += self.class_log_prior_
        return jll


class BernoulliNB(_BaseDiscreteNB):
    """Naive Bayes classifier for multivariate Bernoulli models.

    Like MultinomialNB, this classifier is suitable for discrete data. The
    difference is that while MultinomialNB works with occurrence counts,
    BernoulliNB is designed for binary/boolean features.

    Read more in the :ref:`User Guide <bernoulli_naive_bayes>`.

    Parameters
    ----------
    alpha : float or array-like of shape (n_features,), default=1.0
        Additive (Laplace/Lidstone) smoothing parameter
        (set alpha=0 and force_alpha=True, for no smoothing).

    force_alpha : bool, default=False
        If False and alpha is less than 1e-10, it will set alpha to
        1e-10. If True, alpha will remain unchanged. This may cause
        numerical errors if alpha is too close to 0.

        .. versionadded:: 1.2
        .. deprecated:: 1.2
           The default value of `force_alpha` will change to `True` in v1.4.

    binarize : float or None, default=0.0
        Threshold for binarizing (mapping to booleans) of sample features.
        If None, input is presumed to already consist of binary vectors.

    fit_prior : bool, default=True
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    class_prior : array-like of shape (n_classes,), default=None
        Prior probabilities of the classes. If specified, the priors are not
        adjusted according to the data.

    Attributes
    ----------
    class_count_ : ndarray of shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    class_log_prior_ : ndarray of shape (n_classes,)
        Log probability of each class (smoothed).

    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier

    feature_count_ : ndarray of shape (n_classes, n_features)
        Number of samples encountered for each (class, feature)
        during fitting. This value is weighted by the sample weight when
        provided.

    feature_log_prob_ : ndarray of shape (n_classes, n_features)
        Empirical log probability of features given a class, P(x_i|y).

    n_features_ : int
        Number of features of each sample.

        .. deprecated:: 1.0
            Attribute `n_features_` was deprecated in version 1.0 and will be
            removed in 1.2. Use `n_features_in_` instead.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    CategoricalNB : Naive Bayes classifier for categorical features.
    ComplementNB : The Complement Naive Bayes classifier
        described in Rennie et al. (2003).
    GaussianNB : Gaussian Naive Bayes (GaussianNB).
    MultinomialNB : Naive Bayes classifier for multinomial models.

    References
    ----------
    C.D. Manning, P. Raghavan and H. Schuetze (2008). Introduction to
    Information Retrieval. Cambridge University Press, pp. 234-265.
    https://nlp.stanford.edu/IR-book/html/htmledition/the-bernoulli-model-1.html

    A. McCallum and K. Nigam (1998). A comparison of event models for naive
    Bayes text classification. Proc. AAAI/ICML-98 Workshop on Learning for
    Text Categorization, pp. 41-48.

    V. Metsis, I. Androutsopoulos and G. Paliouras (2006). Spam filtering with
    naive Bayes -- Which naive Bayes? 3rd Conf. on Email and Anti-Spam (CEAS).

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(1)
    >>> X = rng.randint(5, size=(6, 100))
    >>> Y = np.array([1, 2, 3, 4, 4, 5])
    >>> from sklearn.naive_bayes import BernoulliNB
    >>> clf = BernoulliNB(force_alpha=True)
    >>> clf.fit(X, Y)
    BernoulliNB(force_alpha=True)
    >>> print(clf.predict(X[2:3]))
    [3]
    """

    _parameter_constraints: dict = {
        **_BaseDiscreteNB._parameter_constraints,
        "binarize": [None, Interval(Real, 0, None, closed="left")],
    }

    def __init__(
        self,
        *,
        alpha=1.0,
        force_alpha="warn",
        binarize=0.0,
        fit_prior=True,
        class_prior=None,
    ):
        super().__init__(
            alpha=alpha,
            fit_prior=fit_prior,
            class_prior=class_prior,
            force_alpha=force_alpha,
        )
        self.binarize = binarize

    def _check_X(self, X):
        """Validate X, used only in predict* methods."""
        X = super()._check_X(X)
        if self.binarize is not None:
            X = binarize(X, threshold=self.binarize)
        return X

    def _check_X_y(self, X, y, reset=True):
        X, y = super()._check_X_y(X, y, reset=reset)
        if self.binarize is not None:
            X = binarize(X, threshold=self.binarize)
        return X, y

    def _count(self, X, Y):
        """Count and smooth feature occurrences."""
        self.feature_count_ += safe_sparse_dot(Y.T, X)
        self.class_count_ += Y.sum(axis=0)

    def _update_feature_log_prob(self, alpha):
        """Apply smoothing to raw counts and recompute log probabilities"""
        smoothed_fc = self.feature_count_ + alpha
        smoothed_cc = self.class_count_ + alpha * 2

        self.feature_log_prob_ = np.log(smoothed_fc) - np.log(
            smoothed_cc.reshape(-1, 1)
        )

    def _joint_log_likelihood(self, X):
        """Calculate the posterior log probability of the samples X"""
        n_features = self.feature_log_prob_.shape[1]
        n_features_X = X.shape[1]

        if n_features_X != n_features:
            raise ValueError(
                "Expected input with %d features, got %d instead"
                % (n_features, n_features_X)
            )

        neg_prob = np.log(1 - np.exp(self.feature_log_prob_))
        # Compute  neg_prob · (1 - X).T  as  ∑neg_prob - X · neg_prob
        jll = safe_sparse_dot(X, (self.feature_log_prob_ - neg_prob).T)
        jll += self.class_log_prior_ + neg_prob.sum(axis=1)

        return jll


class CategoricalNB(_BaseDiscreteNB):
    """Naive Bayes classifier for categorical features.

    The categorical Naive Bayes classifier is suitable for classification with
    discrete features that are categorically distributed. The categories of
    each feature are drawn from a categorical distribution.

    Read more in the :ref:`User Guide <categorical_naive_bayes>`.

    Parameters
    ----------
    alpha : float, default=1.0
        Additive (Laplace/Lidstone) smoothing parameter
        (set alpha=0 and force_alpha=True, for no smoothing).

    force_alpha : bool, default=False
        If False and alpha is less than 1e-10, it will set alpha to
        1e-10. If True, alpha will remain unchanged. This may cause
        numerical errors if alpha is too close to 0.

        .. versionadded:: 1.2
        .. deprecated:: 1.2
           The default value of `force_alpha` will change to `True` in v1.4.

    fit_prior : bool, default=True
        Whether to learn class prior probabilities or not.
        If false, a uniform prior will be used.

    class_prior : array-like of shape (n_classes,), default=None
        Prior probabilities of the classes. If specified, the priors are not
        adjusted according to the data.

    min_categories : int or array-like of shape (n_features,), default=None
        Minimum number of categories per feature.

        - integer: Sets the minimum number of categories per feature to
          `n_categories` for each features.
        - array-like: shape (n_features,) where `n_categories[i]` holds the
          minimum number of categories for the ith column of the input.
        - None (default): Determines the number of categories automatically
          from the training data.

        .. versionadded:: 0.24

    Attributes
    ----------
    category_count_ : list of arrays of shape (n_features,)
        Holds arrays of shape (n_classes, n_categories of respective feature)
        for each feature. Each array provides the number of samples
        encountered for each class and category of the specific feature.

    class_count_ : ndarray of shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    class_log_prior_ : ndarray of shape (n_classes,)
        Smoothed empirical log probability for each class.

    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier

    feature_log_prob_ : list of arrays of shape (n_features,)
        Holds arrays of shape (n_classes, n_categories of respective feature)
        for each feature. Each array provides the empirical log probability
        of categories given the respective feature and class, ``P(x_i|y)``.

    n_features_ : int
        Number of features of each sample.

        .. deprecated:: 1.0
            Attribute `n_features_` was deprecated in version 1.0 and will be
            removed in 1.2. Use `n_features_in_` instead.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_categories_ : ndarray of shape (n_features,), dtype=np.int64
        Number of categories for each feature. This value is
        inferred from the data or set by the minimum number of categories.

        .. versionadded:: 0.24

    See Also
    --------
    BernoulliNB : Naive Bayes classifier for multivariate Bernoulli models.
    ComplementNB : Complement Naive Bayes classifier.
    GaussianNB : Gaussian Naive Bayes.
    MultinomialNB : Naive Bayes classifier for multinomial models.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(1)
    >>> X = rng.randint(5, size=(6, 100))
    >>> y = np.array([1, 2, 3, 4, 5, 6])
    >>> from sklearn.naive_bayes import CategoricalNB
    >>> clf = CategoricalNB(force_alpha=True)
    >>> clf.fit(X, y)
    CategoricalNB(force_alpha=True)
    >>> print(clf.predict(X[2:3]))
    [3]
    """

    _parameter_constraints: dict = {
        **_BaseDiscreteNB._parameter_constraints,
        "min_categories": [
            None,
            "array-like",
            Interval(Integral, 1, None, closed="left"),
        ],
        "alpha": [Interval(Real, 0, None, closed="left")],
    }

    def __init__(
        self,
        *,
        alpha=1.0,
        force_alpha="warn",
        fit_prior=True,
        class_prior=None,
        min_categories=None,
    ):
        super().__init__(
            alpha=alpha,
            force_alpha=force_alpha,
            fit_prior=fit_prior,
            class_prior=class_prior,
        )
        self.min_categories = min_categories

    def fit(self, X, y, sample_weight=None):
        """Fit Naive Bayes classifier according to X, y.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features. Here, each feature of X is
            assumed to be from a different categorical distribution.
            It is further assumed that all categories of each feature are
            represented by the numbers 0, ..., n - 1, where n refers to the
            total number of categories for the given feature. This can, for
            instance, be achieved with the help of OrdinalEncoder.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        return super().fit(X, y, sample_weight=sample_weight)

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
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features. Here, each feature of X is
            assumed to be from a different categorical distribution.
            It is further assumed that all categories of each feature are
            represented by the numbers 0, ..., n - 1, where n refers to the
            total number of categories for the given feature. This can, for
            instance, be achieved with the help of OrdinalEncoder.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        return super().partial_fit(X, y, classes, sample_weight=sample_weight)

    def _more_tags(self):
        return {"requires_positive_X": True}

    def _check_X(self, X):
        """Validate X, used only in predict* methods."""
        X = self._validate_data(
            X, dtype="int", accept_sparse=False, force_all_finite=True, reset=False
        )
        check_non_negative(X, "CategoricalNB (input X)")
        return X

    def _check_X_y(self, X, y, reset=True):
        X, y = self._validate_data(
            X, y, dtype="int", accept_sparse=False, force_all_finite=True, reset=reset
        )
        check_non_negative(X, "CategoricalNB (input X)")
        return X, y

    def _init_counters(self, n_classes, n_features):
        self.class_count_ = np.zeros(n_classes, dtype=np.float64)
        self.category_count_ = [np.zeros((n_classes, 0)) for _ in range(n_features)]

    @staticmethod
    def _validate_n_categories(X, min_categories):
        # rely on max for n_categories categories are encoded between 0...n-1
        n_categories_X = X.max(axis=0) + 1
        min_categories_ = np.array(min_categories)
        if min_categories is not None:
            if not np.issubdtype(min_categories_.dtype, np.signedinteger):
                raise ValueError(
                    "'min_categories' should have integral type. Got "
                    f"{min_categories_.dtype} instead."
                )
            n_categories_ = np.maximum(n_categories_X, min_categories_, dtype=np.int64)
            if n_categories_.shape != n_categories_X.shape:
                raise ValueError(
                    f"'min_categories' should have shape ({X.shape[1]},"
                    ") when an array-like is provided. Got"
                    f" {min_categories_.shape} instead."
                )
            return n_categories_
        else:
            return n_categories_X

    def _count(self, X, Y):
        def _update_cat_count_dims(cat_count, highest_feature):
            diff = highest_feature + 1 - cat_count.shape[1]
            if diff > 0:
                # we append a column full of zeros for each new category
                return np.pad(cat_count, [(0, 0), (0, diff)], "constant")
            return cat_count

        def _update_cat_count(X_feature, Y, cat_count, n_classes):
            for j in range(n_classes):
                mask = Y[:, j].astype(bool)
                if Y.dtype.type == np.int64:
                    weights = None
                else:
                    weights = Y[mask, j]
                counts = np.bincount(X_feature[mask], weights=weights)
                indices = np.nonzero(counts)[0]
                cat_count[j, indices] += counts[indices]

        self.class_count_ += Y.sum(axis=0)
        self.n_categories_ = self._validate_n_categories(X, self.min_categories)
        for i in range(self.n_features_in_):
            X_feature = X[:, i]
            self.category_count_[i] = _update_cat_count_dims(
                self.category_count_[i], self.n_categories_[i] - 1
            )
            _update_cat_count(
                X_feature, Y, self.category_count_[i], self.class_count_.shape[0]
            )

    def _update_feature_log_prob(self, alpha):
        feature_log_prob = []
        for i in range(self.n_features_in_):
            smoothed_cat_count = self.category_count_[i] + alpha
            smoothed_class_count = smoothed_cat_count.sum(axis=1)
            feature_log_prob.append(
                np.log(smoothed_cat_count) - np.log(smoothed_class_count.reshape(-1, 1))
            )
        self.feature_log_prob_ = feature_log_prob

    def _joint_log_likelihood(self, X):
        self._check_n_features(X, reset=False)
        jll = np.zeros((X.shape[0], self.class_count_.shape[0]))
        for i in range(self.n_features_in_):
            indices = X[:, i]
            jll += self.feature_log_prob_[i][:, indices].T
        total_ll = jll + self.class_log_prior_
        return total_ll


def _fit_one(estimator, X, y, message_clsname="", message=None, **fit_params):
    """Call ``estimator.fit`` and print elapsed time message.

    See :func:`sklearn.pipeline._fit_one`.
    """
    with _print_elapsed_time(message_clsname, message):
        return estimator.fit(X, y, **fit_params)


def _partial_fit_one(estimator, X, y, message_clsname="", message=None, **fit_params):
    """Call ``estimator.partial_fit`` and print elapsed time message.

    See :func:`sklearn.pipeline._fit_one`.
    """
    with _print_elapsed_time(message_clsname, message):
        return estimator.partial_fit(X, y, **fit_params)


def _jll_one(estimator, X):
    """Call ``estimator.predict_joint_log_proba``.

    See :func:`sklearn.pipeline._transform_one`.
    """
    return estimator.predict_joint_log_proba(X)


class ColumnwiseNB(_BaseNB, _BaseComposition):
    """
    Column-wise Naive Bayes meta-estimator.

    This estimator combines various naive Bayes estimators by applying them
    to different column subsets of the input and joining their predictions
    according to the naive Bayes assumption. This is useful when features are
    heterogeneous and follow different kinds of distributions.

    Parameters
    ----------
    nb_estimators : list of tuples
        List of (name, nb_estimator, columns) tuples specifying the naive Bayes
        estimators to be combined into a single naive Bayes meta-estimator.

        name : str
            Name of the naive Bayes estimator. Like in
            :class:`~sklearn.pipeline.Pipeline`,
            :class:`~sklearn.pipeline.FeatureUnion`,
            and :class:`~sklearn.compose.ColumnTransformer`, this allows the
            subestimator and its parameters to be set using :term:`set_params`
            and searched in grid search.
        nb_estimator : estimator
            The estimator must support :term:`fit` or :term:`partial_fit`,
            depending on how the meta-estimator is fitted. In addition, the
            estimator must support ``predict_joint_log_proba`` method, which
            takes :term:`X` of shape (n_samples, n_features) and returns a
            numpy array of shape (n_samples, n_classes) containing joint
            log-probabilities, ``log P(x,y)`` for each sample point and class.
        columns :  str, array-like of str, int, array-like of int, \
                array-like of bool, slice or callable
            Indexes the data on its second axis. Integers are interpreted as
            positional columns, while strings can reference DataFrame columns
            by name.  A scalar string or int should be used where
            ``nb_estimator`` expects X to be a 1d array-like (vector),
            otherwise a 2d array will be passed to the transformer.
            A callable is passed the input data `X` and can return any of the
            above. To select multiple columns by name or dtype, you can use
            :obj:`~sklearn.compose.make_column_selector`.

    priors : array-like of shape (n_classes,) or str, default=None
        Prior probabilities of classes. If unspecified, the priors are
        calculated as relative frequencies of classes in the training data.
        If str, the priors are taken from the estimator with the given name.
        If array-like, the same priors might have to be specified manually in
        each sub-estimator, in order to ensure consistent predictions.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : bool, default=False
        If True, the time elapsed while fitting each estimator will be
        printed as it is completed.

    Attributes
    ----------
    estimators_ : list of tuples
        List of ``(name, fitted_nb_estimator, columns)`` tuples, which follow
        the order of `nb_estimators`. Here, ``fitted_nb_estimator`` is a fitted naive
        Bayes estimator, except when ``columns`` presents an empty selection of
        columns, in which case it is the original unfitted ``nb_estimator``. If
        the original specification of ``columns`` in ``nb_estimators`` was a
        callable, then ``columns`` is converted to a list of column indices.

    named_estimators_ : :class:`~sklearn.utils.Bunch`
        Read-only attribute to access any subestimator by given name.
        Keys are estimator names and values are the fitted estimators, except
        when a subestimator does not require fitting (i.e., when ``columns`` is
        an empty set of indices).

    class_prior_ : ndarray of shape (n_classes,)
        Prior probabilities of classes used in the naive Bayes meta-estimator,
        which are calculated as relative frequencies, extracted from
        subestimators, or provided, according to the value of `priors`
        at initialization.

    class_count_ : ndarray of shape (n_classes,)
        Number of samples encountered for each class during fitting. This
        value is weighted by the sample weight when provided.

    n_classes_ : int
        The number of classes known to the naive Bayes classifier, `n_classes`.

    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if `X` has
        feature names that are all strings.

    See Also
    --------
    BernoulliNB : Naive Bayes classifier for multivariate Bernoulli models.
    CategoricalNB : Naive Bayes classifier for categorical features.
    ComplementNB : Complement Naive Bayes classifier.
    MultinomialNB : Naive Bayes classifier for multinomial models.
    GaussianNB : Gaussian Naive Bayes.
    :class:`~sklearn.compose.ColumnTransformer` : Applies transformers to columns.

    Notes
    -----
    ColumnwiseNB combines multiple naive Bayes estimators by expressing the
    overall joint probability ``P(x,y)`` through ``P(x_i,y)``, the joint
    probabilities of the subestimators::

        ``Log P(x,y) = Log P(x_1,y) + ... + Log P(x_N,y) - (N - 1) Log P(y)``,

    where ``N`` denotes ``n_estimators``, the number of estimators.
    It is implicitly assumed that the class log priors are finite and agree
    between the estimators and the subestimator::

        ``- inf < Log P(y) = Log P(y|1) = ... = Log P(y|N)``.

    The meta-estimators does not check if this condition holds. Meaningless
    results, including ``NaN``, may be produced by ColumnwiseNB if the class
    priors differ or contain a zero probability.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.RandomState(1)
    >>> X = rng.randint(5, size=(6, 100))
    >>> y = np.array([0, 0, 1, 1, 2, 2])
    >>> from sklearn.naive_bayes import MultinomialNB, GaussianNB, ColumnwiseNB
    >>> clf = ColumnwiseNB(nb_estimators=[('mnb1', MultinomialNB(), [0, 1]),
    ...                                   ('mnb2', MultinomialNB(), [3, 4]),
    ...                                   ('gnb1', GaussianNB(), [5])])
    >>> clf.fit(X, y)
    ColumnwiseNB(nb_estimators=[('mnb1', MultinomialNB(), [0, 1]),
                            ('mnb2', MultinomialNB(), [3, 4]),
                            ('gnb1', GaussianNB(), [5])])
    >>> print(clf.predict(X))
    [0 0 1 0 2 2]
    """

    _required_parameters = ["nb_estimators"]

    _parameter_constraints = {
        "nb_estimators": "no_validation",
        "priors": ["array-like", str, None],
        "n_jobs": [Integral, None],
        "verbose": ["verbose"],
    }

    def _log_message(self, name, idx, total):
        if not self.verbose:
            return None
        return "(%d of %d) Processing %s" % (idx, total, name)

    def __init__(self, nb_estimators, *, priors=None, n_jobs=None, verbose=False):
        self.nb_estimators = nb_estimators
        self.priors = priors
        self.n_jobs = n_jobs
        self.verbose = verbose

    def _check_X(self, X):
        """Validate X, used only in predict* methods."""
        # The meta-estimator checks for feature names only. Other checks
        # and conversion to numpy array are performed by subestimators.
        # It is important that X is not modified by the meta-estimator, and X's
        # columns are passed to an estimator as they are. Note that estimators
        # may modify (a copy of) X. E.g., BernoulliNB._check_X binarises the
        # input.
        self._check_feature_names(X, reset=False)
        return X

    def _joint_log_likelihood(self, X):
        """Calculate the meta-estimator's joint log-probability ``log P(x,y)``."""
        estimators = self._iter(fitted=True, replace_strings=True)
        all_jlls = Parallel(n_jobs=self.n_jobs)(
            delayed(_jll_one)(estimator=nb_estimator, X=_safe_indexing(X, cols, axis=1))
            for (_, nb_estimator, cols) in estimators
        )
        n_estimators = len(all_jlls)
        log_prior = np.log(self.class_prior_)
        return np.sum(all_jlls, axis=0) - (n_estimators - 1) * log_prior

    def _validate_estimators(self, check_partial=False):
        # Check if estimators have fit/partial_fit and joint log prob methods
        # Validate estimator names via _BaseComposition._validate_names(self, names)
        if not self.nb_estimators:
            raise ValueError(
                "A list of naive Bayes estimators must be provided "
                "in the form [(name, nb_estimator, columns), ... ]."
            )
        names, estimators, _ = zip(*self.nb_estimators)
        for e in estimators:
            if (not check_partial) and (
                not (hasattr(e, "fit") and hasattr(e, "predict_joint_log_proba"))
            ):
                raise TypeError(
                    "Estimators must be naive Bayes estimators implementing "
                    "`fit` and `predict_joint_log_proba` methods."
                )
            if check_partial and (
                not (
                    hasattr(e, "partial_fit") and hasattr(e, "predict_joint_log_proba")
                )
            ):
                raise TypeError(
                    "Estimators must be Naive Bayes estimators implementing "
                    "`partial_fit` and `predict_joint_log_proba` methods."
                )
        self._validate_names(names)

    def _validate_column_callables(self, X):
        """
        Convert callable column specifications and store into self._columns.

        Empty-set columns do not enjoy any special treatment.
        """
        # Almost a verbatim copy of ColumnTransformer._validate_column_callables().
        # Consider refactoring in the future.
        # Unlike ColumnTransformer, this estimator does not need to output a
        # dataframe or validate a the remainder, so _estimator_to_input_indices
        # is not really needed, but retained for consistency with
        # ColumnTransformer code.
        all_columns = []
        estimator_to_input_indices = {}
        for name, _, columns in self.nb_estimators:
            if callable(columns):
                columns = columns(X)
            all_columns.append(columns)
            estimator_to_input_indices[name] = _get_column_indices(X, columns)
        self._columns = all_columns
        self._estimator_to_input_indices = estimator_to_input_indices

    @property
    def named_estimators_(self):
        """Access the fitted naive Bayes subestimators by name.

        Read-only attribute to access any estimators by given name.
        Keys are estimators names and values are the fitted estimator
        objects.
        """
        # Almost a verbatim copy of ColumnTransformer.named_transformers_
        # Use Bunch object to improve autocomplete
        return Bunch(**{name: e for name, e, _ in self.estimators_})

    def _iter(self, *, fitted=False, replace_strings=False):
        """Generate ``(name, nb_estimator, columns)`` tuples.

        This is a private method, similar to ColumnTransformer._iter.
        Must not be called before _validate_column_callables.

        Parameters
        ----------
        fitted : bool, default=False
            If False, returns tuples from self.estimators (user-specified), but
            callable columns are replaced with a list column names or indices.
            If True, returns tuples from self.estimators_ (fitted), where
            columns are processed as well.

        replace_strings : bool, default=False
            If True, omits the estimators that do not require fitting, i.e those
            with empty-set columns. The name `replace_strings` is a relic of
            ColumnTransformer implementation, where `passthrough` and `drop`
            required replacement and omission, respectively.

        Yields
        ------
        tuple
            of the form ``(name, nb_estimator, columns)``.

        Notes
        -----
        Loop through estimators from this generator with the following
        parameters, depending on the purpose:

        self._iter(fitted=False, replace_strings=True) :
            fit, 1st partial_fit
        self._iter(fitted=True, replace_strings=True) :
            further partial_fit, predict
        self._iter(fitted=False, replace_strings=False) :
            update fitted estimators. Note that special treatment is required
            for unfitted estimators (those with empty-set columns)!
        self._iter(fitted=True, replace_strings=False) :
            not used here. The usecase in ColumnTransformer would be sorting
            out the transformed output and its column names.
        do not use in :
            a Bunch accessor named_estimators_;
            input validation _validate_estimators, _validate_column_callables;
            parameter management: get_params_, set_params_, _estimators.
        """
        if fitted:
            for name, estimator, cols in self.estimators_:
                if replace_strings and _is_empty_column_selection(cols):
                    continue
                else:
                    yield (name, estimator, cols)
        else:  # fitted=False
            for (name, estimator, _), cols in zip(self.nb_estimators, self._columns):
                if replace_strings and _is_empty_column_selection(cols):
                    continue
                else:
                    yield (name, estimator, cols)

    def _update_class_prior(self):
        """Update class prior after most of the fitting as done."""
        if self.priors is None:  # calculcate empirical prior from counts
            priors = self.class_count_ / self.class_count_.sum()
        elif isinstance(self.priors, str):  # extract prior from estimator
            name = self.priors
            e = self.named_estimators_[name]
            if getattr(e, "class_prior_", None) is not None:
                priors = e.class_prior_
            elif getattr(e, "class_log_prior_", None) is not None:
                priors = np.exp(e.class_log_prior_)
            else:
                raise AttributeError(
                    f"Unable to extract class prior from estimator {name}, as "
                    "it does not have class_prior_ or class_log_prior_ "
                    "attributes."
                )
        else:  # check the provided prior
            priors = np.asarray(self.priors)
        # Check the prior in any case.
        if len(priors) != self.n_classes_:
            raise ValueError("Number of priors must match number of classes.")
        if not np.isclose(priors.sum(), 1.0):
            raise ValueError("The sum of the priors should be 1.")
        if (priors < 0).any():
            raise ValueError("Priors must be non-negative.")
        self.class_prior_ = priors

    def _update_fitted_estimators(self, fitted_estimators):
        """Update tuples in self.estimators_ with fitted_estimators provided.

        Callable columns are replaced with sets of actual str or int indices.
        Estimators that don't require fitting are passed as they were,
        without cloning.
        """
        estimators_ = []
        fitted_estimators = iter(fitted_estimators)

        for name, nb_estimator, cols in self._iter():
            if not _is_empty_column_selection(cols):
                updated_nb_estimator = next(fitted_estimators)
            else:  # don't advance fitted_estimators; use original
                updated_nb_estimator = nb_estimator
            estimators_.append((name, updated_nb_estimator, cols))
        self.estimators_ = estimators_

    def fit(self, X, y, sample_weight=None):
        """Fit the naive Bayes meta-estimator.

        Calls `fit` of each subestimator ``nb_estimator``.  Only a corresponding
        subset of columns of `X` is passed to each subestimator; `sample_weight`
        and `y` are passed to the subestimators as they are.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples
            and `n_features` is the number of features.
        y : array-like of shape (n_samples,)
            Target values.
        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._validate_params()
        self._check_feature_names(X, reset=True)
        # TODO: Consider overriding BaseEstimator._check_feature_names
        # Currently, when X has all str feature names, all features are
        # registered in self.feature_names_in no matter if they are used or not.
        self._check_n_features(X, reset=True)
        self._validate_estimators()
        self._validate_column_callables(X)
        # Consistency checks for X, y are delegated to subestimators

        # Subestimators get original sample_weight. This is for class counts:
        if sample_weight is not None:
            weights = _check_sample_weight(sample_weight, X=y, copy=True)

        # We would use sklearn.utils.multiclass.class_distribution, but it does
        # not return class_count, which we want as well.
        if sample_weight is None:
            self.classes_, self.class_count_ = np.unique(
                column_or_1d(y), return_counts=True
            )
        else:
            self.classes_ = np.unique(column_or_1d(y))
            counts = np.zeros(len(self.classes_))
            for i, c in enumerate(self.classes_):
                counts[i] = (weights * (column_or_1d(y) == c)).sum()
            self.class_count_ = counts
        self.n_classes_ = len(self.classes_)

        estimators = list(self._iter(fitted=False, replace_strings=True))
        fitted_estimators = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_one)(
                estimator=clone(nb_estimator),
                X=_safe_indexing(X, cols, axis=1),
                y=y,
                message_clsname="ColumnwiseNB",
                message=self._log_message(name, idx, len(estimators)),
                sample_weight=sample_weight,
            )
            for idx, (name, nb_estimator, cols) in enumerate(estimators, 1)
        )
        self._update_fitted_estimators(fitted_estimators)
        self._update_class_prior()
        return self

    def partial_fit(self, X, y, classes=None, sample_weight=None):
        """Fit incrementally the naive Bayes meta-estimator on a batch of samples.

        Calls `partial_fit` of each subestimator. Only a corresponding
        subset of columns of `X` is passed to each subestimator. `classes`,
        `sample_weight` and 'y' are passed to the subestimators as they are.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like of shape (n_samples,), default=None
            Weights applied to individual samples (1. for unweighted).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        first_call = not hasattr(self, "classes_")
        if first_call:
            self._validate_params()
            self._check_feature_names(X, reset=True)
            self._check_n_features(X, reset=True)
            self._validate_estimators(check_partial=True)
            self._validate_column_callables(X)
        else:
            self._check_feature_names(X, reset=False)
            self._check_n_features(X, reset=False)
        # Consistency checks for X, y are delegated to subestimators

        # Subestimators get original sample_weight. This is for class counts:
        if sample_weight is not None:
            weights = _check_sample_weight(sample_weight, X=y, copy=True)

        # Subestimators should've checked classes. We set classes_ for counts
        # and so that first_call becomes False at next partial_fit call.
        _check_partial_fit_first_call(self, classes)

        # We don't use sklearn.utils.multiclass.class_distribution, because it
        # neither returns class_count, nor is suitable for partial_fit.
        if sample_weight is None:
            counts = np.zeros(len(self.classes_))
            for i, c in enumerate(self.classes_):
                counts[i] = (column_or_1d(y) == c).sum()
        else:
            counts = np.zeros(len(self.classes_))
            for i, c in enumerate(self.classes_):
                counts[i] = (weights * (column_or_1d(y) == c)).sum()

        if first_call:
            self.n_classes_ = len(self.classes_)
            self.class_count_ = counts
        else:
            self.class_count_ += counts

        estimators = list(self._iter(fitted=not first_call, replace_strings=True))
        fitted_estimators = Parallel(n_jobs=self.n_jobs)(
            delayed(_partial_fit_one)(
                estimator=clone(nb_estimator) if first_call else nb_estimator,
                X=_safe_indexing(X, cols, axis=1),
                y=y,
                message_clsname="ColumnwiseNB",
                message=self._log_message(name, idx, len(estimators)),
                classes=classes,
                sample_weight=sample_weight,
            )
            for idx, (name, nb_estimator, cols) in enumerate(estimators, 1)
        )
        self._update_fitted_estimators(fitted_estimators)
        self._update_class_prior()
        return self

    @property
    def _estimators(self):
        """Internal list of subestimators.

        This is for the implementation of get_params via BaseComposition._get_params,
        which expects lists of tuples of len 2.
        """
        # Implemented in the image and likeness of ColumnTranformer._transformers
        try:
            return [(name, e) for name, e, _ in self.nb_estimators]
        except (TypeError, ValueError):
            # This try-except clause is needed to pass the test from test_common.py:
            # test_estimators_do_not_raise_errors_in_init_or_set_params().
            # ColumnTransformer does the same. See PR #21355 for details.
            return self.nb_estimators

    @_estimators.setter
    def _estimators(self, value):
        # Implemented in the image and likeness of ColumnTranformer._transformers
        # TODO: Is renaming or changing the order legal? Swap `name` and `_`?
        self.nb_estimators = [
            (name, e, col)
            for ((name, e), (_, _, col)) in zip(value, self.nb_estimators)
        ]

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Returns the parameters listed in the constructor as well as the
        subestimators contained within the `estimators` of the `ColumnwiseNB`
        instance.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        # Implemented in the image and likeness of ColumnTranformer.get_params
        return self._get_params("_estimators", deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``. Note that you
        can directly set the parameters of the estimators contained in
        `estimators` of `ColumnwiseNB`.

        Parameters
        ----------
        **kwargs : dict
            Estimator parameters.

        Returns
        -------
        self : ColumnwiseNB
            This estimator.
        """
        # Implemented in the image and likeness of ColumnTranformer.set_params
        self._set_params("_estimators", **kwargs)
        return self

    def _sk_visual_block_(self):
        """HTML representation of this estimator."""
        names, estimators, name_details = zip(*self.nb_estimators)
        return _VisualBlock(
            "parallel", estimators, names=names, name_details=name_details
        )
