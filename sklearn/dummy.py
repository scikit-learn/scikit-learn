# Author: Mathieu Blondel <mathieu@mblondel.org>
#         Arnaud Joly <a.joly@ulg.ac.be>
#         Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
# License: BSD 3 clause

import warnings
import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .base import MultiOutputMixin
from .utils import check_random_state
from .utils.validation import _num_samples
from .utils.validation import check_array
from .utils.validation import check_consistent_length
from .utils.validation import check_is_fitted
from .utils.random import random_choice_csc
from .utils.stats import _weighted_percentile
from .utils.multiclass import class_distribution
from .utils import deprecated


class DummyClassifier(MultiOutputMixin, ClassifierMixin, BaseEstimator):
    """
    DummyClassifier is a classifier that makes predictions using simple rules.

    This classifier is useful as a simple baseline to compare with other
    (real) classifiers. Do not use it for real problems.

    Read more in the :ref:`User Guide <dummy_estimators>`.

    Parameters
    ----------
    strategy : str, default="stratified"
        Strategy to use to generate predictions.

        * "stratified": generates predictions by respecting the training
          set's class distribution.
        * "most_frequent": always predicts the most frequent label in the
          training set.
        * "prior": always predicts the class that maximizes the class prior
          (like "most_frequent") and ``predict_proba`` returns the class prior.
        * "uniform": generates predictions uniformly at random.
        * "constant": always predicts a constant label that is provided by
          the user. This is useful for metrics that evaluate a non-majority
          class

          .. versionadded:: 0.17
             Dummy Classifier now supports prior fitting strategy using
             parameter *prior*.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    constant : int or str or array of shape = [n_outputs]
        The explicit constant as predicted by the "constant" strategy. This
        parameter is useful only for the "constant" strategy.

    Attributes
    ----------
    classes_ : array or list of array of shape = [n_classes]
        Class labels for each output.

    n_classes_ : array or list of array of shape = [n_classes]
        Number of label for each output.

    class_prior_ : array or list of array of shape = [n_classes]
        Probability of each class for each output.

    n_outputs_ : int,
        Number of outputs.

    sparse_output_ : bool,
        True if the array returned from predict is to be in sparse CSC format.
        Is automatically set to True if the input y is passed in sparse format.
    """

    def __init__(self, strategy="stratified", random_state=None,
                 constant=None):
        self.strategy = strategy
        self.random_state = random_state
        self.constant = constant

    def fit(self, X, y, sample_weight=None):
        """Fit the random classifier.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        self : object
        """
        allowed_strategies = ("most_frequent", "stratified", "uniform",
                              "constant", "prior")
        if self.strategy not in allowed_strategies:
            raise ValueError("Unknown strategy type: %s, expected one of %s."
                             % (self.strategy, allowed_strategies))

        if self.strategy == "uniform" and sp.issparse(y):
            y = y.toarray()
            warnings.warn('A local copy of the target data has been converted '
                          'to a numpy array. Predicting on sparse target data '
                          'with the uniform strategy would not save memory '
                          'and would be slower.',
                          UserWarning)

        self.sparse_output_ = sp.issparse(y)

        if not self.sparse_output_:
            y = np.atleast_1d(y)

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]

        check_consistent_length(X, y, sample_weight)

        if self.strategy == "constant":
            if self.constant is None:
                raise ValueError("Constant target value has to be specified "
                                 "when the constant strategy is used.")
            else:
                constant = np.reshape(np.atleast_1d(self.constant), (-1, 1))
                if constant.shape[0] != self.n_outputs_:
                    raise ValueError("Constant target value should have "
                                     "shape (%d, 1)." % self.n_outputs_)

        (self.classes_,
         self.n_classes_,
         self.class_prior_) = class_distribution(y, sample_weight)

        if self.strategy == "constant":
            for k in range(self.n_outputs_):
                if not any(constant[k][0] == c for c in self.classes_[k]):
                    # Checking in case of constant strategy if the constant
                    # provided by the user is in y.
                    err_msg = ("The constant target value must be present in "
                               "the training data. You provided constant={}. "
                               "Possible values are: {}."
                               .format(self.constant, list(self.classes_[k])))
                    raise ValueError(err_msg)

        if self.n_outputs_ == 1:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]
            self.class_prior_ = self.class_prior_[0]

        return self

    def predict(self, X):
        """Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        Returns
        -------
        y : array, shape = [n_samples] or [n_samples, n_outputs]
            Predicted target values for X.
        """
        check_is_fitted(self)

        # numpy random_state expects Python int and not long as size argument
        # under Windows
        n_samples = _num_samples(X)
        rs = check_random_state(self.random_state)

        n_classes_ = self.n_classes_
        classes_ = self.classes_
        class_prior_ = self.class_prior_
        constant = self.constant
        if self.n_outputs_ == 1:
            # Get same type even for self.n_outputs_ == 1
            n_classes_ = [n_classes_]
            classes_ = [classes_]
            class_prior_ = [class_prior_]
            constant = [constant]
        # Compute probability only once
        if self.strategy == "stratified":
            proba = self.predict_proba(X)
            if self.n_outputs_ == 1:
                proba = [proba]

        if self.sparse_output_:
            class_prob = None
            if self.strategy in ("most_frequent", "prior"):
                classes_ = [np.array([cp.argmax()]) for cp in class_prior_]

            elif self.strategy == "stratified":
                class_prob = class_prior_

            elif self.strategy == "uniform":
                raise ValueError("Sparse target prediction is not "
                                 "supported with the uniform strategy")

            elif self.strategy == "constant":
                classes_ = [np.array([c]) for c in constant]

            y = random_choice_csc(n_samples, classes_, class_prob,
                                  self.random_state)
        else:
            if self.strategy in ("most_frequent", "prior"):
                y = np.tile([classes_[k][class_prior_[k].argmax()] for
                             k in range(self.n_outputs_)], [n_samples, 1])

            elif self.strategy == "stratified":
                y = np.vstack([classes_[k][proba[k].argmax(axis=1)] for
                               k in range(self.n_outputs_)]).T

            elif self.strategy == "uniform":
                ret = [classes_[k][rs.randint(n_classes_[k], size=n_samples)]
                       for k in range(self.n_outputs_)]
                y = np.vstack(ret).T

            elif self.strategy == "constant":
                y = np.tile(self.constant, (n_samples, 1))

            if self.n_outputs_ == 1:
                y = np.ravel(y)

        return y

    def predict_proba(self, X):
        """
        Return probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        Returns
        -------
        P : array-like or list of array-lke of shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered arithmetically, for each
            output.
        """
        check_is_fitted(self)

        # numpy random_state expects Python int and not long as size argument
        # under Windows
        n_samples = _num_samples(X)
        rs = check_random_state(self.random_state)

        n_classes_ = self.n_classes_
        classes_ = self.classes_
        class_prior_ = self.class_prior_
        constant = self.constant
        if self.n_outputs_ == 1:
            # Get same type even for self.n_outputs_ == 1
            n_classes_ = [n_classes_]
            classes_ = [classes_]
            class_prior_ = [class_prior_]
            constant = [constant]

        P = []
        for k in range(self.n_outputs_):
            if self.strategy == "most_frequent":
                ind = class_prior_[k].argmax()
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0
            elif self.strategy == "prior":
                out = np.ones((n_samples, 1)) * class_prior_[k]

            elif self.strategy == "stratified":
                out = rs.multinomial(1, class_prior_[k], size=n_samples)
                out = out.astype(np.float64)

            elif self.strategy == "uniform":
                out = np.ones((n_samples, n_classes_[k]), dtype=np.float64)
                out /= n_classes_[k]

            elif self.strategy == "constant":
                ind = np.where(classes_[k] == constant[k])
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0

            P.append(out)

        if self.n_outputs_ == 1:
            P = P[0]

        return P

    def predict_log_proba(self, X):
        """
        Return log probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        Returns
        -------
        P : array-like or list of array-like of shape = [n_samples, n_classes]
            Returns the log probability of the sample for each class in
            the model, where classes are ordered arithmetically for each
            output.
        """
        proba = self.predict_proba(X)
        if self.n_outputs_ == 1:
            return np.log(proba)
        else:
            return [np.log(p) for p in proba]

    def _more_tags(self):
        return {'poor_score': True, 'no_validation': True}

    def score(self, X, y, sample_weight=None):
        """Returns the mean accuracy on the given test data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : {array-like, None}
            Test samples with shape = (n_samples, n_features) or
            None. Passing None as test samples gives the same result
            as passing real test samples, since DummyClassifier
            operates independently of the sampled observations.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True labels for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of self.predict(X) wrt. y.

        """
        if X is None:
            X = np.zeros(shape=(len(y), 1))
        return super().score(X, y, sample_weight)

    @deprecated(
        "The outputs_2d_ attribute is deprecated in version 0.22 "
        "and will be removed in version 0.24. It is equivalent to "
        "n_outputs_ > 1."
    )
    @property
    def outputs_2d_(self):
        return self.n_outputs_ != 1


class DummyRegressor(MultiOutputMixin, RegressorMixin, BaseEstimator):
    """
    DummyRegressor is a regressor that makes predictions using
    simple rules.

    This regressor is useful as a simple baseline to compare with other
    (real) regressors. Do not use it for real problems.

    Read more in the :ref:`User Guide <dummy_estimators>`.

    Parameters
    ----------
    strategy : str
        Strategy to use to generate predictions.

        * "mean": always predicts the mean of the training set
        * "median": always predicts the median of the training set
        * "quantile": always predicts a specified quantile of the training set,
          provided with the quantile parameter.
        * "constant": always predicts a constant value that is provided by
          the user.

    constant : int or float or array of shape = [n_outputs]
        The explicit constant as predicted by the "constant" strategy. This
        parameter is useful only for the "constant" strategy.

    quantile : float in [0.0, 1.0]
        The quantile to predict using the "quantile" strategy. A quantile of
        0.5 corresponds to the median, while 0.0 to the minimum and 1.0 to the
        maximum.

    Attributes
    ----------
    constant_ : array, shape (1, n_outputs)
        Mean or median or quantile of the training targets or constant value
        given by the user.

    n_outputs_ : int,
        Number of outputs.
    """

    def __init__(self, strategy="mean", constant=None, quantile=None):
        self.strategy = strategy
        self.constant = constant
        self.quantile = quantile

    def fit(self, X, y, sample_weight=None):
        """Fit the random regressor.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        self : object
        """
        allowed_strategies = ("mean", "median", "quantile", "constant")
        if self.strategy not in allowed_strategies:
            raise ValueError("Unknown strategy type: %s, expected one of %s."
                             % (self.strategy, allowed_strategies))

        y = check_array(y, ensure_2d=False)
        if len(y) == 0:
            raise ValueError("y must not be empty.")

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))
        self.n_outputs_ = y.shape[1]

        check_consistent_length(X, y, sample_weight)

        if self.strategy == "mean":
            self.constant_ = np.average(y, axis=0, weights=sample_weight)

        elif self.strategy == "median":
            if sample_weight is None:
                self.constant_ = np.median(y, axis=0)
            else:
                self.constant_ = [_weighted_percentile(y[:, k], sample_weight,
                                                       percentile=50.)
                                  for k in range(self.n_outputs_)]

        elif self.strategy == "quantile":
            if self.quantile is None or not np.isscalar(self.quantile):
                raise ValueError("Quantile must be a scalar in the range "
                                 "[0.0, 1.0], but got %s." % self.quantile)

            percentile = self.quantile * 100.0
            if sample_weight is None:
                self.constant_ = np.percentile(y, axis=0, q=percentile)
            else:
                self.constant_ = [_weighted_percentile(y[:, k], sample_weight,
                                                       percentile=percentile)
                                  for k in range(self.n_outputs_)]

        elif self.strategy == "constant":
            if self.constant is None:
                raise TypeError("Constant target value has to be specified "
                                "when the constant strategy is used.")

            self.constant = check_array(self.constant,
                                        accept_sparse=['csr', 'csc', 'coo'],
                                        ensure_2d=False, ensure_min_samples=0)

            if self.n_outputs_ != 1 and self.constant.shape[0] != y.shape[1]:
                raise ValueError(
                    "Constant target value should have "
                    "shape (%d, 1)." % y.shape[1])

            self.constant_ = self.constant

        self.constant_ = np.reshape(self.constant_, (1, -1))
        return self

    def predict(self, X, return_std=False):
        """
        Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, object with finite length or shape}
            Training data, requires length = n_samples

        return_std : boolean, optional
            Whether to return the standard deviation of posterior prediction.
            All zeros in this case.

        Returns
        -------
        y : array, shape = [n_samples] or [n_samples, n_outputs]
            Predicted target values for X.

        y_std : array, shape = [n_samples] or [n_samples, n_outputs]
            Standard deviation of predictive distribution of query points.
        """
        check_is_fitted(self)
        n_samples = _num_samples(X)

        y = np.full((n_samples, self.n_outputs_), self.constant_,
                    dtype=np.array(self.constant_).dtype)
        y_std = np.zeros((n_samples, self.n_outputs_))

        if self.n_outputs_ == 1:
            y = np.ravel(y)
            y_std = np.ravel(y_std)

        return (y, y_std) if return_std else y

    def _more_tags(self):
        return {'poor_score': True, 'no_validation': True}

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the residual
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the total
        sum of squares ((y_true - y_true.mean()) ** 2).sum().
        The best possible score is 1.0 and it can be negative (because the
        model can be arbitrarily worse). A constant model that always
        predicts the expected value of y, disregarding the input features,
        would get a R^2 score of 0.0.

        Parameters
        ----------
        X : {array-like, None}
            Test samples with shape = (n_samples, n_features) or None.
            For some estimators this may be a
            precomputed kernel matrix instead, shape = (n_samples,
            n_samples_fitted], where n_samples_fitted is the number of
            samples used in the fitting for the estimator.
            Passing None as test samples gives the same result
            as passing real test samples, since DummyRegressor
            operates independently of the sampled observations.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.
        """
        if X is None:
            X = np.zeros(shape=(len(y), 1))
        return super().score(X, y, sample_weight)

    @deprecated(
        "The outputs_2d_ attribute is deprecated in version 0.22 "
        "and will be removed in version 0.24. It is equivalent to "
        "n_outputs_ > 1."
    )
    @property
    def outputs_2d_(self):
        return self.n_outputs_ != 1
