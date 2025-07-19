"""Dummy estimators that implement simple rules of thumb."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Integral, Real

import numpy as np
import scipy.sparse as sp

from .base import (
    BaseEstimator,
    ClassifierMixin,
    MultiOutputMixin,
    RegressorMixin,
    _fit_context,
)
from .utils import check_random_state
from .utils._param_validation import Interval, StrOptions
from .utils.multiclass import class_distribution
from .utils.random import _random_choice_csc
from .utils.stats import _weighted_percentile
from .utils.validation import (
    _check_sample_weight,
    _num_samples,
    check_array,
    check_consistent_length,
    check_is_fitted,
    validate_data,
)


class DummyClassifier(MultiOutputMixin, ClassifierMixin, BaseEstimator):
    """
DummyClassifier is a simple model that doesn't learn from your data.

It is used to test if your real model is better than just guessing.

This model does NOT look at the input features (X).  
It only uses the output labels (y) to make simple predictions based on the strategy you choose.

Strategies (Ways to Guess)
--------------------------
- "most_frequent": always picks the label that appeared the most in training.
- "prior": same as "most_frequent", but used when you only care about label probabilities.
- "stratified": randomly picks labels, but keeps the same proportions as your data.
- "uniform": randomly picks any label with equal chance.
- "constant": always gives the same label that you choose.

Terms Explained
---------------
X : input values (like age, marks, etc.)
y : actual answers or labels (like pass/fail, spam/not spam)
fit : train the model using your data
predict : guess the label for new data
random_state : use a fixed number to get same random results every time

Example
-------
>>> from sklearn.dummy import DummyClassifier
>>> X = [[1], [2], [3], [4]]
>>> y = [0, 1, 1, 1]  # most common label is 1
>>> model = DummyClassifier(strategy="most_frequent")
>>> model.fit(X, y)
DummyClassifier(strategy='most_frequent')
>>> model.predict([[5], [6]])
array([1, 1])  # always predicts 1
"""

    _parameter_constraints: dict = {
        "strategy": [
            StrOptions({"most_frequent", "prior", "stratified", "uniform", "constant"})
        ],
        "random_state": ["random_state"],
        "constant": [Integral, str, "array-like", None],
    }

    def __init__(self, *, strategy="prior", random_state=None, constant=None):
        self.strategy = strategy
        self.random_state = random_state
        self.constant = constant

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """Fit the baseline classifier.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        validate_data(self, X, skip_check_array=True)

        self._strategy = self.strategy

        if self._strategy == "uniform" and sp.issparse(y):
            y = y.toarray()
            warnings.warn(
                (
                    "A local copy of the target data has been converted "
                    "to a numpy array. Predicting on sparse target data "
                    "with the uniform strategy would not save memory "
                    "and would be slower."
                ),
                UserWarning,
            )

        self.sparse_output_ = sp.issparse(y)

        if not self.sparse_output_:
            y = np.asarray(y)
            y = np.atleast_1d(y)

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]

        check_consistent_length(X, y)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        if self._strategy == "constant":
            if self.constant is None:
                raise ValueError(
                    "Constant target value has to be specified "
                    "when the constant strategy is used."
                )
            else:
                constant = np.reshape(np.atleast_1d(self.constant), (-1, 1))
                if constant.shape[0] != self.n_outputs_:
                    raise ValueError(
                        "Constant target value should have shape (%d, 1)."
                        % self.n_outputs_
                    )

        (self.classes_, self.n_classes_, self.class_prior_) = class_distribution(
            y, sample_weight
        )

        if self._strategy == "constant":
            for k in range(self.n_outputs_):
                if not any(constant[k][0] == c for c in self.classes_[k]):
                    # Checking in case of constant strategy if the constant
                    # provided by the user is in y.
                    err_msg = (
                        "The constant target value must be present in "
                        "the training data. You provided constant={}. "
                        "Possible values are: {}.".format(
                            self.constant, self.classes_[k].tolist()
                        )
                    )
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
        X : array-like of shape (n_samples, n_features)
            Test data.

        Returns
        -------
        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
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
        if self._strategy == "stratified":
            proba = self.predict_proba(X)
            if self.n_outputs_ == 1:
                proba = [proba]

        if self.sparse_output_:
            class_prob = None
            if self._strategy in ("most_frequent", "prior"):
                classes_ = [np.array([cp.argmax()]) for cp in class_prior_]

            elif self._strategy == "stratified":
                class_prob = class_prior_

            elif self._strategy == "uniform":
                raise ValueError(
                    "Sparse target prediction is not "
                    "supported with the uniform strategy"
                )

            elif self._strategy == "constant":
                classes_ = [np.array([c]) for c in constant]

            y = _random_choice_csc(n_samples, classes_, class_prob, self.random_state)
        else:
            if self._strategy in ("most_frequent", "prior"):
                y = np.tile(
                    [
                        classes_[k][class_prior_[k].argmax()]
                        for k in range(self.n_outputs_)
                    ],
                    [n_samples, 1],
                )

            elif self._strategy == "stratified":
                y = np.vstack(
                    [
                        classes_[k][proba[k].argmax(axis=1)]
                        for k in range(self.n_outputs_)
                    ]
                ).T

            elif self._strategy == "uniform":
                ret = [
                    classes_[k][rs.randint(n_classes_[k], size=n_samples)]
                    for k in range(self.n_outputs_)
                ]
                y = np.vstack(ret).T

            elif self._strategy == "constant":
                y = np.tile(self.constant, (n_samples, 1))

            if self.n_outputs_ == 1:
                y = np.ravel(y)

        return y

    def predict_proba(self, X):
        """
        Return probability estimates for the test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test data.

        Returns
        -------
        P : ndarray of shape (n_samples, n_classes) or list of such arrays
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
            if self._strategy == "most_frequent":
                ind = class_prior_[k].argmax()
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0
            elif self._strategy == "prior":
                out = np.ones((n_samples, 1)) * class_prior_[k]

            elif self._strategy == "stratified":
                out = rs.multinomial(1, class_prior_[k], size=n_samples)
                out = out.astype(np.float64)

            elif self._strategy == "uniform":
                out = np.ones((n_samples, n_classes_[k]), dtype=np.float64)
                out /= n_classes_[k]

            elif self._strategy == "constant":
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
            Training data.

        Returns
        -------
        P : ndarray of shape (n_samples, n_classes) or list of such arrays
            Returns the log probability of the sample for each class in
            the model, where classes are ordered arithmetically for each
            output.
        """
        proba = self.predict_proba(X)
        if self.n_outputs_ == 1:
            return np.log(proba)
        else:
            return [np.log(p) for p in proba]

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = True
        tags.classifier_tags.poor_score = True
        tags.no_validation = True
        return tags

    def score(self, X, y, sample_weight=None):
        """Return the mean accuracy on the given test data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : None or array-like of shape (n_samples, n_features)
            Test samples. Passing None as test samples gives the same result
            as passing real test samples, since DummyClassifier
            operates independently of the sampled observations.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            True labels for X.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of self.predict(X) w.r.t. y.
        """
        if X is None:
            X = np.zeros(shape=(len(y), 1))
        return super().score(X, y, sample_weight)


class DummyRegressor(MultiOutputMixin, RegressorMixin, BaseEstimator):
    """Regressor that makes predictions using simple rules.

    This regressor is useful as a simple baseline to compare with other
    (real) regressors. Do not use it for real problems.

    Read more in the :ref:`User Guide <dummy_estimators>`.

    .. versionadded:: 0.13

    Parameters
    ----------
    strategy : {"mean", "median", "quantile", "constant"}, default="mean"
        Strategy to use to generate predictions.

        * "mean": always predicts the mean of the training set
        * "median": always predicts the median of the training set
        * "quantile": always predicts a specified quantile of the training set,
          provided with the quantile parameter.
        * "constant": always predicts a constant value that is provided by
          the user.

    constant : int or float or array-like of shape (n_outputs,), default=None
        The explicit constant as predicted by the "constant" strategy. This
        parameter is useful only for the "constant" strategy.

    quantile : float in [0.0, 1.0], default=None
        The quantile to predict using the "quantile" strategy. A quantile of
        0.5 corresponds to the median, while 0.0 to the minimum and 1.0 to the
        maximum.

    Attributes
    ----------
    constant_ : ndarray of shape (1, n_outputs)
        Mean or median or quantile of the training targets or constant value
        given by the user.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X` has
        feature names that are all strings.

    n_outputs_ : int
        Number of outputs.

    See Also
    --------
    DummyClassifier: Classifier that makes predictions using simple rules.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.dummy import DummyRegressor
    >>> X = np.array([1.0, 2.0, 3.0, 4.0])
    >>> y = np.array([2.0, 3.0, 5.0, 10.0])
    >>> dummy_regr = DummyRegressor(strategy="mean")
    >>> dummy_regr.fit(X, y)
    DummyRegressor()
    >>> dummy_regr.predict(X)
    array([5., 5., 5., 5.])
    >>> dummy_regr.score(X, y)
    0.0
    """

    _parameter_constraints: dict = {
        "strategy": [StrOptions({"mean", "median", "quantile", "constant"})],
        "quantile": [Interval(Real, 0.0, 1.0, closed="both"), None],
        "constant": [
            Interval(Real, None, None, closed="neither"),
            "array-like",
            None,
        ],
    }

    def __init__(self, *, strategy="mean", constant=None, quantile=None):
        self.strategy = strategy
        self.constant = constant
        self.quantile = quantile

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """Fit the baseline regressor.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        validate_data(self, X, skip_check_array=True)

        y = check_array(y, ensure_2d=False, input_name="y")
        if len(y) == 0:
            raise ValueError("y must not be empty.")

        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))
        self.n_outputs_ = y.shape[1]

        check_consistent_length(X, y, sample_weight)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        if self.strategy == "mean":
            self.constant_ = np.average(y, axis=0, weights=sample_weight)

        elif self.strategy == "median":
            if sample_weight is None:
                self.constant_ = np.median(y, axis=0)
            else:
                self.constant_ = [
                    _weighted_percentile(y[:, k], sample_weight, percentile_rank=50.0)
                    for k in range(self.n_outputs_)
                ]

        elif self.strategy == "quantile":
            if self.quantile is None:
                raise ValueError(
                    "When using `strategy='quantile', you have to specify the desired "
                    "quantile in the range [0, 1]."
                )
            percentile_rank = self.quantile * 100.0
            if sample_weight is None:
                self.constant_ = np.percentile(y, axis=0, q=percentile_rank)
            else:
                self.constant_ = [
                    _weighted_percentile(
                        y[:, k], sample_weight, percentile_rank=percentile_rank
                    )
                    for k in range(self.n_outputs_)
                ]

        elif self.strategy == "constant":
            if self.constant is None:
                raise TypeError(
                    "Constant target value has to be specified "
                    "when the constant strategy is used."
                )

            self.constant_ = check_array(
                self.constant,
                accept_sparse=["csr", "csc", "coo"],
                ensure_2d=False,
                ensure_min_samples=0,
            )

            if self.n_outputs_ != 1 and self.constant_.shape[0] != y.shape[1]:
                raise ValueError(
                    "Constant target value should have shape (%d, 1)." % y.shape[1]
                )

        self.constant_ = np.reshape(self.constant_, (1, -1))
        return self

    def predict(self, X, return_std=False):
        """Perform classification on test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test data.

        return_std : bool, default=False
            Whether to return the standard deviation of posterior prediction.
            All zeros in this case.

            .. versionadded:: 0.20

        Returns
        -------
        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Predicted target values for X.

        y_std : array-like of shape (n_samples,) or (n_samples, n_outputs)
            Standard deviation of predictive distribution of query points.
        """
        check_is_fitted(self)
        n_samples = _num_samples(X)

        y = np.full(
            (n_samples, self.n_outputs_),
            self.constant_,
            dtype=np.array(self.constant_).dtype,
        )
        y_std = np.zeros((n_samples, self.n_outputs_))

        if self.n_outputs_ == 1:
            y = np.ravel(y)
            y_std = np.ravel(y_std)

        return (y, y_std) if return_std else y

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = True
        tags.regressor_tags.poor_score = True
        tags.no_validation = True
        return tags

    def score(self, X, y, sample_weight=None):
        """Return the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as `(1 - u/v)`, where `u` is the
        residual sum of squares `((y_true - y_pred) ** 2).sum()` and `v` is the
        total sum of squares `((y_true - y_true.mean()) ** 2).sum()`. The best
        possible score is 1.0 and it can be negative (because the model can be
        arbitrarily worse). A constant model that always predicts the expected
        value of y, disregarding the input features, would get a R^2 score of
        0.0.

        Parameters
        ----------
        X : None or array-like of shape (n_samples, n_features)
            Test samples. Passing None as test samples gives the same result
            as passing real test samples, since `DummyRegressor`
            operates independently of the sampled observations.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        score : float
            R^2 of `self.predict(X)` w.r.t. y.
        """
        if X is None:
            X = np.zeros(shape=(len(y), 1))
        return super().score(X, y, sample_weight)
