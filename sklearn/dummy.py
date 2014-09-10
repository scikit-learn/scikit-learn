# Author: Mathieu Blondel <mathieu@mblondel.org>
#         Arnaud Joly <a.joly@ulg.ac.be>
#         Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>
# License: BSD 3 clause
from __future__ import division

import warnings
import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator, ClassifierMixin, RegressorMixin
from .externals.six.moves import xrange
from .utils import check_random_state
from .utils.validation import check_array
from sklearn.utils import deprecated
from sklearn.utils.random import random_choice_csc
from sklearn.utils.multiclass import class_distribution


class DummyClassifier(BaseEstimator, ClassifierMixin):
    """
    DummyClassifier is a classifier that makes predictions using simple rules.

    This classifier is useful as a simple baseline to compare with other
    (real) classifiers. Do not use it for real problems.

    Parameters
    ----------
    strategy : str
        Strategy to use to generate predictions.

        * "stratified": generates predictions by respecting the training
          set's class distribution.
        * "most_frequent": always predicts the most frequent label in the
          training set.
        * "uniform": generates predictions uniformly at random.
        * "constant": always predicts a constant label that is provided by
          the user. This is useful for metrics that evaluate a non-majority
          class

    random_state : int seed, RandomState instance, or None (default)
        The seed of the pseudo random number generator to use.

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

    outputs_2d_ : bool,
        True if the output at fit is 2d, else false.

    `sparse_output_` : bool,
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
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        sample_weight : array-like of shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        self : object
            Returns self.
        """
        if self.strategy not in ("most_frequent", "stratified", "uniform",
                                 "constant"):
            raise ValueError("Unknown strategy type.")

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

        self.output_2d_ = y.ndim == 2
        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        self.n_outputs_ = y.shape[1]

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
                # Checking in case of constant strategy if the constant
                # provided by the user is in y.
                if constant[k] not in self.classes_[k]:
                    raise ValueError("The constant target value must be "
                                     "present in training data")

        if self.n_outputs_ == 1 and not self.output_2d_:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]
            self.class_prior_ = self.class_prior_[0]

        return self

    def predict(self, X):
        """
        Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array, shape = [n_samples] or [n_samples, n_outputs]
            Predicted target values for X.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])
        # numpy random_state expects Python int and not long as size argument
        # under Windows
        n_samples = int(X.shape[0])
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
            if self.strategy == "most_frequent":
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
            if self.strategy == "most_frequent":
                y = np.tile([classes_[k][class_prior_[k].argmax()] for
                             k in range(self.n_outputs_)], [n_samples, 1])

            elif self.strategy == "stratified":
                y = np.vstack(classes_[k][proba[k].argmax(axis=1)] for
                              k in range(self.n_outputs_)).T

            elif self.strategy == "uniform":
                ret = [classes_[k][rs.randint(n_classes_[k], size=n_samples)]
                       for k in range(self.n_outputs_)]
                y = np.vstack(ret).T

            elif self.strategy == "constant":
                y = np.tile(self.constant, (n_samples, 1))

            if self.n_outputs_ == 1 and not self.output_2d_:
                y = np.ravel(y)

        return y

    def predict_proba(self, X):
        """
        Return probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        P : array-like or list of array-lke of shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered arithmetically, for each
            output.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("DummyClassifier not fitted.")

        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])
        # numpy random_state expects Python int and not long as size argument
        # under Windows
        n_samples = int(X.shape[0])
        rs = check_random_state(self.random_state)

        n_classes_ = self.n_classes_
        classes_ = self.classes_
        class_prior_ = self.class_prior_
        constant = self.constant
        if self.n_outputs_ == 1 and not self.output_2d_:
            # Get same type even for self.n_outputs_ == 1
            n_classes_ = [n_classes_]
            classes_ = [classes_]
            class_prior_ = [class_prior_]
            constant = [constant]

        P = []
        for k in xrange(self.n_outputs_):
            if self.strategy == "most_frequent":
                ind = np.ones(n_samples, dtype=int) * class_prior_[k].argmax()
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0

            elif self.strategy == "stratified":
                out = rs.multinomial(1, class_prior_[k], size=n_samples)

            elif self.strategy == "uniform":
                out = np.ones((n_samples, n_classes_[k]), dtype=np.float64)
                out /= n_classes_[k]

            elif self.strategy == "constant":
                ind = np.where(classes_[k] == constant[k])
                out = np.zeros((n_samples, n_classes_[k]), dtype=np.float64)
                out[:, ind] = 1.0

            P.append(out)

        if self.n_outputs_ == 1 and not self.output_2d_:
            P = P[0]

        return P

    def predict_log_proba(self, X):
        """
        Return log probability estimates for the test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

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


class DummyRegressor(BaseEstimator, RegressorMixin):
    """
    DummyRegressor is a regressor that makes predictions using
    simple rules.

    This regressor is useful as a simple baseline to compare with other
    (real) regressors. Do not use it for real problems.

    Parameters
    ----------
    strategy : str
        Strategy to use to generate predictions.

        * "mean": always predicts the mean of the training set
        * "median": always predicts the median of the training set
        * "constant": always predicts a constant value that is provided by
          the user.

    constant : int or float or array of shape = [n_outputs]
        The explicit constant as predicted by the "constant" strategy. This
        parameter is useful only for the "constant" strategy.

    Attributes
    ----------
    constant_ : float or array of shape [n_outputs]
        Mean or median of the training targets or constant value given the by
        the user.

    n_outputs_ : int,
        Number of outputs.

    outputs_2d_ : bool,
        True if the output at fit is 2d, else false.
    """

    def __init__(self, strategy="mean", constant=None):
        self.strategy = strategy
        self.constant = constant

    @property
    @deprecated('This will be removed in version 0.17')
    def y_mean_(self):
        if self.strategy == 'mean':
            return self.constant_
        raise AttributeError

    def fit(self, X, y):
        """Fit the random regressor.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            Target values.

        Returns
        -------
        self : object
            Returns self.
        """

        if self.strategy not in ("mean", "median", "constant"):
            raise ValueError("Unknown strategy type: %s, "
                             "expected 'mean', 'median' or 'constant'"
                             % self.strategy)

        y = check_array(y, accept_sparse='csr', ensure_2d=False)
        self.output_2d_ = (y.ndim == 2)

        if self.strategy == "mean":
            self.constant_ = np.reshape(np.mean(y, axis=0), (1, -1))

        elif self.strategy == "median":
            self.constant_ = np.reshape(np.median(y, axis=0), (1, -1))

        elif self.strategy == "constant":
            if self.constant is None:
                raise TypeError("Constant target value has to be specified "
                                "when the constant strategy is used.")

            self.constant = check_array(self.constant,
                                        accept_sparse=['csr', 'csc', 'coo'],
                                        ensure_2d=False)

            if self.output_2d_ and self.constant.shape[0] != y.shape[1]:
                raise ValueError(
                    "Constant target value should have "
                    "shape (%d, 1)." % y.shape[1])

            self.constant_ = np.reshape(self.constant, (1, -1))

        self.n_outputs_ = np.size(self.constant_)  # y.shape[1] is not safe
        return self

    def predict(self, X):
        """
        Perform classification on test vectors X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Input vectors, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array, shape = [n_samples]  or [n_samples, n_outputs]
            Predicted target values for X.
        """
        if not hasattr(self, "constant_"):
            raise ValueError("DummyRegressor not fitted.")

        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])
        n_samples = X.shape[0]

        y = np.ones((n_samples, 1)) * self.constant_

        if self.n_outputs_ == 1 and not self.output_2d_:
            y = np.ravel(y)

        return y
