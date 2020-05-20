import numpy as np

from ..base import BaseEstimator, ClassifierMixin
from .validation import _num_samples, check_array


class ArraySlicingWrapper:
    """
    Parameters
    ----------
    array
    """
    def __init__(self, array):
        self.array = array

    def __getitem__(self, aslice):
        return MockDataFrame(self.array[aslice])


class MockDataFrame:
    """
    Parameters
    ----------
    array
    """
    # have shape and length but don't support indexing.
    def __init__(self, array):
        self.array = array
        self.values = array
        self.shape = array.shape
        self.ndim = array.ndim
        # ugly hack to make iloc work.
        self.iloc = ArraySlicingWrapper(array)

    def __len__(self):
        return len(self.array)

    def __array__(self, dtype=None):
        # Pandas data frames also are array-like: we want to make sure that
        # input validation in cross-validation does not try to call that
        # method.
        return self.array

    def __eq__(self, other):
        return MockDataFrame(self.array == other.array)

    def __ne__(self, other):
        return not self == other


class CheckingClassifier(ClassifierMixin, BaseEstimator):
    """Dummy classifier to test pipelining and meta-estimators.

    Checks some property of `X` and `y`in fit / predict.
    This allows testing whether pipelines / cross-validation or metaestimators
    changed the input.

    Parameters
    ----------
    check_y, check_X : callable, default=None
        The callable used to validate `X` and `y`. These callable should return
        a bool where `False` will trigger an `AssertionError`.

    check_y_params, check_X_params : dict, default=None
        The optional parameters to pass to `check_X` and `check_y`.

    foo_param : int, default=0
        A `foo` param. When `foo > 1`, the output of :meth:`score` will be 1
        otherwise it is 0.

    expected_fit_params : list of str, default=None
        A list of the expected parameters given when calling `fit`.

    Attributes
    ----------
    classes_ : int
        The classes seen during `fit`.

    n_features_in_ : int
        The number of features seen during `fit`.
    """
    def __init__(self, *, check_y=None, check_y_params=None,
                 check_X=None, check_X_params=None, foo_param=0,
                 expected_fit_params=None):
        self.check_y = check_y
        self.check_y_params = check_y_params
        self.check_X = check_X
        self.check_X_params = check_X_params
        self.foo_param = foo_param
        self.expected_fit_params = expected_fit_params

    def fit(self, X, y, **fit_params):
        """Fit classifier.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like of shape (n_samples, n_output) or (n_samples,), optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of the estimator

        Returns
        -------
        self
        """
        assert _num_samples(X) == _num_samples(y)
        if self.check_X is not None:
            params = {} if self.check_X_params is None else self.check_X_params
            assert self.check_X(X, **params)
        if self.check_y is not None:
            params = {} if self.check_y_params is None else self.check_y_params
            assert self.check_y(y)
        self.n_features_in_ = np.shape(X)[1]
        self.classes_ = np.unique(
            check_array(y, ensure_2d=False, allow_nd=True)
        )
        if self.expected_fit_params:
            missing = set(self.expected_fit_params) - set(fit_params)
            if missing:
                raise AssertionError(
                    f'Expected fit parameter(s) {list(missing)} not seen.'
                )
            for key, value in fit_params.items():
                if _num_samples(value) != _num_samples(X):
                    raise AssertionError(
                        f'Fit parameter {key} has length {_num_samples(value)}'
                        f'; expected {_num_samples(X)}.'
                    )

        return self

    def predict(self, X):
        """Predict the first class seen in `classes_`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        preds : ndarray of shape (n_samples,)
            Predictions of the first class seens in `classes_`.
        """
        if self.check_X is not None:
            params = {} if self.check_X_params is None else self.check_X_params
            assert self.check_X(X, **params)
        return self.classes_[np.zeros(_num_samples(X), dtype=np.int)]

    def predict_proba(self, X):
        """Predict probabilities for each class.

        Here, the dummy classifier will provide a probability of 1 for the
        first class of `classes_` and 0 otherwise.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        proba : ndarray of shape (n_samples, n_classes)
            The probabilities for each sample and class.
        """
        proba = np.zeros((_num_samples(X), len(self.classes_)))
        proba[:, 0] = 1
        return proba

    def decision_function(self, X):
        """Confidence score.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        decision : ndarray of shape (n_samples,) if n_classes == 2\
                else (n_samples, n_classes)
            Confidence score.
        """
        if len(self.classes_) == 2:
            # for binary classifier, the confidence score is related to
            # classes_[1] and therefore should be null.
            return np.zeros(_num_samples(X))
        else:
            return self.predict_proba(X)

    def score(self, X=None, Y=None):
        """Fake score.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data, where n_samples is the number of samples and
            n_features is the number of features.

        Y : array-like of shape (n_samples, n_output) or (n_samples,)
            Target relative to X for classification or regression;
            None for unsupervised learning.

        Returns
        -------
        score : float
            Either 0 or 1 depending of `foo_param` (i.e. `foo_param > 1 =>
            score=1` otherwise `score=0`).
        """
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score

    def _more_tags(self):
        return {'_skip_test': True, 'X_types': ['1dlabel']}


class NoSampleWeightWrapper(BaseEstimator):
    """Wrap estimator which will not expose `sample_weight`.

    Parameters
    ----------
    est : estimator, default=None
        The estimator to wrap.
    """
    def __init__(self, est=None):
        self.est = est

    def fit(self, X, y):
        return self.est.fit(X, y)

    def predict(self, X):
        return self.est.predict(X)

    def predict_proba(self, X):
        return self.est.predict_proba(X)

    def _more_tags(self):
        return {'_skip_test': True}
