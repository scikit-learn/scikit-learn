import numpy as np

from ..base import BaseEstimator, ClassifierMixin
from .testing import assert_true
from .validation import _num_samples, check_array


class ArraySlicingWrapper(object):
    """
    Parameters
    ----------
    array
    """
    def __init__(self, array):
        self.array = array

    def __getitem__(self, aslice):
        return MockDataFrame(self.array[aslice])


class MockDataFrame(object):
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


class CheckingClassifier(BaseEstimator, ClassifierMixin):
    """Dummy classifier to test pipelining and meta-estimators.

    Checks some property of X and y in fit / predict.
    This allows testing whether pipelines / cross-validation or metaestimators
    changed the input.

    Parameters
    ----------
    check_y
    check_X
    foo_param
    expected_fit_params
    """
    def __init__(self, check_y=None, check_X=None, foo_param=0,
                 expected_fit_params=None):
        self.check_y = check_y
        self.check_X = check_X
        self.foo_param = foo_param
        self.expected_fit_params = expected_fit_params

    def fit(self, X, y, **fit_params):
        """
        Fit classifier

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of the estimator
        """
        assert len(X) == len(y)
        if self.check_X is not None:
            assert self.check_X(X)
        if self.check_y is not None:
            assert self.check_y(y)
        self.classes_ = np.unique(check_array(y, ensure_2d=False,
                                              allow_nd=True))
        if self.expected_fit_params:
            missing = set(self.expected_fit_params) - set(fit_params)
            assert_true(len(missing) == 0, 'Expected fit parameter(s) %s not '
                                           'seen.' % list(missing))
            for key, value in fit_params.items():
                assert_true(len(value) == len(X),
                            'Fit parameter %s has length %d; '
                            'expected %d.' % (key, len(value), len(X)))

        return self

    def predict(self, T):
        """
        Parameters
        -----------
        T : indexable, length n_samples
        """
        if self.check_X is not None:
            assert self.check_X(T)
        return self.classes_[np.zeros(_num_samples(T), dtype=np.int)]

    def score(self, X=None, Y=None):
        """
        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Input data, where n_samples is the number of samples and
            n_features is the number of features.

        Y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.
        """
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score
