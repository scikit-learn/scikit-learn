import scipy.sparse as sp
import numpy as np

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.estimator_checks import check_estimator
from sklearn.linear_model import LogisticRegression
from sklearn.utils.validation import check_X_y, check_array


class BaseBadClassifier(BaseEstimator, ClassifierMixin):
    def fit(self, X, y):
        return self

    def predict(self, X):
        return np.ones(X.shape[0])


class NoCheckinPredict(BaseBadClassifier):
    def fit(self, X, y):
        X, y = check_X_y(X, y)
        return self


class NoSparseClassifier(BaseBadClassifier):
    def fit(self, X, y):
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc'])
        if sp.issparse(X):
            raise ValueError("Nonsensical Error")
        return self

    def predict(self, X):
        X = check_array(X)
        return np.ones(X.shape[0])


def test_check_estimator():
    # check that we have a set_params and can clone
    msg = "it does not implement a 'get_params' methods"
    assert_raises_regex(TypeError, msg, check_estimator, object)
    msg = "object has no attribute 'fit'"
    assert_raises_regex(AttributeError, msg, check_estimator, BaseEstimator)
    msg = "argument must be a string or a number"
    assert_raises_regex(AssertionError, msg, check_estimator, BaseBadClassifier)
    msg = "Estimator doesn't check for NaN and inf in predict"
    assert_raises_regex(AssertionError, msg, check_estimator, NoCheckinPredict)
    msg = "Estimator type doesn't seem to fail gracefully on sparse data"
    assert_raises_regex(TypeError, msg, check_estimator, NoSparseClassifier)

    # doesn't error on actual estimator
    check_estimator(LogisticRegression)
