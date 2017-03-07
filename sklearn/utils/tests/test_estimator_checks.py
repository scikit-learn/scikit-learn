import scipy.sparse as sp
import numpy as np
import sys
from sklearn.externals.six.moves import cStringIO as StringIO

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.testing import assert_raises_regex, assert_true
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.estimator_checks import check_estimators_unfitted
from sklearn.utils.estimator_checks import check_no_fit_attributes_set_in_init
from sklearn.ensemble import AdaBoostClassifier
from sklearn.linear_model import MultiTaskElasticNet
from sklearn.utils.validation import check_X_y, check_array


class CorrectNotFittedError(ValueError):
    """Exception class to raise if estimator is used before fitting.

    Like NotFittedError, it inherits from ValueError, but not from
    AttributeError. Used for testing only.
    """


class BaseBadClassifier(BaseEstimator, ClassifierMixin):
    def fit(self, X, y):
        return self

    def predict(self, X):
        return np.ones(X.shape[0])


class ChangesDict(BaseEstimator):
    def __init__(self):
        self.key = 0

    def fit(self, X, y=None):
        X, y = check_X_y(X, y)
        return self

    def predict(self, X):
        X = check_array(X)
        self.key = 1000
        return np.ones(X.shape[0])


class SetsWrongAttribute(BaseEstimator):
    def __init__(self):
        self.acceptable_key = 0

    def fit(self, X, y=None):
        self.wrong_attribute = 0
        X, y = check_X_y(X, y)
        return self


class ChangesWrongAttribute(BaseEstimator):
    def __init__(self):
        self.wrong_attribute = 0

    def fit(self, X, y=None):
        self.wrong_attribute = 1
        X, y = check_X_y(X, y)
        return self


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


class CorrectNotFittedErrorClassifier(BaseBadClassifier):
    def fit(self, X, y):
        X, y = check_X_y(X, y)
        self.coef_ = np.ones(X.shape[1])
        return self

    def predict(self, X):
        if not hasattr(self, 'coef_'):
            raise CorrectNotFittedError("estimator is not fitted yet")
        X = check_array(X)
        return np.ones(X.shape[0])


class NoSampleWeightPandasSeriesType(BaseEstimator):
    def fit(self, X, y, sample_weight=None):
        # Convert data
        X, y = check_X_y(X, y,
                         accept_sparse=("csr", "csc"),
                         multi_output=True,
                         y_numeric=True)
        # Function is only called after we verify that pandas is installed
        from pandas import Series
        if isinstance(sample_weight, Series):
            raise ValueError("Estimator does not accept 'sample_weight'"
                             "of type pandas.Series")
        return self

    def predict(self, X):
        X = check_array(X)
        return np.ones(X.shape[0])


def test_check_estimator():
    # tests that the estimator actually fails on "bad" estimators.
    # not a complete test of all checks, which are very extensive.

    # check that we have a set_params and can clone
    msg = "it does not implement a 'get_params' methods"
    assert_raises_regex(TypeError, msg, check_estimator, object)
    # check that we have a fit method
    msg = "object has no attribute 'fit'"
    assert_raises_regex(AttributeError, msg, check_estimator, BaseEstimator)
    # check that fit does input validation
    msg = "TypeError not raised"
    assert_raises_regex(AssertionError, msg, check_estimator, BaseBadClassifier)
    # check that sample_weights in fit accepts pandas.Series type
    try:
        from pandas import Series  # noqa
        msg = ("Estimator NoSampleWeightPandasSeriesType raises error if "
               "'sample_weight' parameter is of type pandas.Series")
        assert_raises_regex(
            ValueError, msg, check_estimator, NoSampleWeightPandasSeriesType)
    except ImportError:
        pass
    # check that predict does input validation (doesn't accept dicts in input)
    msg = "Estimator doesn't check for NaN and inf in predict"
    assert_raises_regex(AssertionError, msg, check_estimator, NoCheckinPredict)
    # check that estimator state does not change
    # at transform/predict/predict_proba time
    msg = 'Estimator changes __dict__ during predict'
    assert_raises_regex(AssertionError, msg, check_estimator, ChangesDict)
    # check that `fit` only changes attribures that
    # are private (start with an _ or end with a _).
    msg = ('Estimator changes public attribute\(s\) during the fit method.'
           ' Estimators are only allowed to change attributes started'
           ' or ended with _, but wrong_attribute changed')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, ChangesWrongAttribute)
    # check that `fit` doesn't add any public attribute
    msg = ('Estimator adds public attribute\(s\) during the fit method.'
           ' Estimators are only allowed to add private attributes'
           ' either started with _ or ended'
           ' with _ but wrong_attribute added')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, SetsWrongAttribute)
    # check for sparse matrix input handling
    name = NoSparseClassifier.__name__
    msg = "Estimator " + name + " doesn't seem to fail gracefully on sparse data"
    # the check for sparse input handling prints to the stdout,
    # instead of raising an error, so as not to remove the original traceback.
    # that means we need to jump through some hoops to catch it.
    old_stdout = sys.stdout
    string_buffer = StringIO()
    sys.stdout = string_buffer
    try:
        check_estimator(NoSparseClassifier)
    except:
        pass
    finally:
        sys.stdout = old_stdout
    assert_true(msg in string_buffer.getvalue())

    # doesn't error on actual estimator
    check_estimator(AdaBoostClassifier)
    check_estimator(MultiTaskElasticNet)


def test_check_estimators_unfitted():
    # check that a ValueError/AttributeError is raised when calling predict
    # on an unfitted estimator
    msg = "AttributeError or ValueError not raised by predict"
    assert_raises_regex(AssertionError, msg, check_estimators_unfitted,
                        "estimator", NoSparseClassifier)

    # check that CorrectNotFittedError inherit from either ValueError
    # or AttributeError
    check_estimators_unfitted("estimator", CorrectNotFittedErrorClassifier)


def test_check_no_fit_attributes_set_in_init():
    class NonConformantEstimator(object):
        def __init__(self):
            self.you_should_not_set_this_ = None

    msg = ("By convention, attributes ending with '_'.+"
           'should not be initialized in the constructor.+'
           "Attribute 'you_should_not_set_this_' was found.+"
           'in estimator estimator_name')

    assert_raises_regex(AssertionError, msg,
                        check_no_fit_attributes_set_in_init,
                        'estimator_name',
                        NonConformantEstimator)
