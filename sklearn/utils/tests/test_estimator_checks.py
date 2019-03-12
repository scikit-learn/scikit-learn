import unittest
import sys

import numpy as np
import scipy.sparse as sp

from io import StringIO

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import deprecated
from sklearn.utils import _joblib
from sklearn.utils.testing import (assert_raises_regex,
                                   assert_equal, ignore_warnings,
                                   assert_warns, assert_raises)
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.estimator_checks \
    import check_class_weight_balanced_linear_classifier
from sklearn.utils.estimator_checks import set_random_state
from sklearn.utils.estimator_checks import set_checking_parameters
from sklearn.utils.estimator_checks import check_estimators_unfitted
from sklearn.utils.estimator_checks import check_fit_score_takes_y
from sklearn.utils.estimator_checks import check_no_attributes_set_in_init
from sklearn.utils.estimator_checks import check_outlier_corruption
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.linear_model import LinearRegression, SGDClassifier
from sklearn.mixture import GaussianMixture
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import NMF
from sklearn.linear_model import MultiTaskElasticNet
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsRegressor
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
    def __init__(self, key=0):
        self.key = key

    def fit(self, X, y=None):
        X, y = check_X_y(X, y)
        return self

    def predict(self, X):
        X = check_array(X)
        self.key = 1000
        return np.ones(X.shape[0])


class SetsWrongAttribute(BaseEstimator):
    def __init__(self, acceptable_key=0):
        self.acceptable_key = acceptable_key

    def fit(self, X, y=None):
        self.wrong_attribute = 0
        X, y = check_X_y(X, y)
        return self


class ChangesWrongAttribute(BaseEstimator):
    def __init__(self, wrong_attribute=0):
        self.wrong_attribute = wrong_attribute

    def fit(self, X, y=None):
        self.wrong_attribute = 1
        X, y = check_X_y(X, y)
        return self


class ChangesUnderscoreAttribute(BaseEstimator):
    def fit(self, X, y=None):
        self._good_attribute = 1
        X, y = check_X_y(X, y)
        return self


class RaisesErrorInSetParams(BaseEstimator):
    def __init__(self, p=0):
        self.p = p

    def set_params(self, **kwargs):
        if 'p' in kwargs:
            p = kwargs.pop('p')
            if p < 0:
                raise ValueError("p can't be less than 0")
            self.p = p
        return super().set_params(**kwargs)

    def fit(self, X, y=None):
        X, y = check_X_y(X, y)
        return self


class ModifiesValueInsteadOfRaisingError(BaseEstimator):
    def __init__(self, p=0):
        self.p = p

    def set_params(self, **kwargs):
        if 'p' in kwargs:
            p = kwargs.pop('p')
            if p < 0:
                p = 0
            self.p = p
        return super().set_params(**kwargs)

    def fit(self, X, y=None):
        X, y = check_X_y(X, y)
        return self


class ModifiesAnotherValue(BaseEstimator):
    def __init__(self, a=0, b='method1'):
        self.a = a
        self.b = b

    def set_params(self, **kwargs):
        if 'a' in kwargs:
            a = kwargs.pop('a')
            self.a = a
            if a is None:
                kwargs.pop('b')
                self.b = 'method2'
        return super().set_params(**kwargs)

    def fit(self, X, y=None):
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


class BadBalancedWeightsClassifier(BaseBadClassifier):
    def __init__(self, class_weight=None):
        self.class_weight = class_weight

    def fit(self, X, y):
        from sklearn.preprocessing import LabelEncoder
        from sklearn.utils import compute_class_weight

        label_encoder = LabelEncoder().fit(y)
        classes = label_encoder.classes_
        class_weight = compute_class_weight(self.class_weight, classes, y)

        # Intentionally modify the balanced class_weight
        # to simulate a bug and raise an exception
        if self.class_weight == "balanced":
            class_weight += 1.

        # Simply assigning coef_ to the class_weight
        self.coef_ = class_weight
        return self


class BadTransformerWithoutMixin(BaseEstimator):
    def fit(self, X, y=None):
        X = check_array(X)
        return self

    def transform(self, X):
        X = check_array(X)
        return X


class NotInvariantPredict(BaseEstimator):
    def fit(self, X, y):
        # Convert data
        X, y = check_X_y(X, y,
                         accept_sparse=("csr", "csc"),
                         multi_output=True,
                         y_numeric=True)
        return self

    def predict(self, X):
        # return 1 if X has more than one element else return 0
        X = check_array(X)
        if X.shape[0] > 1:
            return np.ones(X.shape[0])
        return np.zeros(X.shape[0])


class LargeSparseNotSupportedClassifier(BaseEstimator):
    def fit(self, X, y):
        X, y = check_X_y(X, y,
                         accept_sparse=("csr", "csc", "coo"),
                         accept_large_sparse=True,
                         multi_output=True,
                         y_numeric=True)
        if sp.issparse(X):
            if X.getformat() == "coo":
                if X.row.dtype == "int64" or X.col.dtype == "int64":
                    raise ValueError(
                        "Estimator doesn't support 64-bit indices")
            elif X.getformat() in ["csc", "csr"]:
                if X.indices.dtype == "int64" or X.indptr.dtype == "int64":
                    raise ValueError(
                        "Estimator doesn't support 64-bit indices")

        return self


class SparseTransformer(BaseEstimator):
    def fit(self, X, y=None):
        self.X_shape_ = check_array(X).shape
        return self

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X)

    def transform(self, X):
        X = check_array(X)
        if X.shape[1] != self.X_shape_[1]:
            raise ValueError('Bad number of features')
        return sp.csr_matrix(X)


def test_check_fit_score_takes_y_works_on_deprecated_fit():
    # Tests that check_fit_score_takes_y works on a class with
    # a deprecated fit method

    class TestEstimatorWithDeprecatedFitMethod(BaseEstimator):
        @deprecated("Deprecated for the purpose of testing "
                    "check_fit_score_takes_y")
        def fit(self, X, y):
            return self

    check_fit_score_takes_y("test", TestEstimatorWithDeprecatedFitMethod())


def test_check_estimator():
    # tests that the estimator actually fails on "bad" estimators.
    # not a complete test of all checks, which are very extensive.

    # check that we have a set_params and can clone
    msg = "it does not implement a 'get_params' methods"
    assert_raises_regex(TypeError, msg, check_estimator, object)
    assert_raises_regex(TypeError, msg, check_estimator, object())
    # check that values returned by get_params match set_params
    msg = "get_params result does not match what was passed to set_params"
    assert_raises_regex(AssertionError, msg, check_estimator,
                        ModifiesValueInsteadOfRaisingError())
    assert_warns(UserWarning, check_estimator, RaisesErrorInSetParams())
    assert_raises_regex(AssertionError, msg, check_estimator,
                        ModifiesAnotherValue())
    # check that we have a fit method
    msg = "object has no attribute 'fit'"
    assert_raises_regex(AttributeError, msg, check_estimator, BaseEstimator)
    assert_raises_regex(AttributeError, msg, check_estimator, BaseEstimator())
    # check that fit does input validation
    msg = "ValueError not raised"
    assert_raises_regex(AssertionError, msg, check_estimator,
                        BaseBadClassifier)
    assert_raises_regex(AssertionError, msg, check_estimator,
                        BaseBadClassifier())
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
    assert_raises_regex(AssertionError, msg, check_estimator,
                        NoCheckinPredict())
    # check that estimator state does not change
    # at transform/predict/predict_proba time
    msg = 'Estimator changes __dict__ during predict'
    assert_raises_regex(AssertionError, msg, check_estimator, ChangesDict)
    # check that `fit` only changes attribures that
    # are private (start with an _ or end with a _).
    msg = ('Estimator ChangesWrongAttribute should not change or mutate  '
           'the parameter wrong_attribute from 0 to 1 during fit.')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, ChangesWrongAttribute)
    check_estimator(ChangesUnderscoreAttribute)
    # check that `fit` doesn't add any public attribute
    msg = (r'Estimator adds public attribute\(s\) during the fit method.'
           ' Estimators are only allowed to add private attributes'
           ' either started with _ or ended'
           ' with _ but wrong_attribute added')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, SetsWrongAttribute)
    # check for invariant method
    name = NotInvariantPredict.__name__
    method = 'predict'
    msg = ("{method} of {name} is not invariant when applied "
           "to a subset.").format(method=method, name=name)
    assert_raises_regex(AssertionError, msg,
                        check_estimator, NotInvariantPredict)
    # check for sparse matrix input handling
    name = NoSparseClassifier.__name__
    msg = "Estimator %s doesn't seem to fail gracefully on sparse data" % name
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
    assert msg in string_buffer.getvalue()

    # Large indices test on bad estimator
    msg = ('Estimator LargeSparseNotSupportedClassifier doesn\'t seem to '
           r'support \S{3}_64 matrix, and is not failing gracefully.*')
    assert_raises_regex(AssertionError, msg, check_estimator,
                        LargeSparseNotSupportedClassifier)

    # non-regression test for estimators transforming to sparse data
    check_estimator(SparseTransformer())

    # doesn't error on actual estimator
    check_estimator(AdaBoostClassifier)
    check_estimator(AdaBoostClassifier())
    check_estimator(MultiTaskElasticNet)
    check_estimator(MultiTaskElasticNet())


def test_check_outlier_corruption():
    # should raise AssertionError
    decision = np.array([0., 1., 1.5, 2.])
    assert_raises(AssertionError, check_outlier_corruption, 1, 2, decision)
    # should pass
    decision = np.array([0., 1., 1., 2.])
    check_outlier_corruption(1, 2, decision)


def test_check_estimator_transformer_no_mixin():
    # check that TransformerMixin is not required for transformer tests to run
    assert_raises_regex(AttributeError, '.*fit_transform.*',
                        check_estimator, BadTransformerWithoutMixin())


def test_check_estimator_clones():
    # check that check_estimator doesn't modify the estimator it receives
    from sklearn.datasets import load_iris
    iris = load_iris()

    for Estimator in [GaussianMixture, LinearRegression,
                      RandomForestClassifier, NMF, SGDClassifier,
                      MiniBatchKMeans]:
        with ignore_warnings(category=(FutureWarning, DeprecationWarning)):
            # when 'est = SGDClassifier()'
            est = Estimator()
            set_checking_parameters(est)
            set_random_state(est)
            # without fitting
            old_hash = _joblib.hash(est)
            check_estimator(est)
        assert_equal(old_hash, _joblib.hash(est))

        with ignore_warnings(category=(FutureWarning, DeprecationWarning)):
            # when 'est = SGDClassifier()'
            est = Estimator()
            set_checking_parameters(est)
            set_random_state(est)
            # with fitting
            est.fit(iris.data + 10, iris.target)
            old_hash = _joblib.hash(est)
            check_estimator(est)
        assert_equal(old_hash, _joblib.hash(est))


def test_check_estimators_unfitted():
    # check that a ValueError/AttributeError is raised when calling predict
    # on an unfitted estimator
    msg = "AttributeError or ValueError not raised by predict"
    assert_raises_regex(AssertionError, msg, check_estimators_unfitted,
                        "estimator", NoSparseClassifier())

    # check that CorrectNotFittedError inherit from either ValueError
    # or AttributeError
    check_estimators_unfitted("estimator", CorrectNotFittedErrorClassifier())


def test_check_no_attributes_set_in_init():
    class NonConformantEstimatorPrivateSet:
        def __init__(self):
            self.you_should_not_set_this_ = None

    class NonConformantEstimatorNoParamSet:
        def __init__(self, you_should_set_this_=None):
            pass

    assert_raises_regex(AssertionError,
                        "Estimator estimator_name should not set any"
                        " attribute apart from parameters during init."
                        r" Found attributes \['you_should_not_set_this_'\].",
                        check_no_attributes_set_in_init,
                        'estimator_name',
                        NonConformantEstimatorPrivateSet())
    assert_raises_regex(AssertionError,
                        "Estimator estimator_name should store all "
                        "parameters as an attribute during init. "
                        "Did not find attributes "
                        r"\['you_should_set_this_'\].",
                        check_no_attributes_set_in_init,
                        'estimator_name',
                        NonConformantEstimatorNoParamSet())


def test_check_estimator_pairwise():
    # check that check_estimator() works on estimator with _pairwise
    # kernel or  metric

    # test precomputed kernel
    est = SVC(kernel='precomputed')
    check_estimator(est)

    # test precomputed metric
    est = KNeighborsRegressor(metric='precomputed')
    check_estimator(est)


def run_tests_without_pytest():
    """Runs the tests in this file without using pytest.
    """
    main_module = sys.modules['__main__']
    test_functions = [getattr(main_module, name) for name in dir(main_module)
                      if name.startswith('test_')]
    test_cases = [unittest.FunctionTestCase(fn) for fn in test_functions]
    suite = unittest.TestSuite()
    suite.addTests(test_cases)
    runner = unittest.TextTestRunner()
    runner.run(suite)


def test_check_class_weight_balanced_linear_classifier():
    # check that ill-computed balanced weights raises an exception
    assert_raises_regex(AssertionError,
                        "Classifier estimator_name is not computing"
                        " class_weight=balanced properly.",
                        check_class_weight_balanced_linear_classifier,
                        'estimator_name',
                        BadBalancedWeightsClassifier)


if __name__ == '__main__':
    # This module is run as a script to check that we have no dependency on
    # pytest for estimator checks.
    run_tests_without_pytest()
