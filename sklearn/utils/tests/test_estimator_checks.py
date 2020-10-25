import unittest
import sys

import numpy as np
import scipy.sparse as sp
import joblib

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import deprecated
from sklearn.utils._testing import (assert_raises_regex,
                                    ignore_warnings,
                                    assert_warns, assert_raises,
                                    SkipTest)
from sklearn.utils.estimator_checks import check_estimator, _NotAnArray
from sklearn.utils.estimator_checks \
    import check_class_weight_balanced_linear_classifier
from sklearn.utils.estimator_checks import set_random_state
from sklearn.utils.estimator_checks import _set_checking_parameters
from sklearn.utils.estimator_checks import check_estimators_unfitted
from sklearn.utils.estimator_checks import check_fit_score_takes_y
from sklearn.utils.estimator_checks import check_no_attributes_set_in_init
from sklearn.utils.estimator_checks import check_classifier_data_not_an_array
from sklearn.utils.estimator_checks import check_regressor_data_not_an_array
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.estimator_checks import check_outlier_corruption
from sklearn.utils.fixes import np_version, parse_version
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, SGDClassifier
from sklearn.mixture import GaussianMixture
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import NMF
from sklearn.linear_model import MultiTaskElasticNet, LogisticRegression
from sklearn.svm import SVC, NuSVC
from sklearn.neighbors import KNeighborsRegressor
from sklearn.utils.validation import check_array
from sklearn.utils import all_estimators
from sklearn.exceptions import SkipTestWarning


class CorrectNotFittedError(ValueError):
    """Exception class to raise if estimator is used before fitting.

    Like NotFittedError, it inherits from ValueError, but not from
    AttributeError. Used for testing only.
    """


class BaseBadClassifier(ClassifierMixin, BaseEstimator):
    def fit(self, X, y):
        return self

    def predict(self, X):
        return np.ones(X.shape[0])


class ChangesDict(BaseEstimator):
    def __init__(self, key=0):
        self.key = key

    def fit(self, X, y=None):
        X, y = self._validate_data(X, y)
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
        X, y = self._validate_data(X, y)
        return self


class ChangesWrongAttribute(BaseEstimator):
    def __init__(self, wrong_attribute=0):
        self.wrong_attribute = wrong_attribute

    def fit(self, X, y=None):
        self.wrong_attribute = 1
        X, y = self._validate_data(X, y)
        return self


class ChangesUnderscoreAttribute(BaseEstimator):
    def fit(self, X, y=None):
        self._good_attribute = 1
        X, y = self._validate_data(X, y)
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
        X, y = self._validate_data(X, y)
        return self


class HasMutableParameters(BaseEstimator):
    def __init__(self, p=object()):
        self.p = p

    def fit(self, X, y=None):
        X, y = self._validate_data(X, y)
        return self


class HasImmutableParameters(BaseEstimator):
    # Note that object is an uninitialized class, thus immutable.
    def __init__(self, p=42, q=np.int32(42), r=object):
        self.p = p
        self.q = q
        self.r = r

    def fit(self, X, y=None):
        X, y = self._validate_data(X, y)
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
        X, y = self._validate_data(X, y)
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
        X, y = self._validate_data(X, y)
        return self


class NoCheckinPredict(BaseBadClassifier):
    def fit(self, X, y):
        X, y = self._validate_data(X, y)
        return self


class NoSparseClassifierException(LogisticRegression):
    """An estimator that does not accept sparse inputs but
    returns a plain Exception instead of TypeError or ValueError"""

    def fit(self, X, y, *args, **kwargs):
        X = check_array(X, accept_sparse=True)
        if sp.issparse(X):
            raise Exception("Nonsensical Error")
        return super(NoSparseClassifierException, self).fit(X, y)


class NoSparseClassifierTypeError(LogisticRegression):
    """An estimator that refuses sparse inputs that returns a TypeError
    (as allowed), but with the wrong message"""

    def fit(self, X, y, *args, **kwargs):
        X = check_array(X, accept_sparse=True)
        if sp.issparse(X):
            raise TypeError("Nonsensical Error")
        return super(NoSparseClassifierTypeError, self).fit(X, y)


class NoSparseClassifierValueError(LogisticRegression):
    """An estimator that refuses sparse inputs that returns a ValueError
    (as allowed), but with the wrong message"""

    def fit(self, X, y, *args, **kwargs):
        X = check_array(X, accept_sparse=True)
        if sp.issparse(X):
            raise ValueError("Nonsensical Error")
        return super(NoSparseClassifierValueError, self).fit(X, y)


class SparseDenseDifferentPredict(LogisticRegression):
    """An estimator that returns a different result for sparse inputs and
    dense inputs when using predict"""

    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        if sp.issparse(X):
            return super(SparseDenseDifferentPredict, self).predict(X)
        else:
            return - super(SparseDenseDifferentPredict, self).predict(X)


class SparseDenseDifferentPredictProba(LogisticRegression):
    """An estimator that returns a different result for sparse inputs and
    dense inputs when using predict_proba"""

    def predict_proba(self, X):
        X = check_array(X, accept_sparse=True)
        if sp.issparse(X):
            return super(SparseDenseDifferentPredictProba,
                         self).predict_proba(X)
        else:
            return - super(SparseDenseDifferentPredictProba,
                           self).predict_proba(X)


class CorrectNotFittedErrorClassifier(BaseBadClassifier):
    def fit(self, X, y):
        X, y = self._validate_data(X, y)
        self.coef_ = np.ones(X.shape[1])
        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)
        return np.ones(X.shape[0])


class NoSampleWeightPandasSeriesType(BaseEstimator):
    def fit(self, X, y, sample_weight=None):
        # Convert data
        X, y = self._validate_data(
            X, y,
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
        class_weight = compute_class_weight(self.class_weight, classes=classes,
                                            y=y)

        # Intentionally modify the balanced class_weight
        # to simulate a bug and raise an exception
        if self.class_weight == "balanced":
            class_weight += 1.

        # Simply assigning coef_ to the class_weight
        self.coef_ = class_weight
        return self


class BadTransformerWithoutMixin(BaseEstimator):
    def fit(self, X, y=None):
        X = self._validate_data(X)
        return self

    def transform(self, X):
        X = check_array(X)
        return X


class NotInvariantPredict(BaseEstimator):
    def fit(self, X, y):
        # Convert data
        X, y = self._validate_data(
            X, y,
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


class NotInvariantSampleOrder(BaseEstimator):
    def fit(self, X, y):
        X, y = self._validate_data(
            X, y,
            accept_sparse=("csr", "csc"),
            multi_output=True,
            y_numeric=True)
        # store the original X to check for sample order later
        self._X = X
        return self

    def predict(self, X):
        X = check_array(X)
        # if the input contains the same elements but different sample order,
        # then just return zeros.
        if (np.array_equiv(np.sort(X, axis=0), np.sort(self._X, axis=0)) and
           (X != self._X).any()):
            return np.zeros(X.shape[0])
        return X[:, 0]


class LargeSparseNotSupportedClassifier(BaseEstimator):
    def fit(self, X, y):
        X, y = self._validate_data(
            X, y,
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
                assert "int64" not in (X.indices.dtype, X.indptr.dtype),\
                    "Estimator doesn't support 64-bit indices"

        return self


class SparseTransformer(BaseEstimator):
    def fit(self, X, y=None):
        self.X_shape_ = self._validate_data(X).shape
        return self

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X)

    def transform(self, X):
        X = check_array(X)
        if X.shape[1] != self.X_shape_[1]:
            raise ValueError('Bad number of features')
        return sp.csr_matrix(X)


class EstimatorInconsistentForPandas(BaseEstimator):
    def fit(self, X, y):
        try:
            from pandas import DataFrame
            if isinstance(X, DataFrame):
                self.value_ = X.iloc[0, 0]
            else:
                X = check_array(X)
                self.value_ = X[1, 0]
            return self

        except ImportError:
            X = check_array(X)
            self.value_ = X[1, 0]
            return self

    def predict(self, X):
        X = check_array(X)
        return np.array([self.value_] * X.shape[0])


class UntaggedBinaryClassifier(SGDClassifier):
    # Toy classifier that only supports binary classification, will fail tests.
    def fit(self, X, y, coef_init=None, intercept_init=None,
            sample_weight=None):
        super().fit(X, y, coef_init, intercept_init, sample_weight)
        if len(self.classes_) > 2:
            raise ValueError('Only 2 classes are supported')
        return self

    def partial_fit(self, X, y, classes=None, sample_weight=None):
        super().partial_fit(X=X, y=y, classes=classes,
                            sample_weight=sample_weight)
        if len(self.classes_) > 2:
            raise ValueError('Only 2 classes are supported')
        return self


class TaggedBinaryClassifier(UntaggedBinaryClassifier):
    # Toy classifier that only supports binary classification.
    def _more_tags(self):
        return {'binary_only': True}


class RequiresPositiveYRegressor(LinearRegression):

    def fit(self, X, y):
        X, y = self._validate_data(X, y, multi_output=True)
        if (y <= 0).any():
            raise ValueError('negative y values not supported!')
        return super().fit(X, y)

    def _more_tags(self):
        return {"requires_positive_y": True}


class PoorScoreLogisticRegression(LogisticRegression):
    def decision_function(self, X):
        return super().decision_function(X) + 1

    def _more_tags(self):
        return {"poor_score": True}


def test_not_an_array_array_function():
    if np_version < parse_version('1.17'):
        raise SkipTest("array_function protocol not supported in numpy <1.17")
    not_array = _NotAnArray(np.ones(10))
    msg = "Don't want to call array_function sum!"
    assert_raises_regex(TypeError, msg, np.sum, not_array)
    # always returns True
    assert np.may_share_memory(not_array, None)


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
    msg = "Passing a class was deprecated"
    assert_raises_regex(TypeError, msg, check_estimator, object)
    msg = "object has no attribute '_get_tags'"
    assert_raises_regex(AttributeError, msg, check_estimator, object())
    msg = (
        "Parameter 'p' of estimator 'HasMutableParameters' is of type "
        "object which is not allowed"
    )
    # check that the "default_constructible" test checks for mutable parameters
    check_estimator(HasImmutableParameters())  # should pass
    assert_raises_regex(
        AssertionError, msg, check_estimator, HasMutableParameters()
    )
    # check that values returned by get_params match set_params
    msg = "get_params result does not match what was passed to set_params"
    assert_raises_regex(AssertionError, msg, check_estimator,
                        ModifiesValueInsteadOfRaisingError())
    assert_warns(UserWarning, check_estimator, RaisesErrorInSetParams())
    assert_raises_regex(AssertionError, msg, check_estimator,
                        ModifiesAnotherValue())
    # check that we have a fit method
    msg = "object has no attribute 'fit'"
    assert_raises_regex(AttributeError, msg, check_estimator, BaseEstimator())
    # check that fit does input validation
    msg = "Did not raise"
    assert_raises_regex(AssertionError, msg, check_estimator,
                        BaseBadClassifier())
    # check that sample_weights in fit accepts pandas.Series type
    try:
        from pandas import Series  # noqa
        msg = ("Estimator NoSampleWeightPandasSeriesType raises error if "
               "'sample_weight' parameter is of type pandas.Series")
        assert_raises_regex(
            ValueError, msg, check_estimator, NoSampleWeightPandasSeriesType())
    except ImportError:
        pass
    # check that predict does input validation (doesn't accept dicts in input)
    msg = "Estimator doesn't check for NaN and inf in predict"
    assert_raises_regex(AssertionError, msg, check_estimator,
                        NoCheckinPredict())
    # check that estimator state does not change
    # at transform/predict/predict_proba time
    msg = 'Estimator changes __dict__ during predict'
    assert_raises_regex(AssertionError, msg, check_estimator, ChangesDict())
    # check that `fit` only changes attribures that
    # are private (start with an _ or end with a _).
    msg = ('Estimator ChangesWrongAttribute should not change or mutate  '
           'the parameter wrong_attribute from 0 to 1 during fit.')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, ChangesWrongAttribute())
    check_estimator(ChangesUnderscoreAttribute())
    # check that `fit` doesn't add any public attribute
    msg = (r'Estimator adds public attribute\(s\) during the fit method.'
           ' Estimators are only allowed to add private attributes'
           ' either started with _ or ended'
           ' with _ but wrong_attribute added')
    assert_raises_regex(AssertionError, msg,
                        check_estimator, SetsWrongAttribute())
    # check for sample order invariance
    name = NotInvariantSampleOrder.__name__
    method = 'predict'
    msg = ("{method} of {name} is not invariant when applied to a dataset"
           "with different sample order.").format(method=method, name=name)
    assert_raises_regex(AssertionError, msg,
                        check_estimator, NotInvariantSampleOrder())
    # check for invariant method
    name = NotInvariantPredict.__name__
    method = 'predict'
    msg = ("{method} of {name} is not invariant when applied "
           "to a subset.").format(method=method, name=name)
    assert_raises_regex(AssertionError, msg,
                        check_estimator, NotInvariantPredict())
    # check for sparse matrix input handling
    name = NoSparseClassifierException.__name__
    msg = ("Estimator %s doesn't seem to fail gracefully on sparse data: "
           "it should raise a TypeError or a ValueError if sparse input "
           "is explicitly not supported." % name)
    assert_raises_regex(
        AssertionError, msg, check_estimator, NoSparseClassifierException()
    )

    # check for sparse matrix input handling, if the wrong TypeError is thrown
    for clf in [NoSparseClassifierTypeError, NoSparseClassifierValueError]:
        name = clf.__name__
        msg = ("Estimator %s doesn't seem to fail gracefully on sparse data: "
               "error message should state explicitly that sparse input is"
               " not supported if this is not the case." % name)
        assert_raises_regex(
            AssertionError, msg, check_estimator, clf()
        )
    # check that the equality of the result for sparse/dense input is checked
    assert_raises_regex(AssertionError,
                        "Not equal to tolerance",
                        check_estimator, SparseDenseDifferentPredict())
    assert_raises_regex(AssertionError,
                        "Not equal to tolerance",
                        check_estimator,
                        SparseDenseDifferentPredictProba())

    # Large indices test on bad estimator
    msg = ('Estimator LargeSparseNotSupportedClassifier doesn\'t seem to '
           r'support \S{3}_64 matrix, and is not failing gracefully.*')
    assert_raises_regex(AssertionError, msg, check_estimator,
                        LargeSparseNotSupportedClassifier())

    # does error on binary_only untagged estimator
    msg = 'Only 2 classes are supported'
    assert_raises_regex(ValueError, msg, check_estimator,
                        UntaggedBinaryClassifier())

    # non-regression test for estimators transforming to sparse data
    check_estimator(SparseTransformer())

    # doesn't error on actual estimator
    check_estimator(LogisticRegression())
    check_estimator(LogisticRegression(C=0.01))
    check_estimator(MultiTaskElasticNet())

    # doesn't error on binary_only tagged estimator
    check_estimator(TaggedBinaryClassifier())

    # Check regressor with requires_positive_y estimator tag
    msg = 'negative y values not supported!'
    assert_raises_regex(ValueError, msg, check_estimator,
                        RequiresPositiveYRegressor())

    # Does not raise error on classifier with poor_score tag
    check_estimator(PoorScoreLogisticRegression())


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
        with ignore_warnings(category=FutureWarning):
            # when 'est = SGDClassifier()'
            est = Estimator()
            _set_checking_parameters(est)
            set_random_state(est)
            # without fitting
            old_hash = joblib.hash(est)
            check_estimator(est)
        assert old_hash == joblib.hash(est)

        with ignore_warnings(category=FutureWarning):
            # when 'est = SGDClassifier()'
            est = Estimator()
            _set_checking_parameters(est)
            set_random_state(est)
            # with fitting
            est.fit(iris.data + 10, iris.target)
            old_hash = joblib.hash(est)
            check_estimator(est)
        assert old_hash == joblib.hash(est)


def test_check_estimators_unfitted():
    # check that a ValueError/AttributeError is raised when calling predict
    # on an unfitted estimator
    msg = "Did not raise"
    assert_raises_regex(AssertionError, msg, check_estimators_unfitted,
                        "estimator", BaseBadClassifier())

    # check that CorrectNotFittedError inherit from either ValueError
    # or AttributeError
    check_estimators_unfitted("estimator", CorrectNotFittedErrorClassifier())


def test_check_no_attributes_set_in_init():
    class NonConformantEstimatorPrivateSet(BaseEstimator):
        def __init__(self):
            self.you_should_not_set_this_ = None

    class NonConformantEstimatorNoParamSet(BaseEstimator):
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
    # kernel or metric

    # test precomputed kernel
    est = SVC(kernel='precomputed')
    check_estimator(est)

    # test precomputed metric
    est = KNeighborsRegressor(metric='precomputed')
    check_estimator(est)


def test_check_classifier_data_not_an_array():
    assert_raises_regex(AssertionError,
                        'Not equal to tolerance',
                        check_classifier_data_not_an_array,
                        'estimator_name',
                        EstimatorInconsistentForPandas())


def test_check_regressor_data_not_an_array():
    assert_raises_regex(AssertionError,
                        'Not equal to tolerance',
                        check_regressor_data_not_an_array,
                        'estimator_name',
                        EstimatorInconsistentForPandas())


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


def test_all_estimators_all_public():
    # all_estimator should not fail when pytest is not installed and return
    # only public estimators
    estimators = all_estimators()
    for est in estimators:
        assert not est.__class__.__name__.startswith("_")


if __name__ == '__main__':
    # This module is run as a script to check that we have no dependency on
    # pytest for estimator checks.
    run_tests_without_pytest()


def test_xfail_ignored_in_check_estimator():
    # Make sure checks marked as xfail are just ignored and not run by
    # check_estimator(), but still raise a warning.
    assert_warns(SkipTestWarning, check_estimator, NuSVC())
