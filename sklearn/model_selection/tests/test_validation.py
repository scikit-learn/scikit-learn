"""Test the validation module"""
from __future__ import division

import sys
import warnings
import tempfile
import os
from time import sleep

import pytest
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from sklearn.exceptions import FitFailedWarning

from sklearn.model_selection.tests.test_search import FailingClassifier

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

from sklearn.model_selection import cross_val_score, ShuffleSplit
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_validate
from sklearn.model_selection import permutation_test_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import LeavePGroupsOut
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import learning_curve
from sklearn.model_selection import validation_curve
from sklearn.model_selection._validation import _check_is_permutation
from sklearn.model_selection._validation import _fit_and_score
from sklearn.model_selection._validation import _score

from sklearn.datasets import make_regression
from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.datasets import load_digits
from sklearn.metrics import explained_variance_score
from sklearn.metrics import make_scorer
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import precision_score
from sklearn.metrics import r2_score
from sklearn.metrics.scorer import check_scoring

from sklearn.linear_model import Ridge, LogisticRegression, SGDClassifier
from sklearn.linear_model import PassiveAggressiveClassifier, RidgeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.cluster import KMeans

from sklearn.impute import SimpleImputer

from sklearn.preprocessing import LabelEncoder
from sklearn.pipeline import Pipeline

from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.base import BaseEstimator
from sklearn.base import clone
from sklearn.multiclass import OneVsRestClassifier
from sklearn.utils import shuffle
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification

from sklearn.model_selection.tests.common import OneTimeSplitter
from sklearn.model_selection import GridSearchCV


try:
    WindowsError
except NameError:
    WindowsError = None


class MockImprovingEstimator(BaseEstimator):
    """Dummy classifier to test the learning curve"""
    def __init__(self, n_max_train_sizes):
        self.n_max_train_sizes = n_max_train_sizes
        self.train_sizes = 0
        self.X_subset = None

    def fit(self, X_subset, y_subset=None):
        self.X_subset = X_subset
        self.train_sizes = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, Y=None):
        # training score becomes worse (2 -> 1), test error better (0 -> 1)
        if self._is_training_data(X):
            return 2. - float(self.train_sizes) / self.n_max_train_sizes
        else:
            return float(self.train_sizes) / self.n_max_train_sizes

    def _is_training_data(self, X):
        return X is self.X_subset


class MockIncrementalImprovingEstimator(MockImprovingEstimator):
    """Dummy classifier that provides partial_fit"""
    def __init__(self, n_max_train_sizes):
        super(MockIncrementalImprovingEstimator,
              self).__init__(n_max_train_sizes)
        self.x = None

    def _is_training_data(self, X):
        return self.x in X

    def partial_fit(self, X, y=None, **params):
        self.train_sizes += X.shape[0]
        self.x = X[0]


class MockEstimatorWithParameter(BaseEstimator):
    """Dummy classifier to test the validation curve"""
    def __init__(self, param=0.5):
        self.X_subset = None
        self.param = param

    def fit(self, X_subset, y_subset):
        self.X_subset = X_subset
        self.train_sizes = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, y=None):
        return self.param if self._is_training_data(X) else 1 - self.param

    def _is_training_data(self, X):
        return X is self.X_subset


class MockEstimatorWithSingleFitCallAllowed(MockEstimatorWithParameter):
    """Dummy classifier that disallows repeated calls of fit method"""

    def fit(self, X_subset, y_subset):
        assert not hasattr(self, 'fit_called_'), \
                   'fit is called the second time'
        self.fit_called_ = True
        return super(type(self), self).fit(X_subset, y_subset)

    def predict(self, X):
        raise NotImplementedError


class MockClassifier(object):
    """Dummy classifier to test the cross-validation"""

    def __init__(self, a=0, allow_nd=False):
        self.a = a
        self.allow_nd = allow_nd

    def fit(self, X, Y=None, sample_weight=None, class_prior=None,
            sparse_sample_weight=None, sparse_param=None, dummy_int=None,
            dummy_str=None, dummy_obj=None, callback=None):
        """The dummy arguments are to test that this fit function can
        accept non-array arguments through cross-validation, such as:
            - int
            - str (this is actually array-like)
            - object
            - function
        """
        self.dummy_int = dummy_int
        self.dummy_str = dummy_str
        self.dummy_obj = dummy_obj
        if callback is not None:
            callback(self)

        if self.allow_nd:
            X = X.reshape(len(X), -1)
        if X.ndim >= 3 and not self.allow_nd:
            raise ValueError('X cannot be d')
        if sample_weight is not None:
            assert sample_weight.shape[0] == X.shape[0], (
                'MockClassifier extra fit_param ' 
                'sample_weight.shape[0] is {0}, should be {1}'
                .format(sample_weight.shape[0], X.shape[0]))
        if class_prior is not None:
            assert class_prior.shape[0] == len(np.unique(y)), (
                'MockClassifier extra fit_param class_prior.shape[0]'
                ' is {0}, should be {1}'.format(class_prior.shape[0],
                                                len(np.unique(y))))
        if sparse_sample_weight is not None:
            fmt = ('MockClassifier extra fit_param sparse_sample_weight'
                   '.shape[0] is {0}, should be {1}')
            assert sparse_sample_weight.shape[0] == X.shape[0], \
                fmt.format(sparse_sample_weight.shape[0], X.shape[0])
        if sparse_param is not None:
            fmt = ('MockClassifier extra fit_param sparse_param.shape '
                   'is ({0}, {1}), should be ({2}, {3})')
            assert sparse_param.shape == P_sparse.shape, (
                fmt.format(sparse_param.shape[0],
                           sparse_param.shape[1],
                           P_sparse.shape[0], P_sparse.shape[1]))
        return self

    def predict(self, T):
        if self.allow_nd:
            T = T.reshape(len(T), -1)
        return T[:, 0]

    def score(self, X=None, Y=None):
        return 1. / (1 + np.abs(self.a))

    def get_params(self, deep=False):
        return {'a': self.a, 'allow_nd': self.allow_nd}


# XXX: use 2D array, since 1D X is being detected as a single sample in
# check_consistent_length
X = np.ones((10, 2))
X_sparse = coo_matrix(X)
y = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
# The number of samples per class needs to be > n_splits,
# for StratifiedKFold(n_splits=3)
y2 = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3, 3])
P_sparse = coo_matrix(np.eye(5))


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score():
    clf = MockClassifier()

    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        scores = cross_val_score(clf, X, y2)
        assert_array_equal(scores, clf.score(X, y2))

        # test with multioutput y
        multioutput_y = np.column_stack([y2, y2[::-1]])
        scores = cross_val_score(clf, X_sparse, multioutput_y)
        assert_array_equal(scores, clf.score(X_sparse, multioutput_y))

        scores = cross_val_score(clf, X_sparse, y2)
        assert_array_equal(scores, clf.score(X_sparse, y2))

        # test with multioutput y
        scores = cross_val_score(clf, X_sparse, multioutput_y)
        assert_array_equal(scores, clf.score(X_sparse, multioutput_y))

    # test with X and y as list
    list_check = lambda x: isinstance(x, list)
    clf = CheckingClassifier(check_X=list_check)
    scores = cross_val_score(clf, X.tolist(), y2.tolist())

    clf = CheckingClassifier(check_y=list_check)
    scores = cross_val_score(clf, X, y2.tolist())

    assert_raises(ValueError, cross_val_score, clf, X, y2, scoring="sklearn")

    # test with 3d X and
    X_3d = X[:, :, np.newaxis]
    clf = MockClassifier(allow_nd=True)
    scores = cross_val_score(clf, X_3d, y2)

    clf = MockClassifier(allow_nd=False)
    assert_raises(ValueError, cross_val_score, clf, X_3d, y2,
                  error_score='raise')


@pytest.mark.filterwarnings('ignore:You should specify a value for')  # 0.22
def test_cross_validate_many_jobs():
    # regression test for #12154: cv='warn' with n_jobs>1 trigger a copy of
    # the parameters leading to a failure in check_cv due to cv is 'warn'
    # instead of cv == 'warn'.
    X, y = load_iris(return_X_y=True)
    clf = SVC(gamma='auto')
    grid = GridSearchCV(clf, param_grid={'C': [1, 10]})
    cross_validate(grid, X, y, n_jobs=2)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_validate_invalid_scoring_param():
    X, y = make_classification(random_state=0)
    estimator = MockClassifier()

    # Test the errors
    error_message_regexp = ".*must be unique strings.*"

    # List/tuple of callables should raise a message advising users to use
    # dict of names to callables mapping
    assert_raises_regex(ValueError, error_message_regexp,
                        cross_validate, estimator, X, y,
                        scoring=(make_scorer(precision_score),
                                 make_scorer(accuracy_score)))
    assert_raises_regex(ValueError, error_message_regexp,
                        cross_validate, estimator, X, y,
                        scoring=(make_scorer(precision_score),))

    # So should empty lists/tuples
    assert_raises_regex(ValueError, error_message_regexp + "Empty list.*",
                        cross_validate, estimator, X, y, scoring=())

    # So should duplicated entries
    assert_raises_regex(ValueError, error_message_regexp + "Duplicate.*",
                        cross_validate, estimator, X, y,
                        scoring=('f1_micro', 'f1_micro'))

    # Nested Lists should raise a generic error message
    assert_raises_regex(ValueError, error_message_regexp,
                        cross_validate, estimator, X, y,
                        scoring=[[make_scorer(precision_score)]])

    error_message_regexp = (".*should either be.*string or callable.*for "
                            "single.*.*dict.*for multi.*")

    # Empty dict should raise invalid scoring error
    assert_raises_regex(ValueError, "An empty dict",
                        cross_validate, estimator, X, y, scoring=(dict()))

    # And so should any other invalid entry
    assert_raises_regex(ValueError, error_message_regexp,
                        cross_validate, estimator, X, y, scoring=5)

    multiclass_scorer = make_scorer(precision_recall_fscore_support)

    # Multiclass Scorers that return multiple values are not supported yet
    assert_raises_regex(ValueError,
                        "Classification metrics can't handle a mix of "
                        "binary and continuous targets",
                        cross_validate, estimator, X, y,
                        scoring=multiclass_scorer)
    assert_raises_regex(ValueError,
                        "Classification metrics can't handle a mix of "
                        "binary and continuous targets",
                        cross_validate, estimator, X, y,
                        scoring={"foo": multiclass_scorer})

    multivalued_scorer = make_scorer(confusion_matrix)

    # Multiclass Scorers that return multiple values are not supported yet
    assert_raises_regex(ValueError, "scoring must return a number, got",
                        cross_validate, SVC(gamma='scale'), X, y,
                        scoring=multivalued_scorer)
    assert_raises_regex(ValueError, "scoring must return a number, got",
                        cross_validate, SVC(gamma='scale'), X, y,
                        scoring={"foo": multivalued_scorer})

    assert_raises_regex(ValueError, "'mse' is not a valid scoring value.",
                        cross_validate, SVC(), X, y, scoring="mse")


def test_cross_validate():
    # Compute train and test mse/r2 scores
    cv = KFold(n_splits=5)

    # Regression
    X_reg, y_reg = make_regression(n_samples=30, random_state=0)
    reg = Ridge(random_state=0)

    # Classification
    X_clf, y_clf = make_classification(n_samples=30, random_state=0)
    clf = SVC(kernel="linear", random_state=0)

    for X, y, est in ((X_reg, y_reg, reg), (X_clf, y_clf, clf)):
        # It's okay to evaluate regression metrics on classification too
        mse_scorer = check_scoring(est, 'neg_mean_squared_error')
        r2_scorer = check_scoring(est, 'r2')
        train_mse_scores = []
        test_mse_scores = []
        train_r2_scores = []
        test_r2_scores = []
        fitted_estimators = []
        for train, test in cv.split(X, y):
            est = clone(reg).fit(X[train], y[train])
            train_mse_scores.append(mse_scorer(est, X[train], y[train]))
            train_r2_scores.append(r2_scorer(est, X[train], y[train]))
            test_mse_scores.append(mse_scorer(est, X[test], y[test]))
            test_r2_scores.append(r2_scorer(est, X[test], y[test]))
            fitted_estimators.append(est)

        train_mse_scores = np.array(train_mse_scores)
        test_mse_scores = np.array(test_mse_scores)
        train_r2_scores = np.array(train_r2_scores)
        test_r2_scores = np.array(test_r2_scores)
        fitted_estimators = np.array(fitted_estimators)

        scores = (train_mse_scores, test_mse_scores, train_r2_scores,
                  test_r2_scores, fitted_estimators)

        check_cross_validate_single_metric(est, X, y, scores)
        check_cross_validate_multi_metric(est, X, y, scores)


def check_cross_validate_single_metric(clf, X, y, scores):
    (train_mse_scores, test_mse_scores, train_r2_scores,
     test_r2_scores, fitted_estimators) = scores
    # Test single metric evaluation when scoring is string or singleton list
    for (return_train_score, dict_len) in ((True, 4), (False, 3)):
        # Single metric passed as a string
        if return_train_score:
            mse_scores_dict = cross_validate(clf, X, y, cv=5,
                                             scoring='neg_mean_squared_error',
                                             return_train_score=True)
            assert_array_almost_equal(mse_scores_dict['train_score'],
                                      train_mse_scores)
        else:
            mse_scores_dict = cross_validate(clf, X, y, cv=5,
                                             scoring='neg_mean_squared_error',
                                             return_train_score=False)
        assert isinstance(mse_scores_dict, dict)
        assert_equal(len(mse_scores_dict), dict_len)
        assert_array_almost_equal(mse_scores_dict['test_score'],
                                  test_mse_scores)

        # Single metric passed as a list
        if return_train_score:
            # It must be True by default - deprecated
            r2_scores_dict = cross_validate(clf, X, y, cv=5, scoring=['r2'],
                                            return_train_score=True)
            assert_array_almost_equal(r2_scores_dict['train_r2'],
                                      train_r2_scores, True)
        else:
            r2_scores_dict = cross_validate(clf, X, y, cv=5, scoring=['r2'],
                                            return_train_score=False)
        assert isinstance(r2_scores_dict, dict)
        assert_equal(len(r2_scores_dict), dict_len)
        assert_array_almost_equal(r2_scores_dict['test_r2'], test_r2_scores)

    # Test return_estimator option
    mse_scores_dict = cross_validate(clf, X, y, cv=5,
                                     scoring='neg_mean_squared_error',
                                     return_estimator=True)
    for k, est in enumerate(mse_scores_dict['estimator']):
        assert_almost_equal(est.coef_, fitted_estimators[k].coef_)
        assert_almost_equal(est.intercept_, fitted_estimators[k].intercept_)


def check_cross_validate_multi_metric(clf, X, y, scores):
    # Test multimetric evaluation when scoring is a list / dict
    (train_mse_scores, test_mse_scores, train_r2_scores,
     test_r2_scores, fitted_estimators) = scores
    all_scoring = (('r2', 'neg_mean_squared_error'),
                   {'r2': make_scorer(r2_score),
                    'neg_mean_squared_error': 'neg_mean_squared_error'})

    keys_sans_train = set(('test_r2', 'test_neg_mean_squared_error',
                           'fit_time', 'score_time'))
    keys_with_train = keys_sans_train.union(
        set(('train_r2', 'train_neg_mean_squared_error')))

    for return_train_score in (True, False):
        for scoring in all_scoring:
            if return_train_score:
                # return_train_score must be True by default - deprecated
                cv_results = cross_validate(clf, X, y, cv=5, scoring=scoring,
                                            return_train_score=True)
                assert_array_almost_equal(cv_results['train_r2'],
                                          train_r2_scores)
                assert_array_almost_equal(
                    cv_results['train_neg_mean_squared_error'],
                    train_mse_scores)
            else:
                cv_results = cross_validate(clf, X, y, cv=5, scoring=scoring,
                                            return_train_score=False)
            assert isinstance(cv_results, dict)
            assert_equal(set(cv_results.keys()),
                         keys_with_train if return_train_score
                         else keys_sans_train)
            assert_array_almost_equal(cv_results['test_r2'], test_r2_scores)
            assert_array_almost_equal(
                cv_results['test_neg_mean_squared_error'], test_mse_scores)

            # Make sure all the arrays are of np.ndarray type
            assert type(cv_results['test_r2']) == np.ndarray
            assert (type(cv_results['test_neg_mean_squared_error']) ==
                    np.ndarray)
            assert type(cv_results['fit_time']) == np.ndarray
            assert type(cv_results['score_time']) == np.ndarray

            # Ensure all the times are within sane limits
            assert np.all(cv_results['fit_time'] >= 0)
            assert np.all(cv_results['fit_time'] < 10)
            assert np.all(cv_results['score_time'] >= 0)
            assert np.all(cv_results['score_time'] < 10)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score_predict_groups():
    # Check if ValueError (when groups is None) propagates to cross_val_score
    # and cross_val_predict
    # And also check if groups is correctly passed to the cv object
    X, y = make_classification(n_samples=20, n_classes=2, random_state=0)

    clf = SVC(kernel="linear")

    group_cvs = [LeaveOneGroupOut(), LeavePGroupsOut(2), GroupKFold(),
                 GroupShuffleSplit()]
    for cv in group_cvs:
        assert_raise_message(ValueError,
                             "The 'groups' parameter should not be None.",
                             cross_val_score, estimator=clf, X=X, y=y, cv=cv)
        assert_raise_message(ValueError,
                             "The 'groups' parameter should not be None.",
                             cross_val_predict, estimator=clf, X=X, y=y, cv=cv)


@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score_pandas():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((Series, DataFrame))
    except ImportError:
        pass
    for TargetType, InputFeatureType in types:
        # X dataframe, y series
        # 3 fold cross val is used so we need atleast 3 samples per class
        X_df, y_ser = InputFeatureType(X), TargetType(y2)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        cross_val_score(clf, X_df, y_ser)


def test_cross_val_score_mask():
    # test that cross_val_score works with boolean masks
    svm = SVC(kernel="linear")
    iris = load_iris()
    X, y = iris.data, iris.target
    kfold = KFold(5)
    scores_indices = cross_val_score(svm, X, y, cv=kfold)
    kfold = KFold(5)
    cv_masks = []
    for train, test in kfold.split(X, y):
        mask_train = np.zeros(len(y), dtype=np.bool)
        mask_test = np.zeros(len(y), dtype=np.bool)
        mask_train[train] = 1
        mask_test[test] = 1
        cv_masks.append((train, test))
    scores_masks = cross_val_score(svm, X, y, cv=cv_masks)
    assert_array_equal(scores_indices, scores_masks)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score_precomputed():
    # test for svm with precomputed kernel
    svm = SVC(kernel="precomputed")
    iris = load_iris()
    X, y = iris.data, iris.target
    linear_kernel = np.dot(X, X.T)
    score_precomputed = cross_val_score(svm, linear_kernel, y)
    svm = SVC(kernel="linear")
    score_linear = cross_val_score(svm, X, y)
    assert_array_almost_equal(score_precomputed, score_linear)

    # test with callable
    svm = SVC(gamma='scale', kernel=lambda x, y: np.dot(x, y.T))
    score_callable = cross_val_score(svm, X, y)
    assert_array_almost_equal(score_precomputed, score_callable)

    # Error raised for non-square X
    svm = SVC(kernel="precomputed")
    assert_raises(ValueError, cross_val_score, svm, X, y)

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    assert_raises(ValueError, cross_val_score, svm,
                  linear_kernel.tolist(), y)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score_fit_params():
    clf = MockClassifier()
    n_samples = X.shape[0]
    n_classes = len(np.unique(y))

    W_sparse = coo_matrix((np.array([1]), (np.array([1]), np.array([0]))),
                          shape=(10, 1))
    P_sparse = coo_matrix(np.eye(5))

    DUMMY_INT = 42
    DUMMY_STR = '42'
    DUMMY_OBJ = object()

    def assert_fit_params(clf):
        # Function to test that the values are passed correctly to the
        # classifier arguments for non-array type

        assert_equal(clf.dummy_int, DUMMY_INT)
        assert_equal(clf.dummy_str, DUMMY_STR)
        assert_equal(clf.dummy_obj, DUMMY_OBJ)

    fit_params = {'sample_weight': np.ones(n_samples),
                  'class_prior': np.full(n_classes, 1. / n_classes),
                  'sparse_sample_weight': W_sparse,
                  'sparse_param': P_sparse,
                  'dummy_int': DUMMY_INT,
                  'dummy_str': DUMMY_STR,
                  'dummy_obj': DUMMY_OBJ,
                  'callback': assert_fit_params}
    cross_val_score(clf, X, y, fit_params=fit_params)


def test_cross_val_score_score_func():
    clf = MockClassifier()
    _score_func_args = []

    def score_func(y_test, y_predict):
        _score_func_args.append((y_test, y_predict))
        return 1.0

    with warnings.catch_warnings(record=True):
        scoring = make_scorer(score_func)
        score = cross_val_score(clf, X, y, scoring=scoring, cv=3)
    assert_array_equal(score, [1.0, 1.0, 1.0])
    # Test that score function is called only 3 times (for cv=3)
    assert len(_score_func_args) == 3


def test_cross_val_score_errors():
    class BrokenEstimator:
        pass

    assert_raises(TypeError, cross_val_score, BrokenEstimator(), X)


def test_cross_val_score_with_score_func_classification():
    iris = load_iris()
    clf = SVC(kernel='linear')

    # Default score (should be the accuracy score)
    scores = cross_val_score(clf, iris.data, iris.target, cv=5)
    assert_array_almost_equal(scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # Correct classification score (aka. zero / one score) - should be the
    # same as the default estimator score
    zo_scores = cross_val_score(clf, iris.data, iris.target,
                                scoring="accuracy", cv=5)
    assert_array_almost_equal(zo_scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # F1 score (class are balanced so f1_score should be equal to zero/one
    # score
    f1_scores = cross_val_score(clf, iris.data, iris.target,
                                scoring="f1_weighted", cv=5)
    assert_array_almost_equal(f1_scores, [0.97, 1., 0.97, 0.97, 1.], 2)


def test_cross_val_score_with_score_func_regression():
    X, y = make_regression(n_samples=30, n_features=20, n_informative=5,
                           random_state=0)
    reg = Ridge()

    # Default score of the Ridge regression estimator
    scores = cross_val_score(reg, X, y, cv=5)
    assert_array_almost_equal(scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # R2 score (aka. determination coefficient) - should be the
    # same as the default estimator score
    r2_scores = cross_val_score(reg, X, y, scoring="r2", cv=5)
    assert_array_almost_equal(r2_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # Mean squared error; this is a loss function, so "scores" are negative
    neg_mse_scores = cross_val_score(reg, X, y, cv=5,
                                     scoring="neg_mean_squared_error")
    expected_neg_mse = np.array([-763.07, -553.16, -274.38, -273.26, -1681.99])
    assert_array_almost_equal(neg_mse_scores, expected_neg_mse, 2)

    # Explained variance
    scoring = make_scorer(explained_variance_score)
    ev_scores = cross_val_score(reg, X, y, cv=5, scoring=scoring)
    assert_array_almost_equal(ev_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    X_sparse = coo_matrix(X)
    y = iris.target
    svm = SVC(kernel='linear')
    cv = StratifiedKFold(2)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")
    assert_greater(score, 0.9)
    assert_almost_equal(pvalue, 0.0, 1)

    score_group, _, pvalue_group = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy",
        groups=np.ones(y.size), random_state=0)
    assert score_group == score
    assert pvalue_group == pvalue

    # check that we obtain the same results with a sparse representation
    svm_sparse = SVC(kernel='linear')
    cv_sparse = StratifiedKFold(2)
    score_group, _, pvalue_group = permutation_test_score(
        svm_sparse, X_sparse, y, n_permutations=30, cv=cv_sparse,
        scoring="accuracy", groups=np.ones(y.size), random_state=0)

    assert score_group == score
    assert pvalue_group == pvalue

    # test with custom scoring object
    def custom_score(y_true, y_pred):
        return (((y_true == y_pred).sum() - (y_true != y_pred).sum()) /
                y_true.shape[0])

    scorer = make_scorer(custom_score)
    score, _, pvalue = permutation_test_score(
        svm, X, y, n_permutations=100, scoring=scorer, cv=cv, random_state=0)
    assert_almost_equal(score, .93, 2)
    assert_almost_equal(pvalue, 0.01, 3)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")

    assert_less(score, 0.5)
    assert_greater(pvalue, 0.2)


def test_permutation_test_score_allow_nans():
    # Check that permutation_test_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0] / 2)
    p = Pipeline([
        ('imputer', SimpleImputer(strategy='mean', missing_values=np.nan)),
        ('classifier', MockClassifier()),
    ])
    permutation_test_score(p, X, y, cv=5)


def test_cross_val_score_allow_nans():
    # Check that cross_val_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0] / 2)
    p = Pipeline([
        ('imputer', SimpleImputer(strategy='mean', missing_values=np.nan)),
        ('classifier', MockClassifier()),
    ])
    cross_val_score(p, X, y, cv=5)


def test_cross_val_score_multilabel():
    X = np.array([[-3, 4], [2, 4], [3, 3], [0, 2], [-3, 1],
                  [-2, 1], [0, 0], [-2, -1], [-1, -2], [1, -2]])
    y = np.array([[1, 1], [0, 1], [0, 1], [0, 1], [1, 1],
                  [0, 1], [1, 0], [1, 1], [1, 0], [0, 0]])
    clf = KNeighborsClassifier(n_neighbors=1)
    scoring_micro = make_scorer(precision_score, average='micro')
    scoring_macro = make_scorer(precision_score, average='macro')
    scoring_samples = make_scorer(precision_score, average='samples')
    score_micro = cross_val_score(clf, X, y, scoring=scoring_micro, cv=5)
    score_macro = cross_val_score(clf, X, y, scoring=scoring_macro, cv=5)
    score_samples = cross_val_score(clf, X, y, scoring=scoring_samples, cv=5)
    assert_almost_equal(score_micro, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 3])
    assert_almost_equal(score_macro, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 4])
    assert_almost_equal(score_samples, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 4])


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict():
    boston = load_boston()
    X, y = boston.data, boston.target
    cv = KFold()

    est = Ridge()

    # Naive loop (should be same as cross_val_predict):
    preds2 = np.zeros_like(y)
    for train, test in cv.split(X, y):
        est.fit(X[train], y[train])
        preds2[test] = est.predict(X[test])

    preds = cross_val_predict(est, X, y, cv=cv)
    assert_array_almost_equal(preds, preds2)

    preds = cross_val_predict(est, X, y)
    assert_equal(len(preds), len(y))

    cv = LeaveOneOut()
    preds = cross_val_predict(est, X, y, cv=cv)
    assert_equal(len(preds), len(y))

    Xsp = X.copy()
    Xsp *= (Xsp > np.median(Xsp))
    Xsp = coo_matrix(Xsp)
    preds = cross_val_predict(est, Xsp, y)
    assert_array_almost_equal(len(preds), len(y))

    preds = cross_val_predict(KMeans(), X)
    assert_equal(len(preds), len(y))

    class BadCV():
        def split(self, X, y=None, groups=None):
            for i in range(4):
                yield np.array([0, 1, 2, 3]), np.array([4, 5, 6, 7, 8])

    assert_raises(ValueError, cross_val_predict, est, X, y, cv=BadCV())

    X, y = load_iris(return_X_y=True)

    warning_message = ('Number of classes in training fold (2) does '
                       'not match total number of classes (3). '
                       'Results may not be appropriate for your use case.')
    assert_warns_message(RuntimeWarning, warning_message,
                         cross_val_predict, LogisticRegression(),
                         X, y, method='predict_proba', cv=KFold(2))


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_decision_function_shape():
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='decision_function')
    assert_equal(preds.shape, (50,))

    X, y = load_iris(return_X_y=True)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='decision_function')
    assert_equal(preds.shape, (150, 3))

    # This specifically tests imbalanced splits for binary
    # classification with decision_function. This is only
    # applicable to classifiers that can be fit on a single
    # class.
    X = X[:100]
    y = y[:100]
    assert_raise_message(ValueError,
                         'Only 1 class/es in training fold, this'
                         ' is not supported for decision_function'
                         ' with imbalanced folds. To fix '
                         'this, use a cross-validation technique '
                         'resulting in properly stratified folds',
                         cross_val_predict, RidgeClassifier(), X, y,
                         method='decision_function', cv=KFold(2))

    X, y = load_digits(return_X_y=True)
    est = SVC(kernel='linear', decision_function_shape='ovo')

    preds = cross_val_predict(est,
                              X, y,
                              method='decision_function')
    assert_equal(preds.shape, (1797, 45))

    ind = np.argsort(y)
    X, y = X[ind], y[ind]
    assert_raises_regex(ValueError,
                        r'Output shape \(599L?, 21L?\) of decision_function '
                        r'does not match number of classes \(7\) in fold. '
                        'Irregular decision_function .*',
                        cross_val_predict, est, X, y,
                        cv=KFold(n_splits=3), method='decision_function')


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_predict_proba_shape():
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='predict_proba')
    assert_equal(preds.shape, (50, 2))

    X, y = load_iris(return_X_y=True)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='predict_proba')
    assert_equal(preds.shape, (150, 3))


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_predict_log_proba_shape():
    X, y = make_classification(n_classes=2, n_samples=50, random_state=0)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='predict_log_proba')
    assert_equal(preds.shape, (50, 2))

    X, y = load_iris(return_X_y=True)

    preds = cross_val_predict(LogisticRegression(), X, y,
                              method='predict_log_proba')
    assert_equal(preds.shape, (150, 3))


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_input_types():
    iris = load_iris()
    X, y = iris.data, iris.target
    X_sparse = coo_matrix(X)
    multioutput_y = np.column_stack([y, y[::-1]])

    clf = Ridge(fit_intercept=False, random_state=0)
    # 3 fold cv is used --> atleast 3 samples per class
    # Smoke test
    predictions = cross_val_predict(clf, X, y)
    assert_equal(predictions.shape, (150,))

    # test with multioutput y
    predictions = cross_val_predict(clf, X_sparse, multioutput_y)
    assert_equal(predictions.shape, (150, 2))

    predictions = cross_val_predict(clf, X_sparse, y)
    assert_array_equal(predictions.shape, (150,))

    # test with multioutput y
    predictions = cross_val_predict(clf, X_sparse, multioutput_y)
    assert_array_equal(predictions.shape, (150, 2))

    # test with X and y as list
    list_check = lambda x: isinstance(x, list)
    clf = CheckingClassifier(check_X=list_check)
    predictions = cross_val_predict(clf, X.tolist(), y.tolist())

    clf = CheckingClassifier(check_y=list_check)
    predictions = cross_val_predict(clf, X, y.tolist())

    # test with X and y as list and non empty method
    predictions = cross_val_predict(LogisticRegression(), X.tolist(),
                                    y.tolist(), method='decision_function')
    predictions = cross_val_predict(LogisticRegression(), X,
                                    y.tolist(), method='decision_function')

    # test with 3d X and
    X_3d = X[:, :, np.newaxis]
    check_3d = lambda x: x.ndim == 3
    clf = CheckingClassifier(check_X=check_3d)
    predictions = cross_val_predict(clf, X_3d, y)
    assert_array_equal(predictions.shape, (150,))


@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
# python3.7 deprecation warnings in pandas via matplotlib :-/
def test_cross_val_predict_pandas():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((Series, DataFrame))
    except ImportError:
        pass
    for TargetType, InputFeatureType in types:
        # X dataframe, y series
        X_df, y_ser = InputFeatureType(X), TargetType(y2)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        cross_val_predict(clf, X_df, y_ser)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_score_sparse_fit_params():
    iris = load_iris()
    X, y = iris.data, iris.target
    clf = MockClassifier()
    fit_params = {'sparse_sample_weight': coo_matrix(np.eye(X.shape[0]))}
    a = cross_val_score(clf, X, y, fit_params=fit_params)
    assert_array_equal(a, np.ones(3))


def test_learning_curve():
    n_samples = 30
    n_splits = 3
    X, y = make_classification(n_samples=n_samples, n_features=1,
                               n_informative=1, n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(n_samples * ((n_splits - 1) / n_splits))
    for shuffle_train in [False, True]:
        with warnings.catch_warnings(record=True) as w:
            train_sizes, train_scores, test_scores = learning_curve(
                estimator, X, y, cv=KFold(n_splits=n_splits),
                train_sizes=np.linspace(0.1, 1.0, 10),
                shuffle=shuffle_train)
        if len(w) > 0:
            raise RuntimeError("Unexpected warning: %r" % w[0].message)
        assert_equal(train_scores.shape, (10, 3))
        assert_equal(test_scores.shape, (10, 3))
        assert_array_equal(train_sizes, np.linspace(2, 20, 10))
        assert_array_almost_equal(train_scores.mean(axis=1),
                                  np.linspace(1.9, 1.0, 10))
        assert_array_almost_equal(test_scores.mean(axis=1),
                                  np.linspace(0.1, 1.0, 10))

        # Test a custom cv splitter that can iterate only once
        with warnings.catch_warnings(record=True) as w:
            train_sizes2, train_scores2, test_scores2 = learning_curve(
                estimator, X, y,
                cv=OneTimeSplitter(n_splits=n_splits, n_samples=n_samples),
                train_sizes=np.linspace(0.1, 1.0, 10),
                shuffle=shuffle_train)
        if len(w) > 0:
            raise RuntimeError("Unexpected warning: %r" % w[0].message)
        assert_array_almost_equal(train_scores2, train_scores)
        assert_array_almost_equal(test_scores2, test_scores)


def test_learning_curve_unsupervised():
    X, _ = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y=None, cv=3, train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_learning_curve_verbose():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        train_sizes, train_scores, test_scores = \
            learning_curve(estimator, X, y, cv=3, verbose=1)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    assert("[learning_curve]" in out)


def test_learning_curve_incremental_learning_not_possible():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    # The mockup does not have partial_fit()
    estimator = MockImprovingEstimator(1)
    assert_raises(ValueError, learning_curve, estimator, X, y,
                  exploit_incremental_learning=True)


def test_learning_curve_incremental_learning():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingEstimator(20)
    for shuffle_train in [False, True]:
        train_sizes, train_scores, test_scores = learning_curve(
            estimator, X, y, cv=3, exploit_incremental_learning=True,
            train_sizes=np.linspace(0.1, 1.0, 10), shuffle=shuffle_train)
        assert_array_equal(train_sizes, np.linspace(2, 20, 10))
        assert_array_almost_equal(train_scores.mean(axis=1),
                                  np.linspace(1.9, 1.0, 10))
        assert_array_almost_equal(test_scores.mean(axis=1),
                                  np.linspace(0.1, 1.0, 10))


def test_learning_curve_incremental_learning_unsupervised():
    X, _ = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingEstimator(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y=None, cv=3, exploit_incremental_learning=True,
        train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_learning_curve_batch_and_incremental_learning_are_equal():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    train_sizes = np.linspace(0.2, 1.0, 5)
    estimator = PassiveAggressiveClassifier(max_iter=1, tol=None,
                                            shuffle=False)

    train_sizes_inc, train_scores_inc, test_scores_inc = \
        learning_curve(
            estimator, X, y, train_sizes=train_sizes,
            cv=3, exploit_incremental_learning=True)
    train_sizes_batch, train_scores_batch, test_scores_batch = \
        learning_curve(
            estimator, X, y, cv=3, train_sizes=train_sizes,
            exploit_incremental_learning=False)

    assert_array_equal(train_sizes_inc, train_sizes_batch)
    assert_array_almost_equal(train_scores_inc.mean(axis=1),
                              train_scores_batch.mean(axis=1))
    assert_array_almost_equal(test_scores_inc.mean(axis=1),
                              test_scores_batch.mean(axis=1))


def test_learning_curve_n_sample_range_out_of_bounds():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0, 1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.0, 1.0])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.1, 1.1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0, 20])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[1, 21])


def test_learning_curve_remove_duplicate_sample_sizes():
    X, y = make_classification(n_samples=3, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(2)
    train_sizes, _, _ = assert_warns(
        RuntimeWarning, learning_curve, estimator, X, y, cv=3,
        train_sizes=np.linspace(0.33, 1.0, 3))
    assert_array_equal(train_sizes, [1, 2])


def test_learning_curve_with_boolean_indices():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    cv = KFold(n_splits=3)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


# 0.23. warning about tol not having its correct default value.
@pytest.mark.filterwarnings('ignore:max_iter and tol parameters have been')
def test_learning_curve_with_shuffle():
    # Following test case was designed this way to verify the code
    # changes made in pull request: #7506.
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [11, 12], [13, 14], [15, 16],
                 [17, 18], [19, 20], [7, 8], [9, 10], [11, 12], [13, 14],
                 [15, 16], [17, 18]])
    y = np.array([1, 1, 1, 2, 3, 4, 1, 1, 2, 3, 4, 1, 2, 3, 4])
    groups = np.array([1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 4])
    # Splits on these groups fail without shuffle as the first iteration
    # of the learning curve doesn't contain label 4 in the training set.
    estimator = PassiveAggressiveClassifier(max_iter=5, tol=None,
                                            shuffle=False)

    cv = GroupKFold(n_splits=2)
    train_sizes_batch, train_scores_batch, test_scores_batch = learning_curve(
        estimator, X, y, cv=cv, n_jobs=1, train_sizes=np.linspace(0.3, 1.0, 3),
        groups=groups, shuffle=True, random_state=2)
    assert_array_almost_equal(train_scores_batch.mean(axis=1),
                              np.array([0.75, 0.3, 0.36111111]))
    assert_array_almost_equal(test_scores_batch.mean(axis=1),
                              np.array([0.36111111, 0.25, 0.25]))
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=cv, n_jobs=1,
                  train_sizes=np.linspace(0.3, 1.0, 3), groups=groups,
                  error_score='raise')

    train_sizes_inc, train_scores_inc, test_scores_inc = learning_curve(
        estimator, X, y, cv=cv, n_jobs=1, train_sizes=np.linspace(0.3, 1.0, 3),
        groups=groups, shuffle=True, random_state=2,
        exploit_incremental_learning=True)
    assert_array_almost_equal(train_scores_inc.mean(axis=1),
                              train_scores_batch.mean(axis=1))
    assert_array_almost_equal(test_scores_inc.mean(axis=1),
                              test_scores_batch.mean(axis=1))


def test_validation_curve():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    param_range = np.linspace(0, 1, 10)
    with warnings.catch_warnings(record=True) as w:
        train_scores, test_scores = validation_curve(
            MockEstimatorWithParameter(), X, y, param_name="param",
            param_range=param_range, cv=2
        )
    if len(w) > 0:
        raise RuntimeError("Unexpected warning: %r" % w[0].message)

    assert_array_almost_equal(train_scores.mean(axis=1), param_range)
    assert_array_almost_equal(test_scores.mean(axis=1), 1 - param_range)


def test_validation_curve_clone_estimator():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)

    param_range = np.linspace(1, 0, 10)
    _, _ = validation_curve(
        MockEstimatorWithSingleFitCallAllowed(), X, y,
        param_name="param", param_range=param_range, cv=2
    )


def test_validation_curve_cv_splits_consistency():
    n_samples = 100
    n_splits = 5
    X, y = make_classification(n_samples=100, random_state=0)

    scores1 = validation_curve(SVC(kernel='linear', random_state=0), X, y,
                               'C', [0.1, 0.1, 0.2, 0.2],
                               cv=OneTimeSplitter(n_splits=n_splits,
                                                  n_samples=n_samples))
    # The OneTimeSplitter is a non-re-entrant cv splitter. Unless, the
    # `split` is called for each parameter, the following should produce
    # identical results for param setting 1 and param setting 2 as both have
    # the same C value.
    assert_array_almost_equal(*np.vsplit(np.hstack(scores1)[(0, 2, 1, 3), :],
                                         2))

    scores2 = validation_curve(SVC(kernel='linear', random_state=0), X, y,
                               'C', [0.1, 0.1, 0.2, 0.2],
                               cv=KFold(n_splits=n_splits, shuffle=True))

    # For scores2, compare the 1st and 2nd parameter's scores
    # (Since the C value for 1st two param setting is 0.1, they must be
    # consistent unless the train test folds differ between the param settings)
    assert_array_almost_equal(*np.vsplit(np.hstack(scores2)[(0, 2, 1, 3), :],
                                         2))

    scores3 = validation_curve(SVC(kernel='linear', random_state=0), X, y,
                               'C', [0.1, 0.1, 0.2, 0.2],
                               cv=KFold(n_splits=n_splits))

    # OneTimeSplitter is basically unshuffled KFold(n_splits=5). Sanity check.
    assert_array_almost_equal(np.array(scores3), np.array(scores1))


def test_check_is_permutation():
    rng = np.random.RandomState(0)
    p = np.arange(100)
    rng.shuffle(p)
    assert _check_is_permutation(p, 100)
    assert not _check_is_permutation(np.delete(p, 23), 100)

    p[0] = 23
    assert not _check_is_permutation(p, 100)

    # Check if the additional duplicate indices are caught
    assert not _check_is_permutation(np.hstack((p, 0)), 100)


def test_cross_val_predict_sparse_prediction():
    # check that cross_val_predict gives same result for sparse and dense input
    X, y = make_multilabel_classification(n_classes=2, n_labels=1,
                                          allow_unlabeled=False,
                                          return_indicator=True,
                                          random_state=1)
    X_sparse = csr_matrix(X)
    y_sparse = csr_matrix(y)
    classif = OneVsRestClassifier(SVC(kernel='linear'))
    preds = cross_val_predict(classif, X, y, cv=10)
    preds_sparse = cross_val_predict(classif, X_sparse, y_sparse, cv=10)
    preds_sparse = preds_sparse.toarray()
    assert_array_almost_equal(preds_sparse, preds)


def check_cross_val_predict_with_method(est):
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=0)
    classes = len(set(y))

    kfold = KFold()

    methods = ['decision_function', 'predict_proba', 'predict_log_proba']
    for method in methods:
        predictions = cross_val_predict(est, X, y, method=method)
        assert_equal(len(predictions), len(y))

        expected_predictions = np.zeros([len(y), classes])
        func = getattr(est, method)

        # Naive loop (should be same as cross_val_predict):
        for train, test in kfold.split(X, y):
            est.fit(X[train], y[train])
            expected_predictions[test] = func(X[test])

        predictions = cross_val_predict(est, X, y, method=method,
                                        cv=kfold)
        assert_array_almost_equal(expected_predictions, predictions)

        # Test alternative representations of y
        predictions_y1 = cross_val_predict(est, X, y + 1, method=method,
                                           cv=kfold)
        assert_array_equal(predictions, predictions_y1)

        predictions_y2 = cross_val_predict(est, X, y - 2, method=method,
                                           cv=kfold)
        assert_array_equal(predictions, predictions_y2)

        predictions_ystr = cross_val_predict(est, X, y.astype('str'),
                                             method=method, cv=kfold)
        assert_array_equal(predictions, predictions_ystr)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_with_method():
    check_cross_val_predict_with_method(LogisticRegression())


@pytest.mark.filterwarnings('ignore: max_iter and tol parameters')
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_method_checking():
    # Regression test for issue #9639. Tests that cross_val_predict does not
    # check estimator methods (e.g. predict_proba) before fitting
    est = SGDClassifier(loss='log', random_state=2)
    check_cross_val_predict_with_method(est)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: The default of the `iid`')
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_gridsearchcv_cross_val_predict_with_method():
    est = GridSearchCV(LogisticRegression(random_state=42),
                       {'C': [0.1, 1]},
                       cv=2)
    check_cross_val_predict_with_method(est)


def get_expected_predictions(X, y, cv, classes, est, method):

    expected_predictions = np.zeros([len(y), classes])
    func = getattr(est, method)

    for train, test in cv.split(X, y):
        est.fit(X[train], y[train])
        expected_predictions_ = func(X[test])
        # To avoid 2 dimensional indexing
        if method == 'predict_proba':
            exp_pred_test = np.zeros((len(test), classes))
        else:
            exp_pred_test = np.full((len(test), classes),
                                    np.finfo(expected_predictions.dtype).min)
        exp_pred_test[:, est.classes_] = expected_predictions_
        expected_predictions[test] = exp_pred_test

    return expected_predictions


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_cross_val_predict_class_subset():

    X = np.arange(200).reshape(100, 2)
    y = np.array([x // 10 for x in range(100)])
    classes = 10

    kfold3 = KFold(n_splits=3)
    kfold4 = KFold(n_splits=4)

    le = LabelEncoder()

    methods = ['decision_function', 'predict_proba', 'predict_log_proba']
    for method in methods:
        est = LogisticRegression()

        # Test with n_splits=3
        predictions = cross_val_predict(est, X, y, method=method,
                                        cv=kfold3)

        # Runs a naive loop (should be same as cross_val_predict):
        expected_predictions = get_expected_predictions(X, y, kfold3, classes,
                                                        est, method)
        assert_array_almost_equal(expected_predictions, predictions)

        # Test with n_splits=4
        predictions = cross_val_predict(est, X, y, method=method,
                                        cv=kfold4)
        expected_predictions = get_expected_predictions(X, y, kfold4, classes,
                                                        est, method)
        assert_array_almost_equal(expected_predictions, predictions)

        # Testing unordered labels
        y = shuffle(np.repeat(range(10), 10), random_state=0)
        predictions = cross_val_predict(est, X, y, method=method,
                                        cv=kfold3)
        y = le.fit_transform(y)
        expected_predictions = get_expected_predictions(X, y, kfold3, classes,
                                                        est, method)
        assert_array_almost_equal(expected_predictions, predictions)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_score_memmap():
    # Ensure a scalar score of memmap type is accepted
    iris = load_iris()
    X, y = iris.data, iris.target
    clf = MockClassifier()
    tf = tempfile.NamedTemporaryFile(mode='wb', delete=False)
    tf.write(b'Hello world!!!!!')
    tf.close()
    scores = np.memmap(tf.name, dtype=np.float64)
    score = np.memmap(tf.name, shape=(), mode='r', dtype=np.float64)
    try:
        cross_val_score(clf, X, y, scoring=lambda est, X, y: score)
        # non-scalar should still fail
        assert_raises(ValueError, cross_val_score, clf, X, y,
                      scoring=lambda est, X, y: scores)
    finally:
        # Best effort to release the mmap file handles before deleting the
        # backing file under Windows
        scores, score = None, None
        for _ in range(3):
            try:
                os.unlink(tf.name)
                break
            except WindowsError:
                sleep(1.)


@pytest.mark.filterwarnings('ignore: Using or importing the ABCs from')
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_permutation_test_score_pandas():
    # check permutation_test_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((Series, DataFrame))
    except ImportError:
        pass
    for TargetType, InputFeatureType in types:
        # X dataframe, y series
        iris = load_iris()
        X, y = iris.data, iris.target
        X_df, y_ser = InputFeatureType(X), TargetType(y)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        permutation_test_score(clf, X_df, y_ser)


def test_fit_and_score_failing():
    # Create a failing classifier to deliberately fail
    failing_clf = FailingClassifier(FailingClassifier.FAILING_PARAMETER)
    # dummy X data
    X = np.arange(1, 10)
    y = np.ones(9)
    fit_and_score_args = [failing_clf, X, None, dict(), None, None, 0,
                          None, None]
    # passing error score to trigger the warning message
    fit_and_score_kwargs = {'error_score': 0}
    # check if the warning message type is as expected
    assert_warns(FitFailedWarning, _fit_and_score, *fit_and_score_args,
                 **fit_and_score_kwargs)
    # since we're using FailingClassfier, our error will be the following
    error_message = "ValueError: Failing classifier failed as required"
    # the warning message we're expecting to see
    warning_message = ("Estimator fit failed. The score on this train-test "
                       "partition for these parameters will be set to %f. "
                       "Details: \n%s" % (fit_and_score_kwargs['error_score'],
                                          error_message))
    # check if the same warning is triggered
    assert_warns_message(FitFailedWarning, warning_message, _fit_and_score,
                         *fit_and_score_args, **fit_and_score_kwargs)

    # check if warning was raised, with default error_score argument
    warning_message = ("From version 0.22, errors during fit will result "
                       "in a cross validation score of NaN by default. Use "
                       "error_score='raise' if you want an exception "
                       "raised or error_score=np.nan to adopt the "
                       "behavior from version 0.22.")
    with pytest.raises(ValueError):
        assert_warns_message(FutureWarning, warning_message, _fit_and_score,
                             *fit_and_score_args)

    fit_and_score_kwargs = {'error_score': 'raise'}
    # check if exception was raised, with default error_score='raise'
    assert_raise_message(ValueError, "Failing classifier failed as required",
                         _fit_and_score, *fit_and_score_args,
                         **fit_and_score_kwargs)

    # check that functions upstream pass error_score param to _fit_and_score
    error_message = ("error_score must be the string 'raise' or a"
                     " numeric value. (Hint: if using 'raise', please"
                     " make sure that it has been spelled correctly.)")

    assert_raise_message(ValueError, error_message, cross_validate,
                         failing_clf, X, cv=3, error_score='unvalid-string')

    assert_raise_message(ValueError, error_message, cross_val_score,
                         failing_clf, X, cv=3, error_score='unvalid-string')

    assert_raise_message(ValueError, error_message, learning_curve,
                         failing_clf, X, y, cv=3, error_score='unvalid-string')

    assert_raise_message(ValueError, error_message, validation_curve,
                         failing_clf, X, y, 'parameter',
                         [FailingClassifier.FAILING_PARAMETER], cv=3,
                         error_score='unvalid-string')

    assert_equal(failing_clf.score(), 0.)  # FailingClassifier coverage


def test_fit_and_score_working():
    X, y = make_classification(n_samples=30, random_state=0)
    clf = SVC(kernel="linear", random_state=0)
    train, test = next(ShuffleSplit().split(X))
    # Test return_parameters option
    fit_and_score_args = [clf, X, y, dict(), train, test, 0]
    fit_and_score_kwargs = {'parameters': {'max_iter': 100, 'tol': 0.1},
                            'fit_params': None,
                            'return_parameters': True}
    result = _fit_and_score(*fit_and_score_args,
                            **fit_and_score_kwargs)
    assert result[-1] == fit_and_score_kwargs['parameters']


def three_params_scorer(i, j, k):
    return 3.4213


@pytest.mark.parametrize("return_train_score, scorer, expected", [
    (False, three_params_scorer,
     "[CV] .................................... , score=3.421, total=   0.0s"),
    (True, three_params_scorer,
     "[CV] ................ , score=(train=3.421, test=3.421), total=   0.0s"),
    (True, {'sc1': three_params_scorer, 'sc2': three_params_scorer},
     "[CV]  , sc1=(train=3.421, test=3.421)"
     ", sc2=(train=3.421, test=3.421), total=   0.0s")
])
def test_fit_and_score_verbosity(capsys, return_train_score, scorer, expected):
    X, y = make_classification(n_samples=30, random_state=0)
    clf = SVC(kernel="linear", random_state=0)
    train, test = next(ShuffleSplit().split(X))

    # test print without train score
    fit_and_score_args = [clf, X, y, scorer, train, test, 10, None, None]
    fit_and_score_kwargs = {'return_train_score': return_train_score}
    _fit_and_score(*fit_and_score_args, **fit_and_score_kwargs)
    out, _ = capsys.readouterr()
    assert out.split('\n')[1] == expected


def test_score():
    error_message = "scoring must return a number, got None"

    def two_params_scorer(estimator, X_test):
        return None
    fit_and_score_args = [None, None, None, two_params_scorer]
    assert_raise_message(ValueError, error_message,
                         _score, *fit_and_score_args)
