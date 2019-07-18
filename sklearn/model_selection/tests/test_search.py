"""Test the search module"""

from collections.abc import Iterable, Sized
from io import StringIO
from itertools import chain, product
import pickle
import sys
from types import GeneratorType
import re

import numpy as np
import scipy.sparse as sp
import pytest

from sklearn.utils.fixes import sp_version
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

from scipy.stats import bernoulli, expon, uniform

from sklearn.base import BaseEstimator
from sklearn.base import clone
from sklearn.exceptions import NotFittedError
from sklearn.datasets import make_classification
from sklearn.datasets import make_blobs
from sklearn.datasets import make_multilabel_classification

from sklearn.model_selection import fit_grid_point
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import LeavePGroupsOut
from sklearn.model_selection import GroupKFold
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import ParameterSampler
from sklearn.model_selection._search import BaseSearchCV

from sklearn.model_selection._validation import FitFailedWarning

from sklearn.svm import LinearSVC, SVC
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge, SGDClassifier

from sklearn.model_selection.tests.common import OneTimeSplitter


# Neither of the following two estimators inherit from BaseEstimator,
# to test hyperparameter search on user-defined classifiers.
class MockClassifier(object):
    """Dummy classifier to test the parameter search algorithms"""
    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y):
        assert len(X) == len(Y)
        self.classes_ = np.unique(Y)
        return self

    def predict(self, T):
        return T.shape[0]

    def transform(self, X):
        return X + self.foo_param

    def inverse_transform(self, X):
        return X - self.foo_param

    predict_proba = predict
    predict_log_proba = predict
    decision_function = predict

    def score(self, X=None, Y=None):
        if self.foo_param > 1:
            score = 1.
        else:
            score = 0.
        return score

    def get_params(self, deep=False):
        return {'foo_param': self.foo_param}

    def set_params(self, **params):
        self.foo_param = params['foo_param']
        return self


class LinearSVCNoScore(LinearSVC):
    """An LinearSVC classifier that has no score method."""
    @property
    def score(self):
        raise AttributeError


X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
y = np.array([1, 1, 2, 2])


def assert_grid_iter_equals_getitem(grid):
    assert list(grid) == [grid[i] for i in range(len(grid))]


@pytest.mark.parametrize(
    "input, error_type, error_message",
    [(0, TypeError, r'Parameter grid is not a dict or a list \(0\)'),
     ([{'foo': [0]}, 0], TypeError, r'Parameter grid is not a dict \(0\)'),
     ({'foo': 0}, TypeError, "Parameter grid value is not iterable "
      r"\(key='foo', value=0\)")]
)
def test_validate_parameter_grid_input(input, error_type, error_message):
    with pytest.raises(error_type, match=error_message):
        ParameterGrid(input)


def test_parameter_grid():

    # Test basic properties of ParameterGrid.
    params1 = {"foo": [1, 2, 3]}
    grid1 = ParameterGrid(params1)
    assert isinstance(grid1, Iterable)
    assert isinstance(grid1, Sized)
    assert len(grid1) == 3
    assert_grid_iter_equals_getitem(grid1)

    params2 = {"foo": [4, 2],
               "bar": ["ham", "spam", "eggs"]}
    grid2 = ParameterGrid(params2)
    assert len(grid2) == 6

    # loop to assert we can iterate over the grid multiple times
    for i in range(2):
        # tuple + chain transforms {"a": 1, "b": 2} to ("a", 1, "b", 2)
        points = set(tuple(chain(*(sorted(p.items())))) for p in grid2)
        assert (points ==
                     set(("bar", x, "foo", y)
                         for x, y in product(params2["bar"], params2["foo"])))
    assert_grid_iter_equals_getitem(grid2)

    # Special case: empty grid (useful to get default estimator settings)
    empty = ParameterGrid({})
    assert len(empty) == 1
    assert list(empty) == [{}]
    assert_grid_iter_equals_getitem(empty)
    assert_raises(IndexError, lambda: empty[1])

    has_empty = ParameterGrid([{'C': [1, 10]}, {}, {'C': [.5]}])
    assert len(has_empty) == 4
    assert list(has_empty) == [{'C': 1}, {'C': 10}, {}, {'C': .5}]
    assert_grid_iter_equals_getitem(has_empty)


def test_grid_search():
    # Test that the best estimator contains the right value for foo_param
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=3, verbose=3)
    # make sure it selects the smallest parameter in case of ties
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    grid_search.fit(X, y)
    sys.stdout = old_stdout
    assert grid_search.best_estimator_.foo_param == 2

    assert_array_equal(grid_search.cv_results_["param_foo_param"].data,
                       [1, 2, 3])

    # Smoke test the score etc:
    grid_search.score(X, y)
    grid_search.predict_proba(X)
    grid_search.decision_function(X)
    grid_search.transform(X)

    # Test exception handling on scoring
    grid_search.scoring = 'sklearn'
    assert_raises(ValueError, grid_search.fit, X, y)


def check_hyperparameter_searcher_with_fit_params(klass, **klass_kwargs):
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)
    clf = CheckingClassifier(expected_fit_params=['spam', 'eggs'])
    searcher = klass(clf, {'foo_param': [1, 2, 3]}, cv=2, **klass_kwargs)

    # The CheckingClassifier generates an assertion error if
    # a parameter is missing or has length != len(X).
    assert_raise_message(AssertionError,
                         "Expected fit parameter(s) ['eggs'] not seen.",
                         searcher.fit, X, y, spam=np.ones(10))
    assert_raise_message(AssertionError,
                         "Fit parameter spam has length 1; expected 4.",
                         searcher.fit, X, y, spam=np.ones(1),
                         eggs=np.zeros(10))
    searcher.fit(X, y, spam=np.ones(10), eggs=np.zeros(10))


def test_grid_search_with_fit_params():
    check_hyperparameter_searcher_with_fit_params(GridSearchCV,
                                                  error_score='raise')


def test_random_search_with_fit_params():
    check_hyperparameter_searcher_with_fit_params(RandomizedSearchCV, n_iter=1,
                                                  error_score='raise')


@ignore_warnings
def test_grid_search_no_score():
    # Test grid-search on classifier that has no score function.
    clf = LinearSVC(random_state=0)
    X, y = make_blobs(random_state=0, centers=2)
    Cs = [.1, 1, 10]
    clf_no_score = LinearSVCNoScore(random_state=0)
    grid_search = GridSearchCV(clf, {'C': Cs}, scoring='accuracy')
    grid_search.fit(X, y)

    grid_search_no_score = GridSearchCV(clf_no_score, {'C': Cs},
                                        scoring='accuracy')
    # smoketest grid search
    grid_search_no_score.fit(X, y)

    # check that best params are equal
    assert grid_search_no_score.best_params_ == grid_search.best_params_
    # check that we can call score and that it gives the correct result
    assert grid_search.score(X, y) == grid_search_no_score.score(X, y)

    # giving no scoring function raises an error
    grid_search_no_score = GridSearchCV(clf_no_score, {'C': Cs})
    assert_raise_message(TypeError, "no scoring", grid_search_no_score.fit,
                         [[1]])


def test_grid_search_score_method():
    X, y = make_classification(n_samples=100, n_classes=2, flip_y=.2,
                               random_state=0)
    clf = LinearSVC(random_state=0)
    grid = {'C': [.1]}

    search_no_scoring = GridSearchCV(clf, grid, scoring=None).fit(X, y)
    search_accuracy = GridSearchCV(clf, grid, scoring='accuracy').fit(X, y)
    search_no_score_method_auc = GridSearchCV(LinearSVCNoScore(), grid,
                                              scoring='roc_auc'
                                              ).fit(X, y)
    search_auc = GridSearchCV(clf, grid, scoring='roc_auc').fit(X, y)

    # Check warning only occurs in situation where behavior changed:
    # estimator requires score method to compete with scoring parameter
    score_no_scoring = search_no_scoring.score(X, y)
    score_accuracy = search_accuracy.score(X, y)
    score_no_score_auc = search_no_score_method_auc.score(X, y)
    score_auc = search_auc.score(X, y)

    # ensure the test is sane
    assert score_auc < 1.0
    assert score_accuracy < 1.0
    assert score_auc != score_accuracy

    assert_almost_equal(score_accuracy, score_no_scoring)
    assert_almost_equal(score_auc, score_no_score_auc)


def test_grid_search_groups():
    # Check if ValueError (when groups is None) propagates to GridSearchCV
    # And also check if groups is correctly passed to the cv object
    rng = np.random.RandomState(0)

    X, y = make_classification(n_samples=15, n_classes=2, random_state=0)
    groups = rng.randint(0, 3, 15)

    clf = LinearSVC(random_state=0)
    grid = {'C': [1]}

    group_cvs = [LeaveOneGroupOut(), LeavePGroupsOut(2),
                 GroupKFold(n_splits=3), GroupShuffleSplit()]
    for cv in group_cvs:
        gs = GridSearchCV(clf, grid, cv=cv)
        assert_raise_message(ValueError,
                             "The 'groups' parameter should not be None.",
                             gs.fit, X, y)
        gs.fit(X, y, groups=groups)

    non_group_cvs = [StratifiedKFold(), StratifiedShuffleSplit()]
    for cv in non_group_cvs:
        gs = GridSearchCV(clf, grid, cv=cv)
        # Should not raise an error
        gs.fit(X, y)


def test_classes__property():
    # Test that classes_ property matches best_estimator_.classes_
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)
    Cs = [.1, 1, 10]

    grid_search = GridSearchCV(LinearSVC(random_state=0), {'C': Cs})
    grid_search.fit(X, y)
    assert_array_equal(grid_search.best_estimator_.classes_,
                       grid_search.classes_)

    # Test that regressors do not have a classes_ attribute
    grid_search = GridSearchCV(Ridge(), {'alpha': [1.0, 2.0]})
    grid_search.fit(X, y)
    assert not hasattr(grid_search, 'classes_')

    # Test that the grid searcher has no classes_ attribute before it's fit
    grid_search = GridSearchCV(LinearSVC(random_state=0), {'C': Cs})
    assert not hasattr(grid_search, 'classes_')

    # Test that the grid searcher has no classes_ attribute without a refit
    grid_search = GridSearchCV(LinearSVC(random_state=0),
                               {'C': Cs}, refit=False)
    grid_search.fit(X, y)
    assert not hasattr(grid_search, 'classes_')


def test_trivial_cv_results_attr():
    # Test search over a "grid" with only one point.
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1]}, cv=3)
    grid_search.fit(X, y)
    assert hasattr(grid_search, "cv_results_")

    random_search = RandomizedSearchCV(clf, {'foo_param': [0]}, n_iter=1, cv=3)
    random_search.fit(X, y)
    assert hasattr(grid_search, "cv_results_")


def test_no_refit():
    # Test that GSCV can be used for model selection alone without refitting
    clf = MockClassifier()
    for scoring in [None, ['accuracy', 'precision']]:
        grid_search = GridSearchCV(
            clf, {'foo_param': [1, 2, 3]}, refit=False, cv=3
        )
        grid_search.fit(X, y)
        assert not hasattr(grid_search, "best_estimator_") and \
            hasattr(grid_search, "best_index_") and \
            hasattr(grid_search, "best_params_")

        # Make sure the functions predict/transform etc raise meaningful
        # error messages
        for fn_name in ('predict', 'predict_proba', 'predict_log_proba',
                        'transform', 'inverse_transform'):
            assert_raise_message(NotFittedError,
                                 ('refit=False. %s is available only after '
                                  'refitting on the best parameters'
                                  % fn_name), getattr(grid_search, fn_name), X)

    # Test that an invalid refit param raises appropriate error messages
    for refit in ["", 5, True, 'recall', 'accuracy']:
        assert_raise_message(ValueError, "For multi-metric scoring, the "
                             "parameter refit must be set to a scorer key",
                             GridSearchCV(clf, {}, refit=refit,
                                          scoring={'acc': 'accuracy',
                                                   'prec': 'precision'}
                                          ).fit,
                             X, y)


def test_grid_search_error():
    # Test that grid search will capture errors on data with different length
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, X_[:180], y_)


def test_grid_search_one_grid_point():
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)
    param_dict = {"C": [1.0], "kernel": ["rbf"], "gamma": [0.1]}

    clf = SVC(gamma='auto')
    cv = GridSearchCV(clf, param_dict)
    cv.fit(X_, y_)

    clf = SVC(C=1.0, kernel="rbf", gamma=0.1)
    clf.fit(X_, y_)

    assert_array_equal(clf.dual_coef_, cv.best_estimator_.dual_coef_)


def test_grid_search_when_param_grid_includes_range():
    # Test that the best estimator contains the right value for foo_param
    clf = MockClassifier()
    grid_search = None
    grid_search = GridSearchCV(clf, {'foo_param': range(1, 4)}, cv=3)
    grid_search.fit(X, y)
    assert grid_search.best_estimator_.foo_param == 2


def test_grid_search_bad_param_grid():
    param_dict = {"C": 1.0}
    clf = SVC(gamma='auto')
    assert_raise_message(
        ValueError,
        "Parameter values for parameter (C) need to be a sequence"
        "(but not a string) or np.ndarray.",
        GridSearchCV, clf, param_dict)

    param_dict = {"C": []}
    clf = SVC()
    assert_raise_message(
        ValueError,
        "Parameter values for parameter (C) need to be a non-empty sequence.",
        GridSearchCV, clf, param_dict)

    param_dict = {"C": "1,2,3"}
    clf = SVC(gamma='auto')
    assert_raise_message(
        ValueError,
        "Parameter values for parameter (C) need to be a sequence"
        "(but not a string) or np.ndarray.",
        GridSearchCV, clf, param_dict)

    param_dict = {"C": np.ones((3, 2))}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)


def test_grid_search_sparse():
    # Test that grid search works with both dense and sparse matrices
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(X_[:180], y_[:180])
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator_.C

    X_ = sp.csr_matrix(X_)
    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(X_[:180].tocoo(), y_[:180])
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator_.C

    assert np.mean(y_pred == y_pred2) >= .9
    assert C == C2


def test_grid_search_sparse_scoring():
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, scoring="f1")
    cv.fit(X_[:180], y_[:180])
    y_pred = cv.predict(X_[180:])
    C = cv.best_estimator_.C

    X_ = sp.csr_matrix(X_)
    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, scoring="f1")
    cv.fit(X_[:180], y_[:180])
    y_pred2 = cv.predict(X_[180:])
    C2 = cv.best_estimator_.C

    assert_array_equal(y_pred, y_pred2)
    assert C == C2
    # Smoke test the score
    # np.testing.assert_allclose(f1_score(cv.predict(X_[:180]), y[:180]),
    #                            cv.score(X_[:180], y[:180]))

    # test loss where greater is worse
    def f1_loss(y_true_, y_pred_):
        return -f1_score(y_true_, y_pred_)
    F1Loss = make_scorer(f1_loss, greater_is_better=False)
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]}, scoring=F1Loss)
    cv.fit(X_[:180], y_[:180])
    y_pred3 = cv.predict(X_[180:])
    C3 = cv.best_estimator_.C

    assert C == C3
    assert_array_equal(y_pred, y_pred3)


def test_grid_search_precomputed_kernel():
    # Test that grid search works when the input features are given in the
    # form of a precomputed kernel matrix
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    # compute the training kernel matrix corresponding to the linear kernel
    K_train = np.dot(X_[:180], X_[:180].T)
    y_train = y_[:180]

    clf = SVC(kernel='precomputed')
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    cv.fit(K_train, y_train)

    assert cv.best_score_ >= 0

    # compute the test kernel matrix
    K_test = np.dot(X_[180:], X_[:180].T)
    y_test = y_[180:]

    y_pred = cv.predict(K_test)

    assert np.mean(y_pred == y_test) >= 0

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    assert_raises(ValueError, cv.fit, K_train.tolist(), y_train)


def test_grid_search_precomputed_kernel_error_nonsquare():
    # Test that grid search returns an error with a non-square precomputed
    # training kernel matrix
    K_train = np.zeros((10, 20))
    y_train = np.ones((10, ))
    clf = SVC(kernel='precomputed')
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, K_train, y_train)


class BrokenClassifier(BaseEstimator):
    """Broken classifier that cannot be fit twice"""

    def __init__(self, parameter=None):
        self.parameter = parameter

    def fit(self, X, y):
        assert not hasattr(self, 'has_been_fit_')
        self.has_been_fit_ = True

    def predict(self, X):
        return np.zeros(X.shape[0])


@ignore_warnings
def test_refit():
    # Regression test for bug in refitting
    # Simulates re-fitting a broken estimator; this used to break with
    # sparse SVMs.
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = GridSearchCV(BrokenClassifier(), [{'parameter': [0, 1]}],
                       scoring="precision", refit=True)
    clf.fit(X, y)


def test_refit_callable():
    """
    Test refit=callable, which adds flexibility in identifying the
    "best" estimator.
    """
    def refit_callable(cv_results):
        """
        A dummy function tests `refit=callable` interface.
        Return the index of a model that has the least
        `mean_test_score`.
        """
        # Fit a dummy clf with `refit=True` to get a list of keys in
        # clf.cv_results_.
        X, y = make_classification(n_samples=100, n_features=4,
                                   random_state=42)
        clf = GridSearchCV(LinearSVC(random_state=42), {'C': [0.01, 0.1, 1]},
                           scoring='precision', refit=True)
        clf.fit(X, y)
        # Ensure that `best_index_ != 0` for this dummy clf
        assert clf.best_index_ != 0

        # Assert every key matches those in `cv_results`
        for key in clf.cv_results_.keys():
            assert key in cv_results

        return cv_results['mean_test_score'].argmin()

    X, y = make_classification(n_samples=100, n_features=4,
                               random_state=42)
    clf = GridSearchCV(LinearSVC(random_state=42), {'C': [0.01, 0.1, 1]},
                       scoring='precision', refit=refit_callable)
    clf.fit(X, y)

    assert clf.best_index_ == 0
    # Ensure `best_score_` is disabled when using `refit=callable`
    assert not hasattr(clf, 'best_score_')


def test_refit_callable_invalid_type():
    """
    Test implementation catches the errors when 'best_index_' returns an
    invalid result.
    """
    def refit_callable_invalid_type(cv_results):
        """
        A dummy function tests when returned 'best_index_' is not integer.
        """
        return None

    X, y = make_classification(n_samples=100, n_features=4,
                               random_state=42)

    clf = GridSearchCV(LinearSVC(random_state=42), {'C': [0.1, 1]},
                       scoring='precision', refit=refit_callable_invalid_type)
    with pytest.raises(TypeError,
                       match='best_index_ returned is not an integer'):
        clf.fit(X, y)


@pytest.mark.parametrize('out_bound_value', [-1, 2])
@pytest.mark.parametrize('search_cv', [RandomizedSearchCV, GridSearchCV])
def test_refit_callable_out_bound(out_bound_value, search_cv):
    """
    Test implementation catches the errors when 'best_index_' returns an
    out of bound result.
    """
    def refit_callable_out_bound(cv_results):
        """
        A dummy function tests when returned 'best_index_' is out of bounds.
        """
        return out_bound_value

    X, y = make_classification(n_samples=100, n_features=4,
                               random_state=42)

    clf = search_cv(LinearSVC(random_state=42), {'C': [0.1, 1]},
                    scoring='precision', refit=refit_callable_out_bound)
    with pytest.raises(IndexError, match='best_index_ index out of range'):
        clf.fit(X, y)


def test_refit_callable_multi_metric():
    """
    Test refit=callable in multiple metric evaluation setting
    """
    def refit_callable(cv_results):
        """
        A dummy function tests `refit=callable` interface.
        Return the index of a model that has the least
        `mean_test_prec`.
        """
        assert 'mean_test_prec' in cv_results
        return cv_results['mean_test_prec'].argmin()

    X, y = make_classification(n_samples=100, n_features=4,
                               random_state=42)
    scoring = {'Accuracy': make_scorer(accuracy_score), 'prec': 'precision'}
    clf = GridSearchCV(LinearSVC(random_state=42), {'C': [0.01, 0.1, 1]},
                       scoring=scoring, refit=refit_callable)
    clf.fit(X, y)

    assert clf.best_index_ == 0
    # Ensure `best_score_` is disabled when using `refit=callable`
    assert not hasattr(clf, 'best_score_')


def test_gridsearch_nd():
    # Pass X as list in GridSearchCV
    X_4d = np.arange(10 * 5 * 3 * 2).reshape(10, 5, 3, 2)
    y_3d = np.arange(10 * 7 * 11).reshape(10, 7, 11)
    check_X = lambda x: x.shape[1:] == (5, 3, 2)
    check_y = lambda x: x.shape[1:] == (7, 11)
    clf = CheckingClassifier(check_X=check_X, check_y=check_y)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]})
    grid_search.fit(X_4d, y_3d).score(X, y)
    assert hasattr(grid_search, "cv_results_")


def test_X_as_list():
    # Pass X as list in GridSearchCV
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = CheckingClassifier(check_X=lambda x: isinstance(x, list))
    cv = KFold(n_splits=3)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=cv)
    grid_search.fit(X.tolist(), y).score(X, y)
    assert hasattr(grid_search, "cv_results_")


def test_y_as_list():
    # Pass y as list in GridSearchCV
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = CheckingClassifier(check_y=lambda x: isinstance(x, list))
    cv = KFold(n_splits=3)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=cv)
    grid_search.fit(X, y.tolist()).score(X, y)
    assert hasattr(grid_search, "cv_results_")


@ignore_warnings
def test_pandas_input():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((DataFrame, Series))
    except ImportError:
        pass

    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    for InputFeatureType, TargetType in types:
        # X dataframe, y series
        X_df, y_ser = InputFeatureType(X), TargetType(y)

        def check_df(x):
            return isinstance(x, InputFeatureType)

        def check_series(x):
            return isinstance(x, TargetType)

        clf = CheckingClassifier(check_X=check_df, check_y=check_series)

        grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]})
        grid_search.fit(X_df, y_ser).score(X_df, y_ser)
        grid_search.predict(X_df)
        assert hasattr(grid_search, "cv_results_")


def test_unsupervised_grid_search():
    # test grid-search with unsupervised estimator
    X, y = make_blobs(n_samples=50, random_state=0)
    km = KMeans(random_state=0, init="random", n_init=1)

    # Multi-metric evaluation unsupervised
    scoring = ['adjusted_rand_score', 'fowlkes_mallows_score']
    for refit in ['adjusted_rand_score', 'fowlkes_mallows_score']:
        grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]),
                                   scoring=scoring, refit=refit)
        grid_search.fit(X, y)
        # Both ARI and FMS can find the right number :)
        assert grid_search.best_params_["n_clusters"] == 3

    # Single metric evaluation unsupervised
    grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]),
                               scoring='fowlkes_mallows_score')
    grid_search.fit(X, y)
    assert grid_search.best_params_["n_clusters"] == 3

    # Now without a score, and without y
    grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]))
    grid_search.fit(X)
    assert grid_search.best_params_["n_clusters"] == 4


def test_gridsearch_no_predict():
    # test grid-search with an estimator without predict.
    # slight duplication of a test from KDE
    def custom_scoring(estimator, X):
        return 42 if estimator.bandwidth == .1 else 0
    X, _ = make_blobs(cluster_std=.1, random_state=1,
                      centers=[[0, 1], [1, 0], [0, 0]])
    search = GridSearchCV(KernelDensity(),
                          param_grid=dict(bandwidth=[.01, .1, 1]),
                          scoring=custom_scoring)
    search.fit(X)
    assert search.best_params_['bandwidth'] == .1
    assert search.best_score_ == 42


def test_param_sampler():
    # test basic properties of param sampler
    param_distributions = {"kernel": ["rbf", "linear"],
                           "C": uniform(0, 1)}
    sampler = ParameterSampler(param_distributions=param_distributions,
                               n_iter=10, random_state=0)
    samples = [x for x in sampler]
    assert len(samples) == 10
    for sample in samples:
        assert sample["kernel"] in ["rbf", "linear"]
        assert 0 <= sample["C"] <= 1

    # test that repeated calls yield identical parameters
    param_distributions = {"C": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
    sampler = ParameterSampler(param_distributions=param_distributions,
                               n_iter=3, random_state=0)
    assert [x for x in sampler] == [x for x in sampler]

    if sp_version >= (0, 16):
        param_distributions = {"C": uniform(0, 1)}
        sampler = ParameterSampler(param_distributions=param_distributions,
                                   n_iter=10, random_state=0)
        assert [x for x in sampler] == [x for x in sampler]


def check_cv_results_array_types(search, param_keys, score_keys):
    # Check if the search `cv_results`'s array are of correct types
    cv_results = search.cv_results_
    assert all(isinstance(cv_results[param], np.ma.MaskedArray)
               for param in param_keys)
    assert all(cv_results[key].dtype == object for key in param_keys)
    assert not any(isinstance(cv_results[key], np.ma.MaskedArray)
                   for key in score_keys)
    assert all(cv_results[key].dtype == np.float64
               for key in score_keys if not key.startswith('rank'))

    scorer_keys = search.scorer_.keys() if search.multimetric_ else ['score']

    for key in scorer_keys:
        assert cv_results['rank_test_%s' % key].dtype == np.int32


def check_cv_results_keys(cv_results, param_keys, score_keys, n_cand):
    # Test the search.cv_results_ contains all the required results
    assert_array_equal(sorted(cv_results.keys()),
                       sorted(param_keys + score_keys + ('params',)))
    assert all(cv_results[key].shape == (n_cand,)
               for key in param_keys + score_keys)


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def test_grid_search_cv_results():
    X, y = make_classification(n_samples=50, n_features=4,
                               random_state=42)

    n_splits = 3
    n_grid_points = 6
    params = [dict(kernel=['rbf', ], C=[1, 10], gamma=[0.1, 1]),
              dict(kernel=['poly', ], degree=[1, 2])]

    param_keys = ('param_C', 'param_degree', 'param_gamma', 'param_kernel')
    score_keys = ('mean_test_score', 'mean_train_score',
                  'rank_test_score',
                  'split0_test_score', 'split1_test_score',
                  'split2_test_score',
                  'split0_train_score', 'split1_train_score',
                  'split2_train_score',
                  'std_test_score', 'std_train_score',
                  'mean_fit_time', 'std_fit_time',
                  'mean_score_time', 'std_score_time')
    n_candidates = n_grid_points

    for iid in (False, True):
        search = GridSearchCV(SVC(), cv=n_splits, iid=iid,
                              param_grid=params, return_train_score=True)
        search.fit(X, y)
        assert iid == search.iid
        cv_results = search.cv_results_
        # Check if score and timing are reasonable
        assert all(cv_results['rank_test_score'] >= 1)
        assert (all(cv_results[k] >= 0) for k in score_keys
                if k is not 'rank_test_score')
        assert (all(cv_results[k] <= 1) for k in score_keys
                if 'time' not in k and
                k is not 'rank_test_score')
        # Check cv_results structure
        check_cv_results_array_types(search, param_keys, score_keys)
        check_cv_results_keys(cv_results, param_keys, score_keys, n_candidates)
        # Check masking
        cv_results = search.cv_results_
        n_candidates = len(search.cv_results_['params'])
        assert all((cv_results['param_C'].mask[i] and
                    cv_results['param_gamma'].mask[i] and
                    not cv_results['param_degree'].mask[i])
                   for i in range(n_candidates)
                   if cv_results['param_kernel'][i] == 'linear')
        assert all((not cv_results['param_C'].mask[i] and
                    not cv_results['param_gamma'].mask[i] and
                    cv_results['param_degree'].mask[i])
                   for i in range(n_candidates)
                   if cv_results['param_kernel'][i] == 'rbf')


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def test_random_search_cv_results():
    X, y = make_classification(n_samples=50, n_features=4, random_state=42)

    n_splits = 3
    n_search_iter = 30

    params = dict(C=expon(scale=10), gamma=expon(scale=0.1))
    param_keys = ('param_C', 'param_gamma')
    score_keys = ('mean_test_score', 'mean_train_score',
                  'rank_test_score',
                  'split0_test_score', 'split1_test_score',
                  'split2_test_score',
                  'split0_train_score', 'split1_train_score',
                  'split2_train_score',
                  'std_test_score', 'std_train_score',
                  'mean_fit_time', 'std_fit_time',
                  'mean_score_time', 'std_score_time')
    n_cand = n_search_iter

    for iid in (False, True):
        search = RandomizedSearchCV(SVC(), n_iter=n_search_iter,
                                    cv=n_splits, iid=iid,
                                    param_distributions=params,
                                    return_train_score=True)
        search.fit(X, y)
        assert iid == search.iid
        cv_results = search.cv_results_
        # Check results structure
        check_cv_results_array_types(search, param_keys, score_keys)
        check_cv_results_keys(cv_results, param_keys, score_keys, n_cand)
        # For random_search, all the param array vals should be unmasked
        assert not(any(np.ma.getmaskarray(cv_results['param_C'])) or
                   any(np.ma.getmaskarray(cv_results['param_gamma'])))


@pytest.mark.parametrize(
    "SearchCV, specialized_params",
    [(GridSearchCV, {'param_grid': {'C': [1, 10]}}),
     (RandomizedSearchCV,
      {'param_distributions': {'C': [1, 10]}, 'n_iter': 2})]
)
def test_search_default_iid(SearchCV, specialized_params):
    # Test the IID parameter
    # noise-free simple 2d-data
    X, y = make_blobs(centers=[[0, 0], [1, 0], [0, 1], [1, 1]], random_state=0,
                      cluster_std=0.1, shuffle=False, n_samples=80)
    # split dataset into two folds that are not iid
    # first one contains data of all 4 blobs, second only from two.
    mask = np.ones(X.shape[0], dtype=np.bool)
    mask[np.where(y == 1)[0][::2]] = 0
    mask[np.where(y == 2)[0][::2]] = 0
    # this leads to perfect classification on one fold and a score of 1/3 on
    # the other
    # create "cv" for splits
    cv = [[mask, ~mask], [~mask, mask]]

    common_params = {'estimator': SVC(), 'cv': cv,
                     'return_train_score': True}
    search = SearchCV(**common_params, **specialized_params)
    search.fit(X, y)

    test_cv_scores = np.array(
        [search.cv_results_['split%d_test_score' % s][0]
         for s in range(search.n_splits_)]
    )
    test_mean = search.cv_results_['mean_test_score'][0]
    test_std = search.cv_results_['std_test_score'][0]

    train_cv_scores = np.array(
        [search.cv_results_['split%d_train_score' % s][0]
         for s in range(search.n_splits_)]
    )
    train_mean = search.cv_results_['mean_train_score'][0]
    train_std = search.cv_results_['std_train_score'][0]

    assert search.cv_results_['param_C'][0] == 1
    # scores are the same as above
    assert_allclose(test_cv_scores, [1, 1. / 3.])
    assert_allclose(train_cv_scores, [1, 1])
    # Unweighted mean/std is used
    assert test_mean == pytest.approx(np.mean(test_cv_scores))
    assert test_std == pytest.approx(np.std(test_cv_scores))

    # For the train scores, we do not take a weighted mean irrespective of
    # i.i.d. or not
    assert train_mean == pytest.approx(1)
    assert train_std == pytest.approx(0)


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def test_search_iid_param():
    # Test the IID parameter
    # noise-free simple 2d-data
    X, y = make_blobs(centers=[[0, 0], [1, 0], [0, 1], [1, 1]], random_state=0,
                      cluster_std=0.1, shuffle=False, n_samples=80)
    # split dataset into two folds that are not iid
    # first one contains data of all 4 blobs, second only from two.
    mask = np.ones(X.shape[0], dtype=np.bool)
    mask[np.where(y == 1)[0][::2]] = 0
    mask[np.where(y == 2)[0][::2]] = 0
    # this leads to perfect classification on one fold and a score of 1/3 on
    # the other
    # create "cv" for splits
    cv = [[mask, ~mask], [~mask, mask]]
    # once with iid=True (default)
    grid_search = GridSearchCV(SVC(gamma='auto'), param_grid={'C': [1, 10]},
                               cv=cv, return_train_score=True, iid=True)
    random_search = RandomizedSearchCV(SVC(gamma='auto'), n_iter=2,
                                       param_distributions={'C': [1, 10]},
                                       cv=cv, iid=True,
                                       return_train_score=True)
    for search in (grid_search, random_search):
        search.fit(X, y)
        assert search.iid or search.iid is None

        test_cv_scores = np.array(list(search.cv_results_['split%d_test_score'
                                                          % s_i][0]
                                       for s_i in range(search.n_splits_)))
        test_mean = search.cv_results_['mean_test_score'][0]
        test_std = search.cv_results_['std_test_score'][0]

        train_cv_scores = np.array(list(search.cv_results_['split%d_train_'
                                                           'score' % s_i][0]
                                        for s_i in range(search.n_splits_)))
        train_mean = search.cv_results_['mean_train_score'][0]
        train_std = search.cv_results_['std_train_score'][0]

        # Test the first candidate
        assert search.cv_results_['param_C'][0] == 1
        assert_array_almost_equal(test_cv_scores, [1, 1. / 3.])
        assert_array_almost_equal(train_cv_scores, [1, 1])

        # for first split, 1/4 of dataset is in test, for second 3/4.
        # take weighted average and weighted std
        expected_test_mean = 1 * 1. / 4. + 1. / 3. * 3. / 4.
        expected_test_std = np.sqrt(1. / 4 * (expected_test_mean - 1) ** 2 +
                                    3. / 4 * (expected_test_mean - 1. / 3.) **
                                    2)
        assert_almost_equal(test_mean, expected_test_mean)
        assert_almost_equal(test_std, expected_test_std)
        assert_array_almost_equal(test_cv_scores,
                                  cross_val_score(SVC(C=1, gamma='auto'), X,
                                                  y, cv=cv))

        # For the train scores, we do not take a weighted mean irrespective of
        # i.i.d. or not
        assert_almost_equal(train_mean, 1)
        assert_almost_equal(train_std, 0)

    # once with iid=False
    grid_search = GridSearchCV(SVC(gamma='auto'),
                               param_grid={'C': [1, 10]},
                               cv=cv, iid=False, return_train_score=True)
    random_search = RandomizedSearchCV(SVC(gamma='auto'), n_iter=2,
                                       param_distributions={'C': [1, 10]},
                                       cv=cv, iid=False,
                                       return_train_score=True)

    for search in (grid_search, random_search):
        search.fit(X, y)
        assert not search.iid

        test_cv_scores = np.array(list(search.cv_results_['split%d_test_score'
                                                          % s][0]
                                       for s in range(search.n_splits_)))
        test_mean = search.cv_results_['mean_test_score'][0]
        test_std = search.cv_results_['std_test_score'][0]

        train_cv_scores = np.array(list(search.cv_results_['split%d_train_'
                                                           'score' % s][0]
                                        for s in range(search.n_splits_)))
        train_mean = search.cv_results_['mean_train_score'][0]
        train_std = search.cv_results_['std_train_score'][0]

        assert search.cv_results_['param_C'][0] == 1
        # scores are the same as above
        assert_array_almost_equal(test_cv_scores, [1, 1. / 3.])
        # Unweighted mean/std is used
        assert_almost_equal(test_mean, np.mean(test_cv_scores))
        assert_almost_equal(test_std, np.std(test_cv_scores))

        # For the train scores, we do not take a weighted mean irrespective of
        # i.i.d. or not
        assert_almost_equal(train_mean, 1)
        assert_almost_equal(train_std, 0)


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def test_grid_search_cv_results_multimetric():
    X, y = make_classification(n_samples=50, n_features=4, random_state=42)

    n_splits = 3
    params = [dict(kernel=['rbf', ], C=[1, 10], gamma=[0.1, 1]),
              dict(kernel=['poly', ], degree=[1, 2])]

    for iid in (False, True):
        grid_searches = []
        for scoring in ({'accuracy': make_scorer(accuracy_score),
                         'recall': make_scorer(recall_score)},
                        'accuracy', 'recall'):
            grid_search = GridSearchCV(SVC(), cv=n_splits,
                                       iid=iid, param_grid=params,
                                       scoring=scoring, refit=False)
            grid_search.fit(X, y)
            assert grid_search.iid == iid
            grid_searches.append(grid_search)

        compare_cv_results_multimetric_with_single(*grid_searches, iid=iid)


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def test_random_search_cv_results_multimetric():
    X, y = make_classification(n_samples=50, n_features=4, random_state=42)

    n_splits = 3
    n_search_iter = 30
    scoring = ('accuracy', 'recall')

    # Scipy 0.12's stats dists do not accept seed, hence we use param grid
    params = dict(C=np.logspace(-4, 1, 3),
                  gamma=np.logspace(-5, 0, 3, base=0.1))
    for iid in (True, False):
        for refit in (True, False):
            random_searches = []
            for scoring in (('accuracy', 'recall'), 'accuracy', 'recall'):
                # If True, for multi-metric pass refit='accuracy'
                if refit:
                    probability = True
                    refit = 'accuracy' if isinstance(scoring, tuple) else refit
                else:
                    probability = False
                clf = SVC(probability=probability, random_state=42)
                random_search = RandomizedSearchCV(clf, n_iter=n_search_iter,
                                                   cv=n_splits, iid=iid,
                                                   param_distributions=params,
                                                   scoring=scoring,
                                                   refit=refit, random_state=0)
                random_search.fit(X, y)
                random_searches.append(random_search)

            compare_cv_results_multimetric_with_single(*random_searches,
                                                       iid=iid)
            if refit:
                compare_refit_methods_when_refit_with_acc(
                    random_searches[0], random_searches[1], refit)


@pytest.mark.filterwarnings("ignore:The parameter 'iid' is deprecated")  # 0.24
def compare_cv_results_multimetric_with_single(
        search_multi, search_acc, search_rec, iid):
    """Compare multi-metric cv_results with the ensemble of multiple
    single metric cv_results from single metric grid/random search"""

    assert search_multi.iid == iid
    assert search_multi.multimetric_
    assert_array_equal(sorted(search_multi.scorer_),
                       ('accuracy', 'recall'))

    cv_results_multi = search_multi.cv_results_
    cv_results_acc_rec = {re.sub('_score$', '_accuracy', k): v
                          for k, v in search_acc.cv_results_.items()}
    cv_results_acc_rec.update({re.sub('_score$', '_recall', k): v
                               for k, v in search_rec.cv_results_.items()})

    # Check if score and timing are reasonable, also checks if the keys
    # are present
    assert all((np.all(cv_results_multi[k] <= 1) for k in (
                    'mean_score_time', 'std_score_time', 'mean_fit_time',
                    'std_fit_time')))

    # Compare the keys, other than time keys, among multi-metric and
    # single metric grid search results. np.testing.assert_equal performs a
    # deep nested comparison of the two cv_results dicts
    np.testing.assert_equal({k: v for k, v in cv_results_multi.items()
                             if not k.endswith('_time')},
                            {k: v for k, v in cv_results_acc_rec.items()
                             if not k.endswith('_time')})


def compare_refit_methods_when_refit_with_acc(search_multi, search_acc, refit):
    """Compare refit multi-metric search methods with single metric methods"""
    if refit:
        assert search_multi.refit == 'accuracy'
    else:
        assert not search_multi.refit
    assert search_acc.refit == refit

    X, y = make_blobs(n_samples=100, n_features=4, random_state=42)
    for method in ('predict', 'predict_proba', 'predict_log_proba'):
        assert_almost_equal(getattr(search_multi, method)(X),
                            getattr(search_acc, method)(X))
    assert_almost_equal(search_multi.score(X, y), search_acc.score(X, y))
    for key in ('best_index_', 'best_score_', 'best_params_'):
        assert getattr(search_multi, key) == getattr(search_acc, key)


def test_search_cv_results_rank_tie_breaking():
    X, y = make_blobs(n_samples=50, random_state=42)

    # The two C values are close enough to give similar models
    # which would result in a tie of their mean cv-scores
    param_grid = {'C': [1, 1.001, 0.001]}

    grid_search = GridSearchCV(SVC(), param_grid=param_grid,
                               return_train_score=True)
    random_search = RandomizedSearchCV(SVC(), n_iter=3,
                                       param_distributions=param_grid,
                                       return_train_score=True)

    for search in (grid_search, random_search):
        search.fit(X, y)
        cv_results = search.cv_results_
        # Check tie breaking strategy -
        # Check that there is a tie in the mean scores between
        # candidates 1 and 2 alone
        assert_almost_equal(cv_results['mean_test_score'][0],
                            cv_results['mean_test_score'][1])
        assert_almost_equal(cv_results['mean_train_score'][0],
                            cv_results['mean_train_score'][1])
        assert not np.allclose(cv_results['mean_test_score'][1],
                               cv_results['mean_test_score'][2])
        assert not np.allclose(cv_results['mean_train_score'][1],
                               cv_results['mean_train_score'][2])
        # 'min' rank should be assigned to the tied candidates
        assert_almost_equal(search.cv_results_['rank_test_score'], [1, 1, 3])


def test_search_cv_results_none_param():
    X, y = [[1], [2], [3], [4], [5]], [0, 0, 0, 0, 1]
    estimators = (DecisionTreeRegressor(), DecisionTreeClassifier())
    est_parameters = {"random_state": [0, None]}
    cv = KFold(random_state=0)

    for est in estimators:
        grid_search = GridSearchCV(est, est_parameters, cv=cv,
                                   ).fit(X, y)
        assert_array_equal(grid_search.cv_results_['param_random_state'],
                           [0, None])


@ignore_warnings()
def test_search_cv_timing():
    svc = LinearSVC(random_state=0)

    X = [[1, ], [2, ], [3, ], [4, ]]
    y = [0, 1, 1, 0]

    gs = GridSearchCV(svc, {'C': [0, 1]}, cv=2, error_score=0)
    rs = RandomizedSearchCV(svc, {'C': [0, 1]}, cv=2, error_score=0, n_iter=2)

    for search in (gs, rs):
        search.fit(X, y)
        for key in ['mean_fit_time', 'std_fit_time']:
            # NOTE The precision of time.time in windows is not high
            # enough for the fit/score times to be non-zero for trivial X and y
            assert np.all(search.cv_results_[key] >= 0)
            assert np.all(search.cv_results_[key] < 1)

        for key in ['mean_score_time', 'std_score_time']:
            assert search.cv_results_[key][1] >= 0
            assert search.cv_results_[key][0] == 0.0
            assert np.all(search.cv_results_[key] < 1)

        assert hasattr(search, "refit_time_")
        assert isinstance(search.refit_time_, float)
        assert search.refit_time_ >= 0


def test_grid_search_correct_score_results():
    # test that correct scores are used
    n_splits = 3
    clf = LinearSVC(random_state=0)
    X, y = make_blobs(random_state=0, centers=2)
    Cs = [.1, 1, 10]
    for score in ['f1', 'roc_auc']:
        grid_search = GridSearchCV(clf, {'C': Cs}, scoring=score, cv=n_splits)
        cv_results = grid_search.fit(X, y).cv_results_

        # Test scorer names
        result_keys = list(cv_results.keys())
        expected_keys = (("mean_test_score", "rank_test_score") +
                         tuple("split%d_test_score" % cv_i
                               for cv_i in range(n_splits)))
        assert all(np.in1d(expected_keys, result_keys))

        cv = StratifiedKFold(n_splits=n_splits)
        n_splits = grid_search.n_splits_
        for candidate_i, C in enumerate(Cs):
            clf.set_params(C=C)
            cv_scores = np.array(
                list(grid_search.cv_results_['split%d_test_score'
                                             % s][candidate_i]
                     for s in range(n_splits)))
            for i, (train, test) in enumerate(cv.split(X, y)):
                clf.fit(X[train], y[train])
                if score == "f1":
                    correct_score = f1_score(y[test], clf.predict(X[test]))
                elif score == "roc_auc":
                    dec = clf.decision_function(X[test])
                    correct_score = roc_auc_score(y[test], dec)
                assert_almost_equal(correct_score, cv_scores[i])


def test_fit_grid_point():
    X, y = make_classification(random_state=0)
    cv = StratifiedKFold(random_state=0)
    svc = LinearSVC(random_state=0)
    scorer = make_scorer(accuracy_score)

    for params in ({'C': 0.1}, {'C': 0.01}, {'C': 0.001}):
        for train, test in cv.split(X, y):
            this_scores, this_params, n_test_samples = fit_grid_point(
                X, y, clone(svc), params, train, test,
                scorer, verbose=False)

            est = clone(svc).set_params(**params)
            est.fit(X[train], y[train])
            expected_score = scorer(est, X[test], y[test])

            # Test the return values of fit_grid_point
            assert_almost_equal(this_scores, expected_score)
            assert params == this_params
            assert n_test_samples == test.size

    # Should raise an error upon multimetric scorer
    assert_raise_message(ValueError, "For evaluating multiple scores, use "
                         "sklearn.model_selection.cross_validate instead.",
                         fit_grid_point, X, y, svc, params, train, test,
                         {'score': scorer}, verbose=True)


def test_pickle():
    # Test that a fit search can be pickled
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, refit=True, cv=3)
    grid_search.fit(X, y)
    grid_search_pickled = pickle.loads(pickle.dumps(grid_search))
    assert_array_almost_equal(grid_search.predict(X),
                              grid_search_pickled.predict(X))

    random_search = RandomizedSearchCV(clf, {'foo_param': [1, 2, 3]},
                                       refit=True, n_iter=3, cv=3)
    random_search.fit(X, y)
    random_search_pickled = pickle.loads(pickle.dumps(random_search))
    assert_array_almost_equal(random_search.predict(X),
                              random_search_pickled.predict(X))


@pytest.mark.filterwarnings('ignore: The default value of multioutput')  # 0.23
def test_grid_search_with_multioutput_data():
    # Test search with multi-output estimator

    X, y = make_multilabel_classification(return_indicator=True,
                                          random_state=0)

    est_parameters = {"max_depth": [1, 2, 3, 4]}
    cv = KFold(random_state=0)

    estimators = [DecisionTreeRegressor(random_state=0),
                  DecisionTreeClassifier(random_state=0)]

    # Test with grid search cv
    for est in estimators:
        grid_search = GridSearchCV(est, est_parameters, cv=cv)
        grid_search.fit(X, y)
        res_params = grid_search.cv_results_['params']
        for cand_i in range(len(res_params)):
            est.set_params(**res_params[cand_i])

            for i, (train, test) in enumerate(cv.split(X, y)):
                est.fit(X[train], y[train])
                correct_score = est.score(X[test], y[test])
                assert_almost_equal(
                    correct_score,
                    grid_search.cv_results_['split%d_test_score' % i][cand_i])

    # Test with a randomized search
    for est in estimators:
        random_search = RandomizedSearchCV(est, est_parameters,
                                           cv=cv, n_iter=3)
        random_search.fit(X, y)
        res_params = random_search.cv_results_['params']
        for cand_i in range(len(res_params)):
            est.set_params(**res_params[cand_i])

            for i, (train, test) in enumerate(cv.split(X, y)):
                est.fit(X[train], y[train])
                correct_score = est.score(X[test], y[test])
                assert_almost_equal(
                    correct_score,
                    random_search.cv_results_['split%d_test_score'
                                              % i][cand_i])


def test_predict_proba_disabled():
    # Test predict_proba when disabled on estimator.
    X = np.arange(20).reshape(5, -1)
    y = [0, 0, 1, 1, 1]
    clf = SVC(probability=False)
    gs = GridSearchCV(clf, {}, cv=2).fit(X, y)
    assert not hasattr(gs, "predict_proba")


def test_grid_search_allows_nans():
    # Test GridSearchCV with SimpleImputer
    X = np.arange(20, dtype=np.float64).reshape(5, -1)
    X[2, :] = np.nan
    y = [0, 0, 1, 1, 1]
    p = Pipeline([
        ('imputer', SimpleImputer(strategy='mean', missing_values=np.nan)),
        ('classifier', MockClassifier()),
    ])
    GridSearchCV(p, {'classifier__foo_param': [1, 2, 3]}, cv=2).fit(X, y)


class FailingClassifier(BaseEstimator):
    """Classifier that raises a ValueError on fit()"""

    FAILING_PARAMETER = 2

    def __init__(self, parameter=None):
        self.parameter = parameter

    def fit(self, X, y=None):
        if self.parameter == FailingClassifier.FAILING_PARAMETER:
            raise ValueError("Failing classifier failed as required")

    def predict(self, X):
        return np.zeros(X.shape[0])

    def score(self, X=None, Y=None):
        return 0.


def test_grid_search_failing_classifier():
    # GridSearchCV with on_error != 'raise'
    # Ensures that a warning is raised and score reset where appropriate.

    X, y = make_classification(n_samples=20, n_features=10, random_state=0)

    clf = FailingClassifier()

    # refit=False because we only want to check that errors caused by fits
    # to individual folds will be caught and warnings raised instead. If
    # refit was done, then an exception would be raised on refit and not
    # caught by grid_search (expected behavior), and this would cause an
    # error in this test.
    gs = GridSearchCV(clf, [{'parameter': [0, 1, 2]}], scoring='accuracy',
                      refit=False, error_score=0.0)
    assert_warns(FitFailedWarning, gs.fit, X, y)
    n_candidates = len(gs.cv_results_['params'])

    # Ensure that grid scores were set to zero as required for those fits
    # that are expected to fail.
    def get_cand_scores(i):
        return np.array(list(gs.cv_results_['split%d_test_score' % s][i]
                             for s in range(gs.n_splits_)))

    assert all((np.all(get_cand_scores(cand_i) == 0.0)
                for cand_i in range(n_candidates)
                if gs.cv_results_['param_parameter'][cand_i] ==
                FailingClassifier.FAILING_PARAMETER))

    gs = GridSearchCV(clf, [{'parameter': [0, 1, 2]}], scoring='accuracy',
                      refit=False, error_score=float('nan'))
    assert_warns(FitFailedWarning, gs.fit, X, y)
    n_candidates = len(gs.cv_results_['params'])
    assert all(np.all(np.isnan(get_cand_scores(cand_i)))
               for cand_i in range(n_candidates)
               if gs.cv_results_['param_parameter'][cand_i] ==
               FailingClassifier.FAILING_PARAMETER)

    ranks = gs.cv_results_['rank_test_score']

    # Check that succeeded estimators have lower ranks
    assert ranks[0] <= 2 and ranks[1] <= 2
    # Check that failed estimator has the highest rank
    assert ranks[clf.FAILING_PARAMETER] == 3
    assert gs.best_index_ != clf.FAILING_PARAMETER


def test_grid_search_failing_classifier_raise():
    # GridSearchCV with on_error == 'raise' raises the error

    X, y = make_classification(n_samples=20, n_features=10, random_state=0)

    clf = FailingClassifier()

    # refit=False because we want to test the behaviour of the grid search part
    gs = GridSearchCV(clf, [{'parameter': [0, 1, 2]}], scoring='accuracy',
                      refit=False, error_score='raise')

    # FailingClassifier issues a ValueError so this is what we look for.
    assert_raises(ValueError, gs.fit, X, y)


def test_parameters_sampler_replacement():
    # raise warning if n_iter is bigger than total parameter space
    params = {'first': [0, 1], 'second': ['a', 'b', 'c']}
    sampler = ParameterSampler(params, n_iter=7)
    n_iter = 7
    grid_size = 6
    expected_warning = ('The total space of parameters %d is smaller '
                        'than n_iter=%d. Running %d iterations. For '
                        'exhaustive searches, use GridSearchCV.'
                        % (grid_size, n_iter, grid_size))
    assert_warns_message(UserWarning, expected_warning,
                         list, sampler)

    # degenerates to GridSearchCV if n_iter the same as grid_size
    sampler = ParameterSampler(params, n_iter=6)
    samples = list(sampler)
    assert len(samples) == 6
    for values in ParameterGrid(params):
        assert values in samples

    # test sampling without replacement in a large grid
    params = {'a': range(10), 'b': range(10), 'c': range(10)}
    sampler = ParameterSampler(params, n_iter=99, random_state=42)
    samples = list(sampler)
    assert len(samples) == 99
    hashable_samples = ["a%db%dc%d" % (p['a'], p['b'], p['c'])
                        for p in samples]
    assert len(set(hashable_samples)) == 99

    # doesn't go into infinite loops
    params_distribution = {'first': bernoulli(.5), 'second': ['a', 'b', 'c']}
    sampler = ParameterSampler(params_distribution, n_iter=7)
    samples = list(sampler)
    assert len(samples) == 7


def test_stochastic_gradient_loss_param():
    # Make sure the predict_proba works when loss is specified
    # as one of the parameters in the param_grid.
    param_grid = {
        'loss': ['log'],
    }
    X = np.arange(24).reshape(6, -1)
    y = [0, 0, 0, 1, 1, 1]
    clf = GridSearchCV(estimator=SGDClassifier(tol=1e-3, loss='hinge'),
                       param_grid=param_grid, cv=3)

    # When the estimator is not fitted, `predict_proba` is not available as the
    # loss is 'hinge'.
    assert not hasattr(clf, "predict_proba")
    clf.fit(X, y)
    clf.predict_proba(X)
    clf.predict_log_proba(X)

    # Make sure `predict_proba` is not available when setting loss=['hinge']
    # in param_grid
    param_grid = {
        'loss': ['hinge'],
    }
    clf = GridSearchCV(estimator=SGDClassifier(tol=1e-3, loss='hinge'),
                       param_grid=param_grid, cv=3)
    assert not hasattr(clf, "predict_proba")
    clf.fit(X, y)
    assert not hasattr(clf, "predict_proba")


def test_search_train_scores_set_to_false():
    X = np.arange(6).reshape(6, -1)
    y = [0, 0, 0, 1, 1, 1]
    clf = LinearSVC(random_state=0)

    gs = GridSearchCV(clf, param_grid={'C': [0.1, 0.2]}, cv=3)
    gs.fit(X, y)


def test_grid_search_cv_splits_consistency():
    # Check if a one time iterable is accepted as a cv parameter.
    n_samples = 100
    n_splits = 5
    X, y = make_classification(n_samples=n_samples, random_state=0)

    gs = GridSearchCV(LinearSVC(random_state=0),
                      param_grid={'C': [0.1, 0.2, 0.3]},
                      cv=OneTimeSplitter(n_splits=n_splits,
                                         n_samples=n_samples),
                      return_train_score=True)
    gs.fit(X, y)

    gs2 = GridSearchCV(LinearSVC(random_state=0),
                       param_grid={'C': [0.1, 0.2, 0.3]},
                       cv=KFold(n_splits=n_splits), return_train_score=True)
    gs2.fit(X, y)

    # Give generator as a cv parameter
    assert isinstance(KFold(n_splits=n_splits,
                            shuffle=True, random_state=0).split(X, y),
                      GeneratorType)
    gs3 = GridSearchCV(LinearSVC(random_state=0),
                       param_grid={'C': [0.1, 0.2, 0.3]},
                       cv=KFold(n_splits=n_splits, shuffle=True,
                                random_state=0).split(X, y),
                       return_train_score=True)
    gs3.fit(X, y)

    gs4 = GridSearchCV(LinearSVC(random_state=0),
                       param_grid={'C': [0.1, 0.2, 0.3]},
                       cv=KFold(n_splits=n_splits, shuffle=True,
                                random_state=0), return_train_score=True)
    gs4.fit(X, y)

    def _pop_time_keys(cv_results):
        for key in ('mean_fit_time', 'std_fit_time',
                    'mean_score_time', 'std_score_time'):
            cv_results.pop(key)
        return cv_results

    # Check if generators are supported as cv and
    # that the splits are consistent
    np.testing.assert_equal(_pop_time_keys(gs3.cv_results_),
                            _pop_time_keys(gs4.cv_results_))

    # OneTimeSplitter is a non-re-entrant cv where split can be called only
    # once if ``cv.split`` is called once per param setting in GridSearchCV.fit
    # the 2nd and 3rd parameter will not be evaluated as no train/test indices
    # will be generated for the 2nd and subsequent cv.split calls.
    # This is a check to make sure cv.split is not called once per param
    # setting.
    np.testing.assert_equal({k: v for k, v in gs.cv_results_.items()
                             if not k.endswith('_time')},
                            {k: v for k, v in gs2.cv_results_.items()
                             if not k.endswith('_time')})

    # Check consistency of folds across the parameters
    gs = GridSearchCV(LinearSVC(random_state=0),
                      param_grid={'C': [0.1, 0.1, 0.2, 0.2]},
                      cv=KFold(n_splits=n_splits, shuffle=True),
                      return_train_score=True)
    gs.fit(X, y)

    # As the first two param settings (C=0.1) and the next two param
    # settings (C=0.2) are same, the test and train scores must also be
    # same as long as the same train/test indices are generated for all
    # the cv splits, for both param setting
    for score_type in ('train', 'test'):
        per_param_scores = {}
        for param_i in range(4):
            per_param_scores[param_i] = list(
                gs.cv_results_['split%d_%s_score' % (s, score_type)][param_i]
                for s in range(5))

        assert_array_almost_equal(per_param_scores[0],
                                  per_param_scores[1])
        assert_array_almost_equal(per_param_scores[2],
                                  per_param_scores[3])


def test_transform_inverse_transform_round_trip():
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=3, verbose=3)

    grid_search.fit(X, y)
    X_round_trip = grid_search.inverse_transform(grid_search.transform(X))
    assert_array_equal(X, X_round_trip)


def test_custom_run_search():
    def check_results(results, gscv):
        exp_results = gscv.cv_results_
        assert sorted(results.keys()) == sorted(exp_results)
        for k in results:
            if not k.endswith('_time'):
                # XXX: results['params'] is a list :|
                results[k] = np.asanyarray(results[k])
                if results[k].dtype.kind == 'O':
                    assert_array_equal(exp_results[k], results[k],
                                       err_msg='Checking ' + k)
                else:
                    assert_allclose(exp_results[k], results[k],
                                    err_msg='Checking ' + k)

    def fit_grid(param_grid):
        return GridSearchCV(clf, param_grid,
                            return_train_score=True).fit(X, y)

    class CustomSearchCV(BaseSearchCV):
        def __init__(self, estimator, **kwargs):
            super().__init__(estimator, **kwargs)

        def _run_search(self, evaluate):
            results = evaluate([{'max_depth': 1}, {'max_depth': 2}])
            check_results(results, fit_grid({'max_depth': [1, 2]}))
            results = evaluate([{'min_samples_split': 5},
                                {'min_samples_split': 10}])
            check_results(results, fit_grid([{'max_depth': [1, 2]},
                                             {'min_samples_split': [5, 10]}]))

    # Using regressor to make sure each score differs
    clf = DecisionTreeRegressor(random_state=0)
    X, y = make_classification(n_samples=100, n_informative=4,
                               random_state=0)
    mycv = CustomSearchCV(clf, return_train_score=True).fit(X, y)
    gscv = fit_grid([{'max_depth': [1, 2]},
                     {'min_samples_split': [5, 10]}])

    results = mycv.cv_results_
    check_results(results, gscv)
    for attr in dir(gscv):
        if attr[0].islower() and attr[-1:] == '_' and \
           attr not in {'cv_results_', 'best_estimator_',
                        'refit_time_'}:
            assert getattr(gscv, attr) == getattr(mycv, attr), \
                   "Attribute %s not equal" % attr


def test__custom_fit_no_run_search():
    class NoRunSearchSearchCV(BaseSearchCV):
        def __init__(self, estimator, **kwargs):
            super().__init__(estimator, **kwargs)

        def fit(self, X, y=None, groups=None, **fit_params):
            return self

    # this should not raise any exceptions
    NoRunSearchSearchCV(SVC()).fit(X, y)

    class BadSearchCV(BaseSearchCV):
        def __init__(self, estimator, **kwargs):
            super().__init__(estimator, **kwargs)

    with pytest.raises(NotImplementedError,
                       match="_run_search not implemented."):
        # this should raise a NotImplementedError
        BadSearchCV(SVC()).fit(X, y)


@pytest.mark.parametrize("iid", [False, True])
def test_deprecated_grid_search_iid(iid):
    # FIXME: remove in 0.24
    depr_msg = "The parameter 'iid' is deprecated in 0.22 and will be removed"
    X, y = make_blobs(n_samples=54, random_state=0, centers=2)
    grid = GridSearchCV(
        SVC(random_state=0), param_grid={'C': [10]}, cv=3, iid=iid
    )
    with pytest.warns(DeprecationWarning, match=depr_msg):
        grid.fit(X, y)


def test_empty_cv_iterator_error():
    # Use global X, y

    # create cv
    cv = KFold(n_splits=3).split(X)

    # pop all of it, this should cause the expected ValueError
    [u for u in cv]
    # cv is empty now

    train_size = 100
    ridge = RandomizedSearchCV(Ridge(), {'alpha': [1e-3, 1e-2, 1e-1]},
                               cv=cv, n_jobs=4)

    # assert that this raises an error
    with pytest.raises(ValueError,
                       match='No fits were performed. '
                             'Was the CV iterator empty\\? '
                             'Were there no candidates\\?'):
        ridge.fit(X[:train_size], y[:train_size])


def test_random_search_bad_cv():
    # Use global X, y

    class BrokenKFold(KFold):
        def get_n_splits(self, *args, **kw):
            return 1

    # create bad cv
    cv = BrokenKFold(n_splits=3)

    train_size = 100
    ridge = RandomizedSearchCV(Ridge(), {'alpha': [1e-3, 1e-2, 1e-1]},
                               cv=cv, n_jobs=4)

    # assert that this raises an error
    with pytest.raises(ValueError,
                       match='cv.split and cv.get_n_splits returned '
                             'inconsistent results. Expected \\d+ '
                             'splits, got \\d+'):
        ridge.fit(X[:train_size], y[:train_size])
