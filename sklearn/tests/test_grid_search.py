"""
Testing for grid search module (sklearn.grid_search)

"""

from collections import Iterable, Sized
from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.externals.six.moves import xrange
from itertools import chain, product
import pickle
import sys

import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_false, assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

from scipy.stats import bernoulli, expon, uniform

from sklearn.externals.six.moves import zip
from sklearn.base import BaseEstimator
from sklearn.datasets import make_classification
from sklearn.datasets import make_blobs
from sklearn.datasets import make_multilabel_classification
from sklearn.grid_search import (GridSearchCV, RandomizedSearchCV,
                                 ParameterGrid, ParameterSampler,
                                 ChangedBehaviorWarning)
from sklearn.svm import LinearSVC, SVC
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn.cluster import KMeans
from sklearn.neighbors import KernelDensity
from sklearn.metrics import f1_score
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
from sklearn.cross_validation import KFold, StratifiedKFold, FitFailedWarning
from sklearn.preprocessing import Imputer
from sklearn.pipeline import Pipeline


# Neither of the following two estimators inherit from BaseEstimator,
# to test hyperparameter search on user-defined classifiers.
class MockClassifier(object):
    """Dummy classifier to test the cross-validation"""
    def __init__(self, foo_param=0):
        self.foo_param = foo_param

    def fit(self, X, Y):
        assert_true(len(X) == len(Y))
        return self

    def predict(self, T):
        return T.shape[0]

    predict_proba = predict
    decision_function = predict
    transform = predict

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


def test_parameter_grid():
    # Test basic properties of ParameterGrid.
    params1 = {"foo": [1, 2, 3]}
    grid1 = ParameterGrid(params1)
    assert_true(isinstance(grid1, Iterable))
    assert_true(isinstance(grid1, Sized))
    assert_equal(len(grid1), 3)

    params2 = {"foo": [4, 2],
               "bar": ["ham", "spam", "eggs"]}
    grid2 = ParameterGrid(params2)
    assert_equal(len(grid2), 6)

    # loop to assert we can iterate over the grid multiple times
    for i in xrange(2):
        # tuple + chain transforms {"a": 1, "b": 2} to ("a", 1, "b", 2)
        points = set(tuple(chain(*(sorted(p.items())))) for p in grid2)
        assert_equal(points,
                     set(("bar", x, "foo", y)
                         for x, y in product(params2["bar"], params2["foo"])))

    # Special case: empty grid (useful to get default estimator settings)
    empty = ParameterGrid({})
    assert_equal(len(empty), 1)
    assert_equal(list(empty), [{}])

    has_empty = ParameterGrid([{'C': [1, 10]}, {}])
    assert_equal(len(has_empty), 3)
    assert_equal(list(has_empty), [{'C': 1}, {'C': 10}, {}])


def test_grid_search():
    # Test that the best estimator contains the right value for foo_param
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, verbose=3)
    # make sure it selects the smallest parameter in case of ties
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    grid_search.fit(X, y)
    sys.stdout = old_stdout
    assert_equal(grid_search.best_estimator_.foo_param, 2)

    for i, foo_i in enumerate([1, 2, 3]):
        assert_true(grid_search.grid_scores_[i][0]
                    == {'foo_param': foo_i})
    # Smoke test the score etc:
    grid_search.score(X, y)
    grid_search.predict_proba(X)
    grid_search.decision_function(X)
    grid_search.transform(X)

    # Test exception handling on scoring
    grid_search.scoring = 'sklearn'
    assert_raises(ValueError, grid_search.fit, X, y)


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
    assert_equal(grid_search_no_score.best_params_, grid_search.best_params_)
    # check that we can call score and that it gives the correct result
    assert_equal(grid_search.score(X, y), grid_search_no_score.score(X, y))

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
                                              scoring='roc_auc').fit(X, y)
    search_auc = GridSearchCV(clf, grid, scoring='roc_auc').fit(X, y)

    # Check warning only occurs in situation where behavior changed:
    # estimator requires score method to compete with scoring parameter
    score_no_scoring = assert_no_warnings(search_no_scoring.score, X, y)
    score_accuracy = assert_warns(ChangedBehaviorWarning,
                                  search_accuracy.score, X, y)
    score_no_score_auc = assert_no_warnings(search_no_score_method_auc.score,
                                            X, y)
    score_auc = assert_warns(ChangedBehaviorWarning,
                             search_auc.score, X, y)
    # ensure the test is sane
    assert_true(score_auc < 1.0)
    assert_true(score_accuracy < 1.0)
    assert_not_equal(score_auc, score_accuracy)

    assert_almost_equal(score_accuracy, score_no_scoring)
    assert_almost_equal(score_auc, score_no_score_auc)


def test_trivial_grid_scores():
    # Test search over a "grid" with only one point.
    # Non-regression test: grid_scores_ wouldn't be set by GridSearchCV.
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1]})
    grid_search.fit(X, y)
    assert_true(hasattr(grid_search, "grid_scores_"))

    random_search = RandomizedSearchCV(clf, {'foo_param': [0]}, n_iter=1)
    random_search.fit(X, y)
    assert_true(hasattr(random_search, "grid_scores_"))


def test_no_refit():
    # Test that grid search can be used for model selection only
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, refit=False)
    grid_search.fit(X, y)
    assert_true(hasattr(grid_search, "best_params_"))


def test_grid_search_error():
    # Test that grid search will capture errors on data with different
    # length
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)

    clf = LinearSVC()
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, X_[:180], y_)


def test_grid_search_iid():
    # test the iid parameter
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
    svm = SVC(kernel='linear')
    # create "cv" for splits
    cv = [[mask, ~mask], [~mask, mask]]
    # once with iid=True (default)
    grid_search = GridSearchCV(svm, param_grid={'C': [1, 10]}, cv=cv)
    grid_search.fit(X, y)
    first = grid_search.grid_scores_[0]
    assert_equal(first.parameters['C'], 1)
    assert_array_almost_equal(first.cv_validation_scores, [1, 1. / 3.])
    # for first split, 1/4 of dataset is in test, for second 3/4.
    # take weighted average
    assert_almost_equal(first.mean_validation_score,
                        1 * 1. / 4. + 1. / 3. * 3. / 4.)

    # once with iid=False
    grid_search = GridSearchCV(svm, param_grid={'C': [1, 10]}, cv=cv,
                               iid=False)
    grid_search.fit(X, y)
    first = grid_search.grid_scores_[0]
    assert_equal(first.parameters['C'], 1)
    # scores are the same as above
    assert_array_almost_equal(first.cv_validation_scores, [1, 1. / 3.])
    # averaged score is just mean of scores
    assert_almost_equal(first.mean_validation_score,
                        np.mean(first.cv_validation_scores))


def test_grid_search_one_grid_point():
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)
    param_dict = {"C": [1.0], "kernel": ["rbf"], "gamma": [0.1]}

    clf = SVC()
    cv = GridSearchCV(clf, param_dict)
    cv.fit(X_, y_)

    clf = SVC(C=1.0, kernel="rbf", gamma=0.1)
    clf.fit(X_, y_)

    assert_array_equal(clf.dual_coef_, cv.best_estimator_.dual_coef_)


def test_grid_search_bad_param_grid():
    param_dict = {"C": 1.0}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)

    param_dict = {"C": []}
    clf = SVC()
    assert_raises(ValueError, GridSearchCV, clf, param_dict)

    param_dict = {"C": np.ones(6).reshape(3, 2)}
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

    assert_true(np.mean(y_pred == y_pred2) >= .9)
    assert_equal(C, C2)


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
    assert_equal(C, C2)
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

    assert_equal(C, C3)
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

    assert_true(cv.best_score_ >= 0)

    # compute the test kernel matrix
    K_test = np.dot(X_[180:], X_[:180].T)
    y_test = y_[180:]

    y_pred = cv.predict(K_test)

    assert_true(np.mean(y_pred == y_test) >= 0)

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


def test_grid_search_precomputed_kernel_error_kernel_function():
    # Test that grid search returns an error when using a kernel_function
    X_, y_ = make_classification(n_samples=200, n_features=100, random_state=0)
    kernel_function = lambda x1, x2: np.dot(x1, x2.T)
    clf = SVC(kernel=kernel_function)
    cv = GridSearchCV(clf, {'C': [0.1, 1.0]})
    assert_raises(ValueError, cv.fit, X_, y_)


class BrokenClassifier(BaseEstimator):
    """Broken classifier that cannot be fit twice"""

    def __init__(self, parameter=None):
        self.parameter = parameter

    def fit(self, X, y):
        assert_true(not hasattr(self, 'has_been_fit_'))
        self.has_been_fit_ = True

    def predict(self, X):
        return np.zeros(X.shape[0])


def test_refit():
    # Regression test for bug in refitting
    # Simulates re-fitting a broken estimator; this used to break with
    # sparse SVMs.
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = GridSearchCV(BrokenClassifier(), [{'parameter': [0, 1]}],
                       scoring="precision", refit=True)
    clf.fit(X, y)


def test_gridsearch_nd():
    # Pass X as list in GridSearchCV
    X_4d = np.arange(10 * 5 * 3 * 2).reshape(10, 5, 3, 2)
    y_3d = np.arange(10 * 7 * 11).reshape(10, 7, 11)
    check_X = lambda x: x.shape[1:] == (5, 3, 2)
    check_y = lambda x: x.shape[1:] == (7, 11)
    clf = CheckingClassifier(check_X=check_X, check_y=check_y)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]})
    grid_search.fit(X_4d, y_3d).score(X, y)
    assert_true(hasattr(grid_search, "grid_scores_"))


def test_X_as_list():
    # Pass X as list in GridSearchCV
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = CheckingClassifier(check_X=lambda x: isinstance(x, list))
    cv = KFold(n=len(X), n_folds=3)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=cv)
    grid_search.fit(X.tolist(), y).score(X, y)
    assert_true(hasattr(grid_search, "grid_scores_"))


def test_y_as_list():
    # Pass y as list in GridSearchCV
    X = np.arange(100).reshape(10, 10)
    y = np.array([0] * 5 + [1] * 5)

    clf = CheckingClassifier(check_y=lambda x: isinstance(x, list))
    cv = KFold(n=len(X), n_folds=3)
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, cv=cv)
    grid_search.fit(X, y.tolist()).score(X, y)
    assert_true(hasattr(grid_search, "grid_scores_"))


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
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)

        grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]})
        grid_search.fit(X_df, y_ser).score(X_df, y_ser)
        grid_search.predict(X_df)
        assert_true(hasattr(grid_search, "grid_scores_"))


def test_unsupervised_grid_search():
    # test grid-search with unsupervised estimator
    X, y = make_blobs(random_state=0)
    km = KMeans(random_state=0)
    grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]),
                               scoring='adjusted_rand_score')
    grid_search.fit(X, y)
    # ARI can find the right number :)
    assert_equal(grid_search.best_params_["n_clusters"], 3)

    # Now without a score, and without y
    grid_search = GridSearchCV(km, param_grid=dict(n_clusters=[2, 3, 4]))
    grid_search.fit(X)
    assert_equal(grid_search.best_params_["n_clusters"], 4)


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
    assert_equal(search.best_params_['bandwidth'], .1)
    assert_equal(search.best_score_, 42)


def test_param_sampler():
    # test basic properties of param sampler
    param_distributions = {"kernel": ["rbf", "linear"],
                           "C": uniform(0, 1)}
    sampler = ParameterSampler(param_distributions=param_distributions,
                               n_iter=10, random_state=0)
    samples = [x for x in sampler]
    assert_equal(len(samples), 10)
    for sample in samples:
        assert_true(sample["kernel"] in ["rbf", "linear"])
        assert_true(0 <= sample["C"] <= 1)


def test_randomized_search_grid_scores():
    # Make a dataset with a lot of noise to get various kind of prediction
    # errors across CV folds and parameter settings
    X, y = make_classification(n_samples=200, n_features=100, n_informative=3,
                               random_state=0)

    # XXX: as of today (scipy 0.12) it's not possible to set the random seed
    # of scipy.stats distributions: the assertions in this test should thus
    # not depend on the randomization
    params = dict(C=expon(scale=10),
                  gamma=expon(scale=0.1))
    n_cv_iter = 3
    n_search_iter = 30
    search = RandomizedSearchCV(SVC(), n_iter=n_search_iter, cv=n_cv_iter,
                                param_distributions=params, iid=False)
    search.fit(X, y)
    assert_equal(len(search.grid_scores_), n_search_iter)

    # Check consistency of the structure of each cv_score item
    for cv_score in search.grid_scores_:
        assert_equal(len(cv_score.cv_validation_scores), n_cv_iter)
        # Because we set iid to False, the mean_validation score is the
        # mean of the fold mean scores instead of the aggregate sample-wise
        # mean score
        assert_almost_equal(np.mean(cv_score.cv_validation_scores),
                            cv_score.mean_validation_score)
        assert_equal(list(sorted(cv_score.parameters.keys())),
                     list(sorted(params.keys())))

    # Check the consistency with the best_score_ and best_params_ attributes
    sorted_grid_scores = list(sorted(search.grid_scores_,
                              key=lambda x: x.mean_validation_score))
    best_score = sorted_grid_scores[-1].mean_validation_score
    assert_equal(search.best_score_, best_score)

    tied_best_params = [s.parameters for s in sorted_grid_scores
                        if s.mean_validation_score == best_score]
    assert_true(search.best_params_ in tied_best_params,
                "best_params_={0} is not part of the"
                " tied best models: {1}".format(
                    search.best_params_, tied_best_params))


def test_grid_search_score_consistency():
    # test that correct scores are used
    clf = LinearSVC(random_state=0)
    X, y = make_blobs(random_state=0, centers=2)
    Cs = [.1, 1, 10]
    for score in ['f1', 'roc_auc']:
        grid_search = GridSearchCV(clf, {'C': Cs}, scoring=score)
        grid_search.fit(X, y)
        cv = StratifiedKFold(n_folds=3, y=y)
        for C, scores in zip(Cs, grid_search.grid_scores_):
            clf.set_params(C=C)
            scores = scores[2]  # get the separate runs from grid scores
            i = 0
            for train, test in cv:
                clf.fit(X[train], y[train])
                if score == "f1":
                    correct_score = f1_score(y[test], clf.predict(X[test]))
                elif score == "roc_auc":
                    dec = clf.decision_function(X[test])
                    correct_score = roc_auc_score(y[test], dec)
                assert_almost_equal(correct_score, scores[i])
                i += 1


def test_pickle():
    # Test that a fit search can be pickled
    clf = MockClassifier()
    grid_search = GridSearchCV(clf, {'foo_param': [1, 2, 3]}, refit=True)
    grid_search.fit(X, y)
    pickle.dumps(grid_search)  # smoke test

    random_search = RandomizedSearchCV(clf, {'foo_param': [1, 2, 3]},
                                       refit=True, n_iter=3)
    random_search.fit(X, y)
    pickle.dumps(random_search)  # smoke test


def test_grid_search_with_multioutput_data():
    # Test search with multi-output estimator

    X, y = make_multilabel_classification(return_indicator=True,
                                          random_state=0)

    est_parameters = {"max_depth": [1, 2, 3, 4]}
    cv = KFold(y.shape[0], random_state=0)

    estimators = [DecisionTreeRegressor(random_state=0),
                  DecisionTreeClassifier(random_state=0)]

    # Test with grid search cv
    for est in estimators:
        grid_search = GridSearchCV(est, est_parameters, cv=cv)
        grid_search.fit(X, y)
        for parameters, _, cv_validation_scores in grid_search.grid_scores_:
            est.set_params(**parameters)

            for i, (train, test) in enumerate(cv):
                est.fit(X[train], y[train])
                correct_score = est.score(X[test], y[test])
                assert_almost_equal(correct_score,
                                    cv_validation_scores[i])

    # Test with a randomized search
    for est in estimators:
        random_search = RandomizedSearchCV(est, est_parameters,
                                           cv=cv, n_iter=3)
        random_search.fit(X, y)
        for parameters, _, cv_validation_scores in random_search.grid_scores_:
            est.set_params(**parameters)

            for i, (train, test) in enumerate(cv):
                est.fit(X[train], y[train])
                correct_score = est.score(X[test], y[test])
                assert_almost_equal(correct_score,
                                    cv_validation_scores[i])


def test_predict_proba_disabled():
    # Test predict_proba when disabled on estimator.
    X = np.arange(20).reshape(5, -1)
    y = [0, 0, 1, 1, 1]
    clf = SVC(probability=False)
    gs = GridSearchCV(clf, {}, cv=2).fit(X, y)
    assert_false(hasattr(gs, "predict_proba"))


def test_grid_search_allows_nans():
    # Test GridSearchCV with Imputer
    X = np.arange(20, dtype=np.float64).reshape(5, -1)
    X[2, :] = np.nan
    y = [0, 0, 1, 1, 1]
    p = Pipeline([
        ('imputer', Imputer(strategy='mean', missing_values='NaN')),
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

    # Ensure that grid scores were set to zero as required for those fits
    # that are expected to fail.
    assert all(np.all(this_point.cv_validation_scores == 0.0)
               for this_point in gs.grid_scores_
               if this_point.parameters['parameter'] ==
               FailingClassifier.FAILING_PARAMETER)

    gs = GridSearchCV(clf, [{'parameter': [0, 1, 2]}], scoring='accuracy',
                      refit=False, error_score=float('nan'))
    assert_warns(FitFailedWarning, gs.fit, X, y)
    assert all(np.all(np.isnan(this_point.cv_validation_scores))
               for this_point in gs.grid_scores_
               if this_point.parameters['parameter'] ==
               FailingClassifier.FAILING_PARAMETER)


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
    # raise error if n_iter too large
    params = {'first': [0, 1], 'second': ['a', 'b', 'c']}
    sampler = ParameterSampler(params, n_iter=7)
    assert_raises(ValueError, list, sampler)
    # degenerates to GridSearchCV if n_iter the same as grid_size
    sampler = ParameterSampler(params, n_iter=6)
    samples = list(sampler)
    assert_equal(len(samples), 6)
    for values in ParameterGrid(params):
        assert_true(values in samples)

    # test sampling without replacement in a large grid
    params = {'a': range(10), 'b': range(10), 'c': range(10)}
    sampler = ParameterSampler(params, n_iter=99, random_state=42)
    samples = list(sampler)
    assert_equal(len(samples), 99)
    hashable_samples = ["a%db%dc%d" % (p['a'], p['b'], p['c'])
                        for p in samples]
    assert_equal(len(set(hashable_samples)), 99)

    # doesn't go into infinite loops
    params_distribution = {'first': bernoulli(.5), 'second': ['a', 'b', 'c']}
    sampler = ParameterSampler(params_distribution, n_iter=7)
    samples = list(sampler)
    assert_equal(len(samples), 7)
