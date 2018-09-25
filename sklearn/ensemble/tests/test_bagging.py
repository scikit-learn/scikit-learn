"""
Testing for the bagging ensemble module (sklearn.ensemble.bagging).
"""

# Author: Gilles Louppe
# License: BSD 3 clause

import pytest
import numpy as np

from sklearn.base import BaseEstimator

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_raise_message

from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.model_selection import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import Perceptron, LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR
from sklearn.random_projection import SparseRandomProjection
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_boston, load_iris, make_hastie_10_2
from sklearn.utils import check_random_state, hash
from sklearn.preprocessing import FunctionTransformer

from scipy.sparse import csc_matrix, csr_matrix

rng = check_random_state(0)

# also load the iris dataset
# and randomly permute it
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification():
    # Check classification for various parameter settings.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "max_features": [1, 2, 4],
                          "bootstrap": [True, False],
                          "bootstrap_features": [True, False]})

    for base_estimator in [None,
                           DummyClassifier(),
                           Perceptron(tol=1e-3),
                           DecisionTreeClassifier(),
                           KNeighborsClassifier(),
                           SVC(gamma="scale")]:
        for params in grid:
            BaggingClassifier(base_estimator=base_estimator,
                              random_state=rng,
                              **params).fit(X_train, y_train).predict(X_test)


def test_sparse_classification():
    # Check classification for various parameter settings on sparse input.

    class CustomSVC(SVC):
        """SVC variant that records the nature of the training set"""

        def fit(self, X, y):
            super(CustomSVC, self).fit(X, y)
            self.data_type_ = type(X)
            return self

    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)
    parameter_sets = [
        {"max_samples": 0.5,
         "max_features": 2,
         "bootstrap": True,
         "bootstrap_features": True},
        {"max_samples": 1.0,
         "max_features": 4,
         "bootstrap": True,
         "bootstrap_features": True},
        {"max_features": 2,
         "bootstrap": False,
         "bootstrap_features": True},
        {"max_samples": 0.5,
         "bootstrap": True,
         "bootstrap_features": False},
    ]

    for sparse_format in [csc_matrix, csr_matrix]:
        X_train_sparse = sparse_format(X_train)
        X_test_sparse = sparse_format(X_test)
        for params in parameter_sets:
            for f in ['predict', 'predict_proba', 'predict_log_proba', 'decision_function']:
                # Trained on sparse format
                sparse_classifier = BaggingClassifier(
                    base_estimator=CustomSVC(gamma='scale',
                                             decision_function_shape='ovr'),
                    random_state=1,
                    **params
                ).fit(X_train_sparse, y_train)
                sparse_results = getattr(sparse_classifier, f)(X_test_sparse)

                # Trained on dense format
                dense_classifier = BaggingClassifier(
                    base_estimator=CustomSVC(gamma='scale',
                                             decision_function_shape='ovr'),
                    random_state=1,
                    **params
                ).fit(X_train, y_train)
                dense_results = getattr(dense_classifier, f)(X_test)
                assert_array_almost_equal(sparse_results, dense_results)

            sparse_type = type(X_train_sparse)
            types = [i.data_type_ for i in sparse_classifier.estimators_]

            assert all([t == sparse_type for t in types])


def test_regression():
    # Check regression for various parameter settings.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "max_features": [0.5, 1.0],
                          "bootstrap": [True, False],
                          "bootstrap_features": [True, False]})

    for base_estimator in [None,
                           DummyRegressor(),
                           DecisionTreeRegressor(),
                           KNeighborsRegressor(),
                           SVR(gamma='scale')]:
        for params in grid:
            BaggingRegressor(base_estimator=base_estimator,
                             random_state=rng,
                             **params).fit(X_train, y_train).predict(X_test)


def test_sparse_regression():
    # Check regression for various parameter settings on sparse input.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data[:50],
                                                        boston.target[:50],
                                                        random_state=rng)

    class CustomSVR(SVR):
        """SVC variant that records the nature of the training set"""

        def fit(self, X, y):
            super(CustomSVR, self).fit(X, y)
            self.data_type_ = type(X)
            return self

    parameter_sets = [
        {"max_samples": 0.5,
         "max_features": 2,
         "bootstrap": True,
         "bootstrap_features": True},
        {"max_samples": 1.0,
         "max_features": 4,
         "bootstrap": True,
         "bootstrap_features": True},
        {"max_features": 2,
         "bootstrap": False,
         "bootstrap_features": True},
        {"max_samples": 0.5,
         "bootstrap": True,
         "bootstrap_features": False},
    ]

    for sparse_format in [csc_matrix, csr_matrix]:
        X_train_sparse = sparse_format(X_train)
        X_test_sparse = sparse_format(X_test)
        for params in parameter_sets:

            # Trained on sparse format
            sparse_classifier = BaggingRegressor(
                base_estimator=CustomSVR(gamma='scale'),
                random_state=1,
                **params
            ).fit(X_train_sparse, y_train)
            sparse_results = sparse_classifier.predict(X_test_sparse)

            # Trained on dense format
            dense_results = BaggingRegressor(
                base_estimator=CustomSVR(gamma='scale'),
                random_state=1,
                **params
            ).fit(X_train, y_train).predict(X_test)

            sparse_type = type(X_train_sparse)
            types = [i.data_type_ for i in sparse_classifier.estimators_]

            assert_array_almost_equal(sparse_results, dense_results)
            assert all([t == sparse_type for t in types])
            assert_array_almost_equal(sparse_results, dense_results)


class DummySizeEstimator(BaseEstimator):

    def fit(self, X, y):
        self.training_size_ = X.shape[0]
        self.training_hash_ = hash(X)


def test_bootstrap_samples():
    # Test that bootstrapping samples generate non-perfect base estimators.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    base_estimator = DecisionTreeRegressor().fit(X_train, y_train)

    # without bootstrap, all trees are perfect on the training set
    ensemble = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                                max_samples=1.0,
                                bootstrap=False,
                                random_state=rng).fit(X_train, y_train)

    assert_equal(base_estimator.score(X_train, y_train),
                 ensemble.score(X_train, y_train))

    # with bootstrap, trees are no longer perfect on the training set
    ensemble = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                                max_samples=1.0,
                                bootstrap=True,
                                random_state=rng).fit(X_train, y_train)

    assert_greater(base_estimator.score(X_train, y_train),
                   ensemble.score(X_train, y_train))

    # check that each sampling correspond to a complete bootstrap resample.
    # the size of each bootstrap should be the same as the input data but
    # the data should be different (checked using the hash of the data).
    ensemble = BaggingRegressor(base_estimator=DummySizeEstimator(),
                                bootstrap=True).fit(X_train, y_train)
    training_hash = []
    for estimator in ensemble.estimators_:
        assert estimator.training_size_ == X_train.shape[0]
        training_hash.append(estimator.training_hash_)
    assert len(set(training_hash)) == len(training_hash)


def test_bootstrap_features():
    # Test that bootstrapping features may generate duplicate features.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    ensemble = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                                max_features=1.0,
                                bootstrap_features=False,
                                random_state=rng).fit(X_train, y_train)

    for features in ensemble.estimators_features_:
        assert_equal(boston.data.shape[1], np.unique(features).shape[0])

    ensemble = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                                max_features=1.0,
                                bootstrap_features=True,
                                random_state=rng).fit(X_train, y_train)

    for features in ensemble.estimators_features_:
        assert_greater(boston.data.shape[1], np.unique(features).shape[0])


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_probability():
    # Predict probabilities.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)

    with np.errstate(divide="ignore", invalid="ignore"):
        # Normal case
        ensemble = BaggingClassifier(base_estimator=DecisionTreeClassifier(),
                                     random_state=rng).fit(X_train, y_train)

        assert_array_almost_equal(np.sum(ensemble.predict_proba(X_test),
                                         axis=1),
                                  np.ones(len(X_test)))

        assert_array_almost_equal(ensemble.predict_proba(X_test),
                                  np.exp(ensemble.predict_log_proba(X_test)))

        # Degenerate case, where some classes are missing
        ensemble = BaggingClassifier(base_estimator=LogisticRegression(),
                                     random_state=rng,
                                     max_samples=5).fit(X_train, y_train)

        assert_array_almost_equal(np.sum(ensemble.predict_proba(X_test),
                                         axis=1),
                                  np.ones(len(X_test)))

        assert_array_almost_equal(ensemble.predict_proba(X_test),
                                  np.exp(ensemble.predict_log_proba(X_test)))


def test_oob_score_classification():
    # Check that oob prediction is a good estimation of the generalization
    # error.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)

    for base_estimator in [DecisionTreeClassifier(), SVC(gamma="scale")]:
        clf = BaggingClassifier(base_estimator=base_estimator,
                                n_estimators=100,
                                bootstrap=True,
                                oob_score=True,
                                random_state=rng).fit(X_train, y_train)

        test_score = clf.score(X_test, y_test)

        assert_less(abs(test_score - clf.oob_score_), 0.1)

        # Test with few estimators
        assert_warns(UserWarning,
                     BaggingClassifier(base_estimator=base_estimator,
                                       n_estimators=1,
                                       bootstrap=True,
                                       oob_score=True,
                                       random_state=rng).fit,
                     X_train,
                     y_train)


def test_oob_score_regression():
    # Check that oob prediction is a good estimation of the generalization
    # error.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    clf = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                           n_estimators=50,
                           bootstrap=True,
                           oob_score=True,
                           random_state=rng).fit(X_train, y_train)

    test_score = clf.score(X_test, y_test)

    assert_less(abs(test_score - clf.oob_score_), 0.1)

    # Test with few estimators
    assert_warns(UserWarning,
                 BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                                  n_estimators=1,
                                  bootstrap=True,
                                  oob_score=True,
                                  random_state=rng).fit,
                 X_train,
                 y_train)


def test_single_estimator():
    # Check singleton ensembles.
    rng = check_random_state(0)
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    clf1 = BaggingRegressor(base_estimator=KNeighborsRegressor(),
                            n_estimators=1,
                            bootstrap=False,
                            bootstrap_features=False,
                            random_state=rng).fit(X_train, y_train)

    clf2 = KNeighborsRegressor().fit(X_train, y_train)

    assert_array_almost_equal(clf1.predict(X_test), clf2.predict(X_test))


def test_error():
    # Test that it gives proper exception on deficient input.
    X, y = iris.data, iris.target
    base = DecisionTreeClassifier()

    # Test max_samples
    assert_raises(ValueError,
                  BaggingClassifier(base, max_samples=-1).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_samples=0.0).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_samples=2.0).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_samples=1000).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_samples="foobar").fit, X, y)

    # Test max_features
    assert_raises(ValueError,
                  BaggingClassifier(base, max_features=-1).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_features=0.0).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_features=2.0).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_features=5).fit, X, y)
    assert_raises(ValueError,
                  BaggingClassifier(base, max_features="foobar").fit, X, y)

    # Test support of decision_function
    assert_false(hasattr(BaggingClassifier(base).fit(X, y), 'decision_function'))


def test_parallel_classification():
    # Check parallel classification.
    rng = check_random_state(0)

    # Classification
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)

    ensemble = BaggingClassifier(DecisionTreeClassifier(),
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    # predict_proba
    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict_proba(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict_proba(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = BaggingClassifier(DecisionTreeClassifier(),
                                 n_jobs=1,
                                 random_state=0).fit(X_train, y_train)

    y3 = ensemble.predict_proba(X_test)
    assert_array_almost_equal(y1, y3)

    # decision_function
    ensemble = BaggingClassifier(SVC(gamma='scale',
                                     decision_function_shape='ovr'),
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    ensemble.set_params(n_jobs=1)
    decisions1 = ensemble.decision_function(X_test)
    ensemble.set_params(n_jobs=2)
    decisions2 = ensemble.decision_function(X_test)
    assert_array_almost_equal(decisions1, decisions2)

    X_err = np.hstack((X_test, np.zeros((X_test.shape[0], 1))))
    assert_raise_message(ValueError, "Number of features of the model "
                         "must match the input. Model n_features is {0} "
                         "and input n_features is {1} "
                         "".format(X_test.shape[1], X_err.shape[1]),
                         ensemble.decision_function, X_err)

    ensemble = BaggingClassifier(SVC(gamma='scale',
                                     decision_function_shape='ovr'),
                                 n_jobs=1,
                                 random_state=0).fit(X_train, y_train)

    decisions3 = ensemble.decision_function(X_test)
    assert_array_almost_equal(decisions1, decisions3)


def test_parallel_regression():
    # Check parallel regression.
    rng = check_random_state(0)

    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    ensemble = BaggingRegressor(DecisionTreeRegressor(),
                                n_jobs=3,
                                random_state=0).fit(X_train, y_train)

    ensemble.set_params(n_jobs=1)
    y1 = ensemble.predict(X_test)
    ensemble.set_params(n_jobs=2)
    y2 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y2)

    ensemble = BaggingRegressor(DecisionTreeRegressor(),
                                n_jobs=1,
                                random_state=0).fit(X_train, y_train)

    y3 = ensemble.predict(X_test)
    assert_array_almost_equal(y1, y3)


@pytest.mark.filterwarnings('ignore: The default of the `iid`')  # 0.22
@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_gridsearch():
    # Check that bagging ensembles can be grid-searched.
    # Transform iris into a binary classification task
    X, y = iris.data, iris.target
    y[y == 2] = 1

    # Grid search with scoring based on decision_function
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__C': (1, 2)}

    GridSearchCV(BaggingClassifier(SVC(gamma="scale")),
                 parameters,
                 scoring="roc_auc").fit(X, y)


def test_base_estimator():
    # Check base_estimator and its default values.
    rng = check_random_state(0)

    # Classification
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)

    ensemble = BaggingClassifier(None,
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    assert_true(isinstance(ensemble.base_estimator_, DecisionTreeClassifier))

    ensemble = BaggingClassifier(DecisionTreeClassifier(),
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    assert_true(isinstance(ensemble.base_estimator_, DecisionTreeClassifier))

    ensemble = BaggingClassifier(Perceptron(tol=1e-3),
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    assert_true(isinstance(ensemble.base_estimator_, Perceptron))

    # Regression
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    ensemble = BaggingRegressor(None,
                                n_jobs=3,
                                random_state=0).fit(X_train, y_train)

    assert_true(isinstance(ensemble.base_estimator_, DecisionTreeRegressor))

    ensemble = BaggingRegressor(DecisionTreeRegressor(),
                                n_jobs=3,
                                random_state=0).fit(X_train, y_train)

    assert_true(isinstance(ensemble.base_estimator_, DecisionTreeRegressor))

    ensemble = BaggingRegressor(SVR(gamma='scale'),
                                n_jobs=3,
                                random_state=0).fit(X_train, y_train)
    assert_true(isinstance(ensemble.base_estimator_, SVR))


def test_bagging_with_pipeline():
    estimator = BaggingClassifier(make_pipeline(SelectKBest(k=1),
                                                DecisionTreeClassifier()),
                                  max_features=2)
    estimator.fit(iris.data, iris.target)
    assert_true(isinstance(estimator[0].steps[-1][1].random_state,
                           int))


class DummyZeroEstimator(BaseEstimator):

    def fit(self, X, y):
        self.classes_ = np.unique(y)
        return self

    def predict(self, X):
        return self.classes_[np.zeros(X.shape[0], dtype=int)]


def test_bagging_sample_weight_unsupported_but_passed():
    estimator = BaggingClassifier(DummyZeroEstimator())
    rng = check_random_state(0)

    estimator.fit(iris.data, iris.target).predict(iris.data)
    assert_raises(ValueError, estimator.fit, iris.data, iris.target,
                  sample_weight=rng.randint(10, size=(iris.data.shape[0])))


def test_warm_start(random_state=42):
    # Test if fitting incrementally with warm start gives a forest of the
    # right size and the same results as a normal fit.
    X, y = make_hastie_10_2(n_samples=20, random_state=1)

    clf_ws = None
    for n_estimators in [5, 10]:
        if clf_ws is None:
            clf_ws = BaggingClassifier(n_estimators=n_estimators,
                                       random_state=random_state,
                                       warm_start=True)
        else:
            clf_ws.set_params(n_estimators=n_estimators)
        clf_ws.fit(X, y)
        assert_equal(len(clf_ws), n_estimators)

    clf_no_ws = BaggingClassifier(n_estimators=10, random_state=random_state,
                                  warm_start=False)
    clf_no_ws.fit(X, y)

    assert_equal(set([tree.random_state for tree in clf_ws]),
                 set([tree.random_state for tree in clf_no_ws]))


def test_warm_start_smaller_n_estimators():
    # Test if warm start'ed second fit with smaller n_estimators raises error.
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    clf = BaggingClassifier(n_estimators=5, warm_start=True)
    clf.fit(X, y)
    clf.set_params(n_estimators=4)
    assert_raises(ValueError, clf.fit, X, y)


def test_warm_start_equal_n_estimators():
    # Test that nothing happens when fitting without increasing n_estimators
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=43)

    clf = BaggingClassifier(n_estimators=5, warm_start=True, random_state=83)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    # modify X to nonsense values, this should not change anything
    X_train += 1.

    assert_warns_message(UserWarning,
                         "Warm-start fitting without increasing n_estimators does not",
                         clf.fit, X_train, y_train)
    assert_array_equal(y_pred, clf.predict(X_test))


def test_warm_start_equivalence():
    # warm started classifier with 5+5 estimators should be equivalent to
    # one classifier with 10 estimators
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=43)

    clf_ws = BaggingClassifier(n_estimators=5, warm_start=True,
                               random_state=3141)
    clf_ws.fit(X_train, y_train)
    clf_ws.set_params(n_estimators=10)
    clf_ws.fit(X_train, y_train)
    y1 = clf_ws.predict(X_test)

    clf = BaggingClassifier(n_estimators=10, warm_start=False,
                            random_state=3141)
    clf.fit(X_train, y_train)
    y2 = clf.predict(X_test)

    assert_array_almost_equal(y1, y2)


def test_warm_start_with_oob_score_fails():
    # Check using oob_score and warm_start simultaneously fails
    X, y = make_hastie_10_2(n_samples=20, random_state=1)
    clf = BaggingClassifier(n_estimators=5, warm_start=True, oob_score=True)
    assert_raises(ValueError, clf.fit, X, y)


def test_oob_score_removed_on_warm_start():
    X, y = make_hastie_10_2(n_samples=2000, random_state=1)

    clf = BaggingClassifier(n_estimators=50, oob_score=True)
    clf.fit(X, y)

    clf.set_params(warm_start=True, oob_score=False, n_estimators=100)
    clf.fit(X, y)

    assert_raises(AttributeError, getattr, clf, "oob_score_")


def test_oob_score_consistency():
    # Make sure OOB scores are identical when random_state, estimator, and
    # training data are fixed and fitting is done twice
    X, y = make_hastie_10_2(n_samples=200, random_state=1)
    bagging = BaggingClassifier(KNeighborsClassifier(), max_samples=0.5,
                                max_features=0.5, oob_score=True,
                                random_state=1)
    assert_equal(bagging.fit(X, y).oob_score_, bagging.fit(X, y).oob_score_)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_estimators_samples():
    # Check that format of estimators_samples_ is correct and that results
    # generated at fit time can be identically reproduced at a later time
    # using data saved in object attributes.
    X, y = make_hastie_10_2(n_samples=200, random_state=1)
    bagging = BaggingClassifier(LogisticRegression(), max_samples=0.5,
                                max_features=0.5, random_state=1,
                                bootstrap=False)
    bagging.fit(X, y)

    # Get relevant attributes
    estimators_samples = bagging.estimators_samples_
    estimators_features = bagging.estimators_features_
    estimators = bagging.estimators_

    # Test for correct formatting
    assert_equal(len(estimators_samples), len(estimators))
    assert_equal(len(estimators_samples[0]), len(X) // 2)
    assert_equal(estimators_samples[0].dtype.kind, 'i')

    # Re-fit single estimator to test for consistent sampling
    estimator_index = 0
    estimator_samples = estimators_samples[estimator_index]
    estimator_features = estimators_features[estimator_index]
    estimator = estimators[estimator_index]

    X_train = (X[estimator_samples])[:, estimator_features]
    y_train = y[estimator_samples]

    orig_coefs = estimator.coef_
    estimator.fit(X_train, y_train)
    new_coefs = estimator.coef_

    assert_array_almost_equal(orig_coefs, new_coefs)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_estimators_samples_deterministic():
    # This test is a regression test to check that with a random step
    # (e.g. SparseRandomProjection) and a given random state, the results
    # generated at fit time can be identically reproduced at a later time using
    # data saved in object attributes. Check issue #9524 for full discussion.

    iris = load_iris()
    X, y = iris.data, iris.target

    base_pipeline = make_pipeline(SparseRandomProjection(n_components=2),
                                  LogisticRegression())
    clf = BaggingClassifier(base_estimator=base_pipeline,
                            max_samples=0.5,
                            random_state=0)
    clf.fit(X, y)
    pipeline_estimator_coef = clf.estimators_[0].steps[-1][1].coef_.copy()

    estimator = clf.estimators_[0]
    estimator_sample = clf.estimators_samples_[0]
    estimator_feature = clf.estimators_features_[0]

    X_train = (X[estimator_sample])[:, estimator_feature]
    y_train = y[estimator_sample]

    estimator.fit(X_train, y_train)
    assert_array_equal(estimator.steps[-1][1].coef_, pipeline_estimator_coef)


def test_max_samples_consistency():
    # Make sure validated max_samples and original max_samples are identical
    # when valid integer max_samples supplied by user
    max_samples = 100
    X, y = make_hastie_10_2(n_samples=2*max_samples, random_state=1)
    bagging = BaggingClassifier(KNeighborsClassifier(),
                                max_samples=max_samples,
                                max_features=0.5, random_state=1)
    bagging.fit(X, y)
    assert_equal(bagging._max_samples, max_samples)


def test_set_oob_score_label_encoding():
    # Make sure the oob_score doesn't change when the labels change
    # See: https://github.com/scikit-learn/scikit-learn/issues/8933
    random_state = 5
    X = [[-1], [0], [1]] * 5
    Y1 = ['A', 'B', 'C'] * 5
    Y2 = [-1, 0, 1] * 5
    Y3 = [0, 1, 2] * 5
    x1 = BaggingClassifier(oob_score=True,
                           random_state=random_state).fit(X, Y1).oob_score_
    x2 = BaggingClassifier(oob_score=True,
                           random_state=random_state).fit(X, Y2).oob_score_
    x3 = BaggingClassifier(oob_score=True,
                           random_state=random_state).fit(X, Y3).oob_score_
    assert_equal([x1, x2], [x3, x3])


def replace(X):
    X = X.copy().astype('float')
    X[~np.isfinite(X)] = 0
    return X


def test_bagging_regressor_with_missing_inputs():
    # Check that BaggingRegressor can accept X with missing/infinite data
    X = np.array([
        [1, 3, 5],
        [2, None, 6],
        [2, np.nan, 6],
        [2, np.inf, 6],
        [2, np.NINF, 6],
    ])
    y_values = [
        np.array([2, 3, 3, 3, 3]),
        np.array([
            [2, 1, 9],
            [3, 6, 8],
            [3, 6, 8],
            [3, 6, 8],
            [3, 6, 8],
        ])
    ]
    for y in y_values:
        regressor = DecisionTreeRegressor()
        pipeline = make_pipeline(
            FunctionTransformer(replace, validate=False),
            regressor
        )
        pipeline.fit(X, y).predict(X)
        bagging_regressor = BaggingRegressor(pipeline)
        y_hat = bagging_regressor.fit(X, y).predict(X)
        assert_equal(y.shape, y_hat.shape)

        # Verify that exceptions can be raised by wrapper regressor
        regressor = DecisionTreeRegressor()
        pipeline = make_pipeline(regressor)
        assert_raises(ValueError, pipeline.fit, X, y)
        bagging_regressor = BaggingRegressor(pipeline)
        assert_raises(ValueError, bagging_regressor.fit, X, y)


def test_bagging_classifier_with_missing_inputs():
    # Check that BaggingClassifier can accept X with missing/infinite data
    X = np.array([
        [1, 3, 5],
        [2, None, 6],
        [2, np.nan, 6],
        [2, np.inf, 6],
        [2, np.NINF, 6],
    ])
    y = np.array([3, 6, 6, 6, 6])
    classifier = DecisionTreeClassifier()
    pipeline = make_pipeline(
        FunctionTransformer(replace, validate=False),
        classifier
    )
    pipeline.fit(X, y).predict(X)
    bagging_classifier = BaggingClassifier(pipeline)
    bagging_classifier.fit(X, y)
    y_hat = bagging_classifier.predict(X)
    assert_equal(y.shape, y_hat.shape)
    bagging_classifier.predict_log_proba(X)
    bagging_classifier.predict_proba(X)

    # Verify that exceptions can be raised by wrapper classifier
    classifier = DecisionTreeClassifier()
    pipeline = make_pipeline(classifier)
    assert_raises(ValueError, pipeline.fit, X, y)
    bagging_classifier = BaggingClassifier(pipeline)
    assert_raises(ValueError, bagging_classifier.fit, X, y)
