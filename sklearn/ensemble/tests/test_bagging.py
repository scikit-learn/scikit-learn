"""
Testing for the bagging ensemble module (sklearn.ensemble.bagging).
"""

# Author: Gilles Louppe
# License: BSD 3 clause

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

from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.model_selection import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import Perceptron, LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR
from sklearn.pipeline import make_pipeline
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_boston, load_iris, make_hastie_10_2
from sklearn.utils import check_random_state

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
                           Perceptron(),
                           DecisionTreeClassifier(),
                           KNeighborsClassifier(),
                           SVC()]:
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
                    base_estimator=CustomSVC(decision_function_shape='ovr'),
                    random_state=1,
                    **params
                ).fit(X_train_sparse, y_train)
                sparse_results = getattr(sparse_classifier, f)(X_test_sparse)

                # Trained on dense format
                dense_classifier = BaggingClassifier(
                    base_estimator=CustomSVC(decision_function_shape='ovr'),
                    random_state=1,
                    **params
                ).fit(X_train, y_train)
                dense_results = getattr(dense_classifier, f)(X_test)
                assert_array_equal(sparse_results, dense_results)

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
                           SVR()]:
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
                base_estimator=CustomSVR(),
                random_state=1,
                **params
            ).fit(X_train_sparse, y_train)
            sparse_results = sparse_classifier.predict(X_test_sparse)

            # Trained on dense format
            dense_results = BaggingRegressor(
                base_estimator=CustomSVR(),
                random_state=1,
                **params
            ).fit(X_train, y_train).predict(X_test)

            sparse_type = type(X_train_sparse)
            types = [i.data_type_ for i in sparse_classifier.estimators_]

            assert_array_equal(sparse_results, dense_results)
            assert all([t == sparse_type for t in types])
            assert_array_equal(sparse_results, dense_results)


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

    for base_estimator in [DecisionTreeClassifier(), SVC()]:
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

    assert_array_equal(clf1.predict(X_test), clf2.predict(X_test))


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
    ensemble = BaggingClassifier(SVC(decision_function_shape='ovr'),
                                 n_jobs=3,
                                 random_state=0).fit(X_train, y_train)

    ensemble.set_params(n_jobs=1)
    decisions1 = ensemble.decision_function(X_test)
    ensemble.set_params(n_jobs=2)
    decisions2 = ensemble.decision_function(X_test)
    assert_array_almost_equal(decisions1, decisions2)

    ensemble = BaggingClassifier(SVC(decision_function_shape='ovr'),
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


def test_gridsearch():
    # Check that bagging ensembles can be grid-searched.
    # Transform iris into a binary classification task
    X, y = iris.data, iris.target
    y[y == 2] = 1

    # Grid search with scoring based on decision_function
    parameters = {'n_estimators': (1, 2),
                  'base_estimator__C': (1, 2)}

    GridSearchCV(BaggingClassifier(SVC()),
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

    ensemble = BaggingClassifier(Perceptron(),
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

    ensemble = BaggingRegressor(SVR(),
                                n_jobs=3,
                                random_state=0).fit(X_train, y_train)
    assert_true(isinstance(ensemble.base_estimator_, SVR))


def test_bagging_with_pipeline():
    estimator = BaggingClassifier(make_pipeline(SelectKBest(k=1),
                                                DecisionTreeClassifier()),
                                  max_features=2)
    estimator.fit(iris.data, iris.target)


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
