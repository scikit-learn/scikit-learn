"""
Testing for the ensemble module (sklearn.ensemble.ensemble).
"""

# Authors: Gilles Louppe
# License: BSD 3 clause

import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less

from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.grid_search import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import SGDClassifier, Perceptron
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR

from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_circles, load_boston, load_iris
from sklearn.utils import check_random_state

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
    """Check classificationfor various parameter settings."""
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "max_features": [1, 2, 3, 4],
                          "bootstrap": [True, False],
                          "bootstrap_features": [True, False]})

    for base_estimator in [DummyClassifier(),
                           Perceptron(),
                           DecisionTreeClassifier(),
                           KNeighborsClassifier(),
                           SVC()]:
        for params in grid:
            base_estimator.fit(X_train, y_train)
            ensemble = BaggingClassifier(base_estimator=base_estimator,
                                         random_state=rng,
                                         **params).fit(X_train, y_train)


def test_regression():
    """Check regression for various parameter settings."""
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "max_features": [0.5, 1.0],
                          "bootstrap": [True, False],
                          "bootstrap_features": [True, False]})

    for base_estimator in [DummyRegressor(),
                           DecisionTreeRegressor(),
                           KNeighborsRegressor(),
                           SVR()]:
        for params in grid:
            base_estimator.fit(X_train, y_train)
            ensemble = BaggingRegressor(base_estimator=base_estimator,
                                        random_state=rng,
                                        **params).fit(X_train, y_train)


def test_bootstrap_samples():
    """Test that bootstraping samples generate non-perfect base estimators."""
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
    """Test that bootstraping features may generate dupplicate features."""
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
    """Predict probabilities."""
    olderr = np.seterr(divide="ignore")

    ensemble = BaggingClassifier(base_estimator=DecisionTreeClassifier())
    ensemble.fit(iris.data, iris.target)

    assert_array_almost_equal(np.sum(ensemble.predict_proba(iris.data), axis=1),
                              np.ones(iris.data.shape[0]))
    assert_array_almost_equal(ensemble.predict_proba(iris.data),
                              np.exp(ensemble.predict_log_proba(iris.data)))

    np.seterr(**olderr)


def test_oob_score_classification():
    """Check that oob prediction is a good estimation of the generalization
    error."""
    X_train, X_test, y_train, y_test = train_test_split(iris.data,
                                                        iris.target,
                                                        random_state=rng)

    clf = BaggingClassifier(base_estimator=DecisionTreeClassifier(),
                            n_estimators=50,
                            bootstrap=True,
                            oob_score=True,
                            random_state=rng).fit(X_train, y_train)

    test_score = clf.score(X_test, y_test)

    assert_less(abs(test_score - clf.oob_score_), 0.1)


def test_oob_score_regression():
    """Check that oob prediction is pessimistic estimate.
    Not really a good test that prediction is independent."""
    X_train, X_test, y_train, y_test = train_test_split(boston.data,
                                                        boston.target,
                                                        random_state=rng)

    clf = BaggingRegressor(base_estimator=DecisionTreeRegressor(),
                           n_estimators=50,
                           bootstrap=True,
                           oob_score=True,
                           random_state=rng).fit(X_train, y_train)

    test_score = clf.score(X_test, y_test)
    assert_greater(test_score, clf.oob_score_)
    assert_greater(clf.oob_score_, .8)



def test_error():
    """Test that it gives proper exception on deficient input."""
    X, y = iris.data, iris.target
    base = DecisionTreeClassifier()

    assert_raises(ValueError, BaggingClassifier(base, max_samples=-1).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_samples=0.0).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_samples=2.0).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_samples=1000).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_samples="foobar").fit, X, y)

    assert_raises(ValueError, BaggingClassifier(base, max_features=-1).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_features=0.0).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_features=2.0).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_features=5).fit, X, y)
    assert_raises(ValueError, BaggingClassifier(base, max_features="foobar").fit, X, y)


if __name__ == "__main__":
    import nose
    nose.runmodule()
