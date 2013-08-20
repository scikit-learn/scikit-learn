"""
Testing for the ensemble module (sklearn.ensemble.ensemble).
"""

# Authors: Gilles Louppe
# License: BSD 3 clause

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal
from nose.tools import assert_false, assert_true

from sklearn.utils.testing import assert_less, assert_greater

from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.grid_search import GridSearchCV, ParameterGrid
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import SGDClassifier, Perceptron
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR

from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_circles, load_boston
from sklearn.utils import check_random_state

rng = check_random_state(0)


def test_classification():
    """Check classificationfor various parameter settings."""
    X, y = make_circles(n_samples=100, noise=0.1, random_state=rng)
    X_train, X_test, y_train, y_test = train_test_split(X,
                                                        y,
                                                        random_state=rng)
    grid = ParameterGrid({"max_samples": [0.5, 1.0],
                          "max_features": [1, 2],
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

            score_base = base_estimator.score(X_test, y_test)
            score_ensemble = ensemble.score(X_test, y_test)


def test_regression():
    """Check regression for various parameter settings."""
    boston = load_boston()
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

            score_base = base_estimator.score(X_test, y_test)
            score_ensemble = ensemble.score(X_test, y_test)


if __name__ == "__main__":
    import nose
    nose.runmodule()
