"""
Testing for the bagging module (sklearn.ensemble.bagging).
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
from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import BaggingClassifier, BaggingRegressor
from sklearn.linear_model import SGDClassifier, Perceptron
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.svm import SVC, SVR

from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_classification



def test_toy():
    """Check classification on iris."""
    # Test with a classifier without predict_proba
    X, y = make_classification(n_samples=500, n_features=20, n_informative=5)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

    for base_estimator_class in [DummyClassifier, Perceptron, DecisionTreeClassifier, KNeighborsClassifier]:
        base_estimator = base_estimator_class()
        base_estimator.fit(X_train, y_train)

        bagging = BaggingClassifier(base_estimator=base_estimator_class(), n_estimators=1000, oob_score=False, random_state=0)
        bagging.fit(X_train, y_train)

        score_base = base_estimator.score(X_test, y_test)
        score_bagging = bagging.score(X_test, y_test)

        print base_estimator_class, score_base, score_bagging, score_base <= score_bagging



if __name__ == "__main__":
    import nose
    nose.runmodule()
