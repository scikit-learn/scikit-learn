"""
General tests for all estimators in sklearn.
"""

from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.base import clone, ClassifierMixin
from sklearn.datasets import load_iris
from sklearn.metrics import zero_one_score

estimators = all_estimators()


def test_all_estimators():
    for name, E in estimators:
        e = E()
        # test cloning
        clone(e)
        # test __repr__
        print(e)


def test_classifiers():
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    iris = load_iris()
    X, y = iris.data, iris.target
    for clf in classifiers:
        clf.fit(X, y)
        y_pred = clf.predict(X)
        assert_greater(zero_one_score(y, y_pred), 0.7)
