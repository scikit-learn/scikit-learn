"""
General tests for all estimators in sklearn.
"""

from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.base import clone, ClassifierMixin, MetaEstimatorMixin
from sklearn.datasets import load_iris
from sklearn.metrics import zero_one_score
from sklearn.lda import LDA
from sklearn.grid_search import GridSearchCV
from sklearn.decomposition import SparseCoder
from sklearn.pipeline import Pipeline


def test_all_estimators():
    estimators = all_estimators()
    clf = LDA()

    # some can just not be sensibly default constructed
    dont_test = [Pipeline, GridSearchCV, SparseCoder]
    for name, E in estimators:
        if E in dont_test:
            continue
        # test default-constructibility
        if issubclass(E, MetaEstimatorMixin):
            e = E(clf)
        else:
            e = E()
        #test cloning
        clone(e)
        ## test __repr__
        #print(e)


def test_classifiers():
    estimators = all_estimators()
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    iris = load_iris()
    X, y = iris.data, iris.target
    for name, clf in classifiers:
        clf.fit(X, y)
        y_pred = clf.predict(X)
        assert_greater(zero_one_score(y, y_pred), 0.7)
