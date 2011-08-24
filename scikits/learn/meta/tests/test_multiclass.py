
import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import assert_raises

from scikits.learn.meta import OneVsRestClassifier
from scikits.learn.linear_model import SGDClassifier
from scikits.learn.naive_bayes import MultinomialNB
from scikits.learn.grid_search import GridSearchCV
from scikits.learn import datasets

iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

def test_ovr_exceptions():
    ovr = OneVsRestClassifier(SGDClassifier())
    assert_raises(ValueError, ovr.predict, [])


def test_ovr_fit_predict():
    # A classifier which implements decision_function.
    ovr = OneVsRestClassifier(SGDClassifier())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    pred2 = SGDClassifier().fit(iris.data, iris.target).predict(iris.data)
    # Note: this is gonna become a circular test if SGDClassifier
    #       is rewritten to use fit_ovr.
    assert_equal(np.mean(iris.target == pred), np.mean(iris.target == pred2))

    # A classifier which implements predict_proba.
    ovr = OneVsRestClassifier(MultinomialNB())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_true(np.mean(iris.target == pred) >= 0.65)


def test_ovr_gridsearch():
    ovr = OneVsRestClassifier(SGDClassifier())
    alphas = [0.0002, 0.0005, 0.0008]
    cv = GridSearchCV(ovr, {'estimator__alpha': alphas})
    cv.fit(iris.data, iris.target)
    best_alpha = cv.best_estimator.estimators_[0].alpha
    assert_true(best_alpha in alphas)
