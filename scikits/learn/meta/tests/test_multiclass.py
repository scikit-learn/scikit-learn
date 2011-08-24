
import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import assert_raises

from scikits.learn.meta import OneVsRestClassifier
from scikits.learn.meta import OneVsOneClassifier
from scikits.learn.meta import OutputCodeClassifier
from scikits.learn.svm import LinearSVC
from scikits.learn.naive_bayes import MultinomialNB
from scikits.learn.grid_search import GridSearchCV
from scikits.learn import datasets

iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
n_classes = 3


def test_ovr_exceptions():
    ovr = OneVsRestClassifier(LinearSVC())
    assert_raises(ValueError, ovr.predict, [])


def test_ovr_fit_predict():
    # A classifier which implements decision_function.
    ovr = OneVsRestClassifier(LinearSVC())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovr.estimators_), n_classes)

    pred2 = LinearSVC().fit(iris.data, iris.target).predict(iris.data)
    assert_equal(np.mean(iris.target == pred), np.mean(iris.target == pred2))

    # A classifier which implements predict_proba.
    ovr = OneVsRestClassifier(MultinomialNB())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_true(np.mean(iris.target == pred) >= 0.65)


def test_ovr_gridsearch():
    ovr = OneVsRestClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovr, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator.estimators_[0].C
    assert_true(best_C in Cs)


def test_ovo_exceptions():
    ovo = OneVsOneClassifier(LinearSVC())
    assert_raises(ValueError, ovo.predict, [])


def test_ovo_fit_predict():
    # A classifier which implements decision_function.
    ovo = OneVsOneClassifier(LinearSVC())
    pred = ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)

    # A classifier which implements predict_proba.
    ovo = OneVsOneClassifier(MultinomialNB())
    pred = ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)


def test_ovo_gridsearch():
    ovo = OneVsOneClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovo, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator.estimators_[0].C
    assert_true(best_C in Cs)


def test_ecoc_exceptions():
    ecoc = OutputCodeClassifier(LinearSVC())
    assert_raises(ValueError, ecoc.predict, [])


def test_ecoc_fit_predict():
    # A classifier which implements decision_function.
    ecoc = OutputCodeClassifier(LinearSVC(), code_size=2)
    pred = ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)

    # A classifier which implements predict_proba.
    ecoc = OutputCodeClassifier(MultinomialNB(), code_size=2)
    pred = ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)


def test_ecoc_gridsearch():
    ecoc = OutputCodeClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ecoc, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator.estimators_[0].C
    assert_true(best_C in Cs)
