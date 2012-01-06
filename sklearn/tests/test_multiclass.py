
import numpy as np

from numpy.testing import assert_array_equal
from nose.tools import assert_equal
from nose.tools import assert_almost_equal
from nose.tools import assert_true
from nose.tools import assert_raises

from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OutputCodeClassifier
from sklearn.svm import LinearSVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.grid_search import GridSearchCV
from sklearn import datasets

iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
n_classes = 3


# FIXME: - should use sets
#        - should move to metrics module
def multilabel_precision(Y_true, Y_pred):
    n_predictions = 0
    n_correct = 0
    for i in range(len(Y_true)):
        n_predictions += len(Y_pred[i])
        for label in Y_pred[i]:
            if label in Y_true[i]:
                n_correct += 1
    return float(n_correct) / n_predictions


def multilabel_recall(Y_true, Y_pred):
    n_labels = 0
    n_correct = 0
    for i in range(len(Y_true)):
        n_labels += len(Y_true[i])
        for label in Y_pred[i]:
            if label in Y_true[i]:
                n_correct += 1
    return float(n_correct) / n_labels


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


def test_ovr_multilabel():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 4, 5], [0, 5, 0], [3, 3, 3], [4, 0, 6], [6, 0, 0]])
    y = [[1, 2], [1], [0, 1, 2], [0, 2], [0]]
    Y = np.array([[0, 1, 1],
                  [0, 1, 0],
                  [1, 1, 1],
                  [1, 0, 1],
                  [1, 0, 0]])

    for base_clf in (MultinomialNB(), LinearSVC()):
        # test input as lists of tuples
        clf = OneVsRestClassifier(base_clf).fit(X, y)
        y_pred = clf.predict([[0, 4, 4]])[0]
        assert_equal(set(y_pred), set([1, 2]))
        assert_true(clf.multilabel_)

        # test input as label indicator matrix
        clf = OneVsRestClassifier(base_clf).fit(X, Y)
        y_pred = clf.predict([[0, 4, 4]])[0]
        assert_array_equal(y_pred, [0, 1, 1])
        assert_true(clf.multilabel_)


def test_ovr_multilabel_dataset():
    base_clf = MultinomialNB(alpha=1)
    for au, prec, recall in zip((True, False), (0.65, 0.74), (0.72, 0.84)):
        X, Y = datasets.make_multilabel_classification(n_samples=100,
                                                       n_features=20,
                                                       n_classes=5,
                                                       n_labels=2,
                                                       length=50,
                                                       allow_unlabeled=au,
                                                       random_state=0)
        X_train, Y_train = X[:80], Y[:80]
        X_test, Y_test = X[80:], Y[80:]
        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)
        Y_pred = clf.predict(X_test)
        assert_true(clf.multilabel_)
        assert_almost_equal(multilabel_precision(Y_test, Y_pred), prec,
                            places=2)
        assert_almost_equal(multilabel_recall(Y_test, Y_pred), recall,
                            places=2)


def test_ovr_gridsearch():
    ovr = OneVsRestClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovr, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def test_ovo_exceptions():
    ovo = OneVsOneClassifier(LinearSVC())
    assert_raises(ValueError, ovo.predict, [])


def test_ovo_fit_predict():
    # A classifier which implements decision_function.
    ovo = OneVsOneClassifier(LinearSVC())
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)

    # A classifier which implements predict_proba.
    ovo = OneVsOneClassifier(MultinomialNB())
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)


def test_ovo_gridsearch():
    ovo = OneVsOneClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovo, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def test_ecoc_exceptions():
    ecoc = OutputCodeClassifier(LinearSVC())
    assert_raises(ValueError, ecoc.predict, [])


def test_ecoc_fit_predict():
    # A classifier which implements decision_function.
    ecoc = OutputCodeClassifier(LinearSVC(), code_size=2)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)

    # A classifier which implements predict_proba.
    ecoc = OutputCodeClassifier(MultinomialNB(), code_size=2)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)


def test_ecoc_gridsearch():
    ecoc = OutputCodeClassifier(LinearSVC())
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ecoc, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)
