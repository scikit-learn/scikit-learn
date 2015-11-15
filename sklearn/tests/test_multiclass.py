import numpy as np
import scipy.sparse as sp
import sys
import re
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raise_message
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OutputCodeClassifier
from sklearn.utils.multiclass import check_classification_targets, type_of_target
from sklearn.utils.testing import assert_allclose
from sklearn.utilities import convert_probabilities, probs_to_class
from sklearn.utilities import get_random_numbers
from sklearn.multiclass import RakelClassifier, _get_possibility
from sklearn.multiclass import _valid_possibility, _binomialCoeff
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.powerset import Powerset

from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

from sklearn.svm import LinearSVC, SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import (LinearRegression, Lasso, ElasticNet, Ridge,
                                  Perceptron, LogisticRegression)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn import svm
from sklearn import datasets
from sklearn.externals.six.moves import zip

iris = datasets.load_iris()
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
n_classes = 3


def test_ovr_exceptions():
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    assert_raises(ValueError, ovr.predict, [])

    # Fail on multioutput data
    assert_raises(ValueError, OneVsRestClassifier(MultinomialNB()).fit,
                  np.array([[1, 0], [0, 1]]),
                  np.array([[1, 2], [3, 1]]))
    assert_raises(ValueError, OneVsRestClassifier(MultinomialNB()).fit,
                  np.array([[1, 0], [0, 1]]),
                  np.array([[1.5, 2.4], [3.1, 0.8]]))


def test_check_classification_targets():
    # Test that check_classification_target return correct type. #5782
    y = np.array([0.0, 1.1, 2.0, 3.0])
    msg = type_of_target(y)
    assert_raise_message(ValueError, msg, check_classification_targets, y)


def test_ovr_fit_predict():
    # A classifier which implements decision_function.
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovr.estimators_), n_classes)

    clf = LinearSVC(random_state=0)
    pred2 = clf.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(np.mean(iris.target == pred), np.mean(iris.target == pred2))

    # A classifier which implements predict_proba.
    ovr = OneVsRestClassifier(MultinomialNB())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_greater(np.mean(iris.target == pred), 0.65)


def test_ovr_ovo_regressor():
    # test that ovr and ovo work on regressors which don't have a decision_function
    ovr = OneVsRestClassifier(DecisionTreeRegressor())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovr.estimators_), n_classes)
    assert_array_equal(np.unique(pred), [0, 1, 2])
    # we are doing something sensible
    assert_greater(np.mean(pred == iris.target), .9)

    ovr = OneVsOneClassifier(DecisionTreeRegressor())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovr.estimators_), n_classes * (n_classes - 1) / 2)
    assert_array_equal(np.unique(pred), [0, 1, 2])
    # we are doing something sensible
    assert_greater(np.mean(pred == iris.target), .9)


def test_ovr_fit_predict_sparse():
    for sparse in [sp.csr_matrix, sp.csc_matrix, sp.coo_matrix, sp.dok_matrix,
                   sp.lil_matrix]:
        base_clf = MultinomialNB(alpha=1)

        X, Y = datasets.make_multilabel_classification(n_samples=100,
                                                       n_features=20,
                                                       n_classes=5,
                                                       n_labels=3,
                                                       length=50,
                                                       allow_unlabeled=True,
                                                       random_state=0)

        X_train, Y_train = X[:80], Y[:80]
        X_test = X[80:]

        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)
        Y_pred = clf.predict(X_test)

        clf_sprs = OneVsRestClassifier(base_clf).fit(X_train, sparse(Y_train))
        Y_pred_sprs = clf_sprs.predict(X_test)

        assert_true(clf.multilabel_)
        assert_true(sp.issparse(Y_pred_sprs))
        assert_array_equal(Y_pred_sprs.toarray(), Y_pred)

        # Test predict_proba
        Y_proba = clf_sprs.predict_proba(X_test)

        # predict assigns a label if the probability that the
        # sample has the label is greater than 0.5.
        pred = Y_proba > .5
        assert_array_equal(pred, Y_pred_sprs.toarray())

        # Test decision_function
        clf_sprs = OneVsRestClassifier(svm.SVC()).fit(X_train, sparse(Y_train))
        dec_pred = (clf_sprs.decision_function(X_test) > 0).astype(int)
        assert_array_equal(dec_pred, clf_sprs.predict(X_test).toarray())


def test_ovr_always_present():
    # Test that ovr works with classes that are always present or absent.
    # Note: tests is the case where _ConstantPredictor is utilised
    X = np.ones((10, 2))
    X[:5, :] = 0

    # Build an indicator matrix where two features are always on.
    # As list of lists, it would be: [[int(i >= 5), 2, 3] for i in range(10)]
    y = np.zeros((10, 3))
    y[5:, 0] = 1
    y[:, 1] = 1
    y[:, 2] = 1

    ovr = OneVsRestClassifier(LogisticRegression())
    assert_warns(UserWarning, ovr.fit, X, y)
    y_pred = ovr.predict(X)
    assert_array_equal(np.array(y_pred), np.array(y))
    y_pred = ovr.decision_function(X)
    assert_equal(np.unique(y_pred[:, -2:]), 1)
    y_pred = ovr.predict_proba(X)
    assert_array_equal(y_pred[:, -1], np.ones(X.shape[0]))

    # y has a constantly absent label
    y = np.zeros((10, 2))
    y[5:, 0] = 1  # variable label
    ovr = OneVsRestClassifier(LogisticRegression())
    assert_warns(UserWarning, ovr.fit, X, y)
    y_pred = ovr.predict_proba(X)
    assert_array_equal(y_pred[:, -1], np.zeros(X.shape[0]))


def test_ovr_multiclass():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 0, 5], [0, 5, 0], [3, 0, 0], [0, 0, 6], [6, 0, 0]])
    y = ["eggs", "spam", "ham", "eggs", "ham"]
    Y = np.array([[0, 0, 1],
                  [0, 1, 0],
                  [1, 0, 0],
                  [0, 0, 1],
                  [1, 0, 0]])

    classes = set("ham eggs spam".split())

    for base_clf in (MultinomialNB(), LinearSVC(random_state=0),
                     LinearRegression(), Ridge(),
                     ElasticNet()):

        clf = OneVsRestClassifier(base_clf).fit(X, y)
        assert_equal(set(clf.classes_), classes)
        y_pred = clf.predict(np.array([[0, 0, 4]]))[0]
        assert_equal(set(y_pred), set("eggs"))

        # test input as label indicator matrix
        clf = OneVsRestClassifier(base_clf).fit(X, Y)
        y_pred = clf.predict([[0, 0, 4]])[0]
        assert_array_equal(y_pred, [0, 0, 1])


def test_ovr_binary():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 0, 5], [0, 5, 0], [3, 0, 0], [0, 0, 6], [6, 0, 0]])
    y = ["eggs", "spam", "spam", "eggs", "spam"]
    Y = np.array([[0, 1, 1, 0, 1]]).T

    classes = set("eggs spam".split())

    def conduct_test(base_clf, test_predict_proba=False):
        clf = OneVsRestClassifier(base_clf).fit(X, y)
        assert_equal(set(clf.classes_), classes)
        y_pred = clf.predict(np.array([[0, 0, 4]]))[0]
        assert_equal(set(y_pred), set("eggs"))

        if test_predict_proba:
            X_test = np.array([[0, 0, 4]])
            probabilities = clf.predict_proba(X_test)
            assert_equal(2, len(probabilities[0]))
            assert_equal(clf.classes_[np.argmax(probabilities, axis=1)],
                         clf.predict(X_test))

        # test input as label indicator matrix
        clf = OneVsRestClassifier(base_clf).fit(X, Y)
        y_pred = clf.predict([[3, 0, 0]])[0]
        assert_equal(y_pred, 1)

    for base_clf in (LinearSVC(random_state=0), LinearRegression(),
                     Ridge(), ElasticNet()):
        conduct_test(base_clf)

    for base_clf in (MultinomialNB(), SVC(probability=True),
                     LogisticRegression()):
        conduct_test(base_clf, test_predict_proba=True)


def test_ovr_multilabel():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 4, 5], [0, 5, 0], [3, 3, 3], [4, 0, 6], [6, 0, 0]])
    y = np.array([[0, 1, 1],
                  [0, 1, 0],
                  [1, 1, 1],
                  [1, 0, 1],
                  [1, 0, 0]])

    for base_clf in (MultinomialNB(), LinearSVC(random_state=0),
                     LinearRegression(), Ridge(),
                     ElasticNet(), Lasso(alpha=0.5)):
        clf = OneVsRestClassifier(base_clf).fit(X, y)
        y_pred = clf.predict([[0, 4, 4]])[0]
        assert_array_equal(y_pred, [0, 1, 1])
        assert_true(clf.multilabel_)


def test_ovr_fit_predict_svc():
    ovr = OneVsRestClassifier(svm.SVC())
    ovr.fit(iris.data, iris.target)
    assert_equal(len(ovr.estimators_), 3)
    assert_greater(ovr.score(iris.data, iris.target), .9)


def test_ovr_multilabel_dataset():
    base_clf = MultinomialNB(alpha=1)
    for au, prec, recall in zip((True, False), (0.51, 0.66), (0.51, 0.80)):
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
        assert_almost_equal(precision_score(Y_test, Y_pred, average="micro"),
                            prec,
                            decimal=2)
        assert_almost_equal(recall_score(Y_test, Y_pred, average="micro"),
                            recall,
                            decimal=2)


def test_ovr_multilabel_predict_proba():
    base_clf = MultinomialNB(alpha=1)
    for au in (False, True):
        X, Y = datasets.make_multilabel_classification(n_samples=100,
                                                       n_features=20,
                                                       n_classes=5,
                                                       n_labels=3,
                                                       length=50,
                                                       allow_unlabeled=au,
                                                       random_state=0)
        X_train, Y_train = X[:80], Y[:80]
        X_test = X[80:]
        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)

        # decision function only estimator. Fails in current implementation.
        decision_only = OneVsRestClassifier(svm.SVR()).fit(X_train, Y_train)
        assert_raises(AttributeError, decision_only.predict_proba, X_test)

        # Estimator with predict_proba disabled, depending on parameters.
        decision_only = OneVsRestClassifier(svm.SVC(probability=False))
        decision_only.fit(X_train, Y_train)
        assert_raises(AttributeError, decision_only.predict_proba, X_test)

        Y_pred = clf.predict(X_test)
        Y_proba = clf.predict_proba(X_test)

        # predict assigns a label if the probability that the
        # sample has the label is greater than 0.5.
        pred = Y_proba > .5
        assert_array_equal(pred, Y_pred)


def test_ovr_single_label_predict_proba():
    base_clf = MultinomialNB(alpha=1)
    X, Y = iris.data, iris.target
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)

    # decision function only estimator. Fails in current implementation.
    decision_only = OneVsRestClassifier(svm.SVR()).fit(X_train, Y_train)
    assert_raises(AttributeError, decision_only.predict_proba, X_test)

    Y_pred = clf.predict(X_test)
    Y_proba = clf.predict_proba(X_test)

    assert_almost_equal(Y_proba.sum(axis=1), 1.0)
    # predict assigns a label if the probability that the
    # sample has the label is greater than 0.5.
    pred = np.array([l.argmax() for l in Y_proba])
    assert_false((pred - Y_pred).any())


def test_ovr_multilabel_decision_function():
    X, Y = datasets.make_multilabel_classification(n_samples=100,
                                                   n_features=20,
                                                   n_classes=5,
                                                   n_labels=3,
                                                   length=50,
                                                   allow_unlabeled=True,
                                                   random_state=0)
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(svm.SVC()).fit(X_train, Y_train)
    assert_array_equal((clf.decision_function(X_test) > 0).astype(int),
                       clf.predict(X_test))


def test_ovr_single_label_decision_function():
    X, Y = datasets.make_classification(n_samples=100,
                                        n_features=20,
                                        random_state=0)
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(svm.SVC()).fit(X_train, Y_train)
    assert_array_equal(clf.decision_function(X_test).ravel() > 0,
                       clf.predict(X_test))


def test_ovr_gridsearch():
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovr, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def test_ovr_pipeline():
    # Test with pipeline of length one
    # This test is needed because the multiclass estimators may fail to detect
    # the presence of predict_proba or decision_function.
    clf = Pipeline([("tree", DecisionTreeClassifier())])
    ovr_pipe = OneVsRestClassifier(clf)
    ovr_pipe.fit(iris.data, iris.target)
    ovr = OneVsRestClassifier(DecisionTreeClassifier())
    ovr.fit(iris.data, iris.target)
    assert_array_equal(ovr.predict(iris.data), ovr_pipe.predict(iris.data))


def test_ovr_coef_():
    for base_classifier in [SVC(kernel='linear', random_state=0), LinearSVC(random_state=0)]:
        # SVC has sparse coef with sparse input data

        ovr = OneVsRestClassifier(base_classifier)
        for X in [iris.data, sp.csr_matrix(iris.data)]:
            # test with dense and sparse coef
            ovr.fit(X, iris.target)
            shape = ovr.coef_.shape
            assert_equal(shape[0], n_classes)
            assert_equal(shape[1], iris.data.shape[1])
            # don't densify sparse coefficients
            assert_equal(sp.issparse(ovr.estimators_[0].coef_), sp.issparse(ovr.coef_))


def test_ovr_coef_exceptions():
    # Not fitted exception!
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    # lambda is needed because we don't want coef_ to be evaluated right away
    assert_raises(ValueError, lambda x: ovr.coef_, None)

    # Doesn't have coef_ exception!
    ovr = OneVsRestClassifier(DecisionTreeClassifier())
    ovr.fit(iris.data, iris.target)
    assert_raises(AttributeError, lambda x: ovr.coef_, None)


def test_ovo_exceptions():
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    assert_raises(ValueError, ovo.predict, [])


def test_ovo_fit_on_list():
    # Test that OneVsOne fitting works with a list of targets and yields the
    # same output as predict from an array
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    prediction_from_array = ovo.fit(iris.data, iris.target).predict(iris.data)
    prediction_from_list = ovo.fit(iris.data,
                                   list(iris.target)).predict(iris.data)
    assert_array_equal(prediction_from_array, prediction_from_list)


def test_ovo_fit_predict():
    # A classifier which implements decision_function.
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)

    # A classifier which implements predict_proba.
    ovo = OneVsOneClassifier(MultinomialNB())
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)


def test_ovo_decision_function():
    n_samples = iris.data.shape[0]

    ovo_clf = OneVsOneClassifier(LinearSVC(random_state=0))
    ovo_clf.fit(iris.data, iris.target)
    decisions = ovo_clf.decision_function(iris.data)

    assert_equal(decisions.shape, (n_samples, n_classes))
    assert_array_equal(decisions.argmax(axis=1), ovo_clf.predict(iris.data))

    # Compute the votes
    votes = np.zeros((n_samples, n_classes))

    k = 0
    for i in range(n_classes):
        for j in range(i + 1, n_classes):
            pred = ovo_clf.estimators_[k].predict(iris.data)
            votes[pred == 0, i] += 1
            votes[pred == 1, j] += 1
            k += 1

    # Extract votes and verify
    assert_array_equal(votes, np.round(decisions))

    for class_idx in range(n_classes):
        # For each sample and each class, there only 3 possible vote levels
        # because they are only 3 distinct class pairs thus 3 distinct
        # binary classifiers.
        # Therefore, sorting predictions based on votes would yield
        # mostly tied predictions:
        assert_true(set(votes[:, class_idx]).issubset(set([0., 1., 2.])))

        # The OVO decision function on the other hand is able to resolve
        # most of the ties on this data as it combines both the vote counts
        # and the aggregated confidence levels of the binary classifiers
        # to compute the aggregate decision function. The iris dataset
        # has 150 samples with a couple of duplicates. The OvO decisions
        # can resolve most of the ties:
        assert_greater(len(np.unique(decisions[:, class_idx])), 146)


def test_ovo_gridsearch():
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovo, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def test_ovo_ties():
    # Test that ties are broken using the decision function,
    # not defaulting to the smallest label
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y = np.array([2, 0, 1, 2])
    multi_clf = OneVsOneClassifier(Perceptron(shuffle=False))
    ovo_prediction = multi_clf.fit(X, y).predict(X)
    ovo_decision = multi_clf.decision_function(X)

    # Classifiers are in order 0-1, 0-2, 1-2
    # Use decision_function to compute the votes and the normalized
    # sum_of_confidences, which is used to disambiguate when there is a tie in
    # votes.
    votes = np.round(ovo_decision)
    normalized_confidences = ovo_decision - votes

    # For the first point, there is one vote per class
    assert_array_equal(votes[0, :], 1)
    # For the rest, there is no tie and the prediction is the argmax
    assert_array_equal(np.argmax(votes[1:], axis=1), ovo_prediction[1:])
    # For the tie, the prediction is the class with the highest score
    assert_equal(ovo_prediction[0], normalized_confidences[0].argmax())


def test_ovo_ties2():
    # test that ties can not only be won by the first two labels
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y_ref = np.array([2, 0, 1, 2])

    # cycle through labels so that each label wins once
    for i in range(3):
        y = (y_ref + i) % 3
        multi_clf = OneVsOneClassifier(Perceptron(shuffle=False))
        ovo_prediction = multi_clf.fit(X, y).predict(X)
        assert_equal(ovo_prediction[0], i % 3)


def test_ovo_string_y():
    # Test that the OvO doesn't mess up the encoding of string labels
    X = np.eye(4)
    y = np.array(['a', 'b', 'c', 'd'])

    ovo = OneVsOneClassifier(LinearSVC())
    ovo.fit(X, y)
    assert_array_equal(y, ovo.predict(X))


def test_ecoc_exceptions():
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0))
    assert_raises(ValueError, ecoc.predict, [])


def test_ecoc_fit_predict():
    # A classifier which implements decision_function.
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                code_size=2, random_state=0)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)

    # A classifier which implements predict_proba.
    ecoc = OutputCodeClassifier(MultinomialNB(), code_size=2, random_state=0)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 2)


def test_ecoc_gridsearch():
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                random_state=0)
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ecoc, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def count_identical(array, container):
    retv = 0
    for ar in container:
        if (array == ar).all():
            retv = retv + 1
    return retv


def assert_labelsets(rak):
    assert_equal(len(rak.labelsets_), rak.n_estimators)
    assert_equal(type(rak.labelsets_), type([]))
    for l in rak.labelsets_:
        if (rak.k < 1):
            assert_equal(np.sum(l), int(rak.k * rak.n_labels_))
        else:
            assert_equal(np.sum(l), rak.k)
        assert_equal(type(l), type(np.zeros((10,))))
    if rak.no_hole:
        r = np.array([sum(i) for i in zip(*rak.labelsets_)])
        assert_equal(rak.n_labels_, np.count_nonzero(r))
    if rak.uniqueness:
        for l in rak.labelsets_:
            assert_equal(1, count_identical(l, rak.labelsets_))


def test_rakel_binomialCoeff():
    assert_equal(1, _binomialCoeff(0, 0))
    assert_equal(1, _binomialCoeff(1, 1))
    assert_equal(1, _binomialCoeff(1, 0))
    assert_equal(1, _binomialCoeff(5, 10))
    assert_equal((10 * 9 * 8 * 7 * 6) / (5 * 4 * 3 * 2),
                 _binomialCoeff(10, 5))
    assert_equal((10 * 9 * 8 * 7) / (4 * 3 * 2),
                 _binomialCoeff(10, 6))
    assert_equal((10 * 9 * 8 * 7 * 6 * 5) / (6 * 5 * 4 * 3 * 2),
                 _binomialCoeff(10, 4))


def test_rakel_get_possibility():
    array = [['0', '1', '2'],
             ['0', '2', '1'],
             ['1', '0', '2'],  # 2
             ['1', '2', '0'],
             ['2', '0', '1'],  # 4
             ['2', '1', '0']]
    assert_raises(ValueError, _get_possibility, array[3], 10, array)
    del array[4]
    del array[2]
    assert_array_equal(['1', '0', '2'],
                       _get_possibility(array[0], 2, array))
    assert_array_equal(['1', '0', '2'],
                       _get_possibility(array[1], 1, array))
    assert_array_equal(['2', '0', '1'],
                       _get_possibility(array[2], 1, array))
    assert_array_equal(['1', '0', '2'],
                       _get_possibility(array[3], 3, array))
    assert_raises(ValueError, _get_possibility, array[3], 2, array)


def test_rakel_labelsets():
    random_state = 0

    # init
    n_features = 10
    n_labels = 2
    n_samples = 10
    n_classes = 10
    length = 1
    X, y = datasets.make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state, return_indicator=True)
    rakel = RakelClassifier(DecisionTreeClassifier(random_state=random_state),
                            random_state=random_state, no_hole=False,
                            k=1, n_estimators=10, uniqueness=True)
    assert_array_equal([], rakel._get_labelsets(n_labelsets=11,
                                                n_labels=10,
                                                size_labelsets=11))

    rakel.fit(X, y)
    assert_labelsets(rakel)
    rakel.n_estimators = 0
    assert_raises(ValueError, rakel.fit, X, y)
    rakel.n_estimators = -1
    assert_raises(ValueError, rakel.fit, X, y)
    rakel.n_estimators = 10
    rakel.k = 2
    rakel.fit(X, y)
    assert_true(_valid_possibility(rakel.labelsets_[-1],
                                   rakel.labelsets_[:-1]))
    assert_labelsets(rakel)
    rakel.k = 11
    assert_raise_message(ValueError, "Invalid size", rakel.fit, X, y)
    rakel.k = 2
    rakel.n_estimators = 100
    rakel.fit(X, y)
    assert_equal(len(rakel.labelsets_), n_classes * (n_classes - 1) * 0.5)
    assert_raise_message(ValueError, "duplicate",
                         rakel._get_labelsets, rakel.n_estimators,
                         rakel.n_labels_, rakel.k, random_state=None,
                         uniqueness=True)
    assert_raise_message(ValueError, "duplicate",
                         rakel._get_labelsets, 2,
                         rakel.n_labels_, 0, random_state=None,
                         uniqueness=True)
    assert_array_equal(rakel._get_labelsets(4, 3, 0, uniqueness=False),
                       [np.zeros((3,)) for _ in range(4)])
    rakel.k = 5
    rakel.n_estimators = 252
    rakel.fit(X, y)
    assert_labelsets(rakel)
    rakel.n_estimators = 140
    rakel._ACCEPTABLE_COLLISION = 0.99
    rakel.fit(X, y)
    assert_labelsets(rakel)
    rakel._ACCEPTABLE_COLLISION = 0.5
    rakel.k = 2
    rakel.n_estimators = 2
    rakel.fit(X, y)
    assert_labelsets(rakel)
    rakel.no_hole = True
    assert_raise_message(ValueError, "no hole",
                         rakel._force_no_hole,
                         rakel.labelsets_)
    assert_raise_message(ValueError, "no hole", rakel.fit, X, y)
    rakel.n_estimators = 10
    rakel.fit(X, y)
    assert_labelsets(rakel)
    rakel.uniqueness = False
    rakel.n_estimators = 30
    rakel.fit(X, y)
    assert_labelsets(rakel)
    labs = rakel.labelsets_
    rakel.k = float(rakel.k) / rakel.n_labels_
    rakel.fit(X, y)
    assert_labelsets(rakel)
    assert_array_equal(labs, rakel.labelsets_)


def test_rakel_fit():
    # If fit fails, all computed properties should be reset
    random_state = 0

    # init
    n_features = 10
    n_labels = 2
    n_samples = 10
    n_classes = 10
    length = 1
    X, y = datasets.make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state, return_indicator=True)
    Xt, yt = datasets.make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state + n_samples, return_indicator=True)
    rakel = RakelClassifier(DecisionTreeClassifier(random_state=random_state),
                            random_state=random_state, no_hole=False,
                            k=1, n_estimators=10, uniqueness=True)
    rakel.fit(X, y, copy=False)
    assert_array_equal(rakel.X_, X)
    assert_array_equal(rakel.y_, y)
    assert_equal(rakel.n_labels_, n_classes)
    assert_true(rakel.multi_label_)
    assert_equal(len(rakel.labelsets_), 10)
    assert_true(rakel is rakel.fit(X, y, copy=False))
    rakel.k = 11
    assert_raises(ValueError, rakel.fit, X, y, copy=False)
    assert_array_equal(rakel.X_, None)
    assert_array_equal(rakel.y_, None)
    assert_equal(rakel.n_labels_, 0)
    assert_equal(rakel.labelsets_, None)
    assert_raises(ValueError, rakel.fit, X, y, copy=True)
    assert_array_equal(rakel.X_, None)
    assert_array_equal(rakel.y_, None)
    assert_equal(rakel.n_labels_, 0)
    assert_equal(rakel.labelsets_, None)
    rakel.k = 1
    y_ = y[:, 9].ravel()
    rakel.fit(X, y_, copy=False)
    assert_array_equal(rakel.X_, X)
    assert_array_equal(rakel.y_, y_)
    assert_equal(rakel.n_labels_, 1)
    assert_false(rakel.multi_label_)
    assert_equal(len(rakel.labelsets_), 1)
    rakel.fit(X, y, copy=True)
    assert_array_equal(rakel.X_, X)
    assert_array_equal(rakel.y_, y)
    assert_equal(rakel.n_labels_, n_classes)
    assert_true(rakel.multi_label_)
    assert_equal(len(rakel.labelsets_), 10)
    X[0, 0] = 1 - X[0, 0]
    y[0, 0] = 1 - y[0, 0]
    try:
        assert_array_equal(rakel.X_, X)
    except AssertionError:
        pass
    else:
        raise AssertionError("The X copy did not work.")
    try:
        assert_array_equal(rakel.y_, y)
    except AssertionError:
        pass
    else:
        raise AssertionError("The y copy did not work.")


def test_rakel_predicts():
    random_state = 0

    # init
    n_features = 10
    n_labels = 2
    n_samples = 10
    n_classes = 5
    length = 1
    X, y = datasets.make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state, return_indicator=True)
    Xt, yt = datasets.make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state + n_samples, return_indicator=True)
    tree = DecisionTreeClassifier(random_state=get_random_numbers(
                                  random_state=random_state))
    tree.fit(X, y)
    y_predict = tree.predict(Xt)
    y_probas = tree.predict_proba(Xt)
    y_r = y[:, 0]
    p = Powerset()
    p_r = Powerset()
    single = p_r.compressed_powerize(y_r)
    comp = p.compressed_powerize(y)
    tree.fit(X, comp)
    rm1_tree_yt = tree.predict(Xt)
    rm1_tree_ytp = tree.predict_proba(Xt)
    tree.fit(X, y_r)
    rm1r_tree_yt = tree.predict(Xt)
    rm1r_tree_ytp = tree.predict_proba(Xt)
    tree.fit(X, single)
    rm1r_tree_yt_ = tree.predict(Xt)
    rm1r_tree_ytp_ = tree.predict_proba(Xt)
    assert_array_equal(rm1r_tree_yt_, rm1r_tree_yt)
    assert_array_equal(y_predict, probs_to_class(y_probas))
    assert_allclose(rm1r_tree_ytp, rm1r_tree_ytp_, rtol=1e-10)

    rm1 = RakelClassifier(tree, k=1, n_estimators=1,
                          threshold=0.5, no_hole=True,
                          powerset="probabilities",
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(Exception, "initialized", rm1.predict_proba, Xt)
    assert_array_equal(rm1_tree_yt, probs_to_class(rm1_tree_ytp))

    rm1.fit(X, y_r)
    rm1.fill_method = "i am incorrect"
    assert_raises(ValueError, rm1.predict, Xt)
    rm1.fill_method = -0.2
    assert_raises(ValueError, rm1.predict, Xt)
    rm1.fill_method = 1.1
    assert_raises(ValueError, rm1.predict, Xt)
    rm1.fill_method = "most_frequent"
    res1 = rm1.predict(Xt)
    assert_array_equal(p_r.unpowerize(rm1r_tree_yt_), res1)
    assert_allclose(rm1r_tree_ytp_, rm1.predict_proba(Xt), rtol=1e-10)

    rm1.k = n_classes
    rm1.fit(X, y)
    saved_stdout = sys.stdout
    try:
        out = StringIO()
        sys.stdout = out
        rm1.verbose = 5
        res1 = rm1.predict(Xt)
        output = out.getvalue().strip()
        if not re.compile("Number of estimators: \\d+.").search(output):
            raise AssertionError("Did not print number of estimators.")
        if not re.compile("labelset: ").search(output):
            raise AssertionError("Did not print labelsets.")
        rl = RakelClassifier(tree, k=1, n_estimators=1,
                             threshold=0.5, no_hole=False,
                             powerset="probabilities",
                             random_state=random_state, verbose=5)
        rl.fit(X, y)
        out = StringIO()
        sys.stdout = out
        res1 = rl.predict(Xt)
        output = out.getvalue().strip()
        if not re.compile("cannot be estimated").search(output):
            raise AssertionError("Did not print not estimated (predict).")
        out = StringIO()
        sys.stdout = out
        res1 = rl.predict_proba(Xt)
        output = out.getvalue().strip()
        if not re.compile("cannot be estimated").search(output):
            raise AssertionError("Did not print not estimated (probas).")
    finally:
        sys.stdout = saved_stdout
    rm1.verbose = 0
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    trees = ExtraTreesClassifier()
    assert_array_equal(p.unpowerize(rm1_tree_yt), res1)
    assert_allclose(p.probas_unpowerize(rm1_tree_ytp), resp1, rtol=1e-10)
    assert_allclose(
        y_probas, convert_probabilities(y_probas, yt.shape), rtol=1e-10)
    rm1 = RakelClassifier(trees, k=0.05, n_estimators=30, threshold=0.5,
                          powerset="probabilities", no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    resp1 = rm1.predict_proba(Xt)
    rm1 = RakelClassifier(trees, k=0.05, n_estimators=30, threshold=0.5,
                          powerset="probabilities", no_hole=True,
                          random_state=random_state, verbose=0,
                          fill_method=1.0)
    rm1.fit(X, y)
    assert_allclose(resp1, rm1.predict_proba(Xt), rtol=1e-10)
    rm1 = RakelClassifier(trees, k=0.05, n_estimators=30, threshold=0.5,
                          powerset="probabilities", no_hole=True,
                          random_state=random_state, verbose=0,
                          fill_method=0.0)
    rm1.fit(X, y)
    assert_allclose(resp1, rm1.predict_proba(Xt), rtol=1e-10)
    rm1 = RakelClassifier(trees, k=0.05, n_estimators=30, threshold=0.5,
                          powerset="probabilities", no_hole=True,
                          random_state=random_state, verbose=0,
                          fill_method=0.5)
    rm1.fit(X, y)
    assert_allclose(resp1, rm1.predict_proba(Xt), rtol=1e-10)

    rm1 = RakelClassifier(tree, k=n_classes, n_estimators=1,
                          powerset="probabilities",
                          threshold=0.5, no_hole=False,
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(
        Exception, "initialized", rm1.predict_proba, Xt)
    p = Powerset()
    comp = p.compressed_powerize(y)
    rm1_tree_y = tree.fit(X, comp)
    rm1_tree_yt = tree.predict(Xt)
    rm1_tree_ytp = tree.predict_proba(Xt)

    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(p.unpowerize(rm1_tree_yt), res1)
    assert_allclose(p.probas_unpowerize(rm1_tree_ytp), resp1, rtol=1e-10)
    assert_allclose(
        y_probas, convert_probabilities(y_probas, yt.shape), rtol=1e-10)

    rm1 = RakelClassifier(tree, k=0.9999999999999, n_estimators=1,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(
        Exception, "initialized", rm1.predict_proba, Xt)
    p = Powerset()
    comp = p.compressed_powerize(y)
    rm1_tree_y = tree.fit(X, comp)
    rm1_tree_yt = tree.predict(Xt)
    rm1_tree_ytp = tree.predict_proba(Xt)
    del rm1_tree_y

    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(p.unpowerize(rm1_tree_yt), res1)
    assert_allclose(p.probas_unpowerize(rm1_tree_ytp), resp1, rtol=1e-10)
    assert_allclose(
        y_probas, convert_probabilities(y_probas, yt.shape), rtol=1e-10)

    base = DummyClassifier(strategy="most_frequent")
    rm1 = RakelClassifier(base, k=1, n_estimators=1, threshold=0.5,
                          powerset="probabilities",
                          no_hole=False, random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(
        Exception, "initialized", rm1.predict_proba, Xt)
    base.fit(X, y)
    base_predict = base.predict(Xt)
    base_probas = base.predict_proba(Xt)
    rm1.fit(X, y)
    rm1.labelsets_ = []
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(base_predict, res1)
    for label in resp1:
        for i in range(label.shape[0]):
            assert_allclose(label[0, :], label[i, :], rtol=1e-10)

    # Test no error with absence of random_state
    rm1.estimator = OneVsRestClassifier(tree)
    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    # Test with majority unpowerset
    rm1 = RakelClassifier(base, k=1, n_estimators=n_classes, threshold=0.5,
                          powerset="majority",
                          no_hole=True, fill_method=0.5,
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(Exception, "initialized", rm1.predict_proba, Xt)
    base_sum = np.sum(y, axis=0) * 1.0
    result = np.asarray([np.round(base_sum / n_samples)
                         for _ in range(n_samples)])
    base_probas = convert_probabilities(result, yt.shape)
    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(result, res1)
    assert_allclose(
        convert_probabilities(base_probas, yt.shape), resp1, rtol=1e-10)
    #
    rm1 = RakelClassifier(base, k=1, n_estimators=1, threshold=0.5,
                          powerset="probabilities",
                          no_hole=False, fill_method=0.5,
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(
        Exception, "initialized", rm1.predict_proba, Xt)
    base_predict[:] = 0.5
    base_probas = convert_probabilities(base_probas, yt.shape)
    for label in base_probas:
        it = np.nditer(label, flags=['multi_index'])
        while not it.finished:
            label[it.multi_index] = 0.5
            it.iternext()
    rm1.fit(X, y)
    rm1.labelsets_ = []
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(base_predict, res1)
    assert_allclose(
        convert_probabilities(base_probas, yt.shape), resp1, rtol=1e-10)

    rm1 = RakelClassifier(tree, k=1, n_estimators=n_classes,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    assert_raise_message(Exception, "initialized", rm1.predict, Xt)
    assert_raise_message(Exception, "initialized", rm1.predict_proba, Xt)
    base = OneVsRestClassifier(estimator=tree, n_jobs=1)
    base.fit(X, y)
    base_predict = base.predict(Xt)
    base_probas = base.predict_proba(Xt)
    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    assert_array_equal(base_predict, res1)
    assert_allclose(
        convert_probabilities(base_probas, yt.shape), resp1, rtol=1e-10)

    rm1 = RakelClassifier(tree, k=0.2, n_estimators=n_classes, fill_method=0,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)

    rm1 = RakelClassifier(tree, k=0.2, n_estimators=n_classes, fill_method=1,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res2 = rm1.predict(Xt)
    resp2 = rm1.predict_proba(Xt)
    assert_array_equal(res1, res2)
    assert_allclose(resp1, resp2, rtol=1e-10)
    rm1 = RakelClassifier(tree, k=0.2, n_estimators=n_classes, fill_method=0.5,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res2 = rm1.predict(Xt)
    resp2 = rm1.predict_proba(Xt)
    assert_array_equal(res1, res2)
    assert_allclose(resp1, resp2, rtol=1e-10)
    rm1 = RakelClassifier(tree, k=0.2, n_estimators=n_classes,
                          powerset="probabilities",
                          fill_method="most_frequent",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res2 = rm1.predict(Xt)
    resp2 = rm1.predict_proba(Xt)
    assert_array_equal(res1, res2)
    assert_allclose(resp1, resp2, rtol=1e-10)

    rm1 = RakelClassifier(tree, k=0.4, n_estimators=30, fill_method=0,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res1 = rm1.predict(Xt)
    resp1 = rm1.predict_proba(Xt)
    rm1 = RakelClassifier(tree, k=0.4, n_estimators=30, fill_method=1,
                          powerset="probabilities",
                          threshold=0.5, no_hole=True,
                          random_state=random_state, verbose=0)
    rm1.fit(X, y)
    res2 = rm1.predict(Xt)
    resp2 = rm1.predict_proba(Xt)
    assert_array_equal(res1, res2)
    assert_allclose(resp1, resp2, rtol=1e-10)
