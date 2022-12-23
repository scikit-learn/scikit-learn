import warnings

import numpy as np
import scipy.sparse as sp
import pytest
from numpy.testing import assert_allclose

from re import escape

from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_almost_equal
from sklearn.utils._mocking import CheckingClassifier
from sklearn.multiclass import _ConstantPredictor
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OutputCodeClassifier
from sklearn.utils.multiclass import check_classification_targets, type_of_target
from sklearn.utils import (
    check_array,
    shuffle,
)

from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

from sklearn.svm import LinearSVC, SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import (
    LinearRegression,
    Lasso,
    ElasticNet,
    Ridge,
    Perceptron,
    LogisticRegression,
    SGDClassifier,
)
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.impute import SimpleImputer
from sklearn import svm
from sklearn.exceptions import NotFittedError
from sklearn import datasets
from sklearn.datasets import load_breast_cancer, load_digits, make_classification

msg = "The default value for `force_alpha` will change"
pytestmark = pytest.mark.filterwarnings(f"ignore:{msg}:FutureWarning")

iris = datasets.load_iris()
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
n_classes = 3


def test_ovr_exceptions():
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))

    # test predicting without fitting
    with pytest.raises(NotFittedError):
        ovr.predict([])

    # Fail on multioutput data
    msg = "Multioutput target data is not supported with label binarization"
    with pytest.raises(ValueError, match=msg):
        X = np.array([[1, 0], [0, 1]])
        y = np.array([[1, 2], [3, 1]])
        OneVsRestClassifier(MultinomialNB()).fit(X, y)

    with pytest.raises(ValueError, match=msg):
        X = np.array([[1, 0], [0, 1]])
        y = np.array([[1.5, 2.4], [3.1, 0.8]])
        OneVsRestClassifier(MultinomialNB()).fit(X, y)


def test_check_classification_targets():
    # Test that check_classification_target return correct type. #5782
    y = np.array([0.0, 1.1, 2.0, 3.0])
    msg = type_of_target(y)
    with pytest.raises(ValueError, match=msg):
        check_classification_targets(y)


def test_ovr_fit_predict():
    # A classifier which implements decision_function.
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert len(ovr.estimators_) == n_classes

    clf = LinearSVC(random_state=0)
    pred2 = clf.fit(iris.data, iris.target).predict(iris.data)
    assert np.mean(iris.target == pred) == np.mean(iris.target == pred2)

    # A classifier which implements predict_proba.
    ovr = OneVsRestClassifier(MultinomialNB())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert np.mean(iris.target == pred) > 0.65


def test_ovr_partial_fit():
    # Test if partial_fit is working as intended
    X, y = shuffle(iris.data, iris.target, random_state=0)
    ovr = OneVsRestClassifier(MultinomialNB())
    ovr.partial_fit(X[:100], y[:100], np.unique(y))
    ovr.partial_fit(X[100:], y[100:])
    pred = ovr.predict(X)
    ovr2 = OneVsRestClassifier(MultinomialNB())
    pred2 = ovr2.fit(X, y).predict(X)

    assert_almost_equal(pred, pred2)
    assert len(ovr.estimators_) == len(np.unique(y))
    assert np.mean(y == pred) > 0.65

    # Test when mini batches doesn't have all classes
    # with SGDClassifier
    X = np.abs(np.random.randn(14, 2))
    y = [1, 1, 1, 1, 2, 3, 3, 0, 0, 2, 3, 1, 2, 3]

    ovr = OneVsRestClassifier(
        SGDClassifier(max_iter=1, tol=None, shuffle=False, random_state=0)
    )
    ovr.partial_fit(X[:7], y[:7], np.unique(y))
    ovr.partial_fit(X[7:], y[7:])
    pred = ovr.predict(X)
    ovr1 = OneVsRestClassifier(
        SGDClassifier(max_iter=1, tol=None, shuffle=False, random_state=0)
    )
    pred1 = ovr1.fit(X, y).predict(X)
    assert np.mean(pred == y) == np.mean(pred1 == y)

    # test partial_fit only exists if estimator has it:
    ovr = OneVsRestClassifier(SVC())
    assert not hasattr(ovr, "partial_fit")


def test_ovr_partial_fit_exceptions():
    ovr = OneVsRestClassifier(MultinomialNB())
    X = np.abs(np.random.randn(14, 2))
    y = [1, 1, 1, 1, 2, 3, 3, 0, 0, 2, 3, 1, 2, 3]
    ovr.partial_fit(X[:7], y[:7], np.unique(y))
    # If a new class that was not in the first call of partial fit is seen
    # it should raise ValueError
    y1 = [5] + y[7:-1]
    msg = r"Mini-batch contains \[.+\] while classes must be subset of \[.+\]"
    with pytest.raises(ValueError, match=msg):
        ovr.partial_fit(X=X[7:], y=y1)


def test_ovr_ovo_regressor():
    # test that ovr and ovo work on regressors which don't have a decision_
    # function
    ovr = OneVsRestClassifier(DecisionTreeRegressor())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert len(ovr.estimators_) == n_classes
    assert_array_equal(np.unique(pred), [0, 1, 2])
    # we are doing something sensible
    assert np.mean(pred == iris.target) > 0.9

    ovr = OneVsOneClassifier(DecisionTreeRegressor())
    pred = ovr.fit(iris.data, iris.target).predict(iris.data)
    assert len(ovr.estimators_) == n_classes * (n_classes - 1) / 2
    assert_array_equal(np.unique(pred), [0, 1, 2])
    # we are doing something sensible
    assert np.mean(pred == iris.target) > 0.9


def test_ovr_fit_predict_sparse():
    for sparse in [
        sp.csr_matrix,
        sp.csc_matrix,
        sp.coo_matrix,
        sp.dok_matrix,
        sp.lil_matrix,
    ]:
        base_clf = MultinomialNB(alpha=1)

        X, Y = datasets.make_multilabel_classification(
            n_samples=100,
            n_features=20,
            n_classes=5,
            n_labels=3,
            length=50,
            allow_unlabeled=True,
            random_state=0,
        )

        X_train, Y_train = X[:80], Y[:80]
        X_test = X[80:]

        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)
        Y_pred = clf.predict(X_test)

        clf_sprs = OneVsRestClassifier(base_clf).fit(X_train, sparse(Y_train))
        Y_pred_sprs = clf_sprs.predict(X_test)

        assert clf.multilabel_
        assert sp.issparse(Y_pred_sprs)
        assert_array_equal(Y_pred_sprs.toarray(), Y_pred)

        # Test predict_proba
        Y_proba = clf_sprs.predict_proba(X_test)

        # predict assigns a label if the probability that the
        # sample has the label is greater than 0.5.
        pred = Y_proba > 0.5
        assert_array_equal(pred, Y_pred_sprs.toarray())

        # Test decision_function
        clf = svm.SVC()
        clf_sprs = OneVsRestClassifier(clf).fit(X_train, sparse(Y_train))
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
    msg = r"Label .+ is present in all training examples"
    with pytest.warns(UserWarning, match=msg):
        ovr.fit(X, y)
    y_pred = ovr.predict(X)
    assert_array_equal(np.array(y_pred), np.array(y))
    y_pred = ovr.decision_function(X)
    assert np.unique(y_pred[:, -2:]) == 1
    y_pred = ovr.predict_proba(X)
    assert_array_equal(y_pred[:, -1], np.ones(X.shape[0]))

    # y has a constantly absent label
    y = np.zeros((10, 2))
    y[5:, 0] = 1  # variable label
    ovr = OneVsRestClassifier(LogisticRegression())

    msg = r"Label not 1 is present in all training examples"
    with pytest.warns(UserWarning, match=msg):
        ovr.fit(X, y)
    y_pred = ovr.predict_proba(X)
    assert_array_equal(y_pred[:, -1], np.zeros(X.shape[0]))


def test_ovr_multiclass():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 0, 5], [0, 5, 0], [3, 0, 0], [0, 0, 6], [6, 0, 0]])
    y = ["eggs", "spam", "ham", "eggs", "ham"]
    Y = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 1], [1, 0, 0]])

    classes = set("ham eggs spam".split())

    for base_clf in (
        MultinomialNB(),
        LinearSVC(random_state=0),
        LinearRegression(),
        Ridge(),
        ElasticNet(),
    ):
        clf = OneVsRestClassifier(base_clf).fit(X, y)
        assert set(clf.classes_) == classes
        y_pred = clf.predict(np.array([[0, 0, 4]]))[0]
        assert_array_equal(y_pred, ["eggs"])

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
        assert set(clf.classes_) == classes
        y_pred = clf.predict(np.array([[0, 0, 4]]))[0]
        assert_array_equal(y_pred, ["eggs"])
        if hasattr(base_clf, "decision_function"):
            dec = clf.decision_function(X)
            assert dec.shape == (5,)

        if test_predict_proba:
            X_test = np.array([[0, 0, 4]])
            probabilities = clf.predict_proba(X_test)
            assert 2 == len(probabilities[0])
            assert clf.classes_[np.argmax(probabilities, axis=1)] == clf.predict(X_test)

        # test input as label indicator matrix
        clf = OneVsRestClassifier(base_clf).fit(X, Y)
        y_pred = clf.predict([[3, 0, 0]])[0]
        assert y_pred == 1

    for base_clf in (
        LinearSVC(random_state=0),
        LinearRegression(),
        Ridge(),
        ElasticNet(),
    ):
        conduct_test(base_clf)

    for base_clf in (MultinomialNB(), SVC(probability=True), LogisticRegression()):
        conduct_test(base_clf, test_predict_proba=True)


def test_ovr_multilabel():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 4, 5], [0, 5, 0], [3, 3, 3], [4, 0, 6], [6, 0, 0]])
    y = np.array([[0, 1, 1], [0, 1, 0], [1, 1, 1], [1, 0, 1], [1, 0, 0]])

    for base_clf in (
        MultinomialNB(),
        LinearSVC(random_state=0),
        LinearRegression(),
        Ridge(),
        ElasticNet(),
        Lasso(alpha=0.5),
    ):
        clf = OneVsRestClassifier(base_clf).fit(X, y)
        y_pred = clf.predict([[0, 4, 4]])[0]
        assert_array_equal(y_pred, [0, 1, 1])
        assert clf.multilabel_


def test_ovr_fit_predict_svc():
    ovr = OneVsRestClassifier(svm.SVC())
    ovr.fit(iris.data, iris.target)
    assert len(ovr.estimators_) == 3
    assert ovr.score(iris.data, iris.target) > 0.9


def test_ovr_multilabel_dataset():
    base_clf = MultinomialNB(alpha=1)
    for au, prec, recall in zip((True, False), (0.51, 0.66), (0.51, 0.80)):
        X, Y = datasets.make_multilabel_classification(
            n_samples=100,
            n_features=20,
            n_classes=5,
            n_labels=2,
            length=50,
            allow_unlabeled=au,
            random_state=0,
        )
        X_train, Y_train = X[:80], Y[:80]
        X_test, Y_test = X[80:], Y[80:]
        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)
        Y_pred = clf.predict(X_test)

        assert clf.multilabel_
        assert_almost_equal(
            precision_score(Y_test, Y_pred, average="micro"), prec, decimal=2
        )
        assert_almost_equal(
            recall_score(Y_test, Y_pred, average="micro"), recall, decimal=2
        )


def test_ovr_multilabel_predict_proba():
    base_clf = MultinomialNB(alpha=1)
    for au in (False, True):
        X, Y = datasets.make_multilabel_classification(
            n_samples=100,
            n_features=20,
            n_classes=5,
            n_labels=3,
            length=50,
            allow_unlabeled=au,
            random_state=0,
        )
        X_train, Y_train = X[:80], Y[:80]
        X_test = X[80:]
        clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)

        # Decision function only estimator.
        decision_only = OneVsRestClassifier(svm.SVR()).fit(X_train, Y_train)
        assert not hasattr(decision_only, "predict_proba")

        # Estimator with predict_proba disabled, depending on parameters.
        decision_only = OneVsRestClassifier(svm.SVC(probability=False))
        assert not hasattr(decision_only, "predict_proba")
        decision_only.fit(X_train, Y_train)
        assert not hasattr(decision_only, "predict_proba")
        assert hasattr(decision_only, "decision_function")

        # Estimator which can get predict_proba enabled after fitting
        gs = GridSearchCV(
            svm.SVC(probability=False), param_grid={"probability": [True]}
        )
        proba_after_fit = OneVsRestClassifier(gs)
        assert not hasattr(proba_after_fit, "predict_proba")
        proba_after_fit.fit(X_train, Y_train)
        assert hasattr(proba_after_fit, "predict_proba")

        Y_pred = clf.predict(X_test)
        Y_proba = clf.predict_proba(X_test)

        # predict assigns a label if the probability that the
        # sample has the label is greater than 0.5.
        pred = Y_proba > 0.5
        assert_array_equal(pred, Y_pred)


def test_ovr_single_label_predict_proba():
    base_clf = MultinomialNB(alpha=1)
    X, Y = iris.data, iris.target
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(base_clf).fit(X_train, Y_train)

    # Decision function only estimator.
    decision_only = OneVsRestClassifier(svm.SVR()).fit(X_train, Y_train)
    assert not hasattr(decision_only, "predict_proba")

    Y_pred = clf.predict(X_test)
    Y_proba = clf.predict_proba(X_test)

    assert_almost_equal(Y_proba.sum(axis=1), 1.0)
    # predict assigns a label if the probability that the
    # sample has the label with the greatest predictive probability.
    pred = Y_proba.argmax(axis=1)
    assert not (pred - Y_pred).any()


def test_ovr_multilabel_decision_function():
    X, Y = datasets.make_multilabel_classification(
        n_samples=100,
        n_features=20,
        n_classes=5,
        n_labels=3,
        length=50,
        allow_unlabeled=True,
        random_state=0,
    )
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(svm.SVC()).fit(X_train, Y_train)
    assert_array_equal(
        (clf.decision_function(X_test) > 0).astype(int), clf.predict(X_test)
    )


def test_ovr_single_label_decision_function():
    X, Y = datasets.make_classification(n_samples=100, n_features=20, random_state=0)
    X_train, Y_train = X[:80], Y[:80]
    X_test = X[80:]
    clf = OneVsRestClassifier(svm.SVC()).fit(X_train, Y_train)
    assert_array_equal(clf.decision_function(X_test).ravel() > 0, clf.predict(X_test))


def test_ovr_gridsearch():
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovr, {"estimator__C": Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert best_C in Cs


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


def test_ovo_exceptions():
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    with pytest.raises(NotFittedError):
        ovo.predict([])


def test_ovo_fit_on_list():
    # Test that OneVsOne fitting works with a list of targets and yields the
    # same output as predict from an array
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    prediction_from_array = ovo.fit(iris.data, iris.target).predict(iris.data)
    iris_data_list = [list(a) for a in iris.data]
    prediction_from_list = ovo.fit(iris_data_list, list(iris.target)).predict(
        iris_data_list
    )
    assert_array_equal(prediction_from_array, prediction_from_list)


def test_ovo_fit_predict():
    # A classifier which implements decision_function.
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert len(ovo.estimators_) == n_classes * (n_classes - 1) / 2

    # A classifier which implements predict_proba.
    ovo = OneVsOneClassifier(MultinomialNB())
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert len(ovo.estimators_) == n_classes * (n_classes - 1) / 2


def test_ovo_partial_fit_predict():
    temp = datasets.load_iris()
    X, y = temp.data, temp.target
    ovo1 = OneVsOneClassifier(MultinomialNB())
    ovo1.partial_fit(X[:100], y[:100], np.unique(y))
    ovo1.partial_fit(X[100:], y[100:])
    pred1 = ovo1.predict(X)

    ovo2 = OneVsOneClassifier(MultinomialNB())
    ovo2.fit(X, y)
    pred2 = ovo2.predict(X)
    assert len(ovo1.estimators_) == n_classes * (n_classes - 1) / 2
    assert np.mean(y == pred1) > 0.65
    assert_almost_equal(pred1, pred2)

    # Test when mini-batches have binary target classes
    ovo1 = OneVsOneClassifier(MultinomialNB())
    ovo1.partial_fit(X[:60], y[:60], np.unique(y))
    ovo1.partial_fit(X[60:], y[60:])
    pred1 = ovo1.predict(X)
    ovo2 = OneVsOneClassifier(MultinomialNB())
    pred2 = ovo2.fit(X, y).predict(X)

    assert_almost_equal(pred1, pred2)
    assert len(ovo1.estimators_) == len(np.unique(y))
    assert np.mean(y == pred1) > 0.65

    ovo = OneVsOneClassifier(MultinomialNB())
    X = np.random.rand(14, 2)
    y = [1, 1, 2, 3, 3, 0, 0, 4, 4, 4, 4, 4, 2, 2]
    ovo.partial_fit(X[:7], y[:7], [0, 1, 2, 3, 4])
    ovo.partial_fit(X[7:], y[7:])
    pred = ovo.predict(X)
    ovo2 = OneVsOneClassifier(MultinomialNB())
    pred2 = ovo2.fit(X, y).predict(X)
    assert_almost_equal(pred, pred2)

    # raises error when mini-batch does not have classes from all_classes
    ovo = OneVsOneClassifier(MultinomialNB())
    error_y = [0, 1, 2, 3, 4, 5, 2]
    message_re = escape(
        "Mini-batch contains {0} while it must be subset of {1}".format(
            np.unique(error_y), np.unique(y)
        )
    )
    with pytest.raises(ValueError, match=message_re):
        ovo.partial_fit(X[:7], error_y, np.unique(y))

    # test partial_fit only exists if estimator has it:
    ovr = OneVsOneClassifier(SVC())
    assert not hasattr(ovr, "partial_fit")


def test_ovo_decision_function():
    n_samples = iris.data.shape[0]

    ovo_clf = OneVsOneClassifier(LinearSVC(random_state=0))
    # first binary
    ovo_clf.fit(iris.data, iris.target == 0)
    decisions = ovo_clf.decision_function(iris.data)
    assert decisions.shape == (n_samples,)

    # then multi-class
    ovo_clf.fit(iris.data, iris.target)
    decisions = ovo_clf.decision_function(iris.data)

    assert decisions.shape == (n_samples, n_classes)
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
        assert set(votes[:, class_idx]).issubset(set([0.0, 1.0, 2.0]))

        # The OVO decision function on the other hand is able to resolve
        # most of the ties on this data as it combines both the vote counts
        # and the aggregated confidence levels of the binary classifiers
        # to compute the aggregate decision function. The iris dataset
        # has 150 samples with a couple of duplicates. The OvO decisions
        # can resolve most of the ties:
        assert len(np.unique(decisions[:, class_idx])) > 146


def test_ovo_gridsearch():
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovo, {"estimator__C": Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert best_C in Cs


def test_ovo_ties():
    # Test that ties are broken using the decision function,
    # not defaulting to the smallest label
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y = np.array([2, 0, 1, 2])
    multi_clf = OneVsOneClassifier(Perceptron(shuffle=False, max_iter=4, tol=None))
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
    assert ovo_prediction[0] == normalized_confidences[0].argmax()


def test_ovo_ties2():
    # test that ties can not only be won by the first two labels
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y_ref = np.array([2, 0, 1, 2])

    # cycle through labels so that each label wins once
    for i in range(3):
        y = (y_ref + i) % 3
        multi_clf = OneVsOneClassifier(Perceptron(shuffle=False, max_iter=4, tol=None))
        ovo_prediction = multi_clf.fit(X, y).predict(X)
        assert ovo_prediction[0] == i % 3


def test_ovo_string_y():
    # Test that the OvO doesn't mess up the encoding of string labels
    X = np.eye(4)
    y = np.array(["a", "b", "c", "d"])

    ovo = OneVsOneClassifier(LinearSVC())
    ovo.fit(X, y)
    assert_array_equal(y, ovo.predict(X))


def test_ovo_one_class():
    # Test error for OvO with one class
    X = np.eye(4)
    y = np.array(["a"] * 4)

    ovo = OneVsOneClassifier(LinearSVC())
    msg = "when only one class"
    with pytest.raises(ValueError, match=msg):
        ovo.fit(X, y)


def test_ovo_float_y():
    # Test that the OvO errors on float targets
    X = iris.data
    y = iris.data[:, 0]

    ovo = OneVsOneClassifier(LinearSVC())
    msg = "Unknown label type"
    with pytest.raises(ValueError, match=msg):
        ovo.fit(X, y)


@pytest.mark.parametrize("decoding", ["cityblock", "hamming", "loss"])
def test_ecoc_gridsearch(decoding):
    ecoc = OutputCodeClassifier(
        DecisionTreeClassifier(random_state=0), decoding=decoding, random_state=0
    )
    max_depth = [2, 3, 5]
    cv = GridSearchCV(ecoc, {"estimator__max_depth": max_depth})
    cv.fit(iris.data, iris.target)
    best_max_depth = cv.best_estimator_.estimators_[0].max_depth
    assert best_max_depth in max_depth
    y_pred = cv.predict(iris.data)
    assert y_pred.shape == iris.target.shape


def test_ecoc_delegate_sparse_base_estimator():
    # Non-regression test for
    # https://github.com/scikit-learn/scikit-learn/issues/17218
    X, y = iris.data, iris.target
    X_sp = sp.csc_matrix(X)

    # create an estimator that does not support sparse input
    base_estimator = CheckingClassifier(
        check_X=check_array,
        check_X_params={"ensure_2d": True, "accept_sparse": False},
    )
    ecoc = OutputCodeClassifier(base_estimator, random_state=0)

    with pytest.raises(TypeError, match="A sparse matrix was passed"):
        ecoc.fit(X_sp, y)

    ecoc.fit(X, y)
    with pytest.raises(TypeError, match="A sparse matrix was passed"):
        ecoc.predict(X_sp)

    # smoke test to check when sparse input should be supported
    ecoc = OutputCodeClassifier(DecisionTreeClassifier(random_state=0))
    ecoc.fit(X_sp, y).predict(X_sp)
    assert len(ecoc.estimators_) == 4


def test_ecoc_codebook():
    """Check that the generation of the code works as expected."""
    n_classes = 4
    X, y = make_classification(
        n_samples=50, n_classes=n_classes, n_clusters_per_class=1, random_state=0
    )

    # A single column is enough to generate non-redundant codes and no warning should
    # be raised.
    code_size = 0.5
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ecoc = OutputCodeClassifier(
            DecisionTreeClassifier(random_state=0), code_size=code_size, random_state=0
        ).fit(X, y)
    assert ecoc.code_book_.shape == (n_classes, int(n_classes * code_size))
    unique_codes = np.unique(ecoc.code_book_, axis=0)
    assert unique_codes.shape == ecoc.code_book_.shape

    # check the common case where we will have more columns than needed to encode
    # the classes.
    code_size = 2
    ecoc = OutputCodeClassifier(
        DecisionTreeClassifier(random_state=0), code_size=code_size, random_state=0
    ).fit(X, y)
    assert ecoc.code_book_.shape == (n_classes, int(n_classes * code_size))
    unique_codes = np.unique(ecoc.code_book_, axis=0)
    assert unique_codes.shape == ecoc.code_book_.shape

    # error when the code size is too small
    code_size = 0.1
    err_msg = "The code book size is not big enough to encode all classes."
    with pytest.raises(ValueError, match=err_msg):
        OutputCodeClassifier(
            DecisionTreeClassifier(random_state=0), code_size=code_size, random_state=0
        ).fit(X, y)


def test_ecoc_fit():
    """General check for the fit method."""
    n_classes = 4
    X, y = make_classification(
        n_samples=50, n_classes=n_classes, n_clusters_per_class=1, random_state=0
    )
    # use non-trivial classes to make sure that the inner encoding will be correct
    classes = np.array([3, 10, 4, 1])
    y = classes[y]
    assert_array_equal(np.unique(y), np.sort(classes))

    ecoc = OutputCodeClassifier(
        DecisionTreeClassifier(random_state=0), code_size=2.0, random_state=0
    )

    # classes_ is a property and we check it raises an `AttributeError` when the
    # estimator is not fitted
    err_msg = "OutputCodeClassifier object has no attribute 'classes_'"
    with pytest.raises(AttributeError, match=err_msg):
        ecoc.classes_

    ecoc.fit(X, y)
    assert_array_equal(ecoc.classes_, np.sort(classes))

    for estimator in ecoc.estimators_:
        # all estimators but _ConstantPredictor should be fit on a binary problem
        if isinstance(estimator, _ConstantPredictor):
            continue
        assert_array_equal(estimator.classes_, [0, 1])
    assert len(ecoc.estimators_) == ecoc.code_book_.shape[1]


def test_ecoc_requested_prediction_method():
    """Check that we raise an appropriate error message when the estimator does not
    implement the requested prediction method.
    """
    X, y = make_classification(
        n_samples=50, n_classes=4, n_clusters_per_class=1, random_state=0
    )

    err_msg = (
        "The estimator does not have a `predict_proba` method. Thus, the 'cityblock' "
        "decoding strategy is not supported."
    )
    with pytest.raises(ValueError, match=err_msg):
        OutputCodeClassifier(
            LinearSVC(), code_size=2.0, random_state=0, decoding="cityblock"
        ).fit(X, y)

    err_msg = (
        "The estimator does not have either a `predict_proba` or `decision_function` "
        "method. Thus, the 'loss' decoding strategy is not supported."
    )
    with pytest.raises(ValueError, match=err_msg):
        OutputCodeClassifier(
            LinearRegression(), code_size=2.0, decoding="loss", random_state=0
        ).fit(X, y)


@pytest.mark.parametrize(
    "params",
    [
        {"decoding": "cityblock"},
        {"decoding": "hamming"},
        {"decoding": "loss", "loss": "exponential"},
    ],
    ids=["cityblock", "hamming", "loss"],
)
def test_ecoc_outperform_multiclass(params, global_random_seed):
    """Check that the ECOC strategy significantly outperforms a single tree on
    the digits dataset. This would be inline with the results reported in [1]_.

    [1] Dietterich, T. G., & Bakiri, G. (1994).
    Solving multiclass learning problems via error-correcting output codes.
    Journal of artificial intelligence research, 2, 263-286.
    """
    X, y = load_digits(return_X_y=True)

    single_tree = DecisionTreeClassifier(random_state=global_random_seed)
    ecoc = OutputCodeClassifier(single_tree, **params, random_state=0)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=global_random_seed)
    cv_results_single_tree = cross_val_score(single_tree, X, y, cv=cv)
    cv_results_ecoc = cross_val_score(ecoc, X, y, cv=cv)

    assert cv_results_ecoc.mean() > cv_results_single_tree.mean()


# We limit the number of iterations to speed up the test
@pytest.mark.filterwarnings("ignore: Liblinear failed to converge")
@pytest.mark.parametrize(
    "params",
    [
        {"decoding": "loss", "loss": "linear"},
        {"decoding": "loss", "loss": "exponential"},
        {"decoding": "loss", "loss": "hinge"},
        {"decoding": "loss", "loss": "logistic"},
        {"decoding": "loss", "loss": "square"},
    ],
    ids=["linear", "exponential", "hinge", "logistic", "square"],
)
def test_ecoc_loss_outperform_hamming(params, global_random_seed):
    """Check that the ECOC loss strategy significantly outperforms the Hamming approach.
    Even if the loss is not chosing according to the estimator, it will anyway be better
    than the Hamming approach. It is inline with the results reported in [1]_.

    [1] Allwein, E. L., Schapire, R. E., & Singer, Y. (2000).
    Reducing multiclass to binary: A unifying approach for margin classifiers.
    Journal of machine learning research, 1(Dec), 113-141.
    """
    X, y = load_digits(return_X_y=True)

    svc = LinearSVC(max_iter=10, random_state=global_random_seed)
    ecoc_hamming = OutputCodeClassifier(
        svc, decoding="hamming", random_state=global_random_seed
    )
    ecoc_loss = OutputCodeClassifier(svc, **params, random_state=global_random_seed)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=global_random_seed)
    cv_results_hamming = cross_val_score(ecoc_hamming, X, y, cv=cv)
    cv_results_loss = cross_val_score(ecoc_loss, X, y, cv=cv)

    assert cv_results_loss.mean() > cv_results_hamming.mean()


def test_ecoc_large_code_size():
    """Check that we can request a large number of code size.

    A large number of columns should not be an issue with a low number of classes.
    The matrix will only be sparse. However, if the internal implementation does not
    rely on generator, generating the radom codebook could result in a memory error
    because all the potential integer between 0 and 2 ** `n_requested_columns` will
    be stored in memory. A generator will only store the current integer.
    """
    X, y = make_classification(
        n_samples=100,
        n_features=30,
        n_classes=10,
        n_informative=20,
        n_clusters_per_class=1,
        random_state=0,
    )

    # code_size=3 will results in 30 columns in the codebook and code will be sampled
    # from 0 to 2 ** 30 - 1
    ecoc = OutputCodeClassifier(
        DecisionTreeClassifier(max_depth=2, random_state=0),
        code_size=3,
        random_state=0,
    ).fit(X, y)
    assert ecoc.code_book_.shape == (10, 30)
    y_pred = ecoc.predict(X)
    assert y_pred.shape == y.shape


def test_pairwise_indices():
    clf_precomputed = svm.SVC(kernel="precomputed")
    X, y = iris.data, iris.target

    ovr_false = OneVsOneClassifier(clf_precomputed)
    linear_kernel = np.dot(X, X.T)
    ovr_false.fit(linear_kernel, y)

    n_estimators = len(ovr_false.estimators_)
    precomputed_indices = ovr_false.pairwise_indices_

    for idx in precomputed_indices:
        assert (
            idx.shape[0] * n_estimators / (n_estimators - 1) == linear_kernel.shape[0]
        )


def test_pairwise_n_features_in():
    """Check the n_features_in_ attributes of the meta and base estimators

    When the training data is a regular design matrix, everything is intuitive.
    However, when the training data is a precomputed kernel matrix, the
    multiclass strategy can resample the kernel matrix of the underlying base
    estimator both row-wise and column-wise and this has a non-trivial impact
    on the expected value for the n_features_in_ of both the meta and the base
    estimators.
    """
    X, y = iris.data, iris.target

    # Remove the last sample to make the classes not exactly balanced and make
    # the test more interesting.
    assert y[-1] == 0
    X = X[:-1]
    y = y[:-1]

    # Fitting directly on the design matrix:
    assert X.shape == (149, 4)

    clf_notprecomputed = svm.SVC(kernel="linear").fit(X, y)
    assert clf_notprecomputed.n_features_in_ == 4

    ovr_notprecomputed = OneVsRestClassifier(clf_notprecomputed).fit(X, y)
    assert ovr_notprecomputed.n_features_in_ == 4
    for est in ovr_notprecomputed.estimators_:
        assert est.n_features_in_ == 4

    ovo_notprecomputed = OneVsOneClassifier(clf_notprecomputed).fit(X, y)
    assert ovo_notprecomputed.n_features_in_ == 4
    assert ovo_notprecomputed.n_classes_ == 3
    assert len(ovo_notprecomputed.estimators_) == 3
    for est in ovo_notprecomputed.estimators_:
        assert est.n_features_in_ == 4

    # When working with precomputed kernels we have one "feature" per training
    # sample:
    K = X @ X.T
    assert K.shape == (149, 149)

    clf_precomputed = svm.SVC(kernel="precomputed").fit(K, y)
    assert clf_precomputed.n_features_in_ == 149

    ovr_precomputed = OneVsRestClassifier(clf_precomputed).fit(K, y)
    assert ovr_precomputed.n_features_in_ == 149
    assert ovr_precomputed.n_classes_ == 3
    assert len(ovr_precomputed.estimators_) == 3
    for est in ovr_precomputed.estimators_:
        assert est.n_features_in_ == 149

    # This becomes really interesting with OvO and precomputed kernel together:
    # internally, OvO will drop the samples of the classes not part of the pair
    # of classes under consideration for a given binary classifier. Since we
    # use a precomputed kernel, it will also drop the matching columns of the
    # kernel matrix, and therefore we have fewer "features" as result.
    #
    # Since class 0 has 49 samples, and class 1 and 2 have 50 samples each, a
    # single OvO binary classifier works with a sub-kernel matrix of shape
    # either (99, 99) or (100, 100).
    ovo_precomputed = OneVsOneClassifier(clf_precomputed).fit(K, y)
    assert ovo_precomputed.n_features_in_ == 149
    assert ovr_precomputed.n_classes_ == 3
    assert len(ovr_precomputed.estimators_) == 3
    assert ovo_precomputed.estimators_[0].n_features_in_ == 99  # class 0 vs class 1
    assert ovo_precomputed.estimators_[1].n_features_in_ == 99  # class 0 vs class 2
    assert ovo_precomputed.estimators_[2].n_features_in_ == 100  # class 1 vs class 2


@pytest.mark.parametrize(
    "MultiClassClassifier", [OneVsRestClassifier, OneVsOneClassifier]
)
def test_pairwise_tag(MultiClassClassifier):
    clf_precomputed = svm.SVC(kernel="precomputed")
    clf_notprecomputed = svm.SVC()

    ovr_false = MultiClassClassifier(clf_notprecomputed)
    assert not ovr_false._get_tags()["pairwise"]

    ovr_true = MultiClassClassifier(clf_precomputed)
    assert ovr_true._get_tags()["pairwise"]


@pytest.mark.parametrize(
    "MultiClassClassifier", [OneVsRestClassifier, OneVsOneClassifier]
)
def test_pairwise_cross_val_score(MultiClassClassifier):
    clf_precomputed = svm.SVC(kernel="precomputed")
    clf_notprecomputed = svm.SVC(kernel="linear")

    X, y = iris.data, iris.target

    multiclass_clf_notprecomputed = MultiClassClassifier(clf_notprecomputed)
    multiclass_clf_precomputed = MultiClassClassifier(clf_precomputed)

    linear_kernel = np.dot(X, X.T)
    score_not_precomputed = cross_val_score(
        multiclass_clf_notprecomputed, X, y, error_score="raise"
    )
    score_precomputed = cross_val_score(
        multiclass_clf_precomputed, linear_kernel, y, error_score="raise"
    )
    assert_array_equal(score_precomputed, score_not_precomputed)


@pytest.mark.parametrize(
    "MultiClassClassifier", [OneVsRestClassifier, OneVsOneClassifier]
)
# FIXME: we should move this test in `estimator_checks` once we are able
# to construct meta-estimator instances
def test_support_missing_values(MultiClassClassifier):
    # smoke test to check that pipeline OvR and OvO classifiers are letting
    # the validation of missing values to
    # the underlying pipeline or classifiers
    rng = np.random.RandomState(42)
    X, y = iris.data, iris.target
    X = np.copy(X)  # Copy to avoid that the original data is modified
    mask = rng.choice([1, 0], X.shape, p=[0.1, 0.9]).astype(bool)
    X[mask] = np.nan
    lr = make_pipeline(SimpleImputer(), LogisticRegression(random_state=rng))

    MultiClassClassifier(lr).fit(X, y).score(X, y)


@pytest.mark.parametrize("make_y", [np.ones, np.zeros])
def test_constant_int_target(make_y):
    """Check that constant y target does not raise.

    Non-regression test for #21869
    """
    X = np.ones((10, 2))
    y = make_y((10, 1), dtype=np.int32)
    ovr = OneVsRestClassifier(LogisticRegression())

    ovr.fit(X, y)
    y_pred = ovr.predict_proba(X)
    expected = np.zeros((X.shape[0], 2))
    expected[:, 0] = 1
    assert_allclose(y_pred, expected)


def test_ovo_consistent_binary_classification():
    """Check that ovo is consistent with binary classifier.

    Non-regression test for #13617.
    """
    X, y = load_breast_cancer(return_X_y=True)

    clf = KNeighborsClassifier(n_neighbors=8, weights="distance")
    ovo = OneVsOneClassifier(clf)

    clf.fit(X, y)
    ovo.fit(X, y)

    assert_array_equal(clf.predict(X), ovo.predict(X))
