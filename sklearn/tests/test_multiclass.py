import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_greater
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OutputCodeClassifier
from sklearn.multiclass import _random_code_book
from sklearn.multiclass import _max_hamming_code_book

from sklearn.multiclass import fit_ovr
from sklearn.multiclass import fit_ovo
from sklearn.multiclass import fit_ecoc
from sklearn.multiclass import predict_ovr
from sklearn.multiclass import predict_ovo
from sklearn.multiclass import predict_ecoc
from sklearn.multiclass import predict_proba_ovr

from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics.pairwise import pairwise_distances

from sklearn.preprocessing import LabelBinarizer

from sklearn.svm import LinearSVC, SVC
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import (LinearRegression, Lasso, ElasticNet, Ridge,
                                  Perceptron, LogisticRegression)
from sklearn.tree import DecisionTreeClassifier
from sklearn.grid_search import GridSearchCV
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

    with ignore_warnings():
        assert_raises(ValueError, predict_ovr, [LinearSVC(), MultinomialNB()],
                      LabelBinarizer(), [])

    # Fail on multioutput data
    assert_raises(ValueError, OneVsRestClassifier(MultinomialNB()).fit,
                  np.array([[1, 0], [0, 1]]),
                  np.array([[1, 2], [3, 1]]))
    assert_raises(ValueError, OneVsRestClassifier(MultinomialNB()).fit,
                  np.array([[1, 0], [0, 1]]),
                  np.array([[1.5, 2.4], [3.1, 0.8]]))


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
                                                       return_indicator=True,
                                                       random_state=0)

        X_train, Y_train = X[:80], Y[:80]
        X_test, Y_test = X[80:], Y[80:]

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
    """Test that ovr works with classes that are always present or absent."""
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


@ignore_warnings
def test_ovr_multilabel():
    # Toy dataset where features correspond directly to labels.
    X = np.array([[0, 4, 5], [0, 5, 0], [3, 3, 3], [4, 0, 6], [6, 0, 0]])
    y = [["spam", "eggs"], ["spam"], ["ham", "eggs", "spam"],
         ["ham", "eggs"], ["ham"]]
    # y = [[1, 2], [1], [0, 1, 2], [0, 2], [0]]
    Y = np.array([[0, 1, 1],
                  [0, 1, 0],
                  [1, 1, 1],
                  [1, 0, 1],
                  [1, 0, 0]])

    classes = set("ham eggs spam".split())

    for base_clf in (MultinomialNB(), LinearSVC(random_state=0),
                     LinearRegression(), Ridge(),
                     ElasticNet(), Lasso(alpha=0.5)):
        # test input as lists of tuples
        clf = assert_warns(DeprecationWarning,
                           OneVsRestClassifier(base_clf).fit,
                           X, y)
        assert_equal(set(clf.classes_), classes)
        y_pred = clf.predict([[0, 4, 4]])[0]
        assert_equal(set(y_pred), set(["spam", "eggs"]))
        assert_true(clf.multilabel_)

        # test input as label indicator matrix
        clf = OneVsRestClassifier(base_clf).fit(X, Y)
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
                                                       return_indicator=True,
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
                                                       return_indicator=True,
                                                       random_state=0)
        X_train, Y_train = X[:80], Y[:80]
        X_test, Y_test = X[80:], Y[80:]
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
    X_test, Y_test = X[80:], Y[80:]
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
                                                   return_indicator=True,
                                                   random_state=0)
    X_train, Y_train = X[:80], Y[:80]
    X_test, Y_test = X[80:], Y[80:]
    clf = OneVsRestClassifier(svm.SVC()).fit(X_train, Y_train)
    assert_array_equal((clf.decision_function(X_test) > 0).astype(int),
                       clf.predict(X_test))


def test_ovr_single_label_decision_function():
    X, Y = datasets.make_classification(n_samples=100,
                                        n_features=20,
                                        random_state=0)
    X_train, Y_train = X[:80], Y[:80]
    X_test, Y_test = X[80:], Y[80:]
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
    ovr = OneVsRestClassifier(LinearSVC(random_state=0))
    ovr.fit(iris.data, iris.target)
    shape = ovr.coef_.shape
    assert_equal(shape[0], n_classes)
    assert_equal(shape[1], iris.data.shape[1])


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


def test_ovo_fit_predict():
    # A classifier which implements decision_function.
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)

    # A classifier which implements predict_proba.
    ovo = OneVsOneClassifier(MultinomialNB())
    ovo.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ovo.estimators_), n_classes * (n_classes - 1) / 2)


def test_ovo_gridsearch():
    ovo = OneVsOneClassifier(LinearSVC(random_state=0))
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ovo, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)


def test_ovo_ties():
    # test that ties are broken using the decision function, not defaulting to
    # the smallest label
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y = np.array([2, 0, 1, 2])
    multi_clf = OneVsOneClassifier(Perceptron())
    ovo_prediction = multi_clf.fit(X, y).predict(X)

    # recalculate votes to make sure we have a tie
    predictions = np.vstack([clf.predict(X) for clf in multi_clf.estimators_])
    scores = np.vstack([clf.decision_function(X)
                        for clf in multi_clf.estimators_])
    # classifiers are in order 0-1, 0-2, 1-2
    # aggregate votes:
    votes = np.zeros((4, 3))
    votes[np.arange(4), predictions[0]] += 1
    votes[np.arange(4), 2 * predictions[1]] += 1
    votes[np.arange(4), 1 + predictions[2]] += 1
    # for the first point, there is one vote per class
    assert_array_equal(votes[0, :], 1)
    # for the rest, there is no tie and the prediction is the argmax
    assert_array_equal(np.argmax(votes[1:], axis=1), ovo_prediction[1:])
    # for the tie, the prediction is the class with the highest score
    assert_equal(ovo_prediction[0], 0)
    # in the zero-one classifier, the score for 0 is greater than the score for
    # one.
    assert_greater(scores[0][0], scores[0][1])
    # score for one is greater than score for zero
    assert_greater(scores[2, 0] - scores[0, 0], scores[0, 0] + scores[1, 0])
    # score for one is greater than score for two
    assert_greater(scores[2, 0] - scores[0, 0], -scores[1, 0] - scores[2, 0])


def test_ovo_ties2():
    # test that ties can not only be won by the first two labels
    X = np.array([[1, 2], [2, 1], [-2, 1], [-2, -1]])
    y_ref = np.array([2, 0, 1, 2])

    # cycle through labels so that each label wins once
    for i in range(3):
        y = (y_ref + i) % 3
        multi_clf = OneVsOneClassifier(Perceptron())
        ovo_prediction = multi_clf.fit(X, y).predict(X)
        assert_equal(ovo_prediction[0], i % 3)


def test_ovo_string_y():
    "Test that the OvO doesn't screw the encoding of string labels"
    X = np.eye(4)
    y = np.array(['a', 'b', 'c', 'd'])

    svc = LinearSVC()
    ovo = OneVsOneClassifier(svc)
    ovo.fit(X, y)
    assert_array_equal(y, ovo.predict(X))

def test_code_book_functions():
    random_state = np.random
    random_state.seed(0)
    code_book = _random_code_book(3, random_state, 10000)
    proportion_of_1 = np.sum((code_book==1).astype(int)) * 1.0 / 30000
    assert_true(proportion_of_1 > 0.48 and proportion_of_1 < 0.52)
    code_book = _max_hamming_code_book(5, random_state, 15, 10)
    assert_equal(5, code_book.shape[0])
    assert_equal(15, code_book.shape[1])

def test_ecoc_exceptions():
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0))
    assert_raises(ValueError, ecoc.predict, [])
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                strategy="abc")
    assert_raises(ValueError, ecoc.fit, [], [])
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0), code_size=1.5,
                                strategy="max_hamming")
    assert_raises(ValueError, ecoc.fit, [], np.array([0, 1, 2]))
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0), code_size=0.01,
                                strategy="max_hamming")
    assert_raises(ValueError, ecoc.fit, [], np.array([0, 1, 2, 3]))

def test_ecoc_fit_predict():
    # A classifier which implements decision_function.
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                code_size=1, random_state=0)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 1)

    # A classifier which implements predict_proba.
    ecoc = OutputCodeClassifier(MultinomialNB(), code_size=1, random_state=0)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 1)

def test_ecoc_strategy():
    # For irsi dataset, code_size=1.5 will use random_code_book
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                code_size=1.5, random_state=0)
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), int(n_classes * 1.5))

    # Set the strategy to be "random"
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                code_size=1.5, random_state=0,
                                strategy="random")
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), int(n_classes * 1.5))

    # Set the strategy to be "max_hamming"
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                code_size=1.0, random_state=0,
                                strategy="max_hamming")
    ecoc.fit(iris.data, iris.target).predict(iris.data)
    assert_equal(len(ecoc.estimators_), n_classes * 1)

def test_max_hamming_code_book():
    # Test the the code could be improved using larger max_iter 
    random_state = np.random
    random_state.seed(0)
    dist0 = np.sum(pairwise_distances(_max_hamming_code_book(5, random_state,
                                                            10, 1),
                                      metric='hamming'))
    random_state = np.random
    random_state.seed(0)
    dist1 = np.sum(pairwise_distances(_max_hamming_code_book(5, random_state,
                                                            10, 2),
                                      metric='hamming')) 
    assert_true(dist0 >= dist1);

def test_ecoc_gridsearch():
    ecoc = OutputCodeClassifier(LinearSVC(random_state=0),
                                random_state=0)
    Cs = [0.1, 0.5, 0.8]
    cv = GridSearchCV(ecoc, {'estimator__C': Cs})
    cv.fit(iris.data, iris.target)
    best_C = cv.best_estimator_.estimators_[0].C
    assert_true(best_C in Cs)

@ignore_warnings
def test_deprecated():
    base_estimator = DecisionTreeClassifier(random_state=0)
    X, Y = iris.data, iris.target
    X_train, Y_train = X[:80], Y[:80]
    X_test, Y_test = X[80:], Y[80:]

    all_metas = [
        (OneVsRestClassifier, fit_ovr, predict_ovr, predict_proba_ovr),
        (OneVsOneClassifier, fit_ovo, predict_ovo, None),
        (OutputCodeClassifier, fit_ecoc, predict_ecoc, None),
    ]

    for MetaEst, fit_func, predict_func, proba_func in all_metas:
        try:
            meta_est = MetaEst(base_estimator,
                               random_state=0).fit(X_train, Y_train)

            fitted_return = fit_func(base_estimator, X_train, Y_train,
                                     random_state=0)
        except TypeError:
            meta_est = MetaEst(base_estimator).fit(X_train, Y_train)
            fitted_return = fit_func(base_estimator, X_train, Y_train)


        if len(fitted_return) == 2:
            estimators_, classes_or_lb = fitted_return
            assert_almost_equal(predict_func(estimators_, classes_or_lb, X_test),
                                meta_est.predict(X_test))

            if proba_func is not None:
                assert_almost_equal(proba_func(estimators_, X_test,
                                               is_multilabel=False),
                                    meta_est.predict_proba(X_test))

        else:
            estimators_, classes_or_lb, codebook = fitted_return
            assert_almost_equal(predict_func(estimators_, classes_or_lb,
                                             codebook, X_test),
                                meta_est.predict(X_test))


if __name__ == "__main__":
    import nose
    nose.runmodule()
