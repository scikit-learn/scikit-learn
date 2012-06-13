"""
General tests for all estimators in sklearn.
"""
import warnings
import numpy as np
from nose.tools import assert_raises
from numpy.testing import assert_array_equal

from sklearn.utils.testing import all_estimators
from sklearn.utils.testing import assert_greater
from sklearn.base import clone, ClassifierMixin, RegressorMixin
from sklearn.utils import shuffle
from sklearn.preprocessing import Scaler
#from sklearn.datasets import load_digits
from sklearn.datasets import load_iris, load_boston
from sklearn.metrics import zero_one_score
from sklearn.lda import LDA
from sklearn.svm.base import BaseLibSVM, BaseLibLinear

# import "special" estimators
from sklearn.grid_search import GridSearchCV
from sklearn.decomposition import SparseCoder
from sklearn.pipeline import Pipeline
from sklearn.ensemble import BaseEnsemble
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier,\
        OutputCodeClassifier
from sklearn.feature_selection import RFE, RFECV
from sklearn.naive_bayes import MultinomialNB, BernoulliNB
from sklearn.linear_model import RidgeClassifier, RidgeClassifierCV

dont_test = [Pipeline, GridSearchCV, SparseCoder]
meta_estimators = [BaseEnsemble, OneVsOneClassifier, OutputCodeClassifier,
        OneVsRestClassifier, RFE, RFECV]


def test_all_estimators():
    estimators = all_estimators()
    clf = LDA()

    for name, E in estimators:
        # some can just not be sensibly default constructed
        if E in dont_test:
            continue
        # test default-constructibility
        # get rid of deprecation warnings
        with warnings.catch_warnings(record=True) as w:
            if E in meta_estimators:
                e = E(clf)
            else:
                e = E()
            #test cloning
            clone(e)
            # test __repr__
            repr(e)
        print(w)


def test_classifiers():
    estimators = all_estimators()
    classifiers = [(name, E) for name, E in estimators if issubclass(E,
        ClassifierMixin)]
    iris = load_iris()
    X, y = iris.data, iris.target
    X, y = shuffle(X, y, random_state=7)
    #digits = load_digits()
    #X, y = digits.data, digits.target
    X = Scaler().fit_transform(X)
    for name, Clf in classifiers:
        if Clf in dont_test or Clf in meta_estimators:
            continue
        if Clf in [MultinomialNB, BernoulliNB]:
            # TODO also test these!
            continue
        clf = Clf()
        # fit
        clf.fit(X, y)
        y_pred = clf.predict(X)
        # training set performance
        assert_greater(zero_one_score(y, y_pred), 0.78)
        # raises error on malformed input for predict
        if isinstance(clf, BaseLibSVM) or isinstance(clf, BaseLibLinear):
            # TODO: libsvm decision functions, input validation
            continue
        assert_raises(ValueError, clf.predict, X.T)
        if hasattr(clf, "decision_function"):
            try:
                #raises error on malformed input for decision_function
                assert_raises(ValueError, clf.decision_function, X.T)
                #decision_function agrees with predict:
                decision = clf.decision_function(X)
                assert_array_equal(np.argmax(decision, axis=1), y_pred)
            except NotImplementedError:
                pass
        if hasattr(clf, "predict_proba"):
            try:
                assert_raises(ValueError, clf.predict_proba, X.T)
                # decision_function agrees with predict:
                y_prob = clf.predict_proba(X)
                assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
            except NotImplementedError:
                pass


def test_regressors():
    estimators = all_estimators()
    regressors = [(name, E) for name, E in estimators if issubclass(E,
        RegressorMixin)]
    boston = load_boston()
    X, y = boston.data, boston.target
    X, y = shuffle(X, y, random_state=0)
    # TODO: test with intercept
    # TODO: test with multiple responses
    X = Scaler().fit_transform(X)
    y = Scaler().fit_transform(y)
    for name, Reg in regressors:
        if Reg in dont_test or Reg in meta_estimators:
            continue
        if Reg in [RidgeClassifier, RidgeClassifierCV]:
            #TODO this is not a regressor!
            continue
        reg = Reg()
        if hasattr(reg, 'alpha'):
            reg.set_params(alpha=0.01)

        # fit
        reg.fit(X, y)
        reg.predict(X)
        assert_greater(reg.score(X, y), 0.5)
