import numpy as np
import scipy.sparse as sp

from nose.tools import assert_raises, assert_true

from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal

from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import PassiveAggressiveClassifier

iris = datasets.load_iris()


def test_transform_linear_model():
    for clf in (LogisticRegression(C=0.1),
                LinearSVC(C=0.01, dual=False),
                SGDClassifier(alpha=0.001, n_iter=50, shuffle=True,
                              random_state=0)):
        for thresh in (None, ".09*mean", "1e-5 * median"):
            for func in (np.array, sp.csr_matrix):
                X = func(iris.data)
                clf.set_params(penalty="l1")
                clf.fit(X, iris.target)
                X_new = clf.transform(X, thresh)
                if isinstance(clf, SGDClassifier):
                    assert_true(X_new.shape[1] <= X.shape[1])
                else:
                    assert_less(X_new.shape[1], X.shape[1])
                clf.set_params(penalty="l2")
                clf.fit(X_new, iris.target)
                pred = clf.predict(X_new)
                assert_greater(np.mean(pred == iris.target), 0.7)


def test_invalid_input():
    clf = SGDClassifier(alpha=0.1, n_iter=10, shuffle=True, random_state=None)

    clf.fit(iris.data, iris.target)
    assert_raises(ValueError, clf.transform, iris.data, "gobbledigook")
    assert_raises(ValueError, clf.transform, iris.data, ".5 * gobbledigook")


def test_validate_estimator():
    est = RandomForestClassifier()
    transformer = SelectFromModel(estimator=est)
    transformer.fit(iris.data, iris.target)
    assert_equal(transformer.estimator, est)


def test_feature_importances():
    X, y = datasets.make_classification(n_samples=1000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=0)

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    transformer = SelectFromModel(estimator=est)

    transformer.fit(X, y)
    assert_true(hasattr(transformer.estimator_, 'feature_importances_'))

    X_new = transformer.transform(X)
    assert_less(X_new.shape[1], X.shape[1])

    feature_mask = (transformer.estimator_.feature_importances_ >
                    transformer.estimator_.feature_importances_.mean())
    assert_array_almost_equal(X_new, X[:, feature_mask])

    # Check with sample weights
    sample_weight = np.ones(y.shape)
    sample_weight[y == 1] *= 100

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(X, y, sample_weight=sample_weight)
    importances = transformer.estimator_.feature_importances_
    assert_less(importances[1], X.shape[1])

    est = RandomForestClassifier(n_estimators=50, random_state=0)
    transformer = SelectFromModel(estimator=est)
    transformer.fit(X, y, sample_weight=3*sample_weight)
    importances_bis = transformer.estimator_.feature_importances_
    assert_almost_equal(importances, importances_bis)


def test_partial_fit():
    est = PassiveAggressiveClassifier()
    transformer = SelectFromModel(estimator=est)
    transformer.partial_fit(iris.data, iris.target,
                            classes=np.unique(iris.target))
    id_1 = id(transformer.estimator_)
    transformer.partial_fit(iris.data, iris.target,
                            classes=np.unique(iris.target))
    id_2 = id(transformer.estimator_)
    assert_equal(id_1, id_2)


def test_warm_start():
    est = PassiveAggressiveClassifier()
    transformer = SelectFromModel(estimator=est,
                                  warm_start=True)
    transformer.fit(iris.data, iris.target)
    id_1 = id(transformer.estimator_)
    transformer.fit(iris.data, iris.target)
    id_2 = id(transformer.estimator_)
    assert_equal(id_1, id_2)
