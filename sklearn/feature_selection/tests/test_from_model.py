import numpy as np
import scipy.sparse as sp

from nose.tools import assert_raises, assert_true

from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.svm import LinearSVC

iris = load_iris()


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
