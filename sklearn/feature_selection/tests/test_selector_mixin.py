import numpy as np
import scipy.sparse as sp

from nose.tools import assert_true

from sklearn.utils.testing import assert_less

from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.svm import LinearSVC

iris = load_iris()


def test_transform_linear_model():
    for clf in (LogisticRegression(C=0.1),
                LinearSVC(C=0.01, dual=False),
                SGDClassifier(alpha=0.1, n_iter=10, shuffle=True, seed=0)):
        for func in (np.array, sp.csr_matrix):
            X = func(iris.data)
            clf.set_params(penalty="l1")
            clf.fit(X, iris.target)
            X_new = clf.transform(X)
            if isinstance(clf, SGDClassifier):
                assert_true(X_new.shape[1] <= X.shape[1])
            else:
                assert_less(X_new.shape[1], X.shape[1])
            clf.set_params(penalty="l2")
            clf.fit(X_new, iris.target)
            pred = clf.predict(X_new)
            assert_true(np.mean(pred == iris.target) >= 0.7)
