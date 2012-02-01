import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_true

from sklearn.utils import check_random_state
from sklearn.datasets import load_iris
from sklearn.linear_model import Perceptron

iris = load_iris()
random_state = check_random_state(12)
indices = np.arange(iris.data.shape[0])
random_state.shuffle(indices)
X = iris.data[indices]
y = iris.target[indices]
X_csr = sp.csr_matrix(X)
X_csr.sort_indices()


class MyPerceptron(object):

    def __init__(self, n_iter=1):
        self.n_iter = n_iter

    def fit(self, X, y):
        n_samples, n_features = X.shape
        self.w = np.zeros(n_features, dtype=np.float64)
        self.b = 0.0

        for t in range(self.n_iter):
            for i in range(n_samples):
                if self.predict(X[i])[0] != y[i]:
                    self.w += y[i] * X[i]
                    self.b += y[i]

    def project(self, X):
        return np.dot(X, self.w) + self.b

    def predict(self, X):
        X = np.atleast_2d(X)
        return np.sign(self.project(X))


def test_perceptron_accuracy():
    for data in (X, X_csr):
        clf = Perceptron(n_iter=30, shuffle=False, seed=0)
        clf.fit(data, y)
        score = clf.score(data, y)
        assert_true(score >= 0.7)


def test_perceptron_correctness():
    y_bin = y.copy()
    y_bin[y != 1] = -1

    clf1 = MyPerceptron(n_iter=2)
    clf1.fit(X, y_bin)

    clf2 = Perceptron(n_iter=2)
    clf2.fit(X, y_bin)

    assert_array_almost_equal(clf1.w, clf2.coef_.ravel())
