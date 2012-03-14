"""
Testing for the nearest centroid module.
"""

import numpy as np
from scipy import sparse as sp
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.linear_model import NearestCentroid
from sklearn import datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
X_csr = sp.csr_matrix(X)  # Sparse matrix
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    clf = NearestCentroid()
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    # Same test, but with a sparse matrix.
    clf = NearestCentroid()
    clf.fit(X_csr, y)
    assert_array_equal(clf.predict(T), true_result)


def test_iris():
    """Check consistency on dataset iris."""
    for metric in ('euclidean', 'cosine'):
        clf = NearestCentroid().fit(iris.data, iris.target)
        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.9, "Failed with score = " + str(score)


def test_iris_shrinkage():
    """Check consistency on dataset iris, when using shrinkage."""
    for metric in ('euclidean', 'cosine'):
        for shrink_threshold in [None, 0.1, 0.5]:
            clf = NearestCentroid(shrink_threshold=shrink_threshold)
            clf = clf.fit(iris.data, iris.target)
            score = np.mean(clf.predict(iris.data) == iris.target)
            assert score > 0.8, "Failed with score = " + str(score)


def test_iris_shrinkage_sparse():
    """Check quality on iris, when using shrinkage and sparse matrix."""
    iris_sparse = sp.csr_matrix(iris.data)
    for metric in ('euclidean', 'cosine'):
        for shrink_threshold in [None, 0.1, 0.5]:
            clf = NearestCentroid(shrink_threshold=shrink_threshold)
            clf = clf.fit(iris_sparse, iris.target)
            score = np.mean(clf.predict(iris.data) == iris.target)
            assert score > 0.8, "Failed with score = " + str(score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    for metric in ('euclidean', 'cosine'):
        clf = NearestCentroid(metric=metric).fit(iris.data, iris.target)
        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.9, "Failed with score = " + str(score)


def test_pickle():
    import pickle

    # classification
    obj = NearestCentroid()
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert score == score2, "Failed to generate same score " + \
            " after pickling (classification) "


if __name__ == "__main__":
    import nose
    nose.runmodule()
