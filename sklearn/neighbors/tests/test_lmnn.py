import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn import datasets
from sklearn.neighbors import LargeMarginNearestNeighbor as LMNN

rng = np.random.RandomState(0)
# load and shuffle iris dataset
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# load and shuffle digits
digits = datasets.load_digits()
perm = rng.permutation(digits.target.size)
digits.data = digits.data[perm]
digits.target = digits.target[perm]


def test_neighbors_iris():
    # Sanity checks on the iris dataset
    # Puts three points of each label in the plane and performs a
    # nearest neighbor query on points near the decision boundary.

    clf = LMNN(n_neighbors=1)
    clf.fit(iris.data, iris.target)
    assert_array_equal(clf.predict(iris.data), iris.target)

    clf.set_params(n_neighbors=9)
    clf.fit(iris.data, iris.target)

    assert_true(clf.score(iris.data, iris.target) > 0.95)


def test_neighbors_digits():
    # Sanity check on the digits dataset
    # the 'brute' algorithm has been observed to fail if the input
    # dtype is uint8 due to overflow in distance calculations.

    X = digits.data.astype('uint8')
    Y = digits.target
    (n_samples, n_features) = X.shape
    train_test_boundary = int(n_samples * 0.8)
    train = np.arange(0, train_test_boundary)
    test = np.arange(train_test_boundary, n_samples)
    (X_train, Y_train, X_test, Y_test) = X[train], Y[train], X[test], Y[test]

    clf = LMNN(n_neighbors=1)
    score_uint8 = clf.fit(X_train, Y_train).score(X_test, Y_test)
    score_float = clf.fit(X_train.astype(float), Y_train).score(
        X_test.astype(float), Y_test)
    assert_equal(score_uint8, score_float)


def test_params_errors():
    # Test that invalid parameters raise value error
    X = [[3, 2], [1, 6]]
    y = [1, 0]
    clf = LMNN

    assert_raises(ValueError, clf(n_neighbors=-1).fit, X, y)
    assert_raises(ValueError, clf(max_iter=-1).fit, X, y)
    assert_raises(ValueError, clf(verbose='true').fit, X, y)
    assert_raises(ValueError, clf(max_constraints=-1).fit, X, y)
    assert_raises(ValueError, clf(max_corrections=-1).fit, X, y)
    assert_raises(ValueError, clf(iprint=2).fit, X, y)
    assert_raises(ValueError, clf(tol=-0.5).fit, X, y)
    assert_raises(ValueError, clf(L='invalid').fit, X, y)
    assert_raises(ValueError, clf(n_features_out='invalid').fit, X, y)
    assert_raises(ValueError, clf(use_pca=1).fit, X, y)
    assert_raises(ValueError, clf(n_jobs=-0.5).fit, X, y)
    assert_raises(ValueError, clf(warm_start=1).fit, X, y)
    assert_raises(ValueError, clf(use_sparse=-0.5).fit, X, y)


def check_object_arrays(nparray, list_check):
    for ind, ele in enumerate(nparray):
        assert_array_equal(ele, list_check[ind])


def test_same_lmnn_parallel():
    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_redundant=0, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y)

    clf = LMNN(n_neighbors=3)
    clf.fit(X_train, y_train)
    y = clf.predict(X_test)

    clf.set_params(n_jobs=3)
    clf.fit(X_train, y_train)
    y_parallel = clf.predict(X_test)

    assert_array_equal(y, y_parallel)


def test_dtype_convert():
    classifier = LMNN(n_neighbors=1)
    CLASSES = 15
    X = np.eye(CLASSES)
    y = [ch for ch in 'ABCDEFGHIJKLMNOPQRSTU'[:CLASSES]]

    result = classifier.fit(X, y).predict(X)
    assert_array_equal(result, y)
