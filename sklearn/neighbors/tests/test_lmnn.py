import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn import neighbors, datasets

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

    clf = neighbors.LargeMarginNearestNeighbor(n_neighbors=1)
    clf.fit(iris.data, iris.target)
    assert_array_equal(clf.predict(iris.data), iris.target)

    clf.set_params(n_neighbors=9)
    clf.fit(iris.data, iris.target)
    assert_true(np.mean(clf.predict(iris.data) == iris.target) > 0.95)


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

    clf = neighbors.LargeMarginNearestNeighbor(n_neighbors=1)
    score_uint8 = clf.fit(X_train, Y_train).score(X_test, Y_test)
    score_float = clf.fit(X_train.astype(float), Y_train).score(
        X_test.astype(float), Y_test)
    assert_equal(score_uint8, score_float)


def test_neighbors_badargs():
    # Test bad argument values: these should all raise ValueErrors
    cls = neighbors.LargeMarginNearestNeighbor
    assert_raises(ValueError, cls.fit, L='blah')

    X = rng.random_sample((10, 2))
    # Xsparse = csr_matrix(X)
    y = np.ones(10)

    nbrs = cls()
    assert_raises(ValueError, nbrs.predict, X)

    nbrs = cls()
    assert_raises(ValueError, nbrs.fit, np.ones((0, 2)), np.ones(0))
    assert_raises(ValueError, nbrs.fit, X[:, :, None], y)
    nbrs.fit(X, y)
    assert_raises(ValueError, nbrs.predict, [[]])

    # check negative number of neighbors
    nbrs = cls(n_neighbors=-1)
    assert_raises(ValueError, nbrs.fit, X, y)


def check_object_arrays(nparray, list_check):
    for ind, ele in enumerate(nparray):
        assert_array_equal(ele, list_check[ind])


def test_same_lmnn_parallel():
    X, y = datasets.make_classification(n_samples=30, n_features=5,
                                        n_redundant=0, random_state=0)
    X_train, X_test, y_train, y_test = train_test_split(X, y)

    clf = neighbors.LargeMarginNearestNeighbor(n_neighbors=3)
    clf.fit(X_train, y_train)
    y = clf.predict(X_test)

    clf.set_params(n_jobs=3)
    clf.fit(X_train, y_train)
    y_parallel = clf.predict(X_test)

    assert_array_equal(y, y_parallel)


def test_dtype_convert():
    classifier = neighbors.KNeighborsClassifier(n_neighbors=1)
    CLASSES = 15
    X = np.eye(CLASSES)
    y = [ch for ch in 'ABCDEFGHIJKLMNOPQRSTU'[:CLASSES]]

    result = classifier.fit(X, y).predict(X)
    assert_array_equal(result, y)
