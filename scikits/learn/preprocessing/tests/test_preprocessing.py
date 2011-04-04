import numpy as np
import numpy.linalg as la
import scipy.sparse as sp

from numpy.testing import assert_array_almost_equal, assert_array_equal, \
                          assert_almost_equal, assert_equal

from scikits.learn.preprocessing import Scaler, scale, Normalizer, \
                                        LengthNormalizer, Binarizer, \
                                        LabelBinarizer

from scikits.learn.preprocessing.sparse import Normalizer as SparseNormalizer
from scikits.learn.preprocessing.sparse import LengthNormalizer as \
                                               SparseLengthNormalizer
from scikits.learn.preprocessing.sparse import Binarizer as SparseBinarizer

from scikits.learn import datasets
from scikits.learn.linear_model.stochastic_gradient import SGDClassifier

np.random.seed(0)

iris = datasets.load_iris()

def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a

def test_scaler():
    """Test scaling of dataset along all axis
    """
    # First test with 1D data
    X = np.random.randn(5)

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    X = np.random.randn(4, 5)

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert X_scaled is X

    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_array_almost_equal(X_scaled.mean(axis=0), 5*[0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), 5*[1.0])
    # Check that X has not been copied
    assert X_scaled is not X

    X_scaled = scale(X, axis=1, with_std=False)
    assert_array_almost_equal(X_scaled.mean(axis=1), 4*[0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert_array_almost_equal(X_scaled.std(axis=1), 4*[1.0])
    # Check that the data hasn't been modified

def test_normalizer():
    X_ = np.random.randn(4, 5)

    for klass, init in ((Normalizer, np.array),
                        (SparseNormalizer, sp.csr_matrix)):

        X = init(X_.copy())

        normalizer = klass()
        X_norm = normalizer.transform(X, copy=True)
        assert X_norm is not X
        X_norm = toarray(X_norm)
        assert_array_almost_equal(X_norm.sum(axis=1), np.ones(X.shape[0]))

        normalizer = klass()
        X_norm = normalizer.transform(X, copy=False)
        assert X_norm is X
        X_norm = toarray(X_norm)
        assert_array_almost_equal(X_norm.sum(axis=1), np.ones(X.shape[0]))


def test_length_normalizer():
    X_ = np.random.randn(4, 5)

    for klass, init in ((LengthNormalizer, np.array),
                        (SparseLengthNormalizer, sp.csr_matrix)):

        X = init(X_.copy())

        normalizer = klass()
        X_norm1 = normalizer.transform(X, copy=True)
        assert X_norm1 is not X
        X_norm1 = toarray(X_norm1)

        normalizer = klass()
        X_norm2 = normalizer.transform(X, copy=False)
        assert X_norm2 is X
        X_norm2 = toarray(X_norm2)

        for X_norm in (X_norm1, X_norm2):
            for i in xrange(len(X_norm)):
                assert_almost_equal(la.norm(X_norm[i]), 1.0)

def test_binarizer():
    X_ = np.array([[1, 0, 5],
                  [2, 3, 0]])

    for klass, init in ((Binarizer, np.array),
                        (SparseBinarizer, sp.csr_matrix)):

        X = init(X_.copy())

        binarizer = klass(threshold=2.0)
        X_bin = toarray(binarizer.transform(X, copy=True))
        assert_equal(np.sum(X_bin==0), 4)
        assert_equal(np.sum(X_bin==1), 2)

        binarizer = klass()
        X_bin = binarizer.transform(X, copy=True)
        assert X_bin is not X
        X_bin = toarray(X_bin)
        assert_equal(np.sum(X_bin==0), 2)
        assert_equal(np.sum(X_bin==1), 4)

        binarizer = klass()
        X_bin = binarizer.transform(X, copy=False)
        assert X_bin is X
        X_bin = toarray(X_bin)
        assert_equal(np.sum(X_bin==0), 2)
        assert_equal(np.sum(X_bin==1), 4)

def test_label_binarizer():
    lb = LabelBinarizer()

    # two-class case
    inp = np.array([0, 1, 1, 0])
    expected = np.array([[0, 1, 1, 0]]).T
    got = lb.fit_transform(inp)
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)

    # multi-class case
    inp = np.array([3, 2, 1, 2, 0])
    expected = np.array([[0, 0, 0, 1],
                         [0, 0, 1, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [1, 0, 0, 0]])
    got = lb.fit_transform(inp)
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)

def test_label_binarizer_multilabel():
    lb = LabelBinarizer()

    inp = [(2, 3), (1,), (1, 2)]
    expected = np.array([[0, 1, 1],
                         [1, 0, 0],
                         [1, 1, 0]])
    got = lb.fit_transform(inp)
    assert_array_equal(expected, got)
    assert_equal(lb.inverse_transform(got), inp)

def test_label_binarizer_iris():
    lb = LabelBinarizer()
    Y = lb.fit_transform(iris.target)
    clfs = [SGDClassifier().fit(iris.data, Y[:, k])
            for k in range(len(lb.classes_))]
    Y_pred = np.array([clf.decision_function(iris.data) for clf in clfs]).T
    y_pred = lb.inverse_transform(Y_pred)
    accuracy = np.mean(iris.target == y_pred)
    y_pred2 = SGDClassifier().fit(iris.data, iris.target).predict(iris.data)
    accuracy2 = np.mean(iris.target == y_pred2)
    assert_almost_equal(accuracy, accuracy2)

