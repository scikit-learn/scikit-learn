import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false

from sklearn.preprocessing.label import Binarizer
from sklearn.preprocessing.label import LabelBinarizer
from sklearn.preprocessing.label.label import _transform_selected
from sklearn.preprocessing.label import OneHotEncoder
from sklearn.preprocessing.label import LabelEncoder

from sklearn import datasets
from sklearn.linear_model.stochastic_gradient import SGDClassifier

iris = datasets.load_iris()


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def test_label_binarizer():
    lb = LabelBinarizer()

    # two-class case
    inp = ["neg", "pos", "pos", "neg"]
    expected = np.array([[0, 1, 1, 0]]).T
    got = lb.fit_transform(inp)
    assert_false(lb.multilabel_)
    assert_array_equal(lb.classes_, ["neg", "pos"])
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)

    # multi-class case
    inp = ["spam", "ham", "eggs", "ham", "0"]
    expected = np.array([[0, 0, 0, 1],
                         [0, 0, 1, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [1, 0, 0, 0]])
    got = lb.fit_transform(inp)
    assert_array_equal(lb.classes_, ['0', 'eggs', 'ham', 'spam'])
    assert_false(lb.multilabel_)
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)


def test_label_binarizer_set_label_encoding():
    lb = LabelBinarizer(neg_label=-2, pos_label=2)

    # two-class case
    inp = np.array([0, 1, 1, 0])
    expected = np.array([[-2, 2, 2, -2]]).T
    got = lb.fit_transform(inp)
    assert_false(lb.multilabel_)
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)

    # multi-class case
    inp = np.array([3, 2, 1, 2, 0])
    expected = np.array([[-2, -2, -2, +2],
                         [-2, -2, +2, -2],
                         [-2, +2, -2, -2],
                         [-2, -2, +2, -2],
                         [+2, -2, -2, -2]])
    got = lb.fit_transform(inp)
    assert_false(lb.multilabel_)
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)


def test_label_binarizer_multilabel():
    lb = LabelBinarizer()

    # test input as lists of tuples
    inp = [(2, 3), (1,), (1, 2)]
    indicator_mat = np.array([[0, 1, 1],
                              [1, 0, 0],
                              [1, 1, 0]])
    got = lb.fit_transform(inp)
    assert_true(lb.multilabel_)
    assert_array_equal(indicator_mat, got)
    assert_equal(lb.inverse_transform(got), inp)

    # test input as label indicator matrix
    lb.fit(indicator_mat)
    assert_array_equal(indicator_mat,
                       lb.inverse_transform(indicator_mat))

    # regression test for the two-class multilabel case
    lb = LabelBinarizer()
    inp = [[1, 0], [0], [1], [0, 1]]
    expected = np.array([[1, 1],
                         [1, 0],
                         [0, 1],
                         [1, 1]])
    got = lb.fit_transform(inp)
    assert_true(lb.multilabel_)
    assert_array_equal(expected, got)
    assert_equal([set(x) for x in lb.inverse_transform(got)],
                 [set(x) for x in inp])


def test_label_binarizer_errors():
    """Check that invalid arguments yield ValueError"""
    one_class = np.array([0, 0, 0, 0])
    lb = LabelBinarizer().fit(one_class)
    assert_false(lb.multilabel_)

    multi_label = [(2, 3), (0,), (0, 2)]
    assert_raises(ValueError, lb.transform, multi_label)

    lb = LabelBinarizer()
    assert_raises(ValueError, lb.transform, [])
    assert_raises(ValueError, lb.inverse_transform, [])

    assert_raises(ValueError, LabelBinarizer, neg_label=2, pos_label=1)
    assert_raises(ValueError, LabelBinarizer, neg_label=2, pos_label=2)


def test_one_hot_encoder():
    """Test OneHotEncoder's fit and transform."""
    X = [[3, 2, 1], [0, 1, 1]]
    enc = OneHotEncoder()
    # discover max values automatically
    X_trans = enc.fit_transform(X).toarray()
    assert_equal(X_trans.shape, (2, 5))
    assert_array_equal(enc.active_features_,
                       np.where([1, 0, 0, 1, 0, 1, 1, 0, 1])[0])
    assert_array_equal(enc.feature_indices_, [0, 4, 7, 9])

    # check outcome
    assert_array_equal(X_trans,
                       [[0., 1., 0., 1., 1.],
                        [1., 0., 1., 0., 1.]])

    # max value given as 3
    enc = OneHotEncoder(n_values=4)
    X_trans = enc.fit_transform(X)
    assert_equal(X_trans.shape, (2, 4 * 3))
    assert_array_equal(enc.feature_indices_, [0, 4, 8, 12])

    # max value given per feature
    enc = OneHotEncoder(n_values=[3, 2, 2])
    X = [[1, 0, 1], [0, 1, 1]]
    X_trans = enc.fit_transform(X)
    assert_equal(X_trans.shape, (2, 3 + 2 + 2))
    assert_array_equal(enc.n_values_, [3, 2, 2])
    # check that testing with larger feature works:
    X = np.array([[2, 0, 1], [0, 1, 1]])
    enc.transform(X)

    # test that an error is raise when out of bounds:
    X_too_large = [[0, 2, 1], [0, 1, 1]]
    assert_raises(ValueError, enc.transform, X_too_large)

    # test that error is raised when wrong number of features
    assert_raises(ValueError, enc.transform, X[:, :-1])
    # test that error is raised when wrong number of features in fit
    # with prespecified n_values
    assert_raises(ValueError, enc.fit, X[:, :-1])
    # test exception on wrong init param
    assert_raises(TypeError, OneHotEncoder(n_values=np.int).fit, X)

    enc = OneHotEncoder()
    # test negative input to fit
    assert_raises(ValueError, enc.fit, [[0], [-1]])

    # test negative input to transform
    enc.fit([[0], [1]])
    assert_raises(ValueError, enc.transform, [[0], [-1]])


def _check_transform_selected(X, X_expected, sel):
    for M in (X, sparse.csr_matrix(X)):
        Xtr = _transform_selected(M, Binarizer().transform, sel)
        assert_array_equal(toarray(Xtr), X_expected)


def test_transform_selected():
    X = [[3, 2, 1], [0, 1, 1]]

    X_expected = [[1, 2, 1], [0, 1, 1]]
    _check_transform_selected(X, X_expected, [0])
    _check_transform_selected(X, X_expected, [True, False, False])

    X_expected = [[1, 1, 1], [0, 1, 1]]
    _check_transform_selected(X, X_expected, [0, 1, 2])
    _check_transform_selected(X, X_expected, [True, True, True])
    _check_transform_selected(X, X_expected, "all")

    _check_transform_selected(X, X, [])
    _check_transform_selected(X, X, [False, False, False])


def _run_one_hot(X, X2, cat):
    enc = OneHotEncoder(categorical_features=cat)
    Xtr = enc.fit_transform(X)
    X2tr = enc.transform(X2)
    return Xtr, X2tr


def _check_one_hot(X, X2, cat, n_features):
    ind = np.where(cat)[0]
    # With mask
    A, B = _run_one_hot(X, X2, cat)
    # With indices
    C, D = _run_one_hot(X, X2, ind)
    # Check shape
    assert_equal(A.shape, (2, n_features))
    assert_equal(B.shape, (1, n_features))
    assert_equal(C.shape, (2, n_features))
    assert_equal(D.shape, (1, n_features))
    # Check that mask and indices give the same results
    assert_array_equal(toarray(A), toarray(C))
    assert_array_equal(toarray(B), toarray(D))


def test_one_hot_encoder_categorical_features():
    X = np.array([[3, 2, 1], [0, 1, 1]])
    X2 = np.array([[1, 1, 1]])

    cat = [True, False, False]
    _check_one_hot(X, X2, cat, 4)

    # Edge case: all non-categorical
    cat = [False, False, False]
    _check_one_hot(X, X2, cat, 3)

    # Edge case: all categorical
    cat = [True, True, True]
    _check_one_hot(X, X2, cat, 5)


def test_label_encoder():
    """Test LabelEncoder's transform and inverse_transform methods"""
    le = LabelEncoder()
    le.fit([1, 1, 4, 5, -1, 0])
    assert_array_equal(le.classes_, [-1, 0, 1, 4, 5])
    assert_array_equal(le.transform([0, 1, 4, 4, 5, -1, -1]),
                       [1, 2, 3, 3, 4, 0, 0])
    assert_array_equal(le.inverse_transform([1, 2, 3, 3, 4, 0, 0]),
                       [0, 1, 4, 4, 5, -1, -1])
    assert_raises(ValueError, le.transform, [0, 6])


def test_label_encoder_fit_transform():
    """Test fit_transform"""
    le = LabelEncoder()
    ret = le.fit_transform([1, 1, 4, 5, -1, 0])
    assert_array_equal(ret, [2, 2, 3, 4, 0, 1])

    le = LabelEncoder()
    ret = le.fit_transform(["paris", "paris", "tokyo", "amsterdam"])
    assert_array_equal(ret, [1, 1, 2, 0])


def test_label_encoder_string_labels():
    """Test LabelEncoder's transform and inverse_transform methods with
    non-numeric labels"""
    le = LabelEncoder()
    le.fit(["paris", "paris", "tokyo", "amsterdam"])
    assert_array_equal(le.classes_, ["amsterdam", "paris", "tokyo"])
    assert_array_equal(le.transform(["tokyo", "tokyo", "paris"]),
                       [2, 2, 1])
    assert_array_equal(le.inverse_transform([2, 2, 1]),
                       ["tokyo", "tokyo", "paris"])
    assert_raises(ValueError, le.transform, ["london"])


def test_label_encoder_errors():
    """Check that invalid arguments yield ValueError"""
    le = LabelEncoder()
    assert_raises(ValueError, le.transform, [])
    assert_raises(ValueError, le.inverse_transform, [])


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


def test_label_binarizer_multilabel_unlabeled():
    """Check that LabelBinarizer can handle an unlabeled sample"""
    lb = LabelBinarizer()
    y = [[1, 2], [1], []]
    Y = np.array([[1, 1],
                  [1, 0],
                  [0, 0]])
    assert_array_equal(lb.fit_transform(y), Y)
