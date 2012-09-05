import numpy as np
import numpy.linalg as la
import scipy.sparse as sp

from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal

from nose.tools import assert_raises, assert_true, assert_false

from sklearn.utils.sparsefuncs import mean_variance_axis0
from sklearn.preprocessing import Binarizer
from sklearn.preprocessing import KernelCenterer
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import normalize
from sklearn.preprocessing import Scaler
from sklearn.preprocessing import scale

from sklearn import datasets
from sklearn.linear_model.stochastic_gradient import SGDClassifier

iris = datasets.load_iris()


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def test_scaler_1d():
    """Test scaling of dataset along single axis"""
    rng = np.random.RandomState(0)
    X = rng.randn(5)
    X_orig_copy = X.copy()

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    # check inverse transform
    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert_array_almost_equal(X_scaled_back, X_orig_copy)

    # Test with 1D list
    X = [0., 1., 2, 0.4, 1.]
    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    X_scaled = scale(X)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)


def test_scaler_2d_arrays():
    """Test scaling of 2d array along first axis"""
    rng = np.random.RandomState(0)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero

    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    assert_array_almost_equal(X_scaled.mean(axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)

    # check inverse transform
    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert_true(X_scaled_back is not X)
    assert_true(X_scaled_back is not X_scaled)
    assert_array_almost_equal(X_scaled_back, X)

    X_scaled = scale(X, axis=1, with_std=False)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=1), 4 * [0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=1), 4 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=1), 4 * [1.0])
    # Check that the data hasn't been modified
    assert_true(X_scaled is not X)

    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is X)

    X = rng.randn(4, 5)
    X[:, 0] = 1.0  # first feature is a constant, non zero feature
    scaler = Scaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)


def test_scaler_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sp.csr_matrix(X)

    scaler = Scaler(with_mean=False).fit(X)
    X_scaled = scaler.transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    scaler_csr = Scaler(with_mean=False).fit(X_csr)
    X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    assert_equal(scaler.mean_, scaler_csr.mean_)
    assert_array_almost_equal(scaler.std_, scaler_csr.std_)

    assert_array_almost_equal(
        X_scaled.mean(axis=0), [0., -0.01,  2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis0(X_csr_scaled)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))

    # Check that X has not been modified (copy)
    assert_true(X_scaled is not X)
    assert_true(X_csr_scaled is not X_csr)

    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert_true(X_scaled_back is not X)
    assert_true(X_scaled_back is not X_scaled)
    assert_array_almost_equal(X_scaled_back, X)

    X_csr_scaled_back = scaler_csr.inverse_transform(X_csr_scaled)
    assert_true(X_csr_scaled_back is not X_csr)
    assert_true(X_csr_scaled_back is not X_csr_scaled)
    assert_array_almost_equal(X_scaled_back, X)


def test_scaler_without_copy():
    """Check that Scaler.fit does not change input"""
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sp.csr_matrix(X)

    X_copy = X.copy()
    Scaler(copy=False).fit(X)
    assert_array_equal(X, X_copy)

    X_csr_copy = X_csr.copy()
    Scaler(with_mean=False, copy=False).fit(X_csr)
    assert_array_equal(X_csr.toarray(), X_csr_copy.toarray())


def test_scale_sparse_with_mean_raise_exception():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X_csr = sp.csr_matrix(X)

    # check scaling and fit with direct calls on sparse data
    assert_raises(ValueError, scale, X_csr, with_mean=True)
    assert_raises(ValueError, Scaler(with_mean=True).fit, X_csr)

    # check transform and inverse_transform after a fit on a dense array
    scaler = Scaler(with_mean=True).fit(X)
    assert_raises(ValueError, scaler.transform, X_csr)

    X_transformed_csr = sp.csr_matrix(scaler.transform(X))
    assert_raises(ValueError, scaler.inverse_transform, X_transformed_csr)


def test_scale_function_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sp.csr_matrix(X)

    X_scaled = scale(X, with_mean=False)
    assert_false(np.any(np.isnan(X_scaled)))

    X_csr_scaled = scale(X_csr, with_mean=False)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    assert_array_almost_equal(
        X_scaled.mean(axis=0), [0., -0.01,  2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis0(X_csr_scaled)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))


def test_normalizer_l1():
    rng = np.random.RandomState(0)
    X_dense = rng.randn(4, 5)
    X_sparse_unpruned = sp.csr_matrix(X_dense)

    # set the row number 3 to zero
    X_dense[3, :] = 0.0

    # set the row number 3 to zero without pruning (can happen in real life)
    indptr_3 = X_sparse_unpruned.indptr[3]
    indptr_4 = X_sparse_unpruned.indptr[4]
    X_sparse_unpruned.data[indptr_3:indptr_4] = 0.0

    # build the pruned variant using the regular constructor
    X_sparse_pruned = sp.csr_matrix(X_dense)

    # check inputs that support the no-copy optim
    for X in (X_dense, X_sparse_pruned, X_sparse_unpruned):

        normalizer = Normalizer(norm='l1', copy=True)
        X_norm = normalizer.transform(X)
        assert_true(X_norm is not X)
        X_norm1 = toarray(X_norm)

        normalizer = Normalizer(norm='l1', copy=False)
        X_norm = normalizer.transform(X)
        assert_true(X_norm is X)
        X_norm2 = toarray(X_norm)

        for X_norm in (X_norm1, X_norm2):
            row_sums = np.abs(X_norm).sum(axis=1)
            for i in range(3):
                assert_almost_equal(row_sums[i], 1.0)
            assert_almost_equal(row_sums[3], 0.0)

    # check input for which copy=False won't prevent a copy
    for init in (sp.coo_matrix, sp.csc_matrix, sp.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert_true(X_norm is not X)
        assert_true(isinstance(X_norm, sp.csr_matrix))

        X_norm = toarray(X_norm)
        for i in xrange(3):
            assert_almost_equal(row_sums[i], 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalizer_l2():
    rng = np.random.RandomState(0)
    X_dense = rng.randn(4, 5)
    X_sparse_unpruned = sp.csr_matrix(X_dense)

    # set the row number 3 to zero
    X_dense[3, :] = 0.0

    # set the row number 3 to zero without pruning (can happen in real life)
    indptr_3 = X_sparse_unpruned.indptr[3]
    indptr_4 = X_sparse_unpruned.indptr[4]
    X_sparse_unpruned.data[indptr_3:indptr_4] = 0.0

    # build the pruned variant using the regular constructor
    X_sparse_pruned = sp.csr_matrix(X_dense)

    # check inputs that support the no-copy optim
    for X in (X_dense, X_sparse_pruned, X_sparse_unpruned):

        normalizer = Normalizer(norm='l2', copy=True)
        X_norm1 = normalizer.transform(X)
        assert_true(X_norm1 is not X)
        X_norm1 = toarray(X_norm1)

        normalizer = Normalizer(norm='l2', copy=False)
        X_norm2 = normalizer.transform(X)
        assert_true(X_norm2 is X)
        X_norm2 = toarray(X_norm2)

        for X_norm in (X_norm1, X_norm2):
            for i in xrange(3):
                assert_almost_equal(la.norm(X_norm[i]), 1.0)
            assert_almost_equal(la.norm(X_norm[3]), 0.0)

    # check input for which copy=False won't prevent a copy
    for init in (sp.coo_matrix, sp.csc_matrix, sp.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert_true(X_norm is not X)
        assert_true(isinstance(X_norm, sp.csr_matrix))

        X_norm = toarray(X_norm)
        for i in xrange(3):
            assert_almost_equal(la.norm(X_norm[i]), 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalize_errors():
    """Check that invalid arguments yield ValueError"""
    assert_raises(ValueError, normalize, [[0]], axis=2)
    assert_raises(ValueError, normalize, [[0]], norm='l3')


def test_binarizer():
    X_ = np.array([[1, 0, 5], [2, 3, 0]])

    for init in (np.array, sp.csr_matrix):

        X = init(X_.copy())

        binarizer = Binarizer(threshold=2.0, copy=True)
        X_bin = toarray(binarizer.transform(X))
        assert_equal(np.sum(X_bin == 0), 4)
        assert_equal(np.sum(X_bin == 1), 2)

        binarizer = Binarizer(copy=True).fit(X)
        X_bin = toarray(binarizer.transform(X))
        assert_true(X_bin is not X)
        assert_equal(np.sum(X_bin == 0), 2)
        assert_equal(np.sum(X_bin == 1), 4)

        binarizer = Binarizer(copy=True)
        X_bin = binarizer.transform(X)
        assert_true(X_bin is not X)
        X_bin = toarray(X_bin)
        assert_equal(np.sum(X_bin == 0), 2)
        assert_equal(np.sum(X_bin == 1), 4)

        binarizer = Binarizer(copy=False)
        X_bin = binarizer.transform(X)
        assert_true(X_bin is X)
        X_bin = toarray(X_bin)
        assert_equal(np.sum(X_bin == 0), 2)
        assert_equal(np.sum(X_bin == 1), 4)


def test_label_binarizer():
    lb = LabelBinarizer()

    # two-class case
    inp = ["neg", "pos", "pos", "neg"]
    expected = np.array([[0, 1, 1, 0]]).T
    got = lb.fit_transform(inp)
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
    assert_array_equal(expected, got)
    assert_array_equal(lb.inverse_transform(got), inp)


def test_label_binarizer_set_label_encoding():
    lb = LabelBinarizer(neg_label=-2, pos_label=2)

    # two-class case
    inp = np.array([0, 1, 1, 0])
    expected = np.array([[-2, 2, 2, -2]]).T
    got = lb.fit_transform(inp)
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
    assert_array_equal(expected, got)
    assert_equal([set(x) for x in lb.inverse_transform(got)],
                 [set(x) for x in inp])


def test_label_binarizer_errors():
    """Check that invalid arguments yield ValueError"""
    one_class = np.array([0, 0, 0, 0])
    lb = LabelBinarizer().fit(one_class)

    multi_label = [(2, 3), (0,), (0, 2)]
    assert_raises(ValueError, lb.transform, multi_label)

    lb = LabelBinarizer()
    assert_raises(ValueError, lb.transform, [])
    assert_raises(ValueError, lb.inverse_transform, [])

    assert_raises(ValueError, LabelBinarizer, neg_label=2, pos_label=1)
    assert_raises(ValueError, LabelBinarizer, neg_label=2, pos_label=2)


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
    assert_equal(lb.fit_transform(y), Y)


def test_center_kernel():
    """Test that KernelCenterer is equivalent to Scaler in feature space"""
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    scaler = Scaler(with_std=False)
    scaler.fit(X_fit)
    X_fit_centered = scaler.transform(X_fit)
    K_fit = np.dot(X_fit, X_fit.T)

    # center fit time matrix
    centerer = KernelCenterer()
    K_fit_centered = np.dot(X_fit_centered, X_fit_centered.T)
    K_fit_centered2 = centerer.fit_transform(K_fit)
    assert_array_almost_equal(K_fit_centered, K_fit_centered2)

    # center predict time matrix
    X_pred = rng.random_sample((2, 4))
    K_pred = np.dot(X_pred, X_fit.T)
    X_pred_centered = scaler.transform(X_pred)
    K_pred_centered = np.dot(X_pred_centered, X_fit_centered.T)
    K_pred_centered2 = centerer.transform(K_pred)
    assert_array_almost_equal(K_pred_centered, K_pred_centered2)


def test_fit_transform():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    for obj in ((Scaler(), Normalizer(), Binarizer())):
        X_transformed = obj.fit(X).transform(X)
        X_transformed2 = obj.fit_transform(X)
        assert_array_equal(X_transformed, X_transformed2)
