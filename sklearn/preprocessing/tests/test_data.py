import warnings
import numpy as np
import numpy.linalg as la
from scipy import sparse

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_less_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_warns

from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.preprocessing.data import _transform_selected
from sklearn.preprocessing.data import Binarizer
from sklearn.preprocessing.data import KernelCenterer
from sklearn.preprocessing.data import Normalizer
from sklearn.preprocessing.data import normalize
from sklearn.preprocessing.data import OneHotEncoder
from sklearn.preprocessing.data import StandardScaler
from sklearn.preprocessing.data import scale
from sklearn.preprocessing.data import MinMaxScaler
from sklearn.preprocessing.data import add_dummy_feature
from sklearn.preprocessing.data import PolynomialFeatures

from sklearn import datasets

iris = datasets.load_iris()


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def test_polynomial_features():
    """Test Polynomial Features"""
    X1 = np.arange(6)[:, np.newaxis]
    P1 = np.hstack([np.ones_like(X1),
                    X1, X1 ** 2, X1 ** 3])
    deg1 = 3

    X2 = np.arange(6).reshape((3, 2))
    x1 = X2[:, :1]
    x2 = X2[:, 1:]
    P2 = np.hstack([x1 ** 0 * x2 ** 0,
                    x1 ** 1 * x2 ** 0,
                    x1 ** 0 * x2 ** 1,
                    x1 ** 2 * x2 ** 0,
                    x1 ** 1 * x2 ** 1,
                    x1 ** 0 * x2 ** 2])
    deg2 = 2

    for (deg, X, P) in [(deg1, X1, P1), (deg2, X2, P2)]:
        P_test = PolynomialFeatures(deg, include_bias=True).fit_transform(X)
        assert_array_almost_equal(P_test, P)

        P_test = PolynomialFeatures(deg, include_bias=False).fit_transform(X)
        assert_array_almost_equal(P_test, P[:, 1:])

    interact = PolynomialFeatures(2, interaction_only=True, include_bias=True)
    X_poly = interact.fit_transform(X)
    assert_array_almost_equal(X_poly, P2[:, [0, 1, 2, 4]])

    assert_raises(ValueError, interact.transform, X[:, 1:])


def test_scaler_1d():
    """Test scaling of dataset along single axis"""
    rng = np.random.RandomState(0)
    X = rng.randn(5)
    X_orig_copy = X.copy()

    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    # check inverse transform
    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert_array_almost_equal(X_scaled_back, X_orig_copy)

    # Test with 1D list
    X = [0., 1., 2, 0.4, 1.]
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    X_scaled = scale(X)
    assert_array_almost_equal(X_scaled.mean(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.std(axis=0), 1.0)

    X = np.ones(5)
    assert_array_equal(scale(X, with_mean=False), X)


def test_scaler_2d_arrays():
    """Test scaling of 2d array along first axis"""
    rng = np.random.RandomState(0)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero

    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    assert_array_almost_equal(X_scaled.mean(axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has been copied
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
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)


def test_min_max_scaler_iris():
    X = iris.data
    scaler = MinMaxScaler()
    # default params
    X_trans = scaler.fit_transform(X)
    assert_array_almost_equal(X_trans.min(axis=0), 0)
    assert_array_almost_equal(X_trans.min(axis=0), 0)
    assert_array_almost_equal(X_trans.max(axis=0), 1)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    # not default params: min=1, max=2
    scaler = MinMaxScaler(feature_range=(1, 2))
    X_trans = scaler.fit_transform(X)
    assert_array_almost_equal(X_trans.min(axis=0), 1)
    assert_array_almost_equal(X_trans.max(axis=0), 2)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    # min=-.5, max=.6
    scaler = MinMaxScaler(feature_range=(-.5, .6))
    X_trans = scaler.fit_transform(X)
    assert_array_almost_equal(X_trans.min(axis=0), -.5)
    assert_array_almost_equal(X_trans.max(axis=0), .6)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    # raises on invalid range
    scaler = MinMaxScaler(feature_range=(2, 1))
    assert_raises(ValueError, scaler.fit, X)


def test_min_max_scaler_zero_variance_features():
    """Check min max scaler on toy data with zero variance features"""
    X = [[0., 1., +0.5],
         [0., 1., -0.1],
         [0., 1., +1.1]]

    X_new = [[+0., 2., 0.5],
             [-1., 1., 0.0],
             [+0., 1., 1.5]]

    # default params
    scaler = MinMaxScaler()
    X_trans = scaler.fit_transform(X)
    X_expected_0_1 = [[0., 0., 0.5],
                      [0., 0., 0.0],
                      [0., 0., 1.0]]
    assert_array_almost_equal(X_trans, X_expected_0_1)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    X_trans_new = scaler.transform(X_new)
    X_expected_0_1_new = [[+0., 1., 0.500],
                          [-1., 0., 0.083],
                          [+0., 0., 1.333]]
    assert_array_almost_equal(X_trans_new, X_expected_0_1_new, decimal=2)

    # not default params
    scaler = MinMaxScaler(feature_range=(1, 2))
    X_trans = scaler.fit_transform(X)
    X_expected_1_2 = [[1., 1., 1.5],
                      [1., 1., 1.0],
                      [1., 1., 2.0]]
    assert_array_almost_equal(X_trans, X_expected_1_2)


def test_min_max_scaler_1d():
    """Test scaling of dataset along single axis"""
    rng = np.random.RandomState(0)
    X = rng.randn(5)
    X_orig_copy = X.copy()

    scaler = MinMaxScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_array_almost_equal(X_scaled.min(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.max(axis=0), 1.0)

    # check inverse transform
    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert_array_almost_equal(X_scaled_back, X_orig_copy)

    # Test with 1D list
    X = [0., 1., 2, 0.4, 1.]
    scaler = MinMaxScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_array_almost_equal(X_scaled.min(axis=0), 0.0)
    assert_array_almost_equal(X_scaled.max(axis=0), 1.0)

    # Constant feature.
    X = np.zeros(5)
    scaler = MinMaxScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_greater_equal(X_scaled.min(), 0.)
    assert_less_equal(X_scaled.max(), 1.)


def test_scaler_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    assert_raises(ValueError, StandardScaler().fit, X_csr)

    null_transform = StandardScaler(with_mean=False, with_std=False, copy=True)
    X_null = null_transform.fit_transform(X_csr)
    assert_array_equal(X_null.data, X_csr.data)
    X_orig = null_transform.inverse_transform(X_null)
    assert_array_equal(X_orig.data, X_csr.data)

    scaler = StandardScaler(with_mean=False).fit(X)
    X_scaled = scaler.transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    scaler_csr = StandardScaler(with_mean=False).fit(X_csr)
    X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    scaler_csc = StandardScaler(with_mean=False).fit(X_csc)
    X_csc_scaled = scaler_csr.transform(X_csc, copy=True)
    assert_false(np.any(np.isnan(X_csc_scaled.data)))

    assert_equal(scaler.mean_, scaler_csr.mean_)
    assert_array_almost_equal(scaler.std_, scaler_csr.std_)

    assert_equal(scaler.mean_, scaler_csc.mean_)
    assert_array_almost_equal(scaler.std_, scaler_csc.std_)

    assert_array_almost_equal(
        X_scaled.mean(axis=0), [0., -0.01, 2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(X_csr_scaled, 0)
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
    assert_array_almost_equal(X_csr_scaled_back.toarray(), X)

    X_csc_scaled_back = scaler_csr.inverse_transform(X_csc_scaled.tocsc())
    assert_true(X_csc_scaled_back is not X_csc)
    assert_true(X_csc_scaled_back is not X_csc_scaled)
    assert_array_almost_equal(X_csc_scaled_back.toarray(), X)


def test_scaler_int():
    # test that scaler converts integer input to floating
    # for both sparse and dense matrices
    rng = np.random.RandomState(42)
    X = rng.randint(20, size=(4, 5))
    X[:, 0] = 0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    null_transform = StandardScaler(with_mean=False, with_std=False, copy=True)
    with warnings.catch_warnings(record=True):
        X_null = null_transform.fit_transform(X_csr)
    assert_array_equal(X_null.data, X_csr.data)
    X_orig = null_transform.inverse_transform(X_null)
    assert_array_equal(X_orig.data, X_csr.data)

    with warnings.catch_warnings(record=True):
        scaler = StandardScaler(with_mean=False).fit(X)
        X_scaled = scaler.transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    with warnings.catch_warnings(record=True):
        scaler_csr = StandardScaler(with_mean=False).fit(X_csr)
        X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    with warnings.catch_warnings(record=True):
        scaler_csc = StandardScaler(with_mean=False).fit(X_csc)
        X_csc_scaled = scaler_csr.transform(X_csc, copy=True)
    assert_false(np.any(np.isnan(X_csc_scaled.data)))

    assert_equal(scaler.mean_, scaler_csr.mean_)
    assert_array_almost_equal(scaler.std_, scaler_csr.std_)

    assert_equal(scaler.mean_, scaler_csc.mean_)
    assert_array_almost_equal(scaler.std_, scaler_csc.std_)

    assert_array_almost_equal(
        X_scaled.mean(axis=0),
        [0., 1.109, 1.856, 21., 1.559], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(
        X_csr_scaled.astype(np.float), 0)
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
    assert_array_almost_equal(X_csr_scaled_back.toarray(), X)

    X_csc_scaled_back = scaler_csr.inverse_transform(X_csc_scaled.tocsc())
    assert_true(X_csc_scaled_back is not X_csc)
    assert_true(X_csc_scaled_back is not X_csc_scaled)
    assert_array_almost_equal(X_csc_scaled_back.toarray(), X)


def test_scaler_without_copy():
    """Check that StandardScaler.fit does not change input"""
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)

    X_copy = X.copy()
    StandardScaler(copy=False).fit(X)
    assert_array_equal(X, X_copy)

    X_csr_copy = X_csr.copy()
    StandardScaler(with_mean=False, copy=False).fit(X_csr)
    assert_array_equal(X_csr.toarray(), X_csr_copy.toarray())


def test_scale_sparse_with_mean_raise_exception():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X_csr = sparse.csr_matrix(X)

    # check scaling and fit with direct calls on sparse data
    assert_raises(ValueError, scale, X_csr, with_mean=True)
    assert_raises(ValueError, StandardScaler(with_mean=True).fit, X_csr)

    # check transform and inverse_transform after a fit on a dense array
    scaler = StandardScaler(with_mean=True).fit(X)
    assert_raises(ValueError, scaler.transform, X_csr)

    X_transformed_csr = sparse.csr_matrix(scaler.transform(X))
    assert_raises(ValueError, scaler.inverse_transform, X_transformed_csr)


def test_scale_function_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)

    X_scaled = scale(X, with_mean=False)
    assert_false(np.any(np.isnan(X_scaled)))

    X_csr_scaled = scale(X_csr, with_mean=False)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    # test csc has same outcome
    X_csc_scaled = scale(X_csr.tocsc(), with_mean=False)
    assert_array_almost_equal(X_scaled, X_csc_scaled.toarray())

    # raises value error on axis != 0
    assert_raises(ValueError, scale, X_csr, with_mean=False, axis=1)

    assert_array_almost_equal(X_scaled.mean(axis=0),
                              [0., -0.01, 2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(X_csr_scaled, 0)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))


def test_warning_scaling_integers():
    """Check warning when scaling integer data"""
    X = np.array([[1, 2, 0],
                  [0, 0, 0]], dtype=np.uint8)

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        assert_warns(UserWarning, StandardScaler().fit, X)

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        assert_warns(UserWarning, MinMaxScaler().fit, X)


def test_normalizer_l1():
    rng = np.random.RandomState(0)
    X_dense = rng.randn(4, 5)
    X_sparse_unpruned = sparse.csr_matrix(X_dense)

    # set the row number 3 to zero
    X_dense[3, :] = 0.0

    # set the row number 3 to zero without pruning (can happen in real life)
    indptr_3 = X_sparse_unpruned.indptr[3]
    indptr_4 = X_sparse_unpruned.indptr[4]
    X_sparse_unpruned.data[indptr_3:indptr_4] = 0.0

    # build the pruned variant using the regular constructor
    X_sparse_pruned = sparse.csr_matrix(X_dense)

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
    for init in (sparse.coo_matrix, sparse.csc_matrix, sparse.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert_true(X_norm is not X)
        assert_true(isinstance(X_norm, sparse.csr_matrix))

        X_norm = toarray(X_norm)
        for i in range(3):
            assert_almost_equal(row_sums[i], 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalizer_l2():
    rng = np.random.RandomState(0)
    X_dense = rng.randn(4, 5)
    X_sparse_unpruned = sparse.csr_matrix(X_dense)

    # set the row number 3 to zero
    X_dense[3, :] = 0.0

    # set the row number 3 to zero without pruning (can happen in real life)
    indptr_3 = X_sparse_unpruned.indptr[3]
    indptr_4 = X_sparse_unpruned.indptr[4]
    X_sparse_unpruned.data[indptr_3:indptr_4] = 0.0

    # build the pruned variant using the regular constructor
    X_sparse_pruned = sparse.csr_matrix(X_dense)

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
            for i in range(3):
                assert_almost_equal(la.norm(X_norm[i]), 1.0)
            assert_almost_equal(la.norm(X_norm[3]), 0.0)

    # check input for which copy=False won't prevent a copy
    for init in (sparse.coo_matrix, sparse.csc_matrix, sparse.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert_true(X_norm is not X)
        assert_true(isinstance(X_norm, sparse.csr_matrix))

        X_norm = toarray(X_norm)
        for i in range(3):
            assert_almost_equal(la.norm(X_norm[i]), 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalize():
    """Test normalize function"""
    # Only tests functionality not used by the tests for Normalizer.
    X = np.random.RandomState(37).randn(3, 2)
    assert_array_equal(normalize(X, copy=False),
                       normalize(X.T, axis=0, copy=False).T)
    assert_raises(ValueError, normalize, [[0]], axis=2)
    assert_raises(ValueError, normalize, [[0]], norm='l3')


def test_binarizer():
    X_ = np.array([[1, 0, 5], [2, 3, -1]])

    for init in (np.array, list, sparse.csr_matrix, sparse.csc_matrix):

        X = init(X_.copy())

        binarizer = Binarizer(threshold=2.0, copy=True)
        X_bin = toarray(binarizer.transform(X))
        assert_equal(np.sum(X_bin == 0), 4)
        assert_equal(np.sum(X_bin == 1), 2)
        X_bin = binarizer.transform(X)
        assert_equal(sparse.issparse(X), sparse.issparse(X_bin))

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
        if init is not list:
            assert_true(X_bin is X)
        X_bin = toarray(X_bin)
        assert_equal(np.sum(X_bin == 0), 2)
        assert_equal(np.sum(X_bin == 1), 4)

    binarizer = Binarizer(threshold=-0.5, copy=True)
    for init in (np.array, list):
        X = init(X_.copy())

        X_bin = toarray(binarizer.transform(X))
        assert_equal(np.sum(X_bin == 0), 1)
        assert_equal(np.sum(X_bin == 1), 5)
        X_bin = binarizer.transform(X)

    # Cannot use threshold < 0 for sparse
    assert_raises(ValueError, binarizer.transform, sparse.csc_matrix(X))


def test_center_kernel():
    """Test that KernelCenterer is equivalent to StandardScaler
       in feature space"""
    rng = np.random.RandomState(0)
    X_fit = rng.random_sample((5, 4))
    scaler = StandardScaler(with_std=False)
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
    for obj in ((StandardScaler(), Normalizer(), Binarizer())):
        X_transformed = obj.fit(X).transform(X)
        X_transformed2 = obj.fit_transform(X)
        assert_array_equal(X_transformed, X_transformed2)


def test_add_dummy_feature():
    X = [[1, 0], [0, 1], [0, 1]]
    X = add_dummy_feature(X)
    assert_array_equal(X, [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_add_dummy_feature_coo():
    X = sparse.coo_matrix([[1, 0], [0, 1], [0, 1]])
    X = add_dummy_feature(X)
    assert_true(sparse.isspmatrix_coo(X), X)
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_add_dummy_feature_csc():
    X = sparse.csc_matrix([[1, 0], [0, 1], [0, 1]])
    X = add_dummy_feature(X)
    assert_true(sparse.isspmatrix_csc(X), X)
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_add_dummy_feature_csr():
    X = sparse.csr_matrix([[1, 0], [0, 1], [0, 1]])
    X = add_dummy_feature(X)
    assert_true(sparse.isspmatrix_csr(X), X)
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_one_hot_encoder_sparse():
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

    # test that an error is raised when out of bounds:
    X_too_large = [[0, 2, 1], [0, 1, 1]]
    assert_raises(ValueError, enc.transform, X_too_large)
    assert_raises(ValueError, OneHotEncoder(n_values=2).fit_transform, X)

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

def test_one_hot_encoder_dense():
    """check for sparse=False"""
    X = [[3, 2, 1], [0, 1, 1]]
    enc = OneHotEncoder(sparse=False)
    # discover max values automatically
    X_trans = enc.fit_transform(X)
    assert_equal(X_trans.shape, (2, 5))
    assert_array_equal(enc.active_features_,
                       np.where([1, 0, 0, 1, 0, 1, 1, 0, 1])[0])
    assert_array_equal(enc.feature_indices_, [0, 4, 7, 9])

    # check outcome
    assert_array_equal(X_trans,
                       np.array([[0., 1., 0., 1., 1.],
                                 [1., 0., 1., 0., 1.]]))


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
