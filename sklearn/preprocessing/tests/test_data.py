# Authors:
#
#          Giorgio Patrini
#
# License: BSD 3 clause

import warnings
import itertools

import numpy as np
import numpy.linalg as la
from scipy import sparse, stats
from scipy.sparse import random as sparse_random

import pytest

from sklearn.utils import gen_batches

from sklearn.utils._testing import assert_almost_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_less
from sklearn.utils._testing import assert_warns_message
from sklearn.utils._testing import assert_no_warnings
from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import assert_allclose_dense_sparse
from sklearn.utils._testing import skip_if_32bit
from sklearn.utils._testing import _convert_container

from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.preprocessing._data import _handle_zeros_in_scale
from sklearn.preprocessing._data import Binarizer
from sklearn.preprocessing._data import KernelCenterer
from sklearn.preprocessing._data import Normalizer
from sklearn.preprocessing._data import normalize
from sklearn.preprocessing._data import StandardScaler
from sklearn.preprocessing._data import scale
from sklearn.preprocessing._data import MinMaxScaler
from sklearn.preprocessing._data import minmax_scale
from sklearn.preprocessing._data import QuantileTransformer
from sklearn.preprocessing._data import quantile_transform
from sklearn.preprocessing._data import MaxAbsScaler
from sklearn.preprocessing._data import maxabs_scale
from sklearn.preprocessing._data import RobustScaler
from sklearn.preprocessing._data import robust_scale
from sklearn.preprocessing._data import add_dummy_feature
from sklearn.preprocessing._data import PolynomialFeatures
from sklearn.preprocessing._data import PowerTransformer
from sklearn.preprocessing._data import power_transform
from sklearn.preprocessing._data import BOUNDS_THRESHOLD
from sklearn.exceptions import NotFittedError

from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_predict
from sklearn.svm import SVR
from sklearn.utils import shuffle

from sklearn import datasets

iris = datasets.load_iris()

# Make some data to be used many times
rng = np.random.RandomState(0)
n_features = 30
n_samples = 1000
offsets = rng.uniform(-1, 1, size=n_features)
scales = rng.uniform(1, 10, size=n_features)
X_2d = rng.randn(n_samples, n_features) * scales + offsets
X_1row = X_2d[0, :].reshape(1, n_features)
X_1col = X_2d[:, 0].reshape(n_samples, 1)
X_list_1row = X_1row.tolist()
X_list_1col = X_1col.tolist()


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def _check_dim_1axis(a):
    if isinstance(a, list):
        return np.array(a).shape[0]
    return a.shape[0]


def assert_correct_incr(i, batch_start, batch_stop, n, chunk_size,
                        n_samples_seen):
    if batch_stop != n:
        assert (i + 1) * chunk_size == n_samples_seen
    else:
        assert (i * chunk_size + (batch_stop - batch_start) ==
                     n_samples_seen)


def test_polynomial_features():
    # Test Polynomial Features
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

    assert interact.powers_.shape == (interact.n_output_features_,
                                      interact.n_input_features_)


def test_polynomial_feature_names():
    X = np.arange(30).reshape(10, 3)
    poly = PolynomialFeatures(degree=2, include_bias=True).fit(X)
    feature_names = poly.get_feature_names()
    assert_array_equal(['1', 'x0', 'x1', 'x2', 'x0^2', 'x0 x1',
                        'x0 x2', 'x1^2', 'x1 x2', 'x2^2'],
                       feature_names)

    poly = PolynomialFeatures(degree=3, include_bias=False).fit(X)
    feature_names = poly.get_feature_names(["a", "b", "c"])
    assert_array_equal(['a', 'b', 'c', 'a^2', 'a b', 'a c', 'b^2',
                        'b c', 'c^2', 'a^3', 'a^2 b', 'a^2 c',
                        'a b^2', 'a b c', 'a c^2', 'b^3', 'b^2 c',
                        'b c^2', 'c^3'], feature_names)
    # test some unicode
    poly = PolynomialFeatures(degree=1, include_bias=True).fit(X)
    feature_names = poly.get_feature_names(
        ["\u0001F40D", "\u262E", "\u05D0"])
    assert_array_equal(["1", "\u0001F40D", "\u262E", "\u05D0"],
                       feature_names)


def test_polynomial_feature_array_order():
    X = np.arange(10).reshape(5, 2)

    def is_c_contiguous(a):
        return np.isfortran(a.T)

    assert is_c_contiguous(PolynomialFeatures().fit_transform(X))
    assert is_c_contiguous(PolynomialFeatures(order='C').fit_transform(X))
    assert np.isfortran(PolynomialFeatures(order='F').fit_transform(X))


@pytest.mark.parametrize(['deg', 'include_bias', 'interaction_only', 'dtype'],
                         [(1, True, False, int),
                          (2, True, False, int),
                          (2, True, False, np.float32),
                          (2, True, False, np.float64),
                          (3, False, False, np.float64),
                          (3, False, True, np.float64),
                          (4, False, False, np.float64),
                          (4, False, True, np.float64)])
def test_polynomial_features_csc_X(deg, include_bias, interaction_only, dtype):
    rng = np.random.RandomState(0)
    X = rng.randint(0, 2, (100, 2))
    X_csc = sparse.csc_matrix(X)

    est = PolynomialFeatures(deg, include_bias=include_bias,
                             interaction_only=interaction_only)
    Xt_csc = est.fit_transform(X_csc.astype(dtype))
    Xt_dense = est.fit_transform(X.astype(dtype))

    assert isinstance(Xt_csc, sparse.csc_matrix)
    assert Xt_csc.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csc.A, Xt_dense)


@pytest.mark.parametrize(['deg', 'include_bias', 'interaction_only', 'dtype'],
                         [(1, True, False, int),
                          (2, True, False, int),
                          (2, True, False, np.float32),
                          (2, True, False, np.float64),
                          (3, False, False, np.float64),
                          (3, False, True, np.float64)])
def test_polynomial_features_csr_X(deg, include_bias, interaction_only, dtype):
    rng = np.random.RandomState(0)
    X = rng.randint(0, 2, (100, 2))
    X_csr = sparse.csr_matrix(X)

    est = PolynomialFeatures(deg, include_bias=include_bias,
                             interaction_only=interaction_only)
    Xt_csr = est.fit_transform(X_csr.astype(dtype))
    Xt_dense = est.fit_transform(X.astype(dtype, copy=False))

    assert isinstance(Xt_csr, sparse.csr_matrix)
    assert Xt_csr.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csr.A, Xt_dense)


@pytest.mark.parametrize(['deg', 'include_bias', 'interaction_only', 'dtype'],
                         [(2, True, False, np.float32),
                          (2, True, False, np.float64),
                          (3, False, False, np.float64),
                          (3, False, True, np.float64)])
def test_polynomial_features_csr_X_floats(deg, include_bias,
                                          interaction_only, dtype):
    X_csr = sparse_random(1000, 10, 0.5, random_state=0).tocsr()
    X = X_csr.toarray()

    est = PolynomialFeatures(deg, include_bias=include_bias,
                             interaction_only=interaction_only)
    Xt_csr = est.fit_transform(X_csr.astype(dtype))
    Xt_dense = est.fit_transform(X.astype(dtype))

    assert isinstance(Xt_csr, sparse.csr_matrix)
    assert Xt_csr.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csr.A, Xt_dense)


@pytest.mark.parametrize(['zero_row_index', 'deg', 'interaction_only'],
                         [(0, 2, True), (1, 2, True), (2, 2, True),
                          (0, 3, True), (1, 3, True), (2, 3, True),
                          (0, 2, False), (1, 2, False), (2, 2, False),
                          (0, 3, False), (1, 3, False), (2, 3, False)])
def test_polynomial_features_csr_X_zero_row(zero_row_index, deg,
                                            interaction_only):
    X_csr = sparse_random(3, 10, 1.0, random_state=0).tocsr()
    X_csr[zero_row_index, :] = 0.0
    X = X_csr.toarray()

    est = PolynomialFeatures(deg, include_bias=False,
                             interaction_only=interaction_only)
    Xt_csr = est.fit_transform(X_csr)
    Xt_dense = est.fit_transform(X)

    assert isinstance(Xt_csr, sparse.csr_matrix)
    assert Xt_csr.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csr.A, Xt_dense)


# This degree should always be one more than the highest degree supported by
# _csr_expansion.
@pytest.mark.parametrize(['include_bias', 'interaction_only'],
                         [(True, True), (True, False),
                          (False, True), (False, False)])
def test_polynomial_features_csr_X_degree_4(include_bias, interaction_only):
    X_csr = sparse_random(1000, 10, 0.5, random_state=0).tocsr()
    X = X_csr.toarray()

    est = PolynomialFeatures(4, include_bias=include_bias,
                             interaction_only=interaction_only)
    Xt_csr = est.fit_transform(X_csr)
    Xt_dense = est.fit_transform(X)

    assert isinstance(Xt_csr, sparse.csr_matrix)
    assert Xt_csr.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csr.A, Xt_dense)


@pytest.mark.parametrize(['deg', 'dim', 'interaction_only'],
                         [(2, 1, True),
                          (2, 2, True),
                          (3, 1, True),
                          (3, 2, True),
                          (3, 3, True),
                          (2, 1, False),
                          (2, 2, False),
                          (3, 1, False),
                          (3, 2, False),
                          (3, 3, False)])
def test_polynomial_features_csr_X_dim_edges(deg, dim, interaction_only):
    X_csr = sparse_random(1000, dim, 0.5, random_state=0).tocsr()
    X = X_csr.toarray()

    est = PolynomialFeatures(deg, interaction_only=interaction_only)
    Xt_csr = est.fit_transform(X_csr)
    Xt_dense = est.fit_transform(X)

    assert isinstance(Xt_csr, sparse.csr_matrix)
    assert Xt_csr.dtype == Xt_dense.dtype
    assert_array_almost_equal(Xt_csr.A, Xt_dense)


def test_standard_scaler_1d():
    # Test scaling of dataset along single axis
    for X in [X_1row, X_1col, X_list_1row, X_list_1row]:

        scaler = StandardScaler()
        X_scaled = scaler.fit(X).transform(X, copy=True)

        if isinstance(X, list):
            X = np.array(X)  # cast only after scaling done

        if _check_dim_1axis(X) == 1:
            assert_almost_equal(scaler.mean_, X.ravel())
            assert_almost_equal(scaler.scale_, np.ones(n_features))
            assert_array_almost_equal(X_scaled.mean(axis=0),
                                      np.zeros_like(n_features))
            assert_array_almost_equal(X_scaled.std(axis=0),
                                      np.zeros_like(n_features))
        else:
            assert_almost_equal(scaler.mean_, X.mean())
            assert_almost_equal(scaler.scale_, X.std())
            assert_array_almost_equal(X_scaled.mean(axis=0),
                                      np.zeros_like(n_features))
            assert_array_almost_equal(X_scaled.mean(axis=0), .0)
            assert_array_almost_equal(X_scaled.std(axis=0), 1.)
        assert scaler.n_samples_seen_ == X.shape[0]

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones((5, 1))
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_almost_equal(scaler.mean_, 1.)
    assert_almost_equal(scaler.scale_, 1.)
    assert_array_almost_equal(X_scaled.mean(axis=0), .0)
    assert_array_almost_equal(X_scaled.std(axis=0), .0)
    assert scaler.n_samples_seen_ == X.shape[0]


def test_standard_scaler_dtype():
    # Ensure scaling does not affect dtype
    rng = np.random.RandomState(0)
    n_samples = 10
    n_features = 3
    for dtype in [np.float16, np.float32, np.float64]:
        X = rng.randn(n_samples, n_features).astype(dtype)
        scaler = StandardScaler()
        X_scaled = scaler.fit(X).transform(X)
        assert X.dtype == X_scaled.dtype
        assert scaler.mean_.dtype == np.float64
        assert scaler.scale_.dtype == np.float64


def test_scale_1d():
    # 1-d inputs
    X_list = [1., 3., 5., 0.]
    X_arr = np.array(X_list)

    for X in [X_list, X_arr]:
        X_scaled = scale(X)
        assert_array_almost_equal(X_scaled.mean(), 0.0)
        assert_array_almost_equal(X_scaled.std(), 1.0)
        assert_array_equal(scale(X, with_mean=False, with_std=False), X)


@skip_if_32bit
def test_standard_scaler_numerical_stability():
    # Test numerical stability of scaling
    # np.log(1e-5) is taken because of its floating point representation
    # was empirically found to cause numerical problems with np.mean & np.std.

    x = np.full(8, np.log(1e-5), dtype=np.float64)
    # This does not raise a warning as the number of samples is too low
    # to trigger the problem in recent numpy
    x_scaled = assert_no_warnings(scale, x)
    assert_array_almost_equal(scale(x), np.zeros(8))

    # with 2 more samples, the std computation run into numerical issues:
    x = np.full(10, np.log(1e-5), dtype=np.float64)
    w = "standard deviation of the data is probably very close to 0"
    x_scaled = assert_warns_message(UserWarning, w, scale, x)
    assert_array_almost_equal(x_scaled, np.zeros(10))

    x = np.full(10, 1e-100, dtype=np.float64)
    x_small_scaled = assert_no_warnings(scale, x)
    assert_array_almost_equal(x_small_scaled, np.zeros(10))

    # Large values can cause (often recoverable) numerical stability issues:
    x_big = np.full(10, 1e100, dtype=np.float64)
    w = "Dataset may contain too large values"
    x_big_scaled = assert_warns_message(UserWarning, w, scale, x_big)
    assert_array_almost_equal(x_big_scaled, np.zeros(10))
    assert_array_almost_equal(x_big_scaled, x_small_scaled)

    x_big_centered = assert_warns_message(UserWarning, w, scale, x_big,
                                          with_std=False)
    assert_array_almost_equal(x_big_centered, np.zeros(10))
    assert_array_almost_equal(x_big_centered, x_small_scaled)


def test_scaler_2d_arrays():
    # Test scaling of 2d array along first axis
    rng = np.random.RandomState(0)
    n_features = 5
    n_samples = 4
    X = rng.randn(n_samples, n_features)
    X[:, 0] = 0.0  # first feature is always of zero

    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert not np.any(np.isnan(X_scaled))
    assert scaler.n_samples_seen_ == n_samples

    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has been copied
    assert X_scaled is not X

    # check inverse transform
    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert X_scaled_back is not X
    assert X_scaled_back is not X_scaled
    assert_array_almost_equal(X_scaled_back, X)

    X_scaled = scale(X, axis=1, with_std=False)
    assert not np.any(np.isnan(X_scaled))
    assert_array_almost_equal(X_scaled.mean(axis=1), n_samples * [0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert not np.any(np.isnan(X_scaled))
    assert_array_almost_equal(X_scaled.mean(axis=1), n_samples * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=1), n_samples * [1.0])
    # Check that the data hasn't been modified
    assert X_scaled is not X

    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert not np.any(np.isnan(X_scaled))
    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert X_scaled is X

    X = rng.randn(4, 5)
    X[:, 0] = 1.0  # first feature is a constant, non zero feature
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert not np.any(np.isnan(X_scaled))
    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert X_scaled is not X


def test_scaler_float16_overflow():
    # Test if the scaler will not overflow on float16 numpy arrays
    rng = np.random.RandomState(0)
    # float16 has a maximum of 65500.0. On the worst case 5 * 200000 is 100000
    # which is enough to overflow the data type
    X = rng.uniform(5, 10, [200000, 1]).astype(np.float16)

    with np.errstate(over='raise'):
        scaler = StandardScaler().fit(X)
        X_scaled = scaler.transform(X)

    # Calculate the float64 equivalent to verify result
    X_scaled_f64 = StandardScaler().fit_transform(X.astype(np.float64))

    # Overflow calculations may cause -inf, inf, or nan. Since there is no nan
    # input, all of the outputs should be finite. This may be redundant since a
    # FloatingPointError exception will be thrown on overflow above.
    assert np.all(np.isfinite(X_scaled))

    # The normal distribution is very unlikely to go above 4. At 4.0-8.0 the
    # float16 precision is 2^-8 which is around 0.004. Thus only 2 decimals are
    # checked to account for precision differences.
    assert_array_almost_equal(X_scaled, X_scaled_f64, decimal=2)


def test_handle_zeros_in_scale():
    s1 = np.array([0, 1, 2, 3])
    s2 = _handle_zeros_in_scale(s1, copy=True)

    assert not s1[0] == s2[0]
    assert_array_equal(s1, np.array([0, 1, 2, 3]))
    assert_array_equal(s2, np.array([1, 1, 2, 3]))


def test_minmax_scaler_partial_fit():
    # Test if partial_fit run over many batches of size 1 and 50
    # gives the same results as fit
    X = X_2d
    n = X.shape[0]

    for chunk_size in [1, 2, 50, n, n + 42]:
        # Test mean at the end of the process
        scaler_batch = MinMaxScaler().fit(X)

        scaler_incr = MinMaxScaler()
        for batch in gen_batches(n_samples, chunk_size):
            scaler_incr = scaler_incr.partial_fit(X[batch])

        assert_array_almost_equal(scaler_batch.data_min_,
                                  scaler_incr.data_min_)
        assert_array_almost_equal(scaler_batch.data_max_,
                                  scaler_incr.data_max_)
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_
        assert_array_almost_equal(scaler_batch.data_range_,
                                  scaler_incr.data_range_)
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr.scale_)
        assert_array_almost_equal(scaler_batch.min_, scaler_incr.min_)

        # Test std after 1 step
        batch0 = slice(0, chunk_size)
        scaler_batch = MinMaxScaler().fit(X[batch0])
        scaler_incr = MinMaxScaler().partial_fit(X[batch0])

        assert_array_almost_equal(scaler_batch.data_min_,
                                  scaler_incr.data_min_)
        assert_array_almost_equal(scaler_batch.data_max_,
                                  scaler_incr.data_max_)
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_
        assert_array_almost_equal(scaler_batch.data_range_,
                                  scaler_incr.data_range_)
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr.scale_)
        assert_array_almost_equal(scaler_batch.min_, scaler_incr.min_)

        # Test std until the end of partial fits, and
        scaler_batch = MinMaxScaler().fit(X)
        scaler_incr = MinMaxScaler()  # Clean estimator
        for i, batch in enumerate(gen_batches(n_samples, chunk_size)):
            scaler_incr = scaler_incr.partial_fit(X[batch])
            assert_correct_incr(i, batch_start=batch.start,
                                batch_stop=batch.stop, n=n,
                                chunk_size=chunk_size,
                                n_samples_seen=scaler_incr.n_samples_seen_)


def test_standard_scaler_partial_fit():
    # Test if partial_fit run over many batches of size 1 and 50
    # gives the same results as fit
    X = X_2d
    n = X.shape[0]

    for chunk_size in [1, 2, 50, n, n + 42]:
        # Test mean at the end of the process
        scaler_batch = StandardScaler(with_std=False).fit(X)

        scaler_incr = StandardScaler(with_std=False)
        for batch in gen_batches(n_samples, chunk_size):
            scaler_incr = scaler_incr.partial_fit(X[batch])

        assert_array_almost_equal(scaler_batch.mean_, scaler_incr.mean_)
        assert scaler_batch.var_ == scaler_incr.var_  # Nones
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_

        # Test std after 1 step
        batch0 = slice(0, chunk_size)
        scaler_incr = StandardScaler().partial_fit(X[batch0])
        if chunk_size == 1:
            assert_array_almost_equal(np.zeros(n_features, dtype=np.float64),
                                      scaler_incr.var_)
            assert_array_almost_equal(np.ones(n_features, dtype=np.float64),
                                      scaler_incr.scale_)
        else:
            assert_array_almost_equal(np.var(X[batch0], axis=0),
                                      scaler_incr.var_)
            assert_array_almost_equal(np.std(X[batch0], axis=0),
                                      scaler_incr.scale_)  # no constants

        # Test std until the end of partial fits, and
        scaler_batch = StandardScaler().fit(X)
        scaler_incr = StandardScaler()  # Clean estimator
        for i, batch in enumerate(gen_batches(n_samples, chunk_size)):
            scaler_incr = scaler_incr.partial_fit(X[batch])
            assert_correct_incr(i, batch_start=batch.start,
                                batch_stop=batch.stop, n=n,
                                chunk_size=chunk_size,
                                n_samples_seen=scaler_incr.n_samples_seen_)

        assert_array_almost_equal(scaler_batch.var_, scaler_incr.var_)
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_


def test_standard_scaler_partial_fit_numerical_stability():
    # Test if the incremental computation introduces significative errors
    # for large datasets with values of large magniture
    rng = np.random.RandomState(0)
    n_features = 2
    n_samples = 100
    offsets = rng.uniform(-1e15, 1e15, size=n_features)
    scales = rng.uniform(1e3, 1e6, size=n_features)
    X = rng.randn(n_samples, n_features) * scales + offsets

    scaler_batch = StandardScaler().fit(X)
    scaler_incr = StandardScaler()
    for chunk in X:
        scaler_incr = scaler_incr.partial_fit(chunk.reshape(1, n_features))

    # Regardless of abs values, they must not be more diff 6 significant digits
    tol = 10 ** (-6)
    assert_allclose(scaler_incr.mean_, scaler_batch.mean_, rtol=tol)
    assert_allclose(scaler_incr.var_, scaler_batch.var_, rtol=tol)
    assert_allclose(scaler_incr.scale_, scaler_batch.scale_, rtol=tol)
    # NOTE Be aware that for much larger offsets std is very unstable (last
    # assert) while mean is OK.

    # Sparse input
    size = (100, 3)
    scale = 1e20
    X = rng.randint(0, 2, size).astype(np.float64) * scale
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    for X in [X_csr, X_csc]:
        # with_mean=False is required with sparse input
        scaler = StandardScaler(with_mean=False).fit(X)
        scaler_incr = StandardScaler(with_mean=False)

        for chunk in X:
            # chunk = sparse.csr_matrix(data_chunks)
            scaler_incr = scaler_incr.partial_fit(chunk)

        # Regardless of magnitude, they must not differ more than of 6 digits
        tol = 10 ** (-6)
        assert scaler.mean_ is not None
        assert_allclose(scaler_incr.var_, scaler.var_, rtol=tol)
        assert_allclose(scaler_incr.scale_, scaler.scale_, rtol=tol)


def test_partial_fit_sparse_input():
    # Check that sparsity is not destroyed
    X = np.array([[1.], [0.], [0.], [5.]])
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    null_transform = StandardScaler(with_mean=False, with_std=False, copy=True)
    for X in [X_csr, X_csc]:

        X_null = null_transform.partial_fit(X).transform(X)
        assert_array_equal(X_null.data, X.data)
        X_orig = null_transform.inverse_transform(X_null)
        assert_array_equal(X_orig.data, X_null.data)
        assert_array_equal(X_orig.data, X.data)


def test_standard_scaler_trasform_with_partial_fit():
    # Check some postconditions after applying partial_fit and transform
    X = X_2d[:100, :]

    scaler_incr = StandardScaler()
    for i, batch in enumerate(gen_batches(X.shape[0], 1)):

        X_sofar = X[:(i + 1), :]
        chunks_copy = X_sofar.copy()
        scaled_batch = StandardScaler().fit_transform(X_sofar)

        scaler_incr = scaler_incr.partial_fit(X[batch])
        scaled_incr = scaler_incr.transform(X_sofar)

        assert_array_almost_equal(scaled_batch, scaled_incr)
        assert_array_almost_equal(X_sofar, chunks_copy)  # No change
        right_input = scaler_incr.inverse_transform(scaled_incr)
        assert_array_almost_equal(X_sofar, right_input)

        zero = np.zeros(X.shape[1])
        epsilon = np.finfo(float).eps
        assert_array_less(zero, scaler_incr.var_ + epsilon)  # as less or equal
        assert_array_less(zero, scaler_incr.scale_ + epsilon)
        # (i+1) because the Scaler has been already fitted
        assert (i + 1) == scaler_incr.n_samples_seen_


def test_min_max_scaler_iris():
    X = iris.data
    scaler = MinMaxScaler()
    # default params
    X_trans = scaler.fit_transform(X)
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
    with pytest.raises(ValueError):
        scaler.fit(X)


def test_min_max_scaler_zero_variance_features():
    # Check min max scaler on toy data with zero variance features
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

    # function interface
    X_trans = minmax_scale(X)
    assert_array_almost_equal(X_trans, X_expected_0_1)
    X_trans = minmax_scale(X, feature_range=(1, 2))
    assert_array_almost_equal(X_trans, X_expected_1_2)


def test_minmax_scale_axis1():
    X = iris.data
    X_trans = minmax_scale(X, axis=1)
    assert_array_almost_equal(np.min(X_trans, axis=1), 0)
    assert_array_almost_equal(np.max(X_trans, axis=1), 1)


def test_min_max_scaler_1d():
    # Test scaling of dataset along single axis
    for X in [X_1row, X_1col, X_list_1row, X_list_1row]:

        scaler = MinMaxScaler(copy=True)
        X_scaled = scaler.fit(X).transform(X)

        if isinstance(X, list):
            X = np.array(X)  # cast only after scaling done

        if _check_dim_1axis(X) == 1:
            assert_array_almost_equal(X_scaled.min(axis=0),
                                      np.zeros(n_features))
            assert_array_almost_equal(X_scaled.max(axis=0),
                                      np.zeros(n_features))
        else:
            assert_array_almost_equal(X_scaled.min(axis=0), .0)
            assert_array_almost_equal(X_scaled.max(axis=0), 1.)
        assert scaler.n_samples_seen_ == X.shape[0]

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones((5, 1))
    scaler = MinMaxScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert X_scaled.min() >= 0.
    assert X_scaled.max() <= 1.
    assert scaler.n_samples_seen_ == X.shape[0]

    # Function interface
    X_1d = X_1row.ravel()
    min_ = X_1d.min()
    max_ = X_1d.max()
    assert_array_almost_equal((X_1d - min_) / (max_ - min_),
                              minmax_scale(X_1d, copy=True))


def test_scaler_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    with pytest.raises(ValueError):
        StandardScaler().fit(X_csr)
    with pytest.raises(ValueError):
        StandardScaler().fit(X_csc)

    null_transform = StandardScaler(with_mean=False, with_std=False, copy=True)
    X_null = null_transform.fit_transform(X_csr)
    assert_array_equal(X_null.data, X_csr.data)
    X_orig = null_transform.inverse_transform(X_null)
    assert_array_equal(X_orig.data, X_csr.data)

    scaler = StandardScaler(with_mean=False).fit(X)
    X_scaled = scaler.transform(X, copy=True)
    assert not np.any(np.isnan(X_scaled))

    scaler_csr = StandardScaler(with_mean=False).fit(X_csr)
    X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert not np.any(np.isnan(X_csr_scaled.data))

    scaler_csc = StandardScaler(with_mean=False).fit(X_csc)
    X_csc_scaled = scaler_csc.transform(X_csc, copy=True)
    assert not np.any(np.isnan(X_csc_scaled.data))

    assert_array_almost_equal(scaler.mean_, scaler_csr.mean_)
    assert_array_almost_equal(scaler.var_, scaler_csr.var_)
    assert_array_almost_equal(scaler.scale_, scaler_csr.scale_)

    assert_array_almost_equal(scaler.mean_, scaler_csc.mean_)
    assert_array_almost_equal(scaler.var_, scaler_csc.var_)
    assert_array_almost_equal(scaler.scale_, scaler_csc.scale_)

    assert_array_almost_equal(
        X_scaled.mean(axis=0), [0., -0.01, 2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(X_csr_scaled, 0)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))

    # Check that X has not been modified (copy)
    assert X_scaled is not X
    assert X_csr_scaled is not X_csr

    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert X_scaled_back is not X
    assert X_scaled_back is not X_scaled
    assert_array_almost_equal(X_scaled_back, X)

    X_csr_scaled_back = scaler_csr.inverse_transform(X_csr_scaled)
    assert X_csr_scaled_back is not X_csr
    assert X_csr_scaled_back is not X_csr_scaled
    assert_array_almost_equal(X_csr_scaled_back.toarray(), X)

    X_csc_scaled_back = scaler_csr.inverse_transform(X_csc_scaled.tocsc())
    assert X_csc_scaled_back is not X_csc
    assert X_csc_scaled_back is not X_csc_scaled
    assert_array_almost_equal(X_csc_scaled_back.toarray(), X)


@pytest.mark.parametrize("with_mean", [True, False])
@pytest.mark.parametrize("with_std", [True, False])
@pytest.mark.parametrize("array_constructor",
                         [np.asarray, sparse.csc_matrix, sparse.csr_matrix])
def test_scaler_n_samples_seen_with_nan(with_mean, with_std,
                                        array_constructor):
    X = np.array([[0, 1, 3],
                  [np.nan, 6, 10],
                  [5, 4, np.nan],
                  [8, 0, np.nan]],
                 dtype=np.float64)
    X = array_constructor(X)

    if sparse.issparse(X) and with_mean:
        pytest.skip("'with_mean=True' cannot be used with sparse matrix.")

    transformer = StandardScaler(with_mean=with_mean, with_std=with_std)
    transformer.fit(X)

    assert_array_equal(transformer.n_samples_seen_, np.array([3, 4, 2]))


def _check_identity_scalers_attributes(scaler_1, scaler_2):
    assert scaler_1.mean_ is scaler_2.mean_ is None
    assert scaler_1.var_ is scaler_2.var_ is None
    assert scaler_1.scale_ is scaler_2.scale_ is None
    assert scaler_1.n_samples_seen_ == scaler_2.n_samples_seen_


def test_scaler_return_identity():
    # test that the scaler return identity when with_mean and with_std are
    # False
    X_dense = np.array([[0, 1, 3],
                        [5, 6, 0],
                        [8, 0, 10]],
                       dtype=np.float64)
    X_csr = sparse.csr_matrix(X_dense)
    X_csc = X_csr.tocsc()

    transformer_dense = StandardScaler(with_mean=False, with_std=False)
    X_trans_dense = transformer_dense.fit_transform(X_dense)

    transformer_csr = clone(transformer_dense)
    X_trans_csr = transformer_csr.fit_transform(X_csr)

    transformer_csc = clone(transformer_dense)
    X_trans_csc = transformer_csc.fit_transform(X_csc)

    assert_allclose_dense_sparse(X_trans_csr, X_csr)
    assert_allclose_dense_sparse(X_trans_csc, X_csc)
    assert_allclose(X_trans_dense, X_dense)

    for trans_1, trans_2 in itertools.combinations([transformer_dense,
                                                    transformer_csr,
                                                    transformer_csc],
                                                   2):
        _check_identity_scalers_attributes(trans_1, trans_2)

    transformer_dense.partial_fit(X_dense)
    transformer_csr.partial_fit(X_csr)
    transformer_csc.partial_fit(X_csc)

    for trans_1, trans_2 in itertools.combinations([transformer_dense,
                                                    transformer_csr,
                                                    transformer_csc],
                                                   2):
        _check_identity_scalers_attributes(trans_1, trans_2)

    transformer_dense.fit(X_dense)
    transformer_csr.fit(X_csr)
    transformer_csc.fit(X_csc)

    for trans_1, trans_2 in itertools.combinations([transformer_dense,
                                                    transformer_csr,
                                                    transformer_csc],
                                                   2):
        _check_identity_scalers_attributes(trans_1, trans_2)


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
    assert not np.any(np.isnan(X_scaled))

    with warnings.catch_warnings(record=True):
        scaler_csr = StandardScaler(with_mean=False).fit(X_csr)
        X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert not np.any(np.isnan(X_csr_scaled.data))

    with warnings.catch_warnings(record=True):
        scaler_csc = StandardScaler(with_mean=False).fit(X_csc)
        X_csc_scaled = scaler_csc.transform(X_csc, copy=True)
    assert not np.any(np.isnan(X_csc_scaled.data))

    assert_array_almost_equal(scaler.mean_, scaler_csr.mean_)
    assert_array_almost_equal(scaler.var_, scaler_csr.var_)
    assert_array_almost_equal(scaler.scale_, scaler_csr.scale_)

    assert_array_almost_equal(scaler.mean_, scaler_csc.mean_)
    assert_array_almost_equal(scaler.var_, scaler_csc.var_)
    assert_array_almost_equal(scaler.scale_, scaler_csc.scale_)

    assert_array_almost_equal(
        X_scaled.mean(axis=0),
        [0., 1.109, 1.856, 21., 1.559], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(
        X_csr_scaled.astype(np.float), 0)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))

    # Check that X has not been modified (copy)
    assert X_scaled is not X
    assert X_csr_scaled is not X_csr

    X_scaled_back = scaler.inverse_transform(X_scaled)
    assert X_scaled_back is not X
    assert X_scaled_back is not X_scaled
    assert_array_almost_equal(X_scaled_back, X)

    X_csr_scaled_back = scaler_csr.inverse_transform(X_csr_scaled)
    assert X_csr_scaled_back is not X_csr
    assert X_csr_scaled_back is not X_csr_scaled
    assert_array_almost_equal(X_csr_scaled_back.toarray(), X)

    X_csc_scaled_back = scaler_csr.inverse_transform(X_csc_scaled.tocsc())
    assert X_csc_scaled_back is not X_csc
    assert X_csc_scaled_back is not X_csc_scaled
    assert_array_almost_equal(X_csc_scaled_back.toarray(), X)


def test_scaler_without_copy():
    # Check that StandardScaler.fit does not change input
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    X_copy = X.copy()
    StandardScaler(copy=False).fit(X)
    assert_array_equal(X, X_copy)

    X_csr_copy = X_csr.copy()
    StandardScaler(with_mean=False, copy=False).fit(X_csr)
    assert_array_equal(X_csr.toarray(), X_csr_copy.toarray())

    X_csc_copy = X_csc.copy()
    StandardScaler(with_mean=False, copy=False).fit(X_csc)
    assert_array_equal(X_csc.toarray(), X_csc_copy.toarray())


def test_scale_sparse_with_mean_raise_exception():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)

    # check scaling and fit with direct calls on sparse data
    with pytest.raises(ValueError):
        scale(X_csr, with_mean=True)
    with pytest.raises(ValueError):
        StandardScaler(with_mean=True).fit(X_csr)

    with pytest.raises(ValueError):
        scale(X_csc, with_mean=True)
    with pytest.raises(ValueError):
        StandardScaler(with_mean=True).fit(X_csc)

    # check transform and inverse_transform after a fit on a dense array
    scaler = StandardScaler(with_mean=True).fit(X)
    with pytest.raises(ValueError):
        scaler.transform(X_csr)
    with pytest.raises(ValueError):
        scaler.transform(X_csc)

    X_transformed_csr = sparse.csr_matrix(scaler.transform(X))
    with pytest.raises(ValueError):
        scaler.inverse_transform(X_transformed_csr)

    X_transformed_csc = sparse.csc_matrix(scaler.transform(X))
    with pytest.raises(ValueError):
        scaler.inverse_transform(X_transformed_csc)


def test_scale_input_finiteness_validation():
    # Check if non finite inputs raise ValueError
    X = [[np.inf, 5, 6, 7, 8]]
    with pytest.raises(ValueError, match="Input contains infinity "
                       "or a value too large"):
        scale(X)


def test_robust_scaler_error_sparse():
    X_sparse = sparse.rand(1000, 10)
    scaler = RobustScaler(with_centering=True)
    err_msg = "Cannot center sparse matrices"
    with pytest.raises(ValueError, match=err_msg):
        scaler.fit(X_sparse)


@pytest.mark.parametrize("with_centering", [True, False])
@pytest.mark.parametrize("with_scaling", [True, False])
@pytest.mark.parametrize("X", [np.random.randn(10, 3),
                               sparse.rand(10, 3, density=0.5)])
def test_robust_scaler_attributes(X, with_centering, with_scaling):
    # check consistent type of attributes
    if with_centering and sparse.issparse(X):
        pytest.skip("RobustScaler cannot center sparse matrix")

    scaler = RobustScaler(with_centering=with_centering,
                          with_scaling=with_scaling)
    scaler.fit(X)

    if with_centering:
        assert isinstance(scaler.center_, np.ndarray)
    else:
        assert scaler.center_ is None
    if with_scaling:
        assert isinstance(scaler.scale_, np.ndarray)
    else:
        assert scaler.scale_ is None


def test_robust_scaler_col_zero_sparse():
    # check that the scaler is working when there is not data materialized in a
    # column of a sparse matrix
    X = np.random.randn(10, 5)
    X[:, 0] = 0
    X = sparse.csr_matrix(X)

    scaler = RobustScaler(with_centering=False)
    scaler.fit(X)
    assert scaler.scale_[0] == pytest.approx(1)

    X_trans = scaler.transform(X)
    assert_allclose(X[:, 0].toarray(), X_trans[:, 0].toarray())


def test_robust_scaler_2d_arrays():
    # Test robust scaling of 2d array along first axis
    rng = np.random.RandomState(0)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero

    scaler = RobustScaler()
    X_scaled = scaler.fit(X).transform(X)

    assert_array_almost_equal(np.median(X_scaled, axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0)[0], 0)


@pytest.mark.parametrize("density", [0, 0.05, 0.1, 0.5, 1])
@pytest.mark.parametrize("strictly_signed",
                         ['positive', 'negative', 'zeros', None])
def test_robust_scaler_equivalence_dense_sparse(density, strictly_signed):
    # Check the equivalence of the fitting with dense and sparse matrices
    X_sparse = sparse.rand(1000, 5, density=density).tocsc()
    if strictly_signed == 'positive':
        X_sparse.data = np.abs(X_sparse.data)
    elif strictly_signed == 'negative':
        X_sparse.data = - np.abs(X_sparse.data)
    elif strictly_signed == 'zeros':
        X_sparse.data = np.zeros(X_sparse.data.shape, dtype=np.float64)
    X_dense = X_sparse.toarray()

    scaler_sparse = RobustScaler(with_centering=False)
    scaler_dense = RobustScaler(with_centering=False)

    scaler_sparse.fit(X_sparse)
    scaler_dense.fit(X_dense)

    assert_allclose(scaler_sparse.scale_, scaler_dense.scale_)


def test_robust_scaler_transform_one_row_csr():
    # Check RobustScaler on transforming csr matrix with one row
    rng = np.random.RandomState(0)
    X = rng.randn(4, 5)
    single_row = np.array([[0.1, 1., 2., 0., -1.]])
    scaler = RobustScaler(with_centering=False)
    scaler = scaler.fit(X)
    row_trans = scaler.transform(sparse.csr_matrix(single_row))
    row_expected = single_row / scaler.scale_
    assert_array_almost_equal(row_trans.toarray(), row_expected)
    row_scaled_back = scaler.inverse_transform(row_trans)
    assert_array_almost_equal(single_row, row_scaled_back.toarray())


def test_robust_scaler_iris():
    X = iris.data
    scaler = RobustScaler()
    X_trans = scaler.fit_transform(X)
    assert_array_almost_equal(np.median(X_trans, axis=0), 0)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)
    q = np.percentile(X_trans, q=(25, 75), axis=0)
    iqr = q[1] - q[0]
    assert_array_almost_equal(iqr, 1)


def test_robust_scaler_iris_quantiles():
    X = iris.data
    scaler = RobustScaler(quantile_range=(10, 90))
    X_trans = scaler.fit_transform(X)
    assert_array_almost_equal(np.median(X_trans, axis=0), 0)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)
    q = np.percentile(X_trans, q=(10, 90), axis=0)
    q_range = q[1] - q[0]
    assert_array_almost_equal(q_range, 1)


def test_quantile_transform_iris():
    X = iris.data
    # uniform output distribution
    transformer = QuantileTransformer(n_quantiles=30)
    X_trans = transformer.fit_transform(X)
    X_trans_inv = transformer.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)
    # normal output distribution
    transformer = QuantileTransformer(n_quantiles=30,
                                      output_distribution='normal')
    X_trans = transformer.fit_transform(X)
    X_trans_inv = transformer.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)
    # make sure it is possible to take the inverse of a sparse matrix
    # which contain negative value; this is the case in the iris dataset
    X_sparse = sparse.csc_matrix(X)
    X_sparse_tran = transformer.fit_transform(X_sparse)
    X_sparse_tran_inv = transformer.inverse_transform(X_sparse_tran)
    assert_array_almost_equal(X_sparse.A, X_sparse_tran_inv.A)


def test_quantile_transform_check_error():
    X = np.transpose([[0, 25, 50, 0, 0, 0, 75, 0, 0, 100],
                      [2, 4, 0, 0, 6, 8, 0, 10, 0, 0],
                      [0, 0, 2.6, 4.1, 0, 0, 2.3, 0, 9.5, 0.1]])
    X = sparse.csc_matrix(X)
    X_neg = np.transpose([[0, 25, 50, 0, 0, 0, 75, 0, 0, 100],
                          [-2, 4, 0, 0, 6, 8, 0, 10, 0, 0],
                          [0, 0, 2.6, 4.1, 0, 0, 2.3, 0, 9.5, 0.1]])
    X_neg = sparse.csc_matrix(X_neg)

    err_msg = "Invalid value for 'n_quantiles': 0."
    with pytest.raises(ValueError, match=err_msg):
        QuantileTransformer(n_quantiles=0).fit(X)
    err_msg = "Invalid value for 'subsample': 0."
    with pytest.raises(ValueError, match=err_msg):
        QuantileTransformer(subsample=0).fit(X)
    err_msg = ("The number of quantiles cannot be greater than "
               "the number of samples used. Got 1000 quantiles "
               "and 10 samples.")
    with pytest.raises(ValueError, match=err_msg):
        QuantileTransformer(subsample=10).fit(X)

    transformer = QuantileTransformer(n_quantiles=10)
    err_msg = "QuantileTransformer only accepts non-negative sparse matrices."
    with pytest.raises(ValueError, match=err_msg):
        transformer.fit(X_neg)
    transformer.fit(X)
    err_msg = "QuantileTransformer only accepts non-negative sparse matrices."
    with pytest.raises(ValueError, match=err_msg):
        transformer.transform(X_neg)

    X_bad_feat = np.transpose([[0, 25, 50, 0, 0, 0, 75, 0, 0, 100],
                               [0, 0, 2.6, 4.1, 0, 0, 2.3, 0, 9.5, 0.1]])
    err_msg = ("X does not have the same number of features as the previously"
               " fitted " "data. Got 2 instead of 3.")
    with pytest.raises(ValueError, match=err_msg):
        transformer.transform(X_bad_feat)
    err_msg = ("X does not have the same number of features "
               "as the previously fitted data. Got 2 instead of 3.")
    with pytest.raises(ValueError, match=err_msg):
        transformer.inverse_transform(X_bad_feat)

    transformer = QuantileTransformer(n_quantiles=10,
                                      output_distribution='rnd')
    # check that an error is raised at fit time
    err_msg = ("'output_distribution' has to be either 'normal' or "
               "'uniform'. Got 'rnd' instead.")
    with pytest.raises(ValueError, match=err_msg):
        transformer.fit(X)
    # check that an error is raised at transform time
    transformer.output_distribution = 'uniform'
    transformer.fit(X)
    X_tran = transformer.transform(X)
    transformer.output_distribution = 'rnd'
    err_msg = ("'output_distribution' has to be either 'normal' or 'uniform'."
               " Got 'rnd' instead.")
    with pytest.raises(ValueError, match=err_msg):
        transformer.transform(X)
    # check that an error is raised at inverse_transform time
    err_msg = ("'output_distribution' has to be either 'normal' or 'uniform'."
               " Got 'rnd' instead.")
    with pytest.raises(ValueError, match=err_msg):
        transformer.inverse_transform(X_tran)
    # check that an error is raised if input is scalar
    with pytest.raises(ValueError,
                       match='Expected 2D array, got scalar array instead'):
        transformer.transform(10)
    # check that a warning is raised is n_quantiles > n_samples
    transformer = QuantileTransformer(n_quantiles=100)
    warn_msg = "n_quantiles is set to n_samples"
    with pytest.warns(UserWarning, match=warn_msg) as record:
        transformer.fit(X)
    assert len(record) == 1
    assert transformer.n_quantiles_ == X.shape[0]


def test_quantile_transform_sparse_ignore_zeros():
    X = np.array([[0, 1],
                  [0, 0],
                  [0, 2],
                  [0, 2],
                  [0, 1]])
    X_sparse = sparse.csc_matrix(X)
    transformer = QuantileTransformer(ignore_implicit_zeros=True,
                                      n_quantiles=5)

    # dense case -> warning raise
    assert_warns_message(UserWarning, "'ignore_implicit_zeros' takes effect"
                         " only with sparse matrix. This parameter has no"
                         " effect.", transformer.fit, X)

    X_expected = np.array([[0, 0],
                           [0, 0],
                           [0, 1],
                           [0, 1],
                           [0, 0]])
    X_trans = transformer.fit_transform(X_sparse)
    assert_almost_equal(X_expected, X_trans.A)

    # consider the case where sparse entries are missing values and user-given
    # zeros are to be considered
    X_data = np.array([0, 0, 1, 0, 2, 2, 1, 0, 1, 2, 0])
    X_col = np.array([0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    X_row = np.array([0, 4, 0, 1, 2, 3, 4, 5, 6, 7, 8])
    X_sparse = sparse.csc_matrix((X_data, (X_row, X_col)))
    X_trans = transformer.fit_transform(X_sparse)
    X_expected = np.array([[0., 0.5],
                           [0., 0.],
                           [0., 1.],
                           [0., 1.],
                           [0., 0.5],
                           [0., 0.],
                           [0., 0.5],
                           [0., 1.],
                           [0., 0.]])
    assert_almost_equal(X_expected, X_trans.A)

    transformer = QuantileTransformer(ignore_implicit_zeros=True,
                                      n_quantiles=5)
    X_data = np.array([-1, -1, 1, 0, 0, 0, 1, -1, 1])
    X_col = np.array([0, 0, 1, 1, 1, 1, 1, 1, 1])
    X_row = np.array([0, 4, 0, 1, 2, 3, 4, 5, 6])
    X_sparse = sparse.csc_matrix((X_data, (X_row, X_col)))
    X_trans = transformer.fit_transform(X_sparse)
    X_expected = np.array([[0, 1],
                           [0, 0.375],
                           [0, 0.375],
                           [0, 0.375],
                           [0, 1],
                           [0, 0],
                           [0, 1]])
    assert_almost_equal(X_expected, X_trans.A)
    assert_almost_equal(X_sparse.A, transformer.inverse_transform(X_trans).A)

    # check in conjunction with subsampling
    transformer = QuantileTransformer(ignore_implicit_zeros=True,
                                      n_quantiles=5,
                                      subsample=8,
                                      random_state=0)
    X_trans = transformer.fit_transform(X_sparse)
    assert_almost_equal(X_expected, X_trans.A)
    assert_almost_equal(X_sparse.A, transformer.inverse_transform(X_trans).A)


def test_quantile_transform_dense_toy():
    X = np.array([[0, 2, 2.6],
                  [25, 4, 4.1],
                  [50, 6, 2.3],
                  [75, 8, 9.5],
                  [100, 10, 0.1]])

    transformer = QuantileTransformer(n_quantiles=5)
    transformer.fit(X)

    # using the a uniform output, each entry of X should be map between 0 and 1
    # and equally spaced
    X_trans = transformer.fit_transform(X)
    X_expected = np.tile(np.linspace(0, 1, num=5), (3, 1)).T
    assert_almost_equal(np.sort(X_trans, axis=0), X_expected)

    X_test = np.array([
        [-1, 1, 0],
        [101, 11, 10],
    ])
    X_expected = np.array([
        [0, 0, 0],
        [1, 1, 1],
    ])
    assert_array_almost_equal(transformer.transform(X_test), X_expected)

    X_trans_inv = transformer.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)


def test_quantile_transform_subsampling():
    # Test that subsampling the input yield to a consistent results We check
    # that the computed quantiles are almost mapped to a [0, 1] vector where
    # values are equally spaced. The infinite norm is checked to be smaller
    # than a given threshold. This is repeated 5 times.

    # dense support
    n_samples = 1000000
    n_quantiles = 1000
    X = np.sort(np.random.sample((n_samples, 1)), axis=0)
    ROUND = 5
    inf_norm_arr = []
    for random_state in range(ROUND):
        transformer = QuantileTransformer(random_state=random_state,
                                          n_quantiles=n_quantiles,
                                          subsample=n_samples // 10)
        transformer.fit(X)
        diff = (np.linspace(0, 1, n_quantiles) -
                np.ravel(transformer.quantiles_))
        inf_norm = np.max(np.abs(diff))
        assert inf_norm < 1e-2
        inf_norm_arr.append(inf_norm)
    # each random subsampling yield a unique approximation to the expected
    # linspace CDF
    assert len(np.unique(inf_norm_arr)) == len(inf_norm_arr)

    # sparse support

    X = sparse.rand(n_samples, 1, density=.99, format='csc', random_state=0)
    inf_norm_arr = []
    for random_state in range(ROUND):
        transformer = QuantileTransformer(random_state=random_state,
                                          n_quantiles=n_quantiles,
                                          subsample=n_samples // 10)
        transformer.fit(X)
        diff = (np.linspace(0, 1, n_quantiles) -
                np.ravel(transformer.quantiles_))
        inf_norm = np.max(np.abs(diff))
        assert inf_norm < 1e-1
        inf_norm_arr.append(inf_norm)
    # each random subsampling yield a unique approximation to the expected
    # linspace CDF
    assert len(np.unique(inf_norm_arr)) == len(inf_norm_arr)


def test_quantile_transform_sparse_toy():
    X = np.array([[0., 2., 0.],
                  [25., 4., 0.],
                  [50., 0., 2.6],
                  [0., 0., 4.1],
                  [0., 6., 0.],
                  [0., 8., 0.],
                  [75., 0., 2.3],
                  [0., 10., 0.],
                  [0., 0., 9.5],
                  [100., 0., 0.1]])

    X = sparse.csc_matrix(X)

    transformer = QuantileTransformer(n_quantiles=10)
    transformer.fit(X)

    X_trans = transformer.fit_transform(X)
    assert_array_almost_equal(np.min(X_trans.toarray(), axis=0), 0.)
    assert_array_almost_equal(np.max(X_trans.toarray(), axis=0), 1.)

    X_trans_inv = transformer.inverse_transform(X_trans)
    assert_array_almost_equal(X.toarray(), X_trans_inv.toarray())

    transformer_dense = QuantileTransformer(n_quantiles=10).fit(
        X.toarray())

    X_trans = transformer_dense.transform(X)
    assert_array_almost_equal(np.min(X_trans.toarray(), axis=0), 0.)
    assert_array_almost_equal(np.max(X_trans.toarray(), axis=0), 1.)

    X_trans_inv = transformer_dense.inverse_transform(X_trans)
    assert_array_almost_equal(X.toarray(), X_trans_inv.toarray())


@pytest.mark.filterwarnings("ignore: The default value of `copy`")  # 0.23
def test_quantile_transform_axis1():
    X = np.array([[0, 25, 50, 75, 100],
                  [2, 4, 6, 8, 10],
                  [2.6, 4.1, 2.3, 9.5, 0.1]])

    X_trans_a0 = quantile_transform(X.T, axis=0, n_quantiles=5)
    X_trans_a1 = quantile_transform(X, axis=1, n_quantiles=5)
    assert_array_almost_equal(X_trans_a0, X_trans_a1.T)


def test_quantile_transform_bounds():
    # Lower and upper bounds are manually mapped. We checked that in the case
    # of a constant feature and binary feature, the bounds are properly mapped.
    X_dense = np.array([[0, 0],
                        [0, 0],
                        [1, 0]])
    X_sparse = sparse.csc_matrix(X_dense)

    # check sparse and dense are consistent
    X_trans = QuantileTransformer(n_quantiles=3,
                                  random_state=0).fit_transform(X_dense)
    assert_array_almost_equal(X_trans, X_dense)
    X_trans_sp = QuantileTransformer(n_quantiles=3,
                                     random_state=0).fit_transform(X_sparse)
    assert_array_almost_equal(X_trans_sp.A, X_dense)
    assert_array_almost_equal(X_trans, X_trans_sp.A)

    # check the consistency of the bounds by learning on 1 matrix
    # and transforming another
    X = np.array([[0, 1],
                  [0, 0.5],
                  [1, 0]])
    X1 = np.array([[0, 0.1],
                   [0, 0.5],
                   [1, 0.1]])
    transformer = QuantileTransformer(n_quantiles=3).fit(X)
    X_trans = transformer.transform(X1)
    assert_array_almost_equal(X_trans, X1)

    # check that values outside of the range learned will be mapped properly.
    X = np.random.random((1000, 1))
    transformer = QuantileTransformer()
    transformer.fit(X)
    assert (transformer.transform([[-10]]) ==
                 transformer.transform([[np.min(X)]]))
    assert (transformer.transform([[10]]) ==
                 transformer.transform([[np.max(X)]]))
    assert (transformer.inverse_transform([[-10]]) ==
                 transformer.inverse_transform(
                     [[np.min(transformer.references_)]]))
    assert (transformer.inverse_transform([[10]]) ==
                 transformer.inverse_transform(
                     [[np.max(transformer.references_)]]))


def test_quantile_transform_and_inverse():
    X_1 = iris.data
    X_2 = np.array([[0.], [BOUNDS_THRESHOLD / 10], [1.5], [2], [3], [3], [4]])
    for X in [X_1, X_2]:
        transformer = QuantileTransformer(n_quantiles=1000, random_state=0)
        X_trans = transformer.fit_transform(X)
        X_trans_inv = transformer.inverse_transform(X_trans)
        assert_array_almost_equal(X, X_trans_inv, decimal=9)


def test_quantile_transform_nan():
    X = np.array([[np.nan, 0,  0, 1],
                  [np.nan, np.nan, 0, 0.5],
                  [np.nan, 1, 1, 0]])

    transformer = QuantileTransformer(n_quantiles=10, random_state=42)
    transformer.fit_transform(X)

    # check that the quantile of the first column is all NaN
    assert np.isnan(transformer.quantiles_[:, 0]).all()
    # all other column should not contain NaN
    assert not np.isnan(transformer.quantiles_[:, 1:]).any()


def test_deprecated_quantile_transform_copy():
    future_message = ("The default value of `copy` will change from False to "
                      "True in 0.23 in order to make it more consistent with "
                      "the default `copy` values of other functions in "
                      ":mod:`sklearn.preprocessing` and prevent "
                      "unexpected side effects by modifying the value of `X` "
                      "inplace. To avoid inplace modifications of `X`, it is "
                      "recommended to explicitly set `copy=True`")
    assert_warns_message(FutureWarning, future_message, quantile_transform,
                         np.array([[0, 1], [0, 0.5], [1, 0]]))


@pytest.mark.parametrize("array_type", ['array', 'sparse'])
def test_quantile_transformer_sorted_quantiles(array_type):
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/15733
    # Taken from upstream bug report:
    # https://github.com/numpy/numpy/issues/14685
    X = np.array([0, 1, 1, 2, 2, 3, 3, 4, 5, 5, 1, 1, 9, 9, 9, 8, 8, 7] * 10)
    X = 0.1 * X.reshape(-1, 1)
    X = _convert_container(X, array_type)

    n_quantiles = 100
    qt = QuantileTransformer(n_quantiles=n_quantiles).fit(X)

    # Check that the estimated quantile threasholds are monotically
    # increasing:
    quantiles = qt.quantiles_[:, 0]
    assert len(quantiles) == 100
    assert all(np.diff(quantiles) >= 0)


def test_robust_scaler_invalid_range():
    for range_ in [
        (-1, 90),
        (-2, -3),
        (10, 101),
        (100.5, 101),
        (90, 50),
    ]:
        scaler = RobustScaler(quantile_range=range_)

        with pytest.raises(ValueError, match=r'Invalid quantile range: \('):
            scaler.fit(iris.data)


def test_scale_function_without_centering():
    rng = np.random.RandomState(42)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero
    X_csr = sparse.csr_matrix(X)

    X_scaled = scale(X, with_mean=False)
    assert not np.any(np.isnan(X_scaled))

    X_csr_scaled = scale(X_csr, with_mean=False)
    assert not np.any(np.isnan(X_csr_scaled.data))

    # test csc has same outcome
    X_csc_scaled = scale(X_csr.tocsc(), with_mean=False)
    assert_array_almost_equal(X_scaled, X_csc_scaled.toarray())

    # raises value error on axis != 0
    with pytest.raises(ValueError):
        scale(X_csr, with_mean=False, axis=1)

    assert_array_almost_equal(X_scaled.mean(axis=0),
                              [0., -0.01, 2.24, -0.35, -0.78], 2)
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert X_scaled is not X

    X_csr_scaled_mean, X_csr_scaled_std = mean_variance_axis(X_csr_scaled, 0)
    assert_array_almost_equal(X_csr_scaled_mean, X_scaled.mean(axis=0))
    assert_array_almost_equal(X_csr_scaled_std, X_scaled.std(axis=0))

    # null scale
    X_csr_scaled = scale(X_csr, with_mean=False, with_std=False, copy=True)
    assert_array_almost_equal(X_csr.toarray(), X_csr_scaled.toarray())


def test_robust_scale_axis1():
    X = iris.data
    X_trans = robust_scale(X, axis=1)
    assert_array_almost_equal(np.median(X_trans, axis=1), 0)
    q = np.percentile(X_trans, q=(25, 75), axis=1)
    iqr = q[1] - q[0]
    assert_array_almost_equal(iqr, 1)


def test_robust_scale_1d_array():
    X = iris.data[:, 1]
    X_trans = robust_scale(X)
    assert_array_almost_equal(np.median(X_trans), 0)
    q = np.percentile(X_trans, q=(25, 75))
    iqr = q[1] - q[0]
    assert_array_almost_equal(iqr, 1)


def test_robust_scaler_zero_variance_features():
    # Check RobustScaler on toy data with zero variance features
    X = [[0., 1., +0.5],
         [0., 1., -0.1],
         [0., 1., +1.1]]

    scaler = RobustScaler()
    X_trans = scaler.fit_transform(X)

    # NOTE: for such a small sample size, what we expect in the third column
    # depends HEAVILY on the method used to calculate quantiles. The values
    # here were calculated to fit the quantiles produces by np.percentile
    # using numpy 1.9 Calculating quantiles with
    # scipy.stats.mstats.scoreatquantile or scipy.stats.mstats.mquantiles
    # would yield very different results!
    X_expected = [[0., 0., +0.0],
                  [0., 0., -1.0],
                  [0., 0., +1.0]]
    assert_array_almost_equal(X_trans, X_expected)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    # make sure new data gets transformed correctly
    X_new = [[+0., 2., 0.5],
             [-1., 1., 0.0],
             [+0., 1., 1.5]]
    X_trans_new = scaler.transform(X_new)
    X_expected_new = [[+0., 1., +0.],
                      [-1., 0., -0.83333],
                      [+0., 0., +1.66667]]
    assert_array_almost_equal(X_trans_new, X_expected_new, decimal=3)


def test_maxabs_scaler_zero_variance_features():
    # Check MaxAbsScaler on toy data with zero variance features
    X = [[0., 1., +0.5],
         [0., 1., -0.3],
         [0., 1., +1.5],
         [0., 0., +0.0]]

    scaler = MaxAbsScaler()
    X_trans = scaler.fit_transform(X)
    X_expected = [[0., 1., 1.0 / 3.0],
                  [0., 1., -0.2],
                  [0., 1., 1.0],
                  [0., 0., 0.0]]
    assert_array_almost_equal(X_trans, X_expected)
    X_trans_inv = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X, X_trans_inv)

    # make sure new data gets transformed correctly
    X_new = [[+0., 2., 0.5],
             [-1., 1., 0.0],
             [+0., 1., 1.5]]
    X_trans_new = scaler.transform(X_new)
    X_expected_new = [[+0., 2.0, 1.0 / 3.0],
                      [-1., 1.0, 0.0],
                      [+0., 1.0, 1.0]]

    assert_array_almost_equal(X_trans_new, X_expected_new, decimal=2)

    # function interface
    X_trans = maxabs_scale(X)
    assert_array_almost_equal(X_trans, X_expected)

    # sparse data
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)
    X_trans_csr = scaler.fit_transform(X_csr)
    X_trans_csc = scaler.fit_transform(X_csc)
    X_expected = [[0., 1., 1.0 / 3.0],
                  [0., 1., -0.2],
                  [0., 1., 1.0],
                  [0., 0., 0.0]]
    assert_array_almost_equal(X_trans_csr.A, X_expected)
    assert_array_almost_equal(X_trans_csc.A, X_expected)
    X_trans_csr_inv = scaler.inverse_transform(X_trans_csr)
    X_trans_csc_inv = scaler.inverse_transform(X_trans_csc)
    assert_array_almost_equal(X, X_trans_csr_inv.A)
    assert_array_almost_equal(X, X_trans_csc_inv.A)


def test_maxabs_scaler_large_negative_value():
    # Check MaxAbsScaler on toy data with a large negative value
    X = [[0., 1.,   +0.5, -1.0],
         [0., 1.,   -0.3, -0.5],
         [0., 1., -100.0,  0.0],
         [0., 0.,   +0.0, -2.0]]

    scaler = MaxAbsScaler()
    X_trans = scaler.fit_transform(X)
    X_expected = [[0., 1.,  0.005,    -0.5],
                  [0., 1., -0.003,    -0.25],
                  [0., 1., -1.0,       0.0],
                  [0., 0.,  0.0,      -1.0]]
    assert_array_almost_equal(X_trans, X_expected)


def test_maxabs_scaler_transform_one_row_csr():
    # Check MaxAbsScaler on transforming csr matrix with one row
    X = sparse.csr_matrix([[0.5, 1., 1.]])
    scaler = MaxAbsScaler()
    scaler = scaler.fit(X)
    X_trans = scaler.transform(X)
    X_expected = sparse.csr_matrix([[1., 1., 1.]])
    assert_array_almost_equal(X_trans.toarray(), X_expected.toarray())
    X_scaled_back = scaler.inverse_transform(X_trans)
    assert_array_almost_equal(X.toarray(), X_scaled_back.toarray())


def test_maxabs_scaler_1d():
    # Test scaling of dataset along single axis
    for X in [X_1row, X_1col, X_list_1row, X_list_1row]:

        scaler = MaxAbsScaler(copy=True)
        X_scaled = scaler.fit(X).transform(X)

        if isinstance(X, list):
            X = np.array(X)  # cast only after scaling done

        if _check_dim_1axis(X) == 1:
            assert_array_almost_equal(np.abs(X_scaled.max(axis=0)),
                                      np.ones(n_features))
        else:
            assert_array_almost_equal(np.abs(X_scaled.max(axis=0)), 1.)
        assert scaler.n_samples_seen_ == X.shape[0]

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones((5, 1))
    scaler = MaxAbsScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_array_almost_equal(np.abs(X_scaled.max(axis=0)), 1.)
    assert scaler.n_samples_seen_ == X.shape[0]

    # function interface
    X_1d = X_1row.ravel()
    max_abs = np.abs(X_1d).max()
    assert_array_almost_equal(X_1d / max_abs, maxabs_scale(X_1d, copy=True))


def test_maxabs_scaler_partial_fit():
    # Test if partial_fit run over many batches of size 1 and 50
    # gives the same results as fit
    X = X_2d[:100, :]
    n = X.shape[0]

    for chunk_size in [1, 2, 50, n, n + 42]:
        # Test mean at the end of the process
        scaler_batch = MaxAbsScaler().fit(X)

        scaler_incr = MaxAbsScaler()
        scaler_incr_csr = MaxAbsScaler()
        scaler_incr_csc = MaxAbsScaler()
        for batch in gen_batches(n, chunk_size):
            scaler_incr = scaler_incr.partial_fit(X[batch])
            X_csr = sparse.csr_matrix(X[batch])
            scaler_incr_csr = scaler_incr_csr.partial_fit(X_csr)
            X_csc = sparse.csc_matrix(X[batch])
            scaler_incr_csc = scaler_incr_csc.partial_fit(X_csc)

        assert_array_almost_equal(scaler_batch.max_abs_, scaler_incr.max_abs_)
        assert_array_almost_equal(scaler_batch.max_abs_,
                                  scaler_incr_csr.max_abs_)
        assert_array_almost_equal(scaler_batch.max_abs_,
                                  scaler_incr_csc.max_abs_)
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_
        assert (scaler_batch.n_samples_seen_ ==
                     scaler_incr_csr.n_samples_seen_)
        assert (scaler_batch.n_samples_seen_ ==
                     scaler_incr_csc.n_samples_seen_)
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr.scale_)
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr_csr.scale_)
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr_csc.scale_)
        assert_array_almost_equal(scaler_batch.transform(X),
                                  scaler_incr.transform(X))

        # Test std after 1 step
        batch0 = slice(0, chunk_size)
        scaler_batch = MaxAbsScaler().fit(X[batch0])
        scaler_incr = MaxAbsScaler().partial_fit(X[batch0])

        assert_array_almost_equal(scaler_batch.max_abs_, scaler_incr.max_abs_)
        assert scaler_batch.n_samples_seen_ == scaler_incr.n_samples_seen_
        assert_array_almost_equal(scaler_batch.scale_, scaler_incr.scale_)
        assert_array_almost_equal(scaler_batch.transform(X),
                                  scaler_incr.transform(X))

        # Test std until the end of partial fits, and
        scaler_batch = MaxAbsScaler().fit(X)
        scaler_incr = MaxAbsScaler()  # Clean estimator
        for i, batch in enumerate(gen_batches(n, chunk_size)):
            scaler_incr = scaler_incr.partial_fit(X[batch])
            assert_correct_incr(i, batch_start=batch.start,
                                batch_stop=batch.stop, n=n,
                                chunk_size=chunk_size,
                                n_samples_seen=scaler_incr.n_samples_seen_)


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
        assert X_norm is not X
        X_norm1 = toarray(X_norm)

        normalizer = Normalizer(norm='l1', copy=False)
        X_norm = normalizer.transform(X)
        assert X_norm is X
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

        assert X_norm is not X
        assert isinstance(X_norm, sparse.csr_matrix)

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
        assert X_norm1 is not X
        X_norm1 = toarray(X_norm1)

        normalizer = Normalizer(norm='l2', copy=False)
        X_norm2 = normalizer.transform(X)
        assert X_norm2 is X
        X_norm2 = toarray(X_norm2)

        for X_norm in (X_norm1, X_norm2):
            for i in range(3):
                assert_almost_equal(la.norm(X_norm[i]), 1.0)
            assert_almost_equal(la.norm(X_norm[3]), 0.0)

    # check input for which copy=False won't prevent a copy
    for init in (sparse.coo_matrix, sparse.csc_matrix, sparse.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert X_norm is not X
        assert isinstance(X_norm, sparse.csr_matrix)

        X_norm = toarray(X_norm)
        for i in range(3):
            assert_almost_equal(la.norm(X_norm[i]), 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalizer_max():
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

        normalizer = Normalizer(norm='max', copy=True)
        X_norm1 = normalizer.transform(X)
        assert X_norm1 is not X
        X_norm1 = toarray(X_norm1)

        normalizer = Normalizer(norm='max', copy=False)
        X_norm2 = normalizer.transform(X)
        assert X_norm2 is X
        X_norm2 = toarray(X_norm2)

        for X_norm in (X_norm1, X_norm2):
            row_maxs = X_norm.max(axis=1)
            for i in range(3):
                assert_almost_equal(row_maxs[i], 1.0)
            assert_almost_equal(row_maxs[3], 0.0)

    # check input for which copy=False won't prevent a copy
    for init in (sparse.coo_matrix, sparse.csc_matrix, sparse.lil_matrix):
        X = init(X_dense)
        X_norm = normalizer = Normalizer(norm='l2', copy=False).transform(X)

        assert X_norm is not X
        assert isinstance(X_norm, sparse.csr_matrix)

        X_norm = toarray(X_norm)
        for i in range(3):
            assert_almost_equal(row_maxs[i], 1.0)
        assert_almost_equal(la.norm(X_norm[3]), 0.0)


def test_normalize():
    # Test normalize function
    # Only tests functionality not used by the tests for Normalizer.
    X = np.random.RandomState(37).randn(3, 2)
    assert_array_equal(normalize(X, copy=False),
                       normalize(X.T, axis=0, copy=False).T)
    with pytest.raises(ValueError):
        normalize([[0]], axis=2)
    with pytest.raises(ValueError):
        normalize([[0]], norm='l3')

    rs = np.random.RandomState(0)
    X_dense = rs.randn(10, 5)
    X_sparse = sparse.csr_matrix(X_dense)
    ones = np.ones((10))
    for X in (X_dense, X_sparse):
        for dtype in (np.float32, np.float64):
            for norm in ('l1', 'l2'):
                X = X.astype(dtype)
                X_norm = normalize(X, norm=norm)
                assert X_norm.dtype == dtype

                X_norm = toarray(X_norm)
                if norm == 'l1':
                    row_sums = np.abs(X_norm).sum(axis=1)
                else:
                    X_norm_squared = X_norm**2
                    row_sums = X_norm_squared.sum(axis=1)

                assert_array_almost_equal(row_sums, ones)

    # Test return_norm
    X_dense = np.array([[3.0, 0, 4.0], [1.0, 0.0, 0.0], [2.0, 3.0, 0.0]])
    for norm in ('l1', 'l2', 'max'):
        _, norms = normalize(X_dense, norm=norm, return_norm=True)
        if norm == 'l1':
            assert_array_almost_equal(norms, np.array([7.0, 1.0, 5.0]))
        elif norm == 'l2':
            assert_array_almost_equal(norms, np.array([5.0, 1.0, 3.60555127]))
        else:
            assert_array_almost_equal(norms, np.array([4.0, 1.0, 3.0]))

    X_sparse = sparse.csr_matrix(X_dense)
    for norm in ('l1', 'l2'):
        with pytest.raises(NotImplementedError):
            normalize(X_sparse, norm=norm, return_norm=True)
    _, norms = normalize(X_sparse, norm='max', return_norm=True)
    assert_array_almost_equal(norms, np.array([4.0, 1.0, 3.0]))


def test_binarizer():
    X_ = np.array([[1, 0, 5], [2, 3, -1]])

    for init in (np.array, list, sparse.csr_matrix, sparse.csc_matrix):

        X = init(X_.copy())

        binarizer = Binarizer(threshold=2.0, copy=True)
        X_bin = toarray(binarizer.transform(X))
        assert np.sum(X_bin == 0) == 4
        assert np.sum(X_bin == 1) == 2
        X_bin = binarizer.transform(X)
        assert sparse.issparse(X) == sparse.issparse(X_bin)

        binarizer = Binarizer(copy=True).fit(X)
        X_bin = toarray(binarizer.transform(X))
        assert X_bin is not X
        assert np.sum(X_bin == 0) == 2
        assert np.sum(X_bin == 1) == 4

        binarizer = Binarizer(copy=True)
        X_bin = binarizer.transform(X)
        assert X_bin is not X
        X_bin = toarray(X_bin)
        assert np.sum(X_bin == 0) == 2
        assert np.sum(X_bin == 1) == 4

        binarizer = Binarizer(copy=False)
        X_bin = binarizer.transform(X)
        if init is not list:
            assert X_bin is X

        binarizer = Binarizer(copy=False)
        X_float = np.array([[1, 0, 5], [2, 3, -1]], dtype=np.float64)
        X_bin = binarizer.transform(X_float)
        if init is not list:
            assert X_bin is X_float

        X_bin = toarray(X_bin)
        assert np.sum(X_bin == 0) == 2
        assert np.sum(X_bin == 1) == 4

    binarizer = Binarizer(threshold=-0.5, copy=True)
    for init in (np.array, list):
        X = init(X_.copy())

        X_bin = toarray(binarizer.transform(X))
        assert np.sum(X_bin == 0) == 1
        assert np.sum(X_bin == 1) == 5
        X_bin = binarizer.transform(X)

    # Cannot use threshold < 0 for sparse
    with pytest.raises(ValueError):
        binarizer.transform(sparse.csc_matrix(X))


def test_center_kernel():
    # Test that KernelCenterer is equivalent to StandardScaler
    # in feature space
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


def test_cv_pipeline_precomputed():
    # Cross-validate a regression on four coplanar points with the same
    # value. Use precomputed kernel to ensure Pipeline with KernelCenterer
    # is treated as a _pairwise operation.
    X = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 1, 1]])
    y_true = np.ones((4,))
    K = X.dot(X.T)
    kcent = KernelCenterer()
    pipeline = Pipeline([("kernel_centerer", kcent), ("svr", SVR())])

    # did the pipeline set the _pairwise attribute?
    assert pipeline._pairwise

    # test cross-validation, score should be almost perfect
    # NB: this test is pretty vacuous -- it's mainly to test integration
    #     of Pipeline and KernelCenterer
    y_pred = cross_val_predict(pipeline, K, y_true, cv=2)
    assert_array_almost_equal(y_true, y_pred)


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
    assert sparse.isspmatrix_coo(X), X
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_add_dummy_feature_csc():
    X = sparse.csc_matrix([[1, 0], [0, 1], [0, 1]])
    X = add_dummy_feature(X)
    assert sparse.isspmatrix_csc(X), X
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_add_dummy_feature_csr():
    X = sparse.csr_matrix([[1, 0], [0, 1], [0, 1]])
    X = add_dummy_feature(X)
    assert sparse.isspmatrix_csr(X), X
    assert_array_equal(X.toarray(), [[1, 1, 0], [1, 0, 1], [1, 0, 1]])


def test_fit_cold_start():
    X = iris.data
    X_2d = X[:, :2]

    # Scalers that have a partial_fit method
    scalers = [StandardScaler(with_mean=False, with_std=False),
               MinMaxScaler(),
               MaxAbsScaler()]

    for scaler in scalers:
        scaler.fit_transform(X)
        # with a different shape, this may break the scaler unless the internal
        # state is reset
        scaler.fit_transform(X_2d)


@pytest.mark.filterwarnings("ignore: The default value of `copy`")  # 0.23
def test_quantile_transform_valid_axis():
    X = np.array([[0, 25, 50, 75, 100],
                  [2, 4, 6, 8, 10],
                  [2.6, 4.1, 2.3, 9.5, 0.1]])

    with pytest.raises(ValueError, match="axis should be either equal "
                                         "to 0 or 1. Got axis=2"):
        quantile_transform(X.T, axis=2)


@pytest.mark.parametrize("method", ['box-cox', 'yeo-johnson'])
def test_power_transformer_notfitted(method):
    pt = PowerTransformer(method=method)
    X = np.abs(X_1col)
    with pytest.raises(NotFittedError):
        pt.transform(X)
    with pytest.raises(NotFittedError):
        pt.inverse_transform(X)


@pytest.mark.parametrize('method', ['box-cox', 'yeo-johnson'])
@pytest.mark.parametrize('standardize', [True, False])
@pytest.mark.parametrize('X', [X_1col, X_2d])
def test_power_transformer_inverse(method, standardize, X):
    # Make sure we get the original input when applying transform and then
    # inverse transform
    X = np.abs(X) if method == 'box-cox' else X
    pt = PowerTransformer(method=method, standardize=standardize)
    X_trans = pt.fit_transform(X)
    assert_almost_equal(X, pt.inverse_transform(X_trans))


def test_power_transformer_1d():
    X = np.abs(X_1col)

    for standardize in [True, False]:
        pt = PowerTransformer(method='box-cox', standardize=standardize)

        X_trans = pt.fit_transform(X)
        X_trans_func = power_transform(
            X, method='box-cox',
            standardize=standardize
        )

        X_expected, lambda_expected = stats.boxcox(X.flatten())

        if standardize:
            X_expected = scale(X_expected)

        assert_almost_equal(X_expected.reshape(-1, 1), X_trans)
        assert_almost_equal(X_expected.reshape(-1, 1), X_trans_func)

        assert_almost_equal(X, pt.inverse_transform(X_trans))
        assert_almost_equal(lambda_expected, pt.lambdas_[0])

        assert len(pt.lambdas_) == X.shape[1]
        assert isinstance(pt.lambdas_, np.ndarray)


def test_power_transformer_2d():
    X = np.abs(X_2d)

    for standardize in [True, False]:
        pt = PowerTransformer(method='box-cox', standardize=standardize)

        X_trans_class = pt.fit_transform(X)
        X_trans_func = power_transform(
            X, method='box-cox',
            standardize=standardize
        )

        for X_trans in [X_trans_class, X_trans_func]:
            for j in range(X_trans.shape[1]):
                X_expected, lmbda = stats.boxcox(X[:, j].flatten())

                if standardize:
                    X_expected = scale(X_expected)

                assert_almost_equal(X_trans[:, j], X_expected)
                assert_almost_equal(lmbda, pt.lambdas_[j])

            # Test inverse transformation
            X_inv = pt.inverse_transform(X_trans)
            assert_array_almost_equal(X_inv, X)

        assert len(pt.lambdas_) == X.shape[1]
        assert isinstance(pt.lambdas_, np.ndarray)


def test_power_transformer_boxcox_strictly_positive_exception():
    # Exceptions should be raised for negative arrays and zero arrays when
    # method is boxcox

    pt = PowerTransformer(method='box-cox')
    pt.fit(np.abs(X_2d))
    X_with_negatives = X_2d
    not_positive_message = 'strictly positive'

    with pytest.raises(ValueError, match=not_positive_message):
        pt.transform(X_with_negatives)

    with pytest.raises(ValueError, match=not_positive_message):
        pt.fit(X_with_negatives)

    with pytest.raises(ValueError, match=not_positive_message):
        power_transform(X_with_negatives, 'box-cox')

    with pytest.raises(ValueError, match=not_positive_message):
        pt.transform(np.zeros(X_2d.shape))

    with pytest.raises(ValueError, match=not_positive_message):
        pt.fit(np.zeros(X_2d.shape))

    with pytest.raises(ValueError, match=not_positive_message):
        power_transform(np.zeros(X_2d.shape), 'box-cox')


@pytest.mark.parametrize('X', [X_2d, np.abs(X_2d), -np.abs(X_2d),
                               np.zeros(X_2d.shape)])
def test_power_transformer_yeojohnson_any_input(X):
    # Yeo-Johnson method should support any kind of input
    power_transform(X, method='yeo-johnson')


@pytest.mark.parametrize("method", ['box-cox', 'yeo-johnson'])
def test_power_transformer_shape_exception(method):
    pt = PowerTransformer(method=method)
    X = np.abs(X_2d)
    pt.fit(X)

    # Exceptions should be raised for arrays with different num_columns
    # than during fitting
    wrong_shape_message = 'Input data has a different number of features'

    with pytest.raises(ValueError, match=wrong_shape_message):
        pt.transform(X[:, 0:1])

    with pytest.raises(ValueError, match=wrong_shape_message):
        pt.inverse_transform(X[:, 0:1])


def test_power_transformer_method_exception():
    pt = PowerTransformer(method='monty-python')
    X = np.abs(X_2d)

    # An exception should be raised if PowerTransformer.method isn't valid
    bad_method_message = "'method' must be one of"
    with pytest.raises(ValueError, match=bad_method_message):
        pt.fit(X)


def test_power_transformer_lambda_zero():
    pt = PowerTransformer(method='box-cox', standardize=False)
    X = np.abs(X_2d)[:, 0:1]

    # Test the lambda = 0 case
    pt.lambdas_ = np.array([0])
    X_trans = pt.transform(X)
    assert_array_almost_equal(pt.inverse_transform(X_trans), X)


def test_power_transformer_lambda_one():
    # Make sure lambda = 1 corresponds to the identity for yeo-johnson
    pt = PowerTransformer(method='yeo-johnson', standardize=False)
    X = np.abs(X_2d)[:, 0:1]

    pt.lambdas_ = np.array([1])
    X_trans = pt.transform(X)
    assert_array_almost_equal(X_trans, X)


@pytest.mark.parametrize("method, lmbda", [('box-cox', .1),
                                           ('box-cox', .5),
                                           ('yeo-johnson', .1),
                                           ('yeo-johnson', .5),
                                           ('yeo-johnson', 1.),
                                           ])
def test_optimization_power_transformer(method, lmbda):
    # Test the optimization procedure:
    # - set a predefined value for lambda
    # - apply inverse_transform to a normal dist (we get X_inv)
    # - apply fit_transform to X_inv (we get X_inv_trans)
    # - check that X_inv_trans is roughly equal to X

    rng = np.random.RandomState(0)
    n_samples = 20000
    X = rng.normal(loc=0, scale=1, size=(n_samples, 1))

    pt = PowerTransformer(method=method, standardize=False)
    pt.lambdas_ = [lmbda]
    X_inv = pt.inverse_transform(X)

    pt = PowerTransformer(method=method, standardize=False)
    X_inv_trans = pt.fit_transform(X_inv)

    assert_almost_equal(0, np.linalg.norm(X - X_inv_trans) / n_samples,
                        decimal=2)
    assert_almost_equal(0, X_inv_trans.mean(), decimal=1)
    assert_almost_equal(1, X_inv_trans.std(), decimal=1)


def test_yeo_johnson_darwin_example():
    # test from original paper "A new family of power transformations to
    # improve normality or symmetry" by Yeo and Johnson.
    X = [6.1, -8.4, 1.0, 2.0, 0.7, 2.9, 3.5, 5.1, 1.8, 3.6, 7.0, 3.0, 9.3,
         7.5, -6.0]
    X = np.array(X).reshape(-1, 1)
    lmbda = PowerTransformer(method='yeo-johnson').fit(X).lambdas_
    assert np.allclose(lmbda, 1.305, atol=1e-3)


@pytest.mark.parametrize('method', ['box-cox', 'yeo-johnson'])
def test_power_transformer_nans(method):
    # Make sure lambda estimation is not influenced by NaN values
    # and that transform() supports NaN silently

    X = np.abs(X_1col)
    pt = PowerTransformer(method=method)
    pt.fit(X)
    lmbda_no_nans = pt.lambdas_[0]

    # concat nans at the end and check lambda stays the same
    X = np.concatenate([X, np.full_like(X, np.nan)])
    X = shuffle(X, random_state=0)

    pt.fit(X)
    lmbda_nans = pt.lambdas_[0]

    assert_almost_equal(lmbda_no_nans, lmbda_nans, decimal=5)

    X_trans = pt.transform(X)
    assert_array_equal(np.isnan(X_trans), np.isnan(X))


@pytest.mark.parametrize('method', ['box-cox', 'yeo-johnson'])
@pytest.mark.parametrize('standardize', [True, False])
def test_power_transformer_fit_transform(method, standardize):
    # check that fit_transform() and fit().transform() return the same values
    X = X_1col
    if method == 'box-cox':
        X = np.abs(X)

    pt = PowerTransformer(method, standardize)
    assert_array_almost_equal(pt.fit(X).transform(X), pt.fit_transform(X))


@pytest.mark.parametrize('method', ['box-cox', 'yeo-johnson'])
@pytest.mark.parametrize('standardize', [True, False])
def test_power_transformer_copy_True(method, standardize):
    # Check that neither fit, transform, fit_transform nor inverse_transform
    # modify X inplace when copy=True
    X = X_1col
    if method == 'box-cox':
        X = np.abs(X)

    X_original = X.copy()
    assert X is not X_original  # sanity checks
    assert_array_almost_equal(X, X_original)

    pt = PowerTransformer(method, standardize, copy=True)

    pt.fit(X)
    assert_array_almost_equal(X, X_original)
    X_trans = pt.transform(X)
    assert X_trans is not X

    X_trans = pt.fit_transform(X)
    assert_array_almost_equal(X, X_original)
    assert X_trans is not X

    X_inv_trans = pt.inverse_transform(X_trans)
    assert X_trans is not X_inv_trans


@pytest.mark.parametrize('method', ['box-cox', 'yeo-johnson'])
@pytest.mark.parametrize('standardize', [True, False])
def test_power_transformer_copy_False(method, standardize):
    # check that when copy=False fit doesn't change X inplace but transform,
    # fit_transform and inverse_transform do.
    X = X_1col
    if method == 'box-cox':
        X = np.abs(X)

    X_original = X.copy()
    assert X is not X_original  # sanity checks
    assert_array_almost_equal(X, X_original)

    pt = PowerTransformer(method, standardize, copy=False)

    pt.fit(X)
    assert_array_almost_equal(X, X_original)  # fit didn't change X

    X_trans = pt.transform(X)
    assert X_trans is X

    if method == 'box-cox':
        X = np.abs(X)
    X_trans = pt.fit_transform(X)
    assert X_trans is X

    X_inv_trans = pt.inverse_transform(X_trans)
    assert X_trans is X_inv_trans


def test_power_transform_default_method():
    X = np.abs(X_2d)

    future_warning_message = (
        "The default value of 'method' "
        "will change from 'box-cox'"
    )
    assert_warns_message(FutureWarning, future_warning_message,
                         power_transform, X)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        X_trans_default = power_transform(X)

    X_trans_boxcox = power_transform(X, method='box-cox')
    assert_array_equal(X_trans_boxcox, X_trans_default)
