
# Authors:
#
#          Giorgio Patrini
#
# License: BSD 3 clause

import warnings
import numpy as np
import numpy.linalg as la
from scipy import sparse
from distutils.version import LooseVersion

from sklearn.utils import gen_batches

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import clean_warning_registry
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_less
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_less_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import skip_if_32bit

from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.preprocessing.data import _transform_selected
from sklearn.preprocessing.data import _handle_zeros_in_scale
from sklearn.preprocessing.data import Binarizer
from sklearn.preprocessing.data import KernelCenterer
from sklearn.preprocessing.data import Normalizer
from sklearn.preprocessing.data import normalize
from sklearn.preprocessing.data import OneHotEncoder
from sklearn.preprocessing.data import StandardScaler
from sklearn.preprocessing.data import scale
from sklearn.preprocessing.data import MinMaxScaler
from sklearn.preprocessing.data import minmax_scale
from sklearn.preprocessing.data import MaxAbsScaler
from sklearn.preprocessing.data import maxabs_scale
from sklearn.preprocessing.data import RobustScaler
from sklearn.preprocessing.data import robust_scale
from sklearn.preprocessing.data import add_dummy_feature
from sklearn.preprocessing.data import PolynomialFeatures
from sklearn.exceptions import DataConversionWarning

from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_predict
from sklearn.svm import SVR

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
        assert_equal((i + 1) * chunk_size, n_samples_seen)
    else:
        assert_equal(i * chunk_size + (batch_stop - batch_start),
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

    assert_equal(interact.powers_.shape, (interact.n_output_features_,
                 interact.n_input_features_))


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
    feature_names = poly.get_feature_names([u"\u0001F40D", u"\u262E", u"\u05D0"])
    assert_array_equal([u"1", u"\u0001F40D", u"\u262E", u"\u05D0"],
                       feature_names)


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
        assert_equal(scaler.n_samples_seen_, X.shape[0])

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones(5).reshape(5, 1)
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_almost_equal(scaler.mean_, 1.)
    assert_almost_equal(scaler.scale_, 1.)
    assert_array_almost_equal(X_scaled.mean(axis=0), .0)
    assert_array_almost_equal(X_scaled.std(axis=0), .0)
    assert_equal(scaler.n_samples_seen_, X.shape[0])


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

    x = np.zeros(8, dtype=np.float64) + np.log(1e-5, dtype=np.float64)
    if LooseVersion(np.__version__) >= LooseVersion('1.9'):
        # This does not raise a warning as the number of samples is too low
        # to trigger the problem in recent numpy
        x_scaled = assert_no_warnings(scale, x)
        assert_array_almost_equal(scale(x), np.zeros(8))
    else:
        w = "standard deviation of the data is probably very close to 0"
        x_scaled = assert_warns_message(UserWarning, w, scale, x)
        assert_array_almost_equal(x_scaled, np.zeros(8))

    # with 2 more samples, the std computation run into numerical issues:
    x = np.zeros(10, dtype=np.float64) + np.log(1e-5, dtype=np.float64)
    w = "standard deviation of the data is probably very close to 0"
    x_scaled = assert_warns_message(UserWarning, w, scale, x)
    assert_array_almost_equal(x_scaled, np.zeros(10))

    x = np.ones(10, dtype=np.float64) * 1e-100
    x_small_scaled = assert_no_warnings(scale, x)
    assert_array_almost_equal(x_small_scaled, np.zeros(10))

    # Large values can cause (often recoverable) numerical stability issues:
    x_big = np.ones(10, dtype=np.float64) * 1e100
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
    assert_false(np.any(np.isnan(X_scaled)))
    assert_equal(scaler.n_samples_seen_, n_samples)

    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
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
    assert_array_almost_equal(X_scaled.mean(axis=1), n_samples * [0.0])
    X_scaled = scale(X, axis=1, with_std=True)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=1), n_samples * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=1), n_samples * [1.0])
    # Check that the data hasn't been modified
    assert_true(X_scaled is not X)

    X_scaled = scaler.fit(X).transform(X, copy=False)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is X)

    X = rng.randn(4, 5)
    X[:, 0] = 1.0  # first feature is a constant, non zero feature
    scaler = StandardScaler()
    X_scaled = scaler.fit(X).transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))
    assert_array_almost_equal(X_scaled.mean(axis=0), n_features * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0), [0., 1., 1., 1., 1.])
    # Check that X has not been copied
    assert_true(X_scaled is not X)


def test_handle_zeros_in_scale():
    s1 = np.array([0, 1, 2, 3])
    s2 = _handle_zeros_in_scale(s1, copy=True)

    assert_false(s1[0] == s2[0])
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
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)
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
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)
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
        assert_equal(scaler_batch.var_, scaler_incr.var_)  # Nones
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)

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
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)


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
        assert_true(scaler.mean_ is not None)
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
        epsilon = np.nextafter(0, 1)
        assert_array_less(zero, scaler_incr.var_ + epsilon)  # as less or equal
        assert_array_less(zero, scaler_incr.scale_ + epsilon)
        # (i+1) because the Scaler has been already fitted
        assert_equal((i + 1), scaler_incr.n_samples_seen_)


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
    assert_raises(ValueError, scaler.fit, X)


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
        assert_equal(scaler.n_samples_seen_, X.shape[0])

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones(5).reshape(5, 1)
    scaler = MinMaxScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_greater_equal(X_scaled.min(), 0.)
    assert_less_equal(X_scaled.max(), 1.)
    assert_equal(scaler.n_samples_seen_, X.shape[0])

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

    assert_raises(ValueError, StandardScaler().fit, X_csr)
    assert_raises(ValueError, StandardScaler().fit, X_csc)

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
    X_csc_scaled = scaler_csc.transform(X_csc, copy=True)
    assert_false(np.any(np.isnan(X_csc_scaled.data)))

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
    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        X_null = null_transform.fit_transform(X_csr)
    assert_array_equal(X_null.data, X_csr.data)
    X_orig = null_transform.inverse_transform(X_null)
    assert_array_equal(X_orig.data, X_csr.data)

    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        scaler = StandardScaler(with_mean=False).fit(X)
        X_scaled = scaler.transform(X, copy=True)
    assert_false(np.any(np.isnan(X_scaled)))

    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        scaler_csr = StandardScaler(with_mean=False).fit(X_csr)
        X_csr_scaled = scaler_csr.transform(X_csr, copy=True)
    assert_false(np.any(np.isnan(X_csr_scaled.data)))

    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        scaler_csc = StandardScaler(with_mean=False).fit(X_csc)
        X_csc_scaled = scaler_csc.transform(X_csc, copy=True)
    assert_false(np.any(np.isnan(X_csc_scaled.data)))

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
    assert_raises(ValueError, scale, X_csr, with_mean=True)
    assert_raises(ValueError, StandardScaler(with_mean=True).fit, X_csr)

    assert_raises(ValueError, scale, X_csc, with_mean=True)
    assert_raises(ValueError, StandardScaler(with_mean=True).fit, X_csc)

    # check transform and inverse_transform after a fit on a dense array
    scaler = StandardScaler(with_mean=True).fit(X)
    assert_raises(ValueError, scaler.transform, X_csr)
    assert_raises(ValueError, scaler.transform, X_csc)

    X_transformed_csr = sparse.csr_matrix(scaler.transform(X))
    assert_raises(ValueError, scaler.inverse_transform, X_transformed_csr)

    X_transformed_csc = sparse.csc_matrix(scaler.transform(X))
    assert_raises(ValueError, scaler.inverse_transform, X_transformed_csc)


def test_scale_input_finiteness_validation():
    # Check if non finite inputs raise ValueError
    X = [[np.nan, 5, 6, 7, 8]]
    assert_raises_regex(ValueError,
                        "Input contains NaN, infinity or a value too large",
                        scale, X)

    X = [[np.inf, 5, 6, 7, 8]]
    assert_raises_regex(ValueError,
                        "Input contains NaN, infinity or a value too large",
                        scale, X)


def test_robust_scaler_2d_arrays():
    # Test robust scaling of 2d array along first axis
    rng = np.random.RandomState(0)
    X = rng.randn(4, 5)
    X[:, 0] = 0.0  # first feature is always of zero

    scaler = RobustScaler()
    X_scaled = scaler.fit(X).transform(X)

    assert_array_almost_equal(np.median(X_scaled, axis=0), 5 * [0.0])
    assert_array_almost_equal(X_scaled.std(axis=0)[0], 0)


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


def test_robust_scaler_invalid_range():
    for range_ in [
        (-1, 90),
        (-2, -3),
        (10, 101),
        (100.5, 101),
        (90, 50),
    ]:
        scaler = RobustScaler(quantile_range=range_)

        assert_raises_regex(ValueError, 'Invalid quantile range: \(',
                            scaler.fit, iris.data)


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


def test_warning_scaling_integers():
    # Check warning when scaling integer data
    X = np.array([[1, 2, 0],
                  [0, 0, 0]], dtype=np.uint8)

    w = "Data with input dtype uint8 was converted to float64"

    clean_warning_registry()
    assert_warns_message(DataConversionWarning, w, scale, X)
    assert_warns_message(DataConversionWarning, w, StandardScaler().fit, X)
    assert_warns_message(DataConversionWarning, w, MinMaxScaler().fit, X)


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
        assert_equal(scaler.n_samples_seen_, X.shape[0])

        # check inverse transform
        X_scaled_back = scaler.inverse_transform(X_scaled)
        assert_array_almost_equal(X_scaled_back, X)

    # Constant feature
    X = np.ones(5).reshape(5, 1)
    scaler = MaxAbsScaler()
    X_scaled = scaler.fit(X).transform(X)
    assert_array_almost_equal(np.abs(X_scaled.max(axis=0)), 1.)
    assert_equal(scaler.n_samples_seen_, X.shape[0])

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
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)
        assert_equal(scaler_batch.n_samples_seen_,
                     scaler_incr_csr.n_samples_seen_)
        assert_equal(scaler_batch.n_samples_seen_,
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
        assert_equal(scaler_batch.n_samples_seen_, scaler_incr.n_samples_seen_)
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
        assert_true(X_norm1 is not X)
        X_norm1 = toarray(X_norm1)

        normalizer = Normalizer(norm='max', copy=False)
        X_norm2 = normalizer.transform(X)
        assert_true(X_norm2 is X)
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

        assert_true(X_norm is not X)
        assert_true(isinstance(X_norm, sparse.csr_matrix))

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
    assert_raises(ValueError, normalize, [[0]], axis=2)
    assert_raises(ValueError, normalize, [[0]], norm='l3')

    rs = np.random.RandomState(0)
    X_dense = rs.randn(10, 5)
    X_sparse = sparse.csr_matrix(X_dense)
    ones = np.ones((10))
    for X in (X_dense, X_sparse):
        for dtype in (np.float32, np.float64):
            for norm in ('l1', 'l2'):
                X = X.astype(dtype)
                X_norm = normalize(X, norm=norm)
                assert_equal(X_norm.dtype, dtype)

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
        assert_raises(NotImplementedError, normalize, X_sparse,
                      norm=norm, return_norm=True)
    _, norms = normalize(X_sparse, norm='max', return_norm=True)
    assert_array_almost_equal(norms, np.array([4.0, 1.0, 3.0]))


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

        binarizer = Binarizer(copy=False)
        X_float = np.array([[1, 0, 5], [2, 3, -1]], dtype=np.float64)
        X_bin = binarizer.transform(X_float)
        if init is not list:
            assert_true(X_bin is X_float)

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
    assert_true(pipeline._pairwise)

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
    # Test OneHotEncoder's fit and transform.
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
    error_msg = "unknown categorical feature present \[2\] during transform."
    assert_raises_regex(ValueError, error_msg, enc.transform, X_too_large)
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
    # check for sparse=False
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


def test_transform_selected_copy_arg():
    # transformer that alters X
    def _mutating_transformer(X):
        X[0, 0] = X[0, 0] + 1
        return X

    original_X = np.asarray([[1, 2], [3, 4]])
    expected_Xtr = [[2, 2], [3, 4]]

    X = original_X.copy()
    Xtr = _transform_selected(X, _mutating_transformer, copy=True,
                              selected='all')

    assert_array_equal(toarray(X), toarray(original_X))
    assert_array_equal(toarray(Xtr), expected_Xtr)


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


def test_one_hot_encoder_unknown_transform():
    X = np.array([[0, 2, 1], [1, 0, 3], [1, 0, 2]])
    y = np.array([[4, 1, 1]])

    # Test that one hot encoder raises error for unknown features
    # present during transform.
    oh = OneHotEncoder(handle_unknown='error')
    oh.fit(X)
    assert_raises(ValueError, oh.transform, y)

    # Test the ignore option, ignores unknown features.
    oh = OneHotEncoder(handle_unknown='ignore')
    oh.fit(X)
    assert_array_equal(
        oh.transform(y).toarray(),
        np.array([[0.,  0.,  0.,  0.,  1.,  0.,  0.]]))

    # Raise error if handle_unknown is neither ignore or error.
    oh = OneHotEncoder(handle_unknown='42')
    oh.fit(X)
    assert_raises(ValueError, oh.transform, y)


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
