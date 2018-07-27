# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Denis Engemann <denis-alexander.engemann@inria.fr>
#
# License: BSD 3 clause

import numpy as np
from scipy import sparse
from scipy import linalg
from scipy import stats

import pytest

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import skip_if_32bit
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.fixes import np_version

from sklearn.utils.extmath import density
from sklearn.utils.extmath import logsumexp
from sklearn.utils.extmath import norm, squared_norm
from sklearn.utils.extmath import randomized_svd
from sklearn.utils.extmath import row_norms
from sklearn.utils.extmath import weighted_mode
from sklearn.utils.extmath import cartesian
from sklearn.utils.extmath import log_logistic
from sklearn.utils.extmath import svd_flip
from sklearn.utils.extmath import _incremental_mean_and_var
from sklearn.utils.extmath import _deterministic_vector_sign_flip
from sklearn.utils.extmath import softmax
from sklearn.utils.extmath import stable_cumsum
from sklearn.datasets.samples_generator import make_low_rank_matrix


def test_density():
    rng = np.random.RandomState(0)
    X = rng.randint(10, size=(10, 5))
    X[1, 2] = 0
    X[5, 3] = 0
    X_csr = sparse.csr_matrix(X)
    X_csc = sparse.csc_matrix(X)
    X_coo = sparse.coo_matrix(X)
    X_lil = sparse.lil_matrix(X)

    for X_ in (X_csr, X_csc, X_coo, X_lil):
        assert_equal(density(X_), density(X))


def test_uniform_weights():
    # with uniform weights, results should be identical to stats.mode
    rng = np.random.RandomState(0)
    x = rng.randint(10, size=(10, 5))
    weights = np.ones(x.shape)

    for axis in (None, 0, 1):
        mode, score = stats.mode(x, axis)
        mode2, score2 = weighted_mode(x, weights, axis)

        assert_array_equal(mode, mode2)
        assert_array_equal(score, score2)


def test_random_weights():
    # set this up so that each row should have a weighted mode of 6,
    # with a score that is easily reproduced
    mode_result = 6

    rng = np.random.RandomState(0)
    x = rng.randint(mode_result, size=(100, 10))
    w = rng.random_sample(x.shape)

    x[:, :5] = mode_result
    w[:, :5] += 1

    mode, score = weighted_mode(x, w, axis=1)

    assert_array_equal(mode, mode_result)
    assert_array_almost_equal(score.ravel(), w[:, :5].sum(1))


@ignore_warnings  # Test deprecated backport to be removed in 0.21
def test_logsumexp():
    # Try to add some smallish numbers in logspace
    x = np.array([1e-40] * 1000000)
    logx = np.log(x)
    assert_almost_equal(np.exp(logsumexp(logx)), x.sum())

    X = np.vstack([x, x])
    logX = np.vstack([logx, logx])
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=0)), X.sum(axis=0))
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=1)), X.sum(axis=1))


def check_randomized_svd_low_rank(dtype):
    # Check that extmath.randomized_svd is consistent with linalg.svd
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10
    decimal = 5 if dtype == np.float32 else 7
    dtype = np.dtype(dtype)

    # generate a matrix X of approximate effective rank `rank` and no noise
    # component (very structured signal):
    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=0.0,
                             random_state=0).astype(dtype, copy=False)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    U, s, V = linalg.svd(X, full_matrices=False)

    # Convert the singular values to the specific dtype
    U = U.astype(dtype, copy=False)
    s = s.astype(dtype, copy=False)
    V = V.astype(dtype, copy=False)

    for normalizer in ['auto', 'LU', 'QR']:  # 'none' would not be stable
        # compute the singular values of X using the fast approximate method
        Ua, sa, Va = randomized_svd(
            X, k, power_iteration_normalizer=normalizer, random_state=0)

        # If the input dtype is float, then the output dtype is float of the
        # same bit size (f32 is not upcast to f64)
        # But if the input dtype is int, the output dtype is float64
        if dtype.kind == 'f':
            assert Ua.dtype == dtype
            assert sa.dtype == dtype
            assert Va.dtype == dtype
        else:
            assert Ua.dtype == np.float64
            assert sa.dtype == np.float64
            assert Va.dtype == np.float64

        assert_equal(Ua.shape, (n_samples, k))
        assert_equal(sa.shape, (k,))
        assert_equal(Va.shape, (k, n_features))

        # ensure that the singular values of both methods are equal up to the
        # real rank of the matrix
        assert_almost_equal(s[:k], sa, decimal=decimal)

        # check the singular vectors too (while not checking the sign)
        assert_almost_equal(np.dot(U[:, :k], V[:k, :]), np.dot(Ua, Va),
                            decimal=decimal)

        # check the sparse matrix representation
        X = sparse.csr_matrix(X)

        # compute the singular values of X using the fast approximate method
        Ua, sa, Va = \
            randomized_svd(X, k, power_iteration_normalizer=normalizer,
                           random_state=0)
        if dtype.kind == 'f':
            assert Ua.dtype == dtype
            assert sa.dtype == dtype
            assert Va.dtype == dtype
        else:
            assert Ua.dtype.kind == 'f'
            assert sa.dtype.kind == 'f'
            assert Va.dtype.kind == 'f'

        assert_almost_equal(s[:rank], sa[:rank], decimal=decimal)


@pytest.mark.parametrize('dtype',
                         (np.int32, np.int64, np.float32, np.float64))
def test_randomized_svd_low_rank_all_dtypes(dtype):
    check_randomized_svd_low_rank(dtype)


@ignore_warnings  # extmath.norm is deprecated to be removed in 0.21
def test_norm_squared_norm():
    X = np.random.RandomState(42).randn(50, 63)
    X *= 100        # check stability
    X += 200

    assert_almost_equal(np.linalg.norm(X.ravel()), norm(X))
    assert_almost_equal(norm(X) ** 2, squared_norm(X), decimal=6)
    assert_almost_equal(np.linalg.norm(X), np.sqrt(squared_norm(X)), decimal=6)
    # Check the warning with an int array and np.dot potential overflow
    assert_warns_message(
                    UserWarning, 'Array type is integer, np.dot may '
                    'overflow. Data should be float type to avoid this issue',
                    squared_norm, X.astype(int))


@pytest.mark.parametrize('dtype',
                         (np.float32, np.float64))
def test_row_norms(dtype):
    X = np.random.RandomState(42).randn(100, 100)
    if dtype is np.float32:
        precision = 4
    else:
        precision = 5

    X = X.astype(dtype)
    sq_norm = (X ** 2).sum(axis=1)

    assert_array_almost_equal(sq_norm, row_norms(X, squared=True),
                              precision)
    assert_array_almost_equal(np.sqrt(sq_norm), row_norms(X), precision)

    for csr_index_dtype in [np.int32, np.int64]:
        Xcsr = sparse.csr_matrix(X, dtype=dtype)
        # csr_matrix will use int32 indices by default,
        # up-casting those to int64 when necessary
        if csr_index_dtype is np.int64:
            Xcsr.indptr = Xcsr.indptr.astype(csr_index_dtype)
            Xcsr.indices = Xcsr.indices.astype(csr_index_dtype)
        assert Xcsr.indices.dtype == csr_index_dtype
        assert Xcsr.indptr.dtype == csr_index_dtype
        assert_array_almost_equal(sq_norm, row_norms(Xcsr, squared=True),
                                  precision)
        assert_array_almost_equal(np.sqrt(sq_norm), row_norms(Xcsr),
                                  precision)


def test_randomized_svd_low_rank_with_noise():
    # Check that extmath.randomized_svd can handle noisy matrices
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # generate a matrix X wity structure approximate rank `rank` and an
    # important noisy component
    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=0.1,
                             random_state=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    _, s, _ = linalg.svd(X, full_matrices=False)

    for normalizer in ['auto', 'none', 'LU', 'QR']:
        # compute the singular values of X using the fast approximate
        # method without the iterated power method
        _, sa, _ = randomized_svd(X, k, n_iter=0,
                                  power_iteration_normalizer=normalizer,
                                  random_state=0)

        # the approximation does not tolerate the noise:
        assert_greater(np.abs(s[:k] - sa).max(), 0.01)

        # compute the singular values of X using the fast approximate
        # method with iterated power method
        _, sap, _ = randomized_svd(X, k,
                                   power_iteration_normalizer=normalizer,
                                   random_state=0)

        # the iterated power method is helping getting rid of the noise:
        assert_almost_equal(s[:k], sap, decimal=3)


def test_randomized_svd_infinite_rank():
    # Check that extmath.randomized_svd can handle noisy matrices
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # let us try again without 'low_rank component': just regularly but slowly
    # decreasing singular values: the rank of the data matrix is infinite
    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=1.0,
                             random_state=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    _, s, _ = linalg.svd(X, full_matrices=False)
    for normalizer in ['auto', 'none', 'LU', 'QR']:
        # compute the singular values of X using the fast approximate method
        # without the iterated power method
        _, sa, _ = randomized_svd(X, k, n_iter=0,
                                  power_iteration_normalizer=normalizer)

        # the approximation does not tolerate the noise:
        assert_greater(np.abs(s[:k] - sa).max(), 0.1)

        # compute the singular values of X using the fast approximate method
        # with iterated power method
        _, sap, _ = randomized_svd(X, k, n_iter=5,
                                   power_iteration_normalizer=normalizer)

        # the iterated power method is still managing to get most of the
        # structure at the requested rank
        assert_almost_equal(s[:k], sap, decimal=3)


def test_randomized_svd_transpose_consistency():
    # Check that transposing the design matrix has limited impact
    n_samples = 100
    n_features = 500
    rank = 4
    k = 10

    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=0.5,
                             random_state=0)
    assert_equal(X.shape, (n_samples, n_features))

    U1, s1, V1 = randomized_svd(X, k, n_iter=3, transpose=False,
                                random_state=0)
    U2, s2, V2 = randomized_svd(X, k, n_iter=3, transpose=True,
                                random_state=0)
    U3, s3, V3 = randomized_svd(X, k, n_iter=3, transpose='auto',
                                random_state=0)
    U4, s4, V4 = linalg.svd(X, full_matrices=False)

    assert_almost_equal(s1, s4[:k], decimal=3)
    assert_almost_equal(s2, s4[:k], decimal=3)
    assert_almost_equal(s3, s4[:k], decimal=3)

    assert_almost_equal(np.dot(U1, V1), np.dot(U4[:, :k], V4[:k, :]),
                        decimal=2)
    assert_almost_equal(np.dot(U2, V2), np.dot(U4[:, :k], V4[:k, :]),
                        decimal=2)

    # in this case 'auto' is equivalent to transpose
    assert_almost_equal(s2, s3)


def test_randomized_svd_power_iteration_normalizer():
    # randomized_svd with power_iteration_normalized='none' diverges for
    # large number of power iterations on this dataset
    rng = np.random.RandomState(42)
    X = make_low_rank_matrix(100, 500, effective_rank=50, random_state=rng)
    X += 3 * rng.randint(0, 2, size=X.shape)
    n_components = 50

    # Check that it diverges with many (non-normalized) power iterations
    U, s, V = randomized_svd(X, n_components, n_iter=2,
                             power_iteration_normalizer='none')
    A = X - U.dot(np.diag(s).dot(V))
    error_2 = linalg.norm(A, ord='fro')
    U, s, V = randomized_svd(X, n_components, n_iter=20,
                             power_iteration_normalizer='none')
    A = X - U.dot(np.diag(s).dot(V))
    error_20 = linalg.norm(A, ord='fro')
    assert_greater(np.abs(error_2 - error_20), 100)

    for normalizer in ['LU', 'QR', 'auto']:
        U, s, V = randomized_svd(X, n_components, n_iter=2,
                                 power_iteration_normalizer=normalizer,
                                 random_state=0)
        A = X - U.dot(np.diag(s).dot(V))
        error_2 = linalg.norm(A, ord='fro')

        for i in [5, 10, 50]:
            U, s, V = randomized_svd(X, n_components, n_iter=i,
                                     power_iteration_normalizer=normalizer,
                                     random_state=0)
            A = X - U.dot(np.diag(s).dot(V))
            error = linalg.norm(A, ord='fro')
            assert_greater(15, np.abs(error_2 - error))


def test_randomized_svd_sparse_warnings():
    # randomized_svd throws a warning for lil and dok matrix
    rng = np.random.RandomState(42)
    X = make_low_rank_matrix(50, 20, effective_rank=10, random_state=rng)
    n_components = 5
    for cls in (sparse.lil_matrix, sparse.dok_matrix):
        X = cls(X)
        assert_warns_message(
            sparse.SparseEfficiencyWarning,
            "Calculating SVD of a {} is expensive. "
            "csr_matrix is more efficient.".format(cls.__name__),
            randomized_svd, X, n_components, n_iter=1,
            power_iteration_normalizer='none')


def test_svd_flip():
    # Check that svd_flip works in both situations, and reconstructs input.
    rs = np.random.RandomState(1999)
    n_samples = 20
    n_features = 10
    X = rs.randn(n_samples, n_features)

    # Check matrix reconstruction
    U, S, V = linalg.svd(X, full_matrices=False)
    U1, V1 = svd_flip(U, V, u_based_decision=False)
    assert_almost_equal(np.dot(U1 * S, V1), X, decimal=6)

    # Check transposed matrix reconstruction
    XT = X.T
    U, S, V = linalg.svd(XT, full_matrices=False)
    U2, V2 = svd_flip(U, V, u_based_decision=True)
    assert_almost_equal(np.dot(U2 * S, V2), XT, decimal=6)

    # Check that different flip methods are equivalent under reconstruction
    U_flip1, V_flip1 = svd_flip(U, V, u_based_decision=True)
    assert_almost_equal(np.dot(U_flip1 * S, V_flip1), XT, decimal=6)
    U_flip2, V_flip2 = svd_flip(U, V, u_based_decision=False)
    assert_almost_equal(np.dot(U_flip2 * S, V_flip2), XT, decimal=6)


def test_randomized_svd_sign_flip():
    a = np.array([[2.0, 0.0], [0.0, 1.0]])
    u1, s1, v1 = randomized_svd(a, 2, flip_sign=True, random_state=41)
    for seed in range(10):
        u2, s2, v2 = randomized_svd(a, 2, flip_sign=True, random_state=seed)
        assert_almost_equal(u1, u2)
        assert_almost_equal(v1, v2)
        assert_almost_equal(np.dot(u2 * s2, v2), a)
        assert_almost_equal(np.dot(u2.T, u2), np.eye(2))
        assert_almost_equal(np.dot(v2.T, v2), np.eye(2))


def test_randomized_svd_sign_flip_with_transpose():
    # Check if the randomized_svd sign flipping is always done based on u
    # irrespective of transpose.
    # See https://github.com/scikit-learn/scikit-learn/issues/5608
    # for more details.
    def max_loading_is_positive(u, v):
        """
        returns bool tuple indicating if the values maximising np.abs
        are positive across all rows for u and across all columns for v.
        """
        u_based = (np.abs(u).max(axis=0) == u.max(axis=0)).all()
        v_based = (np.abs(v).max(axis=1) == v.max(axis=1)).all()
        return u_based, v_based

    mat = np.arange(10 * 8).reshape(10, -1)

    # Without transpose
    u_flipped, _, v_flipped = randomized_svd(mat, 3, flip_sign=True)
    u_based, v_based = max_loading_is_positive(u_flipped, v_flipped)
    assert_true(u_based)
    assert_false(v_based)

    # With transpose
    u_flipped_with_transpose, _, v_flipped_with_transpose = randomized_svd(
        mat, 3, flip_sign=True, transpose=True)
    u_based, v_based = max_loading_is_positive(
        u_flipped_with_transpose, v_flipped_with_transpose)
    assert_true(u_based)
    assert_false(v_based)


def test_cartesian():
    # Check if cartesian product delivers the right results

    axes = (np.array([1, 2, 3]), np.array([4, 5]), np.array([6, 7]))

    true_out = np.array([[1, 4, 6],
                         [1, 4, 7],
                         [1, 5, 6],
                         [1, 5, 7],
                         [2, 4, 6],
                         [2, 4, 7],
                         [2, 5, 6],
                         [2, 5, 7],
                         [3, 4, 6],
                         [3, 4, 7],
                         [3, 5, 6],
                         [3, 5, 7]])

    out = cartesian(axes)
    assert_array_equal(true_out, out)

    # check single axis
    x = np.arange(3)
    assert_array_equal(x[:, np.newaxis], cartesian((x,)))


def test_logistic_sigmoid():
    # Check correctness and robustness of logistic sigmoid implementation
    def naive_log_logistic(x):
        return np.log(1 / (1 + np.exp(-x)))

    x = np.linspace(-2, 2, 50)
    assert_array_almost_equal(log_logistic(x), naive_log_logistic(x))

    extreme_x = np.array([-100., 100.])
    assert_array_almost_equal(log_logistic(extreme_x), [-100, 0])


def test_incremental_variance_update_formulas():
    # Test Youngs and Cramer incremental variance formulas.
    # Doggie data from http://www.mathsisfun.com/data/standard-deviation.html
    A = np.array([[600, 470, 170, 430, 300],
                  [600, 470, 170, 430, 300],
                  [600, 470, 170, 430, 300],
                  [600, 470, 170, 430, 300]]).T
    idx = 2
    X1 = A[:idx, :]
    X2 = A[idx:, :]

    old_means = X1.mean(axis=0)
    old_variances = X1.var(axis=0)
    old_sample_count = np.full(X1.shape[1], X1.shape[0], dtype=np.int32)
    final_means, final_variances, final_count = \
        _incremental_mean_and_var(X2, old_means, old_variances,
                                  old_sample_count)
    assert_almost_equal(final_means, A.mean(axis=0), 6)
    assert_almost_equal(final_variances, A.var(axis=0), 6)
    assert_almost_equal(final_count, A.shape[0])


def test_incremental_mean_and_variance_ignore_nan():
    old_means = np.array([535., 535., 535., 535.])
    old_variances = np.array([4225., 4225., 4225., 4225.])
    old_sample_count = np.array([2, 2, 2, 2], dtype=np.int32)

    X = np.array([[170, 170, 170, 170],
                  [430, 430, 430, 430],
                  [300, 300, 300, 300]])

    X_nan = np.array([[170, np.nan, 170, 170],
                      [np.nan, 170, 430, 430],
                      [430, 430, np.nan, 300],
                      [300, 300, 300, np.nan]])

    X_means, X_variances, X_count = _incremental_mean_and_var(
        X, old_means, old_variances, old_sample_count)
    X_nan_means, X_nan_variances, X_nan_count = _incremental_mean_and_var(
        X_nan, old_means, old_variances, old_sample_count)

    assert_allclose(X_nan_means, X_means)
    assert_allclose(X_nan_variances, X_variances)
    assert_allclose(X_nan_count, X_count)


@skip_if_32bit
def test_incremental_variance_numerical_stability():
    # Test Youngs and Cramer incremental variance formulas.

    def np_var(A):
        return A.var(axis=0)

    # Naive one pass variance computation - not numerically stable
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    def one_pass_var(X):
        n = X.shape[0]
        exp_x2 = (X ** 2).sum(axis=0) / n
        expx_2 = (X.sum(axis=0) / n) ** 2
        return exp_x2 - expx_2

    # Two-pass algorithm, stable.
    # We use it as a benchmark. It is not an online algorithm
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass_algorithm
    def two_pass_var(X):
        mean = X.mean(axis=0)
        Y = X.copy()
        return np.mean((Y - mean)**2, axis=0)

    # Naive online implementation
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
    # This works only for chunks for size 1
    def naive_mean_variance_update(x, last_mean, last_variance,
                                   last_sample_count):
        updated_sample_count = (last_sample_count + 1)
        samples_ratio = last_sample_count / float(updated_sample_count)
        updated_mean = x / updated_sample_count + last_mean * samples_ratio
        updated_variance = last_variance * samples_ratio + \
            (x - last_mean) * (x - updated_mean) / updated_sample_count
        return updated_mean, updated_variance, updated_sample_count

    # We want to show a case when one_pass_var has error > 1e-3 while
    # _batch_mean_variance_update has less.
    tol = 200
    n_features = 2
    n_samples = 10000
    x1 = np.array(1e8, dtype=np.float64)
    x2 = np.log(1e-5, dtype=np.float64)
    A0 = np.full((n_samples // 2, n_features), x1, dtype=np.float64)
    A1 = np.full((n_samples // 2, n_features), x2, dtype=np.float64)
    A = np.vstack((A0, A1))

    # Older versions of numpy have different precision
    # In some old version, np.var is not stable
    if np.abs(np_var(A) - two_pass_var(A)).max() < 1e-6:
        stable_var = np_var
    else:
        stable_var = two_pass_var

    # Naive one pass var: >tol (=1063)
    assert_greater(np.abs(stable_var(A) - one_pass_var(A)).max(), tol)

    # Starting point for online algorithms: after A0

    # Naive implementation: >tol (436)
    mean, var, n = A0[0, :], np.zeros(n_features), n_samples // 2
    for i in range(A1.shape[0]):
        mean, var, n = \
            naive_mean_variance_update(A1[i, :], mean, var, n)
    assert_equal(n, A.shape[0])
    # the mean is also slightly unstable
    assert_greater(np.abs(A.mean(axis=0) - mean).max(), 1e-6)
    assert_greater(np.abs(stable_var(A) - var).max(), tol)

    # Robust implementation: <tol (177)
    mean, var = A0[0, :], np.zeros(n_features)
    n = np.full(n_features, n_samples // 2, dtype=np.int32)
    for i in range(A1.shape[0]):
        mean, var, n = \
            _incremental_mean_and_var(A1[i, :].reshape((1, A1.shape[1])),
                                      mean, var, n)
    assert_array_equal(n, A.shape[0])
    assert_array_almost_equal(A.mean(axis=0), mean)
    assert_greater(tol, np.abs(stable_var(A) - var).max())


def test_incremental_variance_ddof():
    # Test that degrees of freedom parameter for calculations are correct.
    rng = np.random.RandomState(1999)
    X = rng.randn(50, 10)
    n_samples, n_features = X.shape
    for batch_size in [11, 20, 37]:
        steps = np.arange(0, X.shape[0], batch_size)
        if steps[-1] != X.shape[0]:
            steps = np.hstack([steps, n_samples])

        for i, j in zip(steps[:-1], steps[1:]):
            batch = X[i:j, :]
            if i == 0:
                incremental_means = batch.mean(axis=0)
                incremental_variances = batch.var(axis=0)
                # Assign this twice so that the test logic is consistent
                incremental_count = batch.shape[0]
                sample_count = np.full(batch.shape[1], batch.shape[0],
                                       dtype=np.int32)
            else:
                result = _incremental_mean_and_var(
                    batch, incremental_means, incremental_variances,
                    sample_count)
                (incremental_means, incremental_variances,
                 incremental_count) = result
                sample_count += batch.shape[0]

            calculated_means = np.mean(X[:j], axis=0)
            calculated_variances = np.var(X[:j], axis=0)
            assert_almost_equal(incremental_means, calculated_means, 6)
            assert_almost_equal(incremental_variances,
                                calculated_variances, 6)
            assert_array_equal(incremental_count, sample_count)


def test_vector_sign_flip():
    # Testing that sign flip is working & largest value has positive sign
    data = np.random.RandomState(36).randn(5, 5)
    max_abs_rows = np.argmax(np.abs(data), axis=1)
    data_flipped = _deterministic_vector_sign_flip(data)
    max_rows = np.argmax(data_flipped, axis=1)
    assert_array_equal(max_abs_rows, max_rows)
    signs = np.sign(data[range(data.shape[0]), max_abs_rows])
    assert_array_equal(data, data_flipped * signs[:, np.newaxis])


def test_softmax():
    rng = np.random.RandomState(0)
    X = rng.randn(3, 5)
    exp_X = np.exp(X)
    sum_exp_X = np.sum(exp_X, axis=1).reshape((-1, 1))
    assert_array_almost_equal(softmax(X), exp_X / sum_exp_X)


def test_stable_cumsum():
    if np_version < (1, 9):
        raise SkipTest("Sum is as unstable as cumsum for numpy < 1.9")
    assert_array_equal(stable_cumsum([1, 2, 3]), np.cumsum([1, 2, 3]))
    r = np.random.RandomState(0).rand(100000)
    assert_warns(RuntimeWarning, stable_cumsum, r, rtol=0, atol=0)

    # test axis parameter
    A = np.random.RandomState(36).randint(1000, size=(5, 5, 5))
    assert_array_equal(stable_cumsum(A, axis=0), np.cumsum(A, axis=0))
    assert_array_equal(stable_cumsum(A, axis=1), np.cumsum(A, axis=1))
    assert_array_equal(stable_cumsum(A, axis=2), np.cumsum(A, axis=2))
