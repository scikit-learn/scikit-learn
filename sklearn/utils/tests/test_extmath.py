# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Denis Engemann <d.engemann@fz-juelich.de>
#
# License: BSD 3 clause
import warnings

import numpy as np
from scipy import sparse
from scipy import linalg
from scipy import stats

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises

from sklearn.utils.extmath import density
from sklearn.utils.extmath import logsumexp
from sklearn.utils.extmath import norm, squared_norm
from sklearn.utils.extmath import randomized_svd
from sklearn.utils.extmath import row_norms
from sklearn.utils.extmath import weighted_mode
from sklearn.utils.extmath import cartesian
from sklearn.utils.extmath import log_logistic, logistic_sigmoid
from sklearn.utils.extmath import fast_dot, _fast_dot
from sklearn.utils.extmath import svd_flip
from sklearn.utils.extmath import _batch_mean_variance_update
from sklearn.utils.extmath import _deterministic_vector_sign_flip
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

        assert_true(np.all(mode == mode2))
        assert_true(np.all(score == score2))


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


def test_logsumexp():
    # Try to add some smallish numbers in logspace
    x = np.array([1e-40] * 1000000)
    logx = np.log(x)
    assert_almost_equal(np.exp(logsumexp(logx)), x.sum())

    X = np.vstack([x, x])
    logX = np.vstack([logx, logx])
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=0)), X.sum(axis=0))
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=1)), X.sum(axis=1))


def test_randomized_svd_low_rank():
    # Check that extmath.randomized_svd is consistent with linalg.svd
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # generate a matrix X of approximate effective rank `rank` and no noise
    # component (very structured signal):
    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=0.0,
                             random_state=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    U, s, V = linalg.svd(X, full_matrices=False)

    # compute the singular values of X using the fast approximate method
    Ua, sa, Va = randomized_svd(X, k)
    assert_equal(Ua.shape, (n_samples, k))
    assert_equal(sa.shape, (k,))
    assert_equal(Va.shape, (k, n_features))

    # ensure that the singular values of both methods are equal up to the real
    # rank of the matrix
    assert_almost_equal(s[:k], sa)

    # check the singular vectors too (while not checking the sign)
    assert_almost_equal(np.dot(U[:, :k], V[:k, :]), np.dot(Ua, Va))

    # check the sparse matrix representation
    X = sparse.csr_matrix(X)

    # compute the singular values of X using the fast approximate method
    Ua, sa, Va = randomized_svd(X, k)
    assert_almost_equal(s[:rank], sa[:rank])


def test_norm_squared_norm():
    X = np.random.RandomState(42).randn(50, 63)
    X *= 100        # check stability
    X += 200

    assert_almost_equal(np.linalg.norm(X.ravel()), norm(X))
    assert_almost_equal(norm(X) ** 2, squared_norm(X), decimal=6)
    assert_almost_equal(np.linalg.norm(X), np.sqrt(squared_norm(X)), decimal=6)


def test_row_norms():
    X = np.random.RandomState(42).randn(100, 100)
    sq_norm = (X ** 2).sum(axis=1)

    assert_array_almost_equal(sq_norm, row_norms(X, squared=True), 5)
    assert_array_almost_equal(np.sqrt(sq_norm), row_norms(X))

    Xcsr = sparse.csr_matrix(X, dtype=np.float32)
    assert_array_almost_equal(sq_norm, row_norms(Xcsr, squared=True), 5)
    assert_array_almost_equal(np.sqrt(sq_norm), row_norms(Xcsr))


def test_randomized_svd_low_rank_with_noise():
    # Check that extmath.randomized_svd can handle noisy matrices
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # generate a matrix X wity structure approximate rank `rank` and an
    # important noisy component
    X = make_low_rank_matrix(n_samples=n_samples, n_features=n_features,
                             effective_rank=rank, tail_strength=0.5,
                             random_state=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    _, s, _ = linalg.svd(X, full_matrices=False)

    # compute the singular values of X using the fast approximate method
    # without the iterated power method
    _, sa, _ = randomized_svd(X, k, n_iter=0)

    # the approximation does not tolerate the noise:
    assert_greater(np.abs(s[:k] - sa).max(), 0.05)

    # compute the singular values of X using the fast approximate method with
    # iterated power method
    _, sap, _ = randomized_svd(X, k, n_iter=5)

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

    # compute the singular values of X using the fast approximate method
    # without the iterated power method
    _, sa, _ = randomized_svd(X, k, n_iter=0)

    # the approximation does not tolerate the noise:
    assert_greater(np.abs(s[:k] - sa).max(), 0.1)

    # compute the singular values of X using the fast approximate method with
    # iterated power method
    _, sap, _ = randomized_svd(X, k, n_iter=5)

    # the iterated power method is still managing to get most of the structure
    # at the requested rank
    assert_almost_equal(s[:k], sap, decimal=3)


def test_randomized_svd_transpose_consistency():
    # Check that transposing the design matrix has limit impact
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
    naive_logistic = lambda x: 1 / (1 + np.exp(-x))
    naive_log_logistic = lambda x: np.log(naive_logistic(x))

    x = np.linspace(-2, 2, 50)
    with warnings.catch_warnings(record=True):
        assert_array_almost_equal(logistic_sigmoid(x), naive_logistic(x))
    assert_array_almost_equal(log_logistic(x), naive_log_logistic(x))

    extreme_x = np.array([-100., 100.])
    assert_array_almost_equal(log_logistic(extreme_x), [-100, 0])


def test_fast_dot():
    # Check fast dot blas wrapper function
    if fast_dot is np.dot:
        return

    rng = np.random.RandomState(42)
    A = rng.random_sample([2, 10])
    B = rng.random_sample([2, 10])

    try:
        linalg.get_blas_funcs(['gemm'])[0]
        has_blas = True
    except (AttributeError, ValueError):
        has_blas = False

    if has_blas:
        # Test _fast_dot for invalid input.

        # Maltyped data.
        for dt1, dt2 in [['f8', 'f4'], ['i4', 'i4']]:
            assert_raises(ValueError, _fast_dot, A.astype(dt1),
                          B.astype(dt2).T)

        # Malformed data.

        ## ndim == 0
        E = np.empty(0)
        assert_raises(ValueError, _fast_dot, E, E)

        ## ndim == 1
        assert_raises(ValueError, _fast_dot, A, A[0])

        ## ndim > 2
        assert_raises(ValueError, _fast_dot, A.T, np.array([A, A]))

        ## min(shape) == 1
        assert_raises(ValueError, _fast_dot, A, A[0, :][None, :])

        # test for matrix mismatch error
        assert_raises(ValueError, _fast_dot, A, A)

    # Test cov-like use case + dtypes.
    for dtype in ['f8', 'f4']:
        A = A.astype(dtype)
        B = B.astype(dtype)

        #  col < row
        C = np.dot(A.T, A)
        C_ = fast_dot(A.T, A)
        assert_almost_equal(C, C_, decimal=5)

        C = np.dot(A.T, B)
        C_ = fast_dot(A.T, B)
        assert_almost_equal(C, C_, decimal=5)

        C = np.dot(A, B.T)
        C_ = fast_dot(A, B.T)
        assert_almost_equal(C, C_, decimal=5)

    # Test square matrix * rectangular use case.
    A = rng.random_sample([2, 2])
    for dtype in ['f8', 'f4']:
        A = A.astype(dtype)
        B = B.astype(dtype)

        C = np.dot(A, B)
        C_ = fast_dot(A, B)
        assert_almost_equal(C, C_, decimal=5)

        C = np.dot(A.T, B)
        C_ = fast_dot(A.T, B)
        assert_almost_equal(C, C_, decimal=5)

    if has_blas:
        for x in [np.array([[d] * 10] * 2) for d in [np.inf, np.nan]]:
            assert_raises(ValueError, _fast_dot, x, x.T)


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
    old_sample_count = X1.shape[0]
    final_means, final_variances, final_count = _batch_mean_variance_update(
        X2, old_means, old_variances, old_sample_count)
    assert_almost_equal(final_means, A.mean(axis=0), 6)
    assert_almost_equal(final_variances, A.var(axis=0), 6)
    assert_almost_equal(final_count, A.shape[0])


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
                sample_count = batch.shape[0]
            else:
                result = _batch_mean_variance_update(batch, incremental_means,
                                                    incremental_variances,
                                                    sample_count)
                (incremental_means, incremental_variances,
                 incremental_count) = result
                sample_count += batch.shape[0]

            calculated_means = np.mean(X[:j], axis=0)
            calculated_variances = np.var(X[:j], axis=0)
            assert_almost_equal(incremental_means, calculated_means, 6)
            assert_almost_equal(incremental_variances,
                                calculated_variances, 6)
            assert_equal(incremental_count, sample_count)


def test_vector_sign_flip():
    # Testing that sign flip is working & largest value has positive sign
    data = np.random.RandomState(36).randn(5, 5)
    max_abs_rows = np.argmax(np.abs(data), axis=1)
    data_flipped = _deterministic_vector_sign_flip(data)
    max_rows = np.argmax(data_flipped, axis=1)
    assert_array_equal(max_abs_rows, max_rows)
    signs = np.sign(data[range(data.shape[0]), max_abs_rows])
    assert_array_equal(data, data_flipped * signs[:, np.newaxis])


if __name__ == '__main__':
    import nose
    nose.runmodule()
