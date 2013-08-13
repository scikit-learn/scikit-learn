# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD 3 clause
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

from sklearn.utils.extmath import density
from sklearn.utils.extmath import logsumexp
from sklearn.utils.extmath import randomized_svd
from sklearn.utils.extmath import weighted_mode
from sklearn.utils.extmath import cartesian
from sklearn.utils.extmath import logistic_sigmoid
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

    assert_true(np.all(mode == mode_result))
    assert_true(np.all(score.ravel() == w[:, :5].sum(1)))


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
    """Check that extmath.randomized_svd is consistent with linalg.svd"""
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


def test_randomized_svd_low_rank_with_noise():
    """Check that extmath.randomized_svd can handle noisy matrices"""
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
    """Check that extmath.randomized_svd can handle noisy matrices"""
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
    """Check that transposing the design matrix has limit impact"""
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
    """Check if cartesian product delivers the right results"""

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
    """Check correctness and robustness of logistic sigmoid implementation"""
    naive_logsig = lambda x: 1 / (1 + np.exp(-x))
    naive_log_logsig = lambda x: np.log(naive_logsig(x))

    # Simulate the previous Cython implementations of logistic_sigmoid based on
    #http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression
    def stable_logsig(x):
        out = np.zeros_like(x)
        positive = x > 0
        negative = x <= 0
        out[positive] = 1. / (1 + np.exp(-x[positive]))
        out[negative] = np.exp(x[negative]) / (1. + np.exp(x[negative]))
        return out

    x = np.linspace(-2, 2, 50)
    assert_array_almost_equal(logistic_sigmoid(x), naive_logsig(x))
    assert_array_almost_equal(logistic_sigmoid(x, log=True),
                              naive_log_logsig(x))
    assert_array_almost_equal(logistic_sigmoid(x), stable_logsig(x),
                              decimal=16)

    extreme_x = np.array([-100, 100], dtype=np.float)
    assert_array_almost_equal(logistic_sigmoid(extreme_x), [0, 1])
    assert_array_almost_equal(logistic_sigmoid(extreme_x, log=True), [-100, 0])
    assert_array_almost_equal(logistic_sigmoid(extreme_x),
                              stable_logsig(extreme_x),
                              decimal=16)
