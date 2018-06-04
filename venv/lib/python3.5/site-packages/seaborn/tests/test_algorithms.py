import numpy as np
from scipy import stats
from ..external.six.moves import range

import numpy.testing as npt
from numpy.testing import assert_array_equal
import nose.tools
from nose.tools import assert_equal, raises

from .. import algorithms as algo

rs = np.random.RandomState(sum(map(ord, "test_algorithms")))
a_norm = rs.randn(100)


def test_bootstrap():
    """Test that bootstrapping gives the right answer in dumb cases."""
    a_ones = np.ones(10)
    n_boot = 5
    out1 = algo.bootstrap(a_ones, n_boot=n_boot)
    assert_array_equal(out1, np.ones(n_boot))
    out2 = algo.bootstrap(a_ones, n_boot=n_boot, func=np.median)
    assert_array_equal(out2, np.ones(n_boot))


def test_bootstrap_length():
    """Test that we get a bootstrap array of the right shape."""
    out = algo.bootstrap(a_norm)
    assert_equal(len(out), 10000)

    n_boot = 100
    out = algo.bootstrap(a_norm, n_boot=n_boot)
    assert_equal(len(out), n_boot)


def test_bootstrap_range():
    """Test that boostrapping a random array stays within the right range."""
    min, max = a_norm.min(), a_norm.max()
    out = algo.bootstrap(a_norm)
    nose.tools.assert_less(min, out.min())
    nose.tools.assert_greater_equal(max, out.max())


def test_bootstrap_multiarg():
    """Test that bootstrap works with multiple input arrays."""
    x = np.vstack([[1, 10] for i in range(10)])
    y = np.vstack([[5, 5] for i in range(10)])

    test_func = lambda x, y: np.vstack((x, y)).max(axis=0)
    out_actual = algo.bootstrap(x, y, n_boot=2, func=test_func)
    out_wanted = np.array([[5, 10], [5, 10]])
    assert_array_equal(out_actual, out_wanted)


def test_bootstrap_axis():
    """Test axis kwarg to bootstrap function."""
    x = rs.randn(10, 20)
    n_boot = 100
    out_default = algo.bootstrap(x, n_boot=n_boot)
    assert_equal(out_default.shape, (n_boot,))
    out_axis = algo.bootstrap(x, n_boot=n_boot, axis=0)
    assert_equal(out_axis.shape, (n_boot, 20))


def test_bootstrap_random_seed():
    """Test that we can get reproducible resamples by seeding the RNG."""
    data = rs.randn(50)
    seed = 42
    boots1 = algo.bootstrap(data, random_seed=seed)
    boots2 = algo.bootstrap(data, random_seed=seed)
    assert_array_equal(boots1, boots2)


def test_smooth_bootstrap():
    """Test smooth bootstrap."""
    x = rs.randn(15)
    n_boot = 100
    out_smooth = algo.bootstrap(x, n_boot=n_boot,
                                smooth=True, func=np.median)
    assert(not np.median(out_smooth) in x)


def test_bootstrap_ols():
    """Test bootstrap of OLS model fit."""
    ols_fit = lambda X, y: np.dot(np.dot(np.linalg.inv(
                                  np.dot(X.T, X)), X.T), y)
    X = np.column_stack((rs.randn(50, 4), np.ones(50)))
    w = [2, 4, 0, 3, 5]
    y_noisy = np.dot(X, w) + rs.randn(50) * 20
    y_lownoise = np.dot(X, w) + rs.randn(50)

    n_boot = 500
    w_boot_noisy = algo.bootstrap(X, y_noisy,
                                  n_boot=n_boot,
                                  func=ols_fit)
    w_boot_lownoise = algo.bootstrap(X, y_lownoise,
                                     n_boot=n_boot,
                                     func=ols_fit)

    assert_equal(w_boot_noisy.shape, (n_boot, 5))
    assert_equal(w_boot_lownoise.shape, (n_boot, 5))
    nose.tools.assert_greater(w_boot_noisy.std(),
                              w_boot_lownoise.std())


def test_bootstrap_units():
    """Test that results make sense when passing unit IDs to bootstrap."""
    data = rs.randn(50)
    ids = np.repeat(range(10), 5)
    bwerr = rs.normal(0, 2, 10)
    bwerr = bwerr[ids]
    data_rm = data + bwerr
    seed = 77

    boots_orig = algo.bootstrap(data_rm, random_seed=seed)
    boots_rm = algo.bootstrap(data_rm, units=ids, random_seed=seed)
    nose.tools.assert_greater(boots_rm.std(), boots_orig.std())


@raises(ValueError)
def test_bootstrap_arglength():
    """Test that different length args raise ValueError."""
    algo.bootstrap(np.arange(5), np.arange(10))


@raises(TypeError)
def test_bootstrap_noncallable():
    """Test that we get a TypeError with noncallable algo.unc."""
    non_func = "mean"
    algo.bootstrap(a_norm, 100, non_func)


def test_randomize_corrmat():
    """Test the correctness of the correlation matrix p values."""
    a = rs.randn(30)
    b = a + rs.rand(30) * 3
    c = rs.randn(30)
    d = [a, b, c]

    p_mat, dist = algo.randomize_corrmat(d, tail="upper", corrected=False,
                                         return_dist=True)
    nose.tools.assert_greater(p_mat[2, 0], p_mat[1, 0])

    corrmat = np.corrcoef(d)
    pctile = 100 - stats.percentileofscore(dist[2, 1], corrmat[2, 1])
    nose.tools.assert_almost_equal(p_mat[2, 1] * 100, pctile)

    d[1] = -a + rs.rand(30)
    p_mat = algo.randomize_corrmat(d)
    nose.tools.assert_greater(0.05, p_mat[1, 0])


def test_randomize_corrmat_dist():
    """Test that the distribution looks right."""
    a = rs.randn(3, 20)
    for n_i in [5, 10]:
        p_mat, dist = algo.randomize_corrmat(a, n_iter=n_i, return_dist=True)
        assert_equal(n_i, dist.shape[-1])

    p_mat, dist = algo.randomize_corrmat(a, n_iter=10000, return_dist=True)

    diag_mean = dist[0, 0].mean()
    assert_equal(diag_mean, 1)

    off_diag_mean = dist[0, 1].mean()
    nose.tools.assert_greater(0.05, off_diag_mean)


def test_randomize_corrmat_correction():
    """Test that FWE correction works."""
    a = rs.randn(3, 20)
    p_mat = algo.randomize_corrmat(a, "upper", False)
    p_mat_corr = algo.randomize_corrmat(a, "upper", True)
    triu = np.triu_indices(3, 1)
    npt.assert_array_less(p_mat[triu], p_mat_corr[triu])


def test_randimoize_corrmat_tails():
    """Test that the tail argument works."""
    a = rs.randn(30)
    b = a + rs.rand(30) * 8
    c = rs.randn(30)
    d = [a, b, c]

    p_mat_b = algo.randomize_corrmat(d, "both", False, random_seed=0)
    p_mat_u = algo.randomize_corrmat(d, "upper", False, random_seed=0)
    p_mat_l = algo.randomize_corrmat(d, "lower", False, random_seed=0)
    assert_equal(p_mat_b[0, 1], p_mat_u[0, 1] * 2)
    assert_equal(p_mat_l[0, 1], 1 - p_mat_u[0, 1])


def test_randomise_corrmat_seed():
    """Test that we can seed the corrmat randomization."""
    a = rs.randn(3, 20)
    _, dist1 = algo.randomize_corrmat(a, random_seed=0, return_dist=True)
    _, dist2 = algo.randomize_corrmat(a, random_seed=0, return_dist=True)
    assert_array_equal(dist1, dist2)


@raises(ValueError)
def test_randomize_corrmat_tail_error():
    """Test that we are strict about tail paramete."""
    a = rs.randn(3, 30)
    algo.randomize_corrmat(a, "hello")
