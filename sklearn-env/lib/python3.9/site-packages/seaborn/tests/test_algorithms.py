import numpy as np
import numpy.random as npr

import pytest
from numpy.testing import assert_array_equal
from distutils.version import LooseVersion

from .. import algorithms as algo


@pytest.fixture
def random():
    np.random.seed(sum(map(ord, "test_algorithms")))


def test_bootstrap(random):
    """Test that bootstrapping gives the right answer in dumb cases."""
    a_ones = np.ones(10)
    n_boot = 5
    out1 = algo.bootstrap(a_ones, n_boot=n_boot)
    assert_array_equal(out1, np.ones(n_boot))
    out2 = algo.bootstrap(a_ones, n_boot=n_boot, func=np.median)
    assert_array_equal(out2, np.ones(n_boot))


def test_bootstrap_length(random):
    """Test that we get a bootstrap array of the right shape."""
    a_norm = np.random.randn(1000)
    out = algo.bootstrap(a_norm)
    assert len(out) == 10000

    n_boot = 100
    out = algo.bootstrap(a_norm, n_boot=n_boot)
    assert len(out) == n_boot


def test_bootstrap_range(random):
    """Test that boostrapping a random array stays within the right range."""
    a_norm = np.random.randn(1000)
    amin, amax = a_norm.min(), a_norm.max()
    out = algo.bootstrap(a_norm)
    assert amin <= out.min()
    assert amax >= out.max()


def test_bootstrap_multiarg(random):
    """Test that bootstrap works with multiple input arrays."""
    x = np.vstack([[1, 10] for i in range(10)])
    y = np.vstack([[5, 5] for i in range(10)])

    def f(x, y):
        return np.vstack((x, y)).max(axis=0)

    out_actual = algo.bootstrap(x, y, n_boot=2, func=f)
    out_wanted = np.array([[5, 10], [5, 10]])
    assert_array_equal(out_actual, out_wanted)


def test_bootstrap_axis(random):
    """Test axis kwarg to bootstrap function."""
    x = np.random.randn(10, 20)
    n_boot = 100

    out_default = algo.bootstrap(x, n_boot=n_boot)
    assert out_default.shape == (n_boot,)

    out_axis = algo.bootstrap(x, n_boot=n_boot, axis=0)
    assert out_axis.shape, (n_boot, x.shape[1])


def test_bootstrap_seed(random):
    """Test that we can get reproducible resamples by seeding the RNG."""
    data = np.random.randn(50)
    seed = 42
    boots1 = algo.bootstrap(data, seed=seed)
    boots2 = algo.bootstrap(data, seed=seed)
    assert_array_equal(boots1, boots2)


def test_bootstrap_ols(random):
    """Test bootstrap of OLS model fit."""
    def ols_fit(X, y):
        XtXinv = np.linalg.inv(np.dot(X.T, X))
        return XtXinv.dot(X.T).dot(y)

    X = np.column_stack((np.random.randn(50, 4), np.ones(50)))
    w = [2, 4, 0, 3, 5]
    y_noisy = np.dot(X, w) + np.random.randn(50) * 20
    y_lownoise = np.dot(X, w) + np.random.randn(50)

    n_boot = 500
    w_boot_noisy = algo.bootstrap(X, y_noisy,
                                  n_boot=n_boot,
                                  func=ols_fit)
    w_boot_lownoise = algo.bootstrap(X, y_lownoise,
                                     n_boot=n_boot,
                                     func=ols_fit)

    assert w_boot_noisy.shape == (n_boot, 5)
    assert w_boot_lownoise.shape == (n_boot, 5)
    assert w_boot_noisy.std() > w_boot_lownoise.std()


def test_bootstrap_units(random):
    """Test that results make sense when passing unit IDs to bootstrap."""
    data = np.random.randn(50)
    ids = np.repeat(range(10), 5)
    bwerr = np.random.normal(0, 2, 10)
    bwerr = bwerr[ids]
    data_rm = data + bwerr
    seed = 77

    boots_orig = algo.bootstrap(data_rm, seed=seed)
    boots_rm = algo.bootstrap(data_rm, units=ids, seed=seed)
    assert boots_rm.std() > boots_orig.std()


def test_bootstrap_arglength():
    """Test that different length args raise ValueError."""
    with pytest.raises(ValueError):
        algo.bootstrap(np.arange(5), np.arange(10))


def test_bootstrap_string_func():
    """Test that named numpy methods are the same as the numpy function."""
    x = np.random.randn(100)

    res_a = algo.bootstrap(x, func="mean", seed=0)
    res_b = algo.bootstrap(x, func=np.mean, seed=0)
    assert np.array_equal(res_a, res_b)

    res_a = algo.bootstrap(x, func="std", seed=0)
    res_b = algo.bootstrap(x, func=np.std, seed=0)
    assert np.array_equal(res_a, res_b)

    with pytest.raises(AttributeError):
        algo.bootstrap(x, func="not_a_method_name")


def test_bootstrap_reproducibility(random):
    """Test that bootstrapping uses the internal random state."""
    data = np.random.randn(50)
    boots1 = algo.bootstrap(data, seed=100)
    boots2 = algo.bootstrap(data, seed=100)
    assert_array_equal(boots1, boots2)

    with pytest.warns(UserWarning):
        # Deprecatd, remove when removing random_seed
        boots1 = algo.bootstrap(data, random_seed=100)
        boots2 = algo.bootstrap(data, random_seed=100)
        assert_array_equal(boots1, boots2)


@pytest.mark.skipif(LooseVersion(np.__version__) < "1.17",
                    reason="Tests new numpy random functionality")
def test_seed_new():

    # Can't use pytest parametrize because tests will fail where the new
    # Generator object and related function are not defined

    test_bank = [
        (None, None, npr.Generator, False),
        (npr.RandomState(0), npr.RandomState(0), npr.RandomState, True),
        (npr.RandomState(0), npr.RandomState(1), npr.RandomState, False),
        (npr.default_rng(1), npr.default_rng(1), npr.Generator, True),
        (npr.default_rng(1), npr.default_rng(2), npr.Generator, False),
        (npr.SeedSequence(10), npr.SeedSequence(10), npr.Generator, True),
        (npr.SeedSequence(10), npr.SeedSequence(20), npr.Generator, False),
        (100, 100, npr.Generator, True),
        (100, 200, npr.Generator, False),
    ]
    for seed1, seed2, rng_class, match in test_bank:
        rng1 = algo._handle_random_seed(seed1)
        rng2 = algo._handle_random_seed(seed2)
        assert isinstance(rng1, rng_class)
        assert isinstance(rng2, rng_class)
        assert (rng1.uniform() == rng2.uniform()) == match


@pytest.mark.skipif(LooseVersion(np.__version__) >= "1.17",
                    reason="Tests old numpy random functionality")
@pytest.mark.parametrize("seed1, seed2, match", [
    (None, None, False),
    (npr.RandomState(0), npr.RandomState(0), True),
    (npr.RandomState(0), npr.RandomState(1), False),
    (100, 100, True),
    (100, 200, False),
])
def test_seed_old(seed1, seed2, match):
    rng1 = algo._handle_random_seed(seed1)
    rng2 = algo._handle_random_seed(seed2)
    assert isinstance(rng1, np.random.RandomState)
    assert isinstance(rng2, np.random.RandomState)
    assert (rng1.uniform() == rng2.uniform()) == match


@pytest.mark.skipif(LooseVersion(np.__version__) >= "1.17",
                    reason="Tests old numpy random functionality")
def test_bad_seed_old():

    with pytest.raises(ValueError):
        algo._handle_random_seed("not_a_random_seed")
