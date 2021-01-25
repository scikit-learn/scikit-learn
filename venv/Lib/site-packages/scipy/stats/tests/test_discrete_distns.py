from scipy.stats import betabinom, hypergeom, nhypergeom, bernoulli, boltzmann
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal, assert_allclose


def test_hypergeom_logpmf():
    # symmetries test
    # f(k,N,K,n) = f(n-k,N,N-K,n) = f(K-k,N,K,N-n) = f(k,N,n,K)
    k = 5
    N = 50
    K = 10
    n = 5
    logpmf1 = hypergeom.logpmf(k, N, K, n)
    logpmf2 = hypergeom.logpmf(n - k, N, N - K, n)
    logpmf3 = hypergeom.logpmf(K - k, N, K, N - n)
    logpmf4 = hypergeom.logpmf(k, N, n, K)
    assert_almost_equal(logpmf1, logpmf2, decimal=12)
    assert_almost_equal(logpmf1, logpmf3, decimal=12)
    assert_almost_equal(logpmf1, logpmf4, decimal=12)

    # test related distribution
    # Bernoulli distribution if n = 1
    k = 1
    N = 10
    K = 7
    n = 1
    hypergeom_logpmf = hypergeom.logpmf(k, N, K, n)
    bernoulli_logpmf = bernoulli.logpmf(k, K/N)
    assert_almost_equal(hypergeom_logpmf, bernoulli_logpmf, decimal=12)


def test_nhypergeom_pmf():
    # test with hypergeom
    M, n, r = 45, 13, 8
    k = 6
    NHG = nhypergeom.pmf(k, M, n, r)
    HG = hypergeom.pmf(k, M, n, k+r-1) * (M - n - (r-1)) / (M - (k+r-1))
    assert_allclose(HG, NHG, rtol=1e-10)


def test_nhypergeom_pmfcdf():
    # test pmf and cdf with arbitrary values.
    M = 8
    n = 3
    r = 4
    support = np.arange(n+1)
    pmf = nhypergeom.pmf(support, M, n, r)
    cdf = nhypergeom.cdf(support, M, n, r)
    assert_allclose(pmf, [1/14, 3/14, 5/14, 5/14], rtol=1e-13)
    assert_allclose(cdf, [1/14, 4/14, 9/14, 1.0], rtol=1e-13)


def test_nhypergeom_r0():
    # test with `r = 0`.
    M = 10
    n = 3
    r = 0
    pmf = nhypergeom.pmf([[0, 1, 2, 0], [1, 2, 0, 3]], M, n, r)
    assert_allclose(pmf, [[1, 0, 0, 1], [0, 0, 1, 0]], rtol=1e-13)


def test_boltzmann_upper_bound():
    k = np.arange(-3, 5)

    N = 1
    p = boltzmann.pmf(k, 0.123, N)
    expected = k == 0
    assert_equal(p, expected)

    lam = np.log(2)
    N = 3
    p = boltzmann.pmf(k, lam, N)
    expected = [0, 0, 0, 4/7, 2/7, 1/7, 0, 0]
    assert_allclose(p, expected, rtol=1e-13)

    c = boltzmann.cdf(k, lam, N)
    expected = [0, 0, 0, 4/7, 6/7, 1, 1, 1]
    assert_allclose(c, expected, rtol=1e-13)


def test_betabinom_a_and_b_unity():
    # test limiting case that betabinom(n, 1, 1) is a discrete uniform
    # distribution from 0 to n
    n = 20
    k = np.arange(n + 1)
    p = betabinom(n, 1, 1).pmf(k)
    expected = np.repeat(1 / (n + 1), n + 1)
    assert_almost_equal(p, expected)


def test_betabinom_bernoulli():
    # test limiting case that betabinom(1, a, b) = bernoulli(a / (a + b))
    a = 2.3
    b = 0.63
    k = np.arange(2)
    p = betabinom(1, a, b).pmf(k)
    expected = bernoulli(a / (a + b)).pmf(k)
    assert_almost_equal(p, expected)
