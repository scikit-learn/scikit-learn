from __future__ import division, print_function, absolute_import 

from scipy.stats import hypergeom, bernoulli
import numpy as np
from numpy.testing import assert_almost_equal

def test_hypergeom_logpmf():
    # symmetries test
    # f(k,N,K,n) = f(n-k,N,N-K,n) = f(K-k,N,K,N-n) = f(k,N,n,K)
    k = 5
    N = 50
    K = 10
    n = 5
    logpmf1 = hypergeom.logpmf(k,N,K,n)
    logpmf2 = hypergeom.logpmf(n-k,N,N-K,n)
    logpmf3 = hypergeom.logpmf(K-k,N,K,N-n)
    logpmf4 = hypergeom.logpmf(k,N,n,K)
    assert_almost_equal(logpmf1, logpmf2, decimal=12)
    assert_almost_equal(logpmf1, logpmf3, decimal=12)
    assert_almost_equal(logpmf1, logpmf4, decimal=12)

    # test related distribution
    # Bernoulli distribution if n = 1
    k = 1
    N = 10
    K = 7
    n = 1
    hypergeom_logpmf = hypergeom.logpmf(k,N,K,n)
    bernoulli_logpmf = bernoulli.logpmf(k,K/N)
    assert_almost_equal(hypergeom_logpmf, bernoulli_logpmf, decimal=12)

