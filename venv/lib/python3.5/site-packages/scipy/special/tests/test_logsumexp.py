from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_almost_equal, assert_equal,
                           assert_array_almost_equal, assert_)

from scipy.special import logsumexp


def test_logsumexp():
    # Test whether logsumexp() function correctly handles large inputs.
    a = np.arange(200)
    desired = np.log(np.sum(np.exp(a)))
    assert_almost_equal(logsumexp(a), desired)

    # Now test with large numbers
    b = [1000, 1000]
    desired = 1000.0 + np.log(2.0)
    assert_almost_equal(logsumexp(b), desired)

    n = 1000
    b = np.ones(n) * 10000
    desired = 10000.0 + np.log(n)
    assert_almost_equal(logsumexp(b), desired)

    x = np.array([1e-40] * 1000000)
    logx = np.log(x)

    X = np.vstack([x, x])
    logX = np.vstack([logx, logx])
    assert_array_almost_equal(np.exp(logsumexp(logX)), X.sum())
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=0)), X.sum(axis=0))
    assert_array_almost_equal(np.exp(logsumexp(logX, axis=1)), X.sum(axis=1))

    # Handling special values properly
    assert_equal(logsumexp(np.inf), np.inf)
    assert_equal(logsumexp(-np.inf), -np.inf)
    assert_equal(logsumexp(np.nan), np.nan)
    assert_equal(logsumexp([-np.inf, -np.inf]), -np.inf)

    # Handling an array with different magnitudes on the axes
    assert_array_almost_equal(logsumexp([[1e10, 1e-10],
                                         [-1e10, -np.inf]], axis=-1),
                              [1e10, -1e10])

    # Test keeping dimensions
    assert_array_almost_equal(logsumexp([[1e10, 1e-10],
                                         [-1e10, -np.inf]],
                                        axis=-1,
                                        keepdims=True),
                              [[1e10], [-1e10]])

    # Test multiple axes
    assert_array_almost_equal(logsumexp([[1e10, 1e-10],
                                         [-1e10, -np.inf]],
                                        axis=(-1,-2)),
                              1e10)


def test_logsumexp_b():
    a = np.arange(200)
    b = np.arange(200, 0, -1)
    desired = np.log(np.sum(b*np.exp(a)))
    assert_almost_equal(logsumexp(a, b=b), desired)

    a = [1000, 1000]
    b = [1.2, 1.2]
    desired = 1000 + np.log(2 * 1.2)
    assert_almost_equal(logsumexp(a, b=b), desired)

    x = np.array([1e-40] * 100000)
    b = np.linspace(1, 1000, 100000)
    logx = np.log(x)

    X = np.vstack((x, x))
    logX = np.vstack((logx, logx))
    B = np.vstack((b, b))
    assert_array_almost_equal(np.exp(logsumexp(logX, b=B)), (B * X).sum())
    assert_array_almost_equal(np.exp(logsumexp(logX, b=B, axis=0)),
                                (B * X).sum(axis=0))
    assert_array_almost_equal(np.exp(logsumexp(logX, b=B, axis=1)),
                                (B * X).sum(axis=1))


def test_logsumexp_sign():
    a = [1,1,1]
    b = [1,-1,-1]

    r, s = logsumexp(a, b=b, return_sign=True)
    assert_almost_equal(r,1)
    assert_equal(s,-1)


def test_logsumexp_sign_zero():
    a = [1,1]
    b = [1,-1]

    r, s = logsumexp(a, b=b, return_sign=True)
    assert_(not np.isfinite(r))
    assert_(not np.isnan(r))
    assert_(r < 0)
    assert_equal(s,0)


def test_logsumexp_sign_shape():
    a = np.ones((1,2,3,4))
    b = np.ones_like(a)

    r, s = logsumexp(a, axis=2, b=b, return_sign=True)

    assert_equal(r.shape, s.shape)
    assert_equal(r.shape, (1,2,4))

    r, s = logsumexp(a, axis=(1,3), b=b, return_sign=True)

    assert_equal(r.shape, s.shape)
    assert_equal(r.shape, (1,3))


def test_logsumexp_shape():
    a = np.ones((1, 2, 3, 4))
    b = np.ones_like(a)

    r = logsumexp(a, axis=2, b=b)
    assert_equal(r.shape, (1, 2, 4))

    r = logsumexp(a, axis=(1, 3), b=b)
    assert_equal(r.shape, (1, 3))


def test_logsumexp_b_zero():
    a = [1,10000]
    b = [1,0]

    assert_almost_equal(logsumexp(a, b=b), 1)


def test_logsumexp_b_shape():
    a = np.zeros((4,1,2,1))
    b = np.ones((3,1,5))

    logsumexp(a, b=b)

