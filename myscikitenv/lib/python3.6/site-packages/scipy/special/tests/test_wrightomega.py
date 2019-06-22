from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_, assert_equal

import scipy.special as sc


def test_wrightomega_nan():
    pts = [complex(np.nan, 0),
           complex(0, np.nan),
           complex(np.nan, np.nan),
           complex(np.nan, 1),
           complex(1, np.nan)]
    for p in pts:
        res = sc.wrightomega(p)
        assert_(np.isnan(res.real))
        assert_(np.isnan(res.imag))


def test_wrightomega_inf_branch():
    pts = [complex(-np.inf, np.pi/4),
           complex(-np.inf, -np.pi/4),
           complex(-np.inf, 3*np.pi/4),
           complex(-np.inf, -3*np.pi/4)]
    expected_results = [complex(0.0, 0.0),
                        complex(0.0, -0.0),
                        complex(-0.0, 0.0),
                        complex(-0.0, -0.0)]
    for p, expected in zip(pts, expected_results):
        res = sc.wrightomega(p)
        # We can't use assert_equal(res, expected) because in older versions of
        # numpy, assert_equal doesn't check the sign of the real and imaginary
        # parts when comparing complex zeros. It does check the sign when the
        # arguments are *real* scalars.
        assert_equal(res.real, expected.real)
        assert_equal(res.imag, expected.imag)


def test_wrightomega_inf():
    pts = [complex(np.inf, 10),
           complex(-np.inf, 10),
           complex(10, np.inf),
           complex(10, -np.inf)]
    for p in pts:
        assert_equal(sc.wrightomega(p), p)


def test_wrightomega_singular():
    pts = [complex(-1.0, np.pi),
           complex(-1.0, -np.pi)]
    for p in pts:
        res = sc.wrightomega(p)
        assert_equal(res, -1.0)
        assert_(np.signbit(res.imag) == False)
