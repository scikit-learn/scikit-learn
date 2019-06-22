from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose

import scipy.special as sc
from scipy.special._testutils import FuncData


def test_line():
    # Test on the line a = x where a simpler asymptotic expansion
    # (analog of DLMF 8.12.15) is available.

    def gammainc_line(x):
        c = np.array([-1/3, -1/540, 25/6048, 101/155520,
                      -3184811/3695155200, -2745493/8151736420])
        res = 0
        xfac = 1
        for ck in c:
            res -= ck*xfac
            xfac /= x
        res /= np.sqrt(2*np.pi*x)
        res += 0.5
        return res

    x = np.logspace(np.log10(25), 300, 500)
    a = x.copy()
    dataset = np.vstack((a, x, gammainc_line(x))).T

    FuncData(sc.gammainc, dataset, (0, 1), 2, rtol=1e-11).check()


def test_gammainc_roundtrip():
    a = np.logspace(-5, 10, 100)
    x = np.logspace(-5, 10, 100)

    y = sc.gammaincinv(a, sc.gammainc(a, x))
    assert_allclose(x, y, rtol=1e-10)


def test_gammaincc_roundtrip():
    a = np.logspace(-5, 10, 100)
    x = np.logspace(-5, 10, 100)

    y = sc.gammainccinv(a, sc.gammaincc(a, x))
    assert_allclose(x, y, rtol=1e-14)
