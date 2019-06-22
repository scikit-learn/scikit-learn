from __future__ import division, print_function, absolute_import

import scipy.special as sc
import numpy as np
from numpy.testing import assert_, assert_equal, assert_allclose


def test_zeta():
    assert_allclose(sc.zeta(2,2), np.pi**2/6 - 1, rtol=1e-12)


def test_zeta_1arg():
    assert_allclose(sc.zeta(2), np.pi**2/6, rtol=1e-12)
    assert_allclose(sc.zeta(4), np.pi**4/90, rtol=1e-12)


def test_zetac():
    assert_equal(sc.zetac(0), -1.5)
    assert_equal(sc.zetac(1.0), np.inf)
    # Expected values in the following were computed using
    # Wolfram Alpha `Zeta[x] - 1`:
    rtol = 1e-12
    assert_allclose(sc.zetac(-2.1), -0.9972705002153750, rtol=rtol)
    assert_allclose(sc.zetac(0.8), -5.437538415895550, rtol=rtol)
    assert_allclose(sc.zetac(0.9999), -10000.42279161673, rtol=rtol)
    assert_allclose(sc.zetac(9), 0.002008392826082214, rtol=rtol)
    assert_allclose(sc.zetac(50), 8.881784210930816e-16, rtol=rtol)
    assert_allclose(sc.zetac(75), 2.646977960169853e-23, rtol=rtol)


def test_zetac_negative_even():
    pts = [-2, -50, -100]
    for p in pts:
        assert_equal(sc.zetac(p), -1)


def test_zetac_inf():
    assert_equal(sc.zetac(np.inf), 0.0)
    assert_(np.isnan(sc.zetac(-np.inf)))
