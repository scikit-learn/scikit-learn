# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause

import numpy as np

from nose.tools import assert_equal
from numpy.testing import assert_almost_equal, assert_array_equal

from ..fixes import _in1d, _copysign, divide, expit


def test_in1d():
    a = np.arange(10)
    b = a[a % 2 == 0]
    assert_equal(_in1d(a, b).sum(), 5)


def test_expit():
    """Check numerical stability of expit (logistic function)."""

    # Simulate our previous Cython implementation, based on
    #http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression
    assert_almost_equal(expit(100.), 1. / (1. + np.exp(-100.)), decimal=16)
    assert_almost_equal(expit(-100.), np.exp(-100.) / (1. + np.exp(-100.)),
                        decimal=16)


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)


def test_copysign():
    a = np.array([-1, 1, -1])
    b = np.array([1, -1, 1])

    assert_array_equal(_copysign(a, b), b)
    assert_array_equal(_copysign(b, a), a)
