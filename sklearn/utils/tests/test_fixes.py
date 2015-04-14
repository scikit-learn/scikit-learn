# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import numpy as np

from nose.tools import assert_equal
from numpy.testing import (assert_almost_equal,
                           assert_array_almost_equal)

from sklearn.utils.fixes import divide, expit


def test_expit():
    # Check numerical stability of expit (logistic function).

    # Simulate our previous Cython implementation, based on
    #http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression
    assert_almost_equal(expit(1000.), 1. / (1. + np.exp(-1000.)), decimal=16)
    assert_almost_equal(expit(-1000.), np.exp(-1000.) / (1. + np.exp(-1000.)),
                        decimal=16)

    x = np.arange(10)
    out = np.zeros_like(x, dtype=np.float32)
    assert_array_almost_equal(expit(x), expit(x, out=out))


def test_divide():
    assert_equal(divide(.6, 1), .600000000000)
