# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np

from nose.tools import assert_equal
from numpy.testing import assert_array_equal

from ..fixes import _in1d, _copysign


def test_in1d():
    a = np.arange(10)
    b = a[a % 2 == 0]
    assert_equal(_in1d(a, b).sum(), 5)


def test_copysign():
    a = np.array([-1, 1, -1])
    b = np.array([1, -1, 1])

    assert_array_equal(_copysign(a, b), b)
    assert_array_equal(_copysign(b, a), a)
