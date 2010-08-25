# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np

import nose

from ..fixes import _in1d

def test_in1d():
    a = np.arange(10)
    b = a[a%2 == 0]
    nose.tools.assert_equal(_in1d(a, b).sum(), 5)



