# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id: test_cd.py 450 2010-03-03 14:21:06Z twigster $

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_raises

from ..ridge import Ridge

def test_ridge():
    alpha = 1.0

    # With more samples than features
    nsamples, nfeatures = 10, 5
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)

    # With more features than samples
    nsamples, nfeatures = 5, 10
    np.random.seed(0)
    y = np.random.randn(nsamples)
    X = np.random.randn(nsamples, nfeatures)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
