# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

# $Id$

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_raises

from ..cd import lasso_coordinate_descent_slow
from ..cd import lasso_coordinate_descent_fast
from ..cd import Lasso


def test_lasso_cd_python_cython_sanity():
    n_samples, n_features, maxit = 100, 100, 40
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    model_slow = Lasso(alpha=1)
    model_slow.learner = lasso_coordinate_descent_slow
    model_slow.fit(X, y, maxit=maxit)

    # check the convergence using the KKT condition
    assert_array_almost_equal(model_slow.gap, 0, 1e-6)

    model_fast = Lasso(alpha=1)
    model_fast.learner = lasso_coordinate_descent_fast
    model_fast.fit(X, y, maxit=maxit)

    # check t convergence using the KKT condition
    assert_array_almost_equal(model_slow.gap, 0, 1e-6)

    assert_array_almost_equal(model_slow.w, model_fast.w)
    assert_array_almost_equal(model_slow.E, model_fast.E)

