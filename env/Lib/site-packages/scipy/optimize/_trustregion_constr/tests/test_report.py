import numpy as np
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_raises, assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns,
                           dec)
from scipy.optimize import minimize, Bounds

def test_gh10880():
    # checks that verbose reporting works with trust-constr
    bnds = Bounds(1, 2)
    opts = {'maxiter': 1000, 'verbose': 2}
    res = minimize(lambda x: x**2, x0=2., method='trust-constr', bounds=bnds, options=opts)

    opts = {'maxiter': 1000, 'verbose': 3}
    res = minimize(lambda x: x**2, x0=2., method='trust-constr', bounds=bnds, options=opts)
