""" Test the cross_val module
"""

import numpy as np

import nose

from .. import cross_val

def test_kfold():
    # Check that errors are raise if there is not enough samples
    nose.tools.assert_raises(AssertionError, cross_val.KFold, 3, 3)
    y = [0, 0, 1, 1, 2]
    nose.tools.assert_raises(AssertionError, cross_val.StratifiedKFold, y, 3)

