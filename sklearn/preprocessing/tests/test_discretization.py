import numpy as np

from sklearn.utils import testing
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.utils.testing import assert_equal, assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal

X = [[-2, 1, -4,   -1], \
     [-1, 2, -3, -0.5], \
     [ 0, 3, -2,  0.5], \
     [ 1, 4, -1,    2]]

def _setup_dis():
    dis = KBinsDiscretizer(n_bins=3)
    return dis.fit(X)

def test_cut_points():
    dis = _setup_dis()

    expected = [[-1.,  2., -3.,  0.],
                [ 0.,  3., -2.,  1.]]

    assert_array_equal(expected, dis.cut_points_)

