import numpy as np

from sklearn.preprocessing import KBinsDiscretizer
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises

X = np.array([[-2, 1, -4,   -1],
              [-1, 2, -3, -0.5],
              [ 0, 3, -2,  0.5],
              [ 1, 4, -1,    2]])

def setup_discretizer():
    dis = KBinsDiscretizer(n_bins=3)
    return dis.fit(X)

def test_cut_points():
    dis = setup_discretizer()

    expected = [[-1., 2., -3., 0.],
                [ 0., 3., -2., 1.]]

    assert_array_equal(expected, dis.cut_points_)

def test_transform():
    dis = setup_discretizer()

    expected = [[0, 0, 0, 0],
                [1, 1, 1, 0],
                [2, 2, 2, 1],
                [2, 2, 2, 2]]
    assert_array_equal(expected, dis.transform(X))

def test_invalid_n_bins():
    dis = KBinsDiscretizer(n_bins=1)
    assert_raises(ValueError, dis.fit_transform, X)

def test_invalid_n_features():
    dis = setup_discretizer()
    bad_X = np.arange(25).reshape(5, -1)
    assert_raises(ValueError, dis.transform, bad_X)

def setup_discretizer_categorical():
    categorical = [1]
    dis = KBinsDiscretizer(n_bins=3, categorical_features=categorical)
    return dis.fit(X)

def test_categorical_cut_points():
    dis = setup_discretizer_categorical()
    expected = [[-1., np.nan, -3., 0.],
                [ 0., np.nan, -2., 1.]]

    assert_array_equal(expected, dis.cut_points_)

def test_categorical_transform():
    # Feature at col_idx=1 should not change

    dis = setup_discretizer_categorical()

    expected = [[0., 1, 0., 0.],
                [1., 2, 1., 0.],
                [2., 3, 2., 1.],
                [2., 4, 2., 2.]]

    assert_array_equal(expected, dis.transform(X))

def test_categorical_invalid():
    invalid_categorical = [
        [1, 1],  # Duplicate column
        [-1],  # Invalid index
        [4],  # Invalid index
        ['a'], # Not an integer index
        [4.5], # Not an integer index
        [[1, 2], [3, 4]]  # Invalid shape
    ]

    for invalid in invalid_categorical:
        dis = KBinsDiscretizer(categorical_features=invalid)
        assert_raises(ValueError, dis.fit, X)

def test_transform_1d_behavior():
    pass  # To be determined