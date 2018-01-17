import numpy as np

from sklearn.utils.testing import *

from sklearn.preprocessing.smooth import rectangular_smoother
from sklearn.preprocessing.smooth import RectangularSmoother


def test_rectangular_smoothing_function():
    sample_array = np.array([[1], [2], [3], [4], [5], [6], [7]])
    expected_result_for_constant_0 = np.array([[1.2], [2],    [3], [4], [5], [4.4], [3.6]])
    expected_result_for_nearest = np.array([[1.6], [2.2], [3], [4], [5], [5.8], [6.4]])
    expected_result_for_wrap = np.array([[3.8], [3.4],  [3], [4], [5], [4.6], [4.2]])
    expected_result_for_mirror = np.array([[2.2], [2.4],  [3], [4], [5], [5.6], [5.8]])
    expected_result_for_interp = np.array([[3],   [4],    [5]])

    assert_array_almost_equal(rectangular_smoother(sample_array, size=5, mode='constant', cval=0.0),
                               expected_result_for_constant_0)

    assert_array_almost_equal(rectangular_smoother(sample_array, size=5, mode='nearest'),
                              expected_result_for_nearest)

    assert_array_almost_equal(rectangular_smoother(sample_array, size=5, mode='wrap'),
                              expected_result_for_wrap)

    assert_array_almost_equal(rectangular_smoother(sample_array, size=5, mode='mirror'),
                              expected_result_for_mirror)

    assert_array_almost_equal(rectangular_smoother(sample_array, size=5, mode='interp'),
                              expected_result_for_interp)

def test_rectangular_smoother_class():
    sample_array = np.array([[1], [2], [3], [4], [5], [6], [7]])

    smoother = RectangularSmoother(size=5, mode='constant', cval=0.0)
    smoother.fit(sample_array)
    assert_almost_equal(rectangular_smoother(sample_array, size=5, mode='constant', cval=0.0),
                        smoother.transform(sample_array))

    smoother = RectangularSmoother(size=5, mode='nearest')
    smoother.fit(sample_array)
    assert_almost_equal(rectangular_smoother(sample_array, size=5, mode='nearest'),
                        smoother.transform(sample_array))

    smoother = RectangularSmoother(size=5, mode='wrap')
    smoother.fit(sample_array)
    assert_almost_equal(rectangular_smoother(sample_array, size=5, mode='wrap'),
                        smoother.transform(sample_array))

    smoother = RectangularSmoother(size=5, mode='mirror')
    smoother.fit(sample_array)
    assert_almost_equal(rectangular_smoother(sample_array, size=5, mode='mirror'),
                        smoother.transform(sample_array))

    smoother = RectangularSmoother(size=5, mode='interp')
    smoother.fit(sample_array)
    assert_almost_equal(rectangular_smoother(sample_array, size=5, mode='interp'),
                        smoother.transform(sample_array))