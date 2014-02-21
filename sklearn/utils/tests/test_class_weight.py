import numpy as np

from sklearn.utils.class_weight import compute_class_weight
from sklearn.utils.fixes import unique

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal


def test_compute_class_weight():
    """Test (and demo) compute_class_weight."""
    classes, y = unique(np.asarray([2, 2, 2, 3, 3, 4]), return_inverse=True)
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_true(cw[0] < cw[1] < cw[2])


def test_compute_class_weight_not_present():
    """Test compute_class_weight in case y doesn't contain all classes."""
    classes = np.arange(4)
    y = np.asarray([0, 0, 0, 1, 1, 2])
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_equal(len(cw), len(classes))
    assert_true(cw[0] < cw[1] < cw[2] <= cw[3])


def test_compute_class_weight_auto_negative():
    """Test compute_class_weight when labels are negative"""
    classes = -np.arange(3)
    y = np.asarray([-1, -1, 0, 0, -2, -2])
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_equal(len(cw), len(classes))
    assert_array_almost_equal(cw, np.array([1., 1., 1.]))
