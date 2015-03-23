import numpy as np

from sklearn.utils.class_weight import compute_class_weight
from sklearn.utils.class_weight import compute_sample_weight

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal


def test_compute_class_weight():
    # Test (and demo) compute_class_weight.
    y = np.asarray([2, 2, 2, 3, 3, 4])
    classes = np.unique(y)
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_true(cw[0] < cw[1] < cw[2])


def test_compute_class_weight_not_present():
    # Raise error when y does not contain all class labels
    classes = np.arange(4)
    y = np.asarray([0, 0, 0, 1, 1, 2])
    assert_raises(ValueError, compute_class_weight, "auto", classes, y)


def test_compute_class_weight_auto_negative():
    # Test compute_class_weight when labels are negative
    # Test with balanced class labels.
    classes = np.array([-2, -1, 0])
    y = np.asarray([-1, -1, 0, 0, -2, -2])
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_equal(len(cw), len(classes))
    assert_array_almost_equal(cw, np.array([1., 1., 1.]))

    # Test with unbalanced class labels.
    y = np.asarray([-1, 0, 0, -2, -2, -2])
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_equal(len(cw), len(classes))
    assert_array_almost_equal(cw, np.array([0.545, 1.636, 0.818]), decimal=3)


def test_compute_class_weight_auto_unordered():
    # Test compute_class_weight when classes are unordered
    classes = np.array([1, 0, 3])
    y = np.asarray([1, 0, 0, 3, 3, 3])
    cw = compute_class_weight("auto", classes, y)
    assert_almost_equal(cw.sum(), classes.shape)
    assert_equal(len(cw), len(classes))
    assert_array_almost_equal(cw, np.array([1.636, 0.818, 0.545]), decimal=3)


def test_compute_sample_weight():
    # Test (and demo) compute_sample_weight.
    # Test with balanced classes
    y = np.asarray([1, 1, 1, 2, 2, 2])
    sample_weight = compute_sample_weight("auto", y)
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1.])

    # Test with user-defined weights
    sample_weight = compute_sample_weight({1: 2, 2: 1}, y)
    assert_array_almost_equal(sample_weight, [2., 2., 2., 1., 1., 1.])

    # Test with column vector of balanced classes
    y = np.asarray([[1], [1], [1], [2], [2], [2]])
    sample_weight = compute_sample_weight("auto", y)
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1.])

    # Test with unbalanced classes
    y = np.asarray([1, 1, 1, 2, 2, 2, 3])
    sample_weight = compute_sample_weight("auto", y)
    expected = np.asarray([.6, .6, .6, .6, .6, .6, 1.8])
    assert_array_almost_equal(sample_weight, expected)

    # Test with `None` weights
    sample_weight = compute_sample_weight(None, y)
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1., 1.])

    # Test with multi-output of balanced classes
    y = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1]])
    sample_weight = compute_sample_weight("auto", y)
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1.])

    # Test with multi-output with user-defined weights
    y = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1]])
    sample_weight = compute_sample_weight([{1: 2, 2: 1}, {0: 1, 1: 2}], y)
    assert_array_almost_equal(sample_weight, [2., 2., 2., 2., 2., 2.])

    # Test with multi-output of unbalanced classes
    y = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1], [3, -1]])
    sample_weight = compute_sample_weight("auto", y)
    assert_array_almost_equal(sample_weight, expected ** 2)


def test_compute_sample_weight_with_subsample():
    # Test compute_sample_weight with subsamples specified.
    # Test with balanced classes and all samples present
    y = np.asarray([1, 1, 1, 2, 2, 2])
    sample_weight = compute_sample_weight("auto", y, range(6))
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1.])

    # Test with column vector of balanced classes and all samples present
    y = np.asarray([[1], [1], [1], [2], [2], [2]])
    sample_weight = compute_sample_weight("auto", y, range(6))
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1.])

    # Test with a subsample
    y = np.asarray([1, 1, 1, 2, 2, 2])
    sample_weight = compute_sample_weight("auto", y, range(4))
    assert_array_almost_equal(sample_weight, [.5, .5, .5, 1.5, 1.5, 1.5])

    # Test with a bootstrap subsample
    y = np.asarray([1, 1, 1, 2, 2, 2])
    sample_weight = compute_sample_weight("auto", y, [0, 1, 1, 2, 2, 3])
    expected = np.asarray([1/3., 1/3., 1/3., 5/3., 5/3., 5/3.])
    assert_array_almost_equal(sample_weight, expected)

    # Test with a bootstrap subsample for multi-output
    y = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1]])
    sample_weight = compute_sample_weight("auto", y, [0, 1, 1, 2, 2, 3])
    assert_array_almost_equal(sample_weight, expected ** 2)

    # Test with a missing class
    y = np.asarray([1, 1, 1, 2, 2, 2, 3])
    sample_weight = compute_sample_weight("auto", y, range(6))
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1., 0.])

    # Test with a missing class for multi-output
    y = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1], [2, 2]])
    sample_weight = compute_sample_weight("auto", y, range(6))
    assert_array_almost_equal(sample_weight, [1., 1., 1., 1., 1., 1., 0.])


def test_compute_sample_weight_errors():
    # Test compute_sample_weight raises errors expected.
    # Invalid preset string
    y = np.asarray([1, 1, 1, 2, 2, 2])
    y_ = np.asarray([[1, 0], [1, 0], [1, 0], [2, 1], [2, 1], [2, 1]])
    assert_raises(ValueError, compute_sample_weight, "ni", y)
    assert_raises(ValueError, compute_sample_weight, "ni", y, range(4))
    assert_raises(ValueError, compute_sample_weight, "ni", y_)
    assert_raises(ValueError, compute_sample_weight, "ni", y_, range(4))

    # Not "auto" for subsample
    assert_raises(ValueError,
                  compute_sample_weight, {1: 2, 2: 1}, y, range(4))

    # Not a list or preset for multi-output
    assert_raises(ValueError, compute_sample_weight, {1: 2, 2: 1}, y_)

    # Incorrect length list for multi-output
    assert_raises(ValueError, compute_sample_weight, [{1: 2, 2: 1}], y_)
