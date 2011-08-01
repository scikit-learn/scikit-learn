from scikits.learn.metrics.cluster import homogeneity_score
from scikits.learn.metrics.cluster import completeness_score
from scikits.learn.metrics.cluster import v_measure_score
from scikits.learn.metrics.cluster import v_measure_score
from scikits.learn.metrics.cluster import homogeneity_completeness_v_measure

from nose.tools import assert_almost_equal


def assert_raise_message(exception, message, callable, *args, **kwargs):
    """Helper function to test error messages in exceptions"""
    try:
        callable(*args, **kwargs)
        raise AssertionError("Should have raised %r" % exception(message))
    except exception, e:
        assert e.message == message


def test_error_messages_on_wrong_input():
    expected = ('labels_true and labels_pred must have same size,'
                ' got 2 and 3')
    assert_raise_message(ValueError, expected,
                         v_measure_score, [0, 1], [1, 1, 1])

    expected = "labels_true must be 1D: shape is (2, 2)"
    assert_raise_message(ValueError, expected,
                         homogeneity_score, [[0, 1], [1, 0]], [1, 1, 1])

    expected = "labels_pred must be 1D: shape is (2, 2)"
    assert_raise_message(ValueError, expected,
                         completeness_score, [0, 1, 0], [[1, 1], [0, 0]])


def test_homogeneous_but_not_complete_labeling():
    # homogeneous but not complete clustering
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 2, 2])
    assert_almost_equal(h, 1.00, 2)
    assert_almost_equal(c, 0.69, 2)
    assert_almost_equal(v, 0.81, 2)


def test_complete_but_not_homogeneous_labeling():
    # complete but not homogeneous clustering
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 1, 1, 2, 2],
        [0, 0, 1, 1, 1, 1])
    assert_almost_equal(h, 0.58, 2)
    assert_almost_equal(c, 1.00, 2)
    assert_almost_equal(v, 0.73, 2)


def test_not_complete_and_not_homogeneous_labeling():
    # neither complete nor homogeneous but not so bad either
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)


def test_non_consicutive_labels():
    # regression tests for labels with gaps
    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 2, 2, 2],
        [0, 1, 0, 1, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)

    h, c, v = homogeneity_completeness_v_measure(
        [0, 0, 0, 1, 1, 1],
        [0, 4, 0, 4, 2, 2])
    assert_almost_equal(h, 0.67, 2)
    assert_almost_equal(c, 0.42, 2)
    assert_almost_equal(v, 0.52, 2)
