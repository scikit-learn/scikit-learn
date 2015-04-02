import warnings
import unittest
import sys
import numpy as np

from nose.tools import assert_raises
from sklearn.utils.testing import (
    _assert_less,
    _assert_greater,
    assert_less_equal,
    assert_greater_equal,
    assert_warns,
    assert_no_warnings,
    assert_equal,
    set_random_state,
    assert_raise_message,
    assert_same_model,
    assert_not_same_model,
    assert_array_not_equal)
from sklearn.tree import DecisionTreeClassifier
from sklearn.lda import LDA
from sklearn.datasets import make_blobs
from sklearn.svm import LinearSVC
from sklearn.cluster import KMeans

try:
    from nose.tools import assert_less

    def test_assert_less():
        # Check that the nose implementation of assert_less gives the
        # same thing as the scikit's
        assert_less(0, 1)
        _assert_less(0, 1)
        assert_raises(AssertionError, assert_less, 1, 0)
        assert_raises(AssertionError, _assert_less, 1, 0)

except ImportError:
    pass

try:
    from nose.tools import assert_greater

    def test_assert_greater():
        # Check that the nose implementation of assert_less gives the
        # same thing as the scikit's
        assert_greater(1, 0)
        _assert_greater(1, 0)
        assert_raises(AssertionError, assert_greater, 0, 1)
        assert_raises(AssertionError, _assert_greater, 0, 1)

except ImportError:
    pass


def test_assert_less_equal():
    assert_less_equal(0, 1)
    assert_less_equal(1, 1)
    assert_raises(AssertionError, assert_less_equal, 1, 0)


def test_assert_greater_equal():
    assert_greater_equal(1, 0)
    assert_greater_equal(1, 1)
    assert_raises(AssertionError, assert_greater_equal, 0, 1)


def test_set_random_state():
    lda = LDA()
    tree = DecisionTreeClassifier()
    # LDA doesn't have random state: smoke test
    set_random_state(lda, 3)
    set_random_state(tree, 3)
    assert_equal(tree.random_state, 3)


def test_assert_raise_message():
    def _raise_ValueError(message):
        raise ValueError(message)

    assert_raise_message(ValueError, "test",
                         _raise_ValueError, "test")

    assert_raises(AssertionError,
                  assert_raise_message, ValueError, "something else",
                  _raise_ValueError, "test")

    assert_raises(ValueError,
                  assert_raise_message, TypeError, "something else",
                  _raise_ValueError, "test")


def test_assert_same_not_same_model():
    X1, y1 = make_blobs(n_samples=200, n_features=5, center_box=(-200, -150),
                        centers=2, random_state=0)
    X2, y2 = make_blobs(n_samples=100, n_features=5, center_box=(-1, 1),
                        centers=3, random_state=1)
    X3, y3 = make_blobs(n_samples=50, n_features=5, center_box=(-100, -50),
                        centers=4, random_state=2)

    for Estimator in (LinearSVC, KMeans):
        assert_same_model(X3, Estimator(random_state=0).fit(X1, y1),
                          Estimator(random_state=0).fit(X1, y1))

        assert_raises(AssertionError, assert_not_same_model, X3,
                      Estimator(random_state=0).fit(X1, y1),
                      Estimator(random_state=0).fit(X1, y1))

        assert_raises(AssertionError, assert_same_model, X3,
                      Estimator().fit(X1, y1), Estimator().fit(X2, y2))

        assert_not_same_model(X3, Estimator().fit(X1, y1),
                              Estimator().fit(X2, y2))


def test_array_not_equal():
    # Rank 1
    assert_array_not_equal(np.array([1, 2]), np.array([2, 2]))
    # Rank 2
    assert_array_not_equal(np.array([[1, 2], [3, 4]]),
                           np.array([[1, 2], [3, 5]]))
    # Different shapes
    assert_array_not_equal(np.array([1, 2]), np.array([[1, 2], [1, 2]]))
    # Different likes
    assert_array_not_equal([1, 2, 3], (1, 2, 5))
    # nan in array
    assert_array_not_equal(np.array([1, 2, 3]), np.array([1, 2, np.nan]))
    assert_array_not_equal(np.array(['floupi', 'floupa']),
                           np.array(['floupipi', 'floupa']))
    # Test different types
    for t in tuple('?bhilqpBHILQPfdgFDG') + ('S1', 'U1'):
        a = np.empty((4, 2, 3), t)
        a.fill(1)
        b = a.copy()
        b.fill(0)
        assert_array_not_equal(a, b)


# This class is inspired from numpy 1.7 with an alteration to check
# the reset warning filters after calls to assert_warns.
# This assert_warns behavior is specific to scikit-learn because
# `clean_warning_registry()` is called internally by assert_warns
# and clears all previous filters.
class TestWarns(unittest.TestCase):
    def test_warn(self):
        def f():
            warnings.warn("yo")
            return 3

        # Test that assert_warns is not impacted by externally set
        # filters and is reset internally.
        # This is because `clean_warning_registry()` is called internally by
        # assert_warns and clears all previous filters.
        warnings.simplefilter("ignore", UserWarning)
        assert_equal(assert_warns(UserWarning, f), 3)

        # Test that the warning registry is empty after assert_warns
        assert_equal(sys.modules['warnings'].filters, [])

        assert_raises(AssertionError, assert_no_warnings, f)
        assert_equal(assert_no_warnings(lambda x: x, 1), 1)

    def test_warn_wrong_warning(self):
        def f():
            warnings.warn("yo", DeprecationWarning)

        failed = False
        filters = sys.modules['warnings'].filters[:]
        try:
            try:
                # Should raise an AssertionError
                assert_warns(UserWarning, f)
                failed = True
            except AssertionError:
                pass
        finally:
            sys.modules['warnings'].filters = filters

        if failed:
            raise AssertionError("wrong warning caught by assert_warn")
