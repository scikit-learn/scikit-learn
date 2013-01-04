from nose.tools import assert_raises

from sklearn.utils.testing import (
    _assert_less,
    _assert_greater,
    assert_equal,
    set_random_state,
    assert_raise_message)

from sklearn.tree import DecisionTreeClassifier
from sklearn.lda import LDA

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
