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
    assert_not_equal,
    set_random_state,
    assert_raise_message,
    assert_same_model,
    assert_not_same_model,
    assert_array_equal,
    assert_array_not_equal)
from sklearn.tree import DecisionTreeClassifier
from sklearn.lda import LDA
from sklearn.datasets import make_blobs
from sklearn.svm import LinearSVC

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


class MockClassifier():
    def __init__(self, same_predict=True, same_transform=True,
                 same_decision_function=True, same_predict_proba=True,
                 same_attributes=True, n_outputs=1):
        self._same_predict = same_predict
        self._same_transform = same_transform
        self._same_decision_function = same_decision_function
        self._same_predict_proba = same_predict_proba
        self._n_outputs = n_outputs
        self._same_attributes = same_attributes

    def fit(self, X=None, y=None):
        X, y = np.array(X), np.array(y)
        self._n_features = X.shape[1] if len(X.shape) == 2 else 1
        self._n_classes = np.unique(y).size
        if self._same_attributes:
            random_base_module = np.random.RandomState(0)  # Seeded RNG
        else:
            random_base_module = np.random  # Unseeded RNG
        # A float attribute
        self.attribute1_ = random_base_module.random_sample()
        # A list attribute
        self.attribute2_ = random_base_module.random_sample((10, 10))
        return self

    def _output(self, X, same=False, integer=True, shape=None):
        X = np.array(X)
        n_samples = X.shape[0] if len(X.shape) == 2 else X.size
        if shape is None:
            shape = (n_samples, self._n_outputs)
        if same:
            random_base_module = np.random.RandomState(0)
        else:
            random_base_module = np.random  # Truly random
        if integer:
            return random_base_module.randint(0, self._n_classes, size=shape)
        else:
            return random_base_module.random_sample(shape)

    def predict(self, X):
        return self._output(X, self._same_predict, integer=True)

    def transform(self, X):
        return self._output(X, self._same_transform, integer=False)

    def predict_proba(self, X):
        return self._output(X, self._same_predict_proba, integer=False)

    def decision_function(self, X):
        return self._output(X, self._same_decision_function, integer=False)


def test_assert_xsame_model_():
    X1, y1 = make_blobs(n_samples=200, n_features=5, center_box=(-200, -150),
                        centers=2, random_state=0)
    X2, y2 = make_blobs(n_samples=100, n_features=5, center_box=(-1, 1),
                        centers=3, random_state=1)
    X3, y3 = make_blobs(n_samples=50, n_features=5, center_box=(-100, -50),
                        centers=4, random_state=2)

    # ----> Test using a regular estimator <----
    assert_same_model(X3, LinearSVC(random_state=0).fit(X1, y1),
                      LinearSVC(random_state=0).fit(X1, y1))

    assert_raises(AssertionError, assert_not_same_model, X3,
                  LinearSVC(random_state=0).fit(X1, y1),
                  LinearSVC(random_state=0).fit(X1, y1))

    assert_raises(AssertionError, assert_same_model, X3,
                  LinearSVC().fit(X1, y1), LinearSVC().fit(X2, y2))

    assert_not_same_model(X3, LinearSVC().fit(X1, y1), LinearSVC().fit(X2, y2))

    # ---->  Test using the TestClassifier  <----
    for o in (1, 2, 3):
        t1 = MockClassifier(n_outputs=o).fit(X1, y1)
        t2 = MockClassifier(n_outputs=o).fit(X1, y1)

        # predict alone different
        t1._same_predict = False
        assert_array_not_equal(t1.predict(X3), t2.predict(X3))
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # decision_function alone different
        t1._same_predict = True
        t1._same_decision_function = False
        assert_array_not_equal(t1.decision_function(X3),
                               t2.decision_function(X3))
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # predict_proba alone different
        t1._same_decision_function = True
        t1._same_predict_proba = False
        assert_array_not_equal(t1.predict_proba(X3), t2.predict_proba(X3))
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # transform alone different
        t1._same_predict_proba = True
        t1._same_transform = False
        assert_array_not_equal(t1.transform(X3), t2.transform(X3))
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # Attributes alone differ
        t1._same_transform = True
        t1._same_attributes = False
        t1.fit(X1, y1)  # Attributes are set inside fit only
        assert_not_equal(t1.attribute1_, t2.attribute1_)
        assert_array_not_equal(t1.attribute2_, t2.attribute2_)
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # All set to True
        t1._same_attributes = True
        t1.fit(X1, y1)
        assert_equal(t1.attribute1_, t2.attribute1_)
        assert_array_equal(t1.attribute2_, t2.attribute2_)
        assert_array_equal(t1.predict(X3), t2.predict(X3))
        assert_array_equal(t1.transform(X3), t2.transform(X3))
        assert_array_equal(t1.predict_proba(X3), t2.predict_proba(X3))
        assert_array_equal(t1.decision_function(X3),
                           t2.decision_function(X3))

        # Any two set to differ
        t1._same_predict_proba = False
        t1._same_predict = False
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)

        # Any one set to False along with differing attributes
        t1._same_predict_proba = True
        t1._same_transform = False
        t1._same_attributes = False
        t1.fit(X1, y1)  # Attributes are set inside fit only
        assert_raises(AssertionError, assert_same_model, X3, t1, t1)
        assert_not_same_model(X3, t1, t2)

        # All set to False
        t1._same_decision_function = False
        t1._same_predict_proba = False
        t1._same_predict = False
        assert_raises(AssertionError, assert_same_model, X3, t1, t2)
        assert_not_same_model(X3, t1, t2)


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
