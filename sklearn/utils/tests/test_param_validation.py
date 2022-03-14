from numbers import Integral, Real

import numpy as np
from scipy.sparse import csr_matrix
import pytest

from sklearn.utils._param_validation import Interval
from sklearn.utils._param_validation import StrOptions
from sklearn.utils._param_validation import _ArrayLikes
from sklearn.utils._param_validation import _Callables
from sklearn.utils._param_validation import _InstancesOf
from sklearn.utils._param_validation import _NoneConstraint
from sklearn.utils._param_validation import _RandomStates
from sklearn.utils._param_validation import _SparseMatrices
from sklearn.utils._param_validation import make_constraint


@pytest.mark.parametrize("interval_type", [Integral, Real])
def test_interval_range(interval_type):
    """Check the range of values depending on closed."""
    interval = Interval(interval_type, -2, 2, "left")
    assert -2 in interval and 2 not in interval

    interval = Interval(interval_type, -2, 2, "right")
    assert -2 not in interval and 2 in interval

    interval = Interval(interval_type, -2, 2, "both")
    assert -2 in interval and 2 in interval

    interval = Interval(interval_type, -2, 2, "neither")
    assert -2 not in interval and 2 not in interval


def test_interval_inf_in_bounds():
    """Check that inf is included if a bound is closed and set to None.

    Only valid for real intervals.
    """
    interval = Interval(Real, 0, None, "right")
    assert np.inf in interval

    interval = Interval(Real, None, 0, "left")
    assert -np.inf in interval


@pytest.mark.parametrize(
    "params, error, match",
    [
        (
            [Integral, 1.0, 2, "both"],
            TypeError,
            r"Expecting left to be an int for an interval over the integers",
        ),
        (
            [Integral, 1, 2.0, "neither"],
            TypeError,
            r"Expecting right to be an int for an interval over the integers",
        ),
        (
            [Integral, None, 0, "left"],
            ValueError,
            r"left can't be None when closed == left",
        ),
        (
            [Integral, 0, None, "right"],
            ValueError,
            r"right can't be None when closed == right",
        ),
        (
            [Integral, 1, -1, "both"],
            ValueError,
            r"right can't be less than left",
        ),
    ],
)
def test_interval_errors(params, error, match):
    """Check that informative errors are raised for invalid combination of parameters"""
    with pytest.raises(error, match=match):
        Interval(*params)


def test_stroptions():
    """Sanity check for the StrOptions constraint"""
    options = StrOptions({"a", "b", "c"}, deprecated={"c"})
    assert options.is_satisfied_by("a")
    assert options.is_satisfied_by("c")
    assert not options.is_satisfied_by("d")

    assert "(deprecated)" in repr(options)


@pytest.mark.parametrize(
    "constraint",
    [
        Interval(Real, None, 0, "left"),
        Interval(Real, 0, None, "left"),
        StrOptions({"a", "b", "c"}),
    ],
)
def test_generate_invalid_param_val(constraint):
    """Check that the value generated does not satisfy the constraint"""
    bad_value = constraint.generate_invalid_param_val()
    assert not constraint.is_satisfied_by(bad_value)


# a class to test the _InstancesOf constraint
class _SomeClass:
    pass


@pytest.mark.parametrize(
    "constraint_declaration, value",
    [
        (Interval(Real, 0, 1, "both"), 0.42),
        (Interval(Integral, 0, None, "neither"), 42),
        (StrOptions({"a", "b", "c"}), "b"),
        (callable, lambda x: x + 1),
        (None, None),
        ("array-like", [[1, 2], [3, 4]]),
        ("array-like", np.array([[1, 2], [3, 4]])),
        ("sparse matrix", csr_matrix([[1, 2], [3, 4]])),
        ("random_state", 0),
        ("random_state", np.random.RandomState(0)),
        ("random_state", None),
        (_SomeClass, _SomeClass()),
        (int, 1),
        (Real, 0.5),
    ],
)
def test_is_satisified_by(constraint_declaration, value):
    """Sanity check for the is_satisfied_by method"""
    constraint = make_constraint(constraint_declaration)
    assert constraint.is_satisfied_by(value)


@pytest.mark.parametrize(
    "constraint_declaration, expected_constraint_class",
    [
        (Interval(Real, 0, 1, "both"), Interval),
        (StrOptions({"option1", "option2"}), StrOptions),
        ("array-like", _ArrayLikes),
        ("sparse matrix", _SparseMatrices),
        ("random_state", _RandomStates),
        (None, _NoneConstraint),
        (callable, _Callables),
        (int, _InstancesOf),
    ],
)
def test_make_constraint(constraint_declaration, expected_constraint_class):
    """Check that make_constraint dispaches to the appropriate constraint class"""
    constraint = make_constraint(constraint_declaration)
    assert constraint.__class__ is expected_constraint_class


def test_make_constraint_unknown():
    with pytest.raises(ValueError, match="Unknown constraint"):
        make_constraint("not a valid constraint")
