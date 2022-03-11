from numbers import Integral, Real

import numpy as np
from scipy.sparse import csr_matrix
import pytest

from sklearn.utils._param_validation import Interval
from sklearn.utils._param_validation import StrOptions
from sklearn.utils._param_validation import _ArrayLikes
from sklearn.utils._param_validation import _SparseMatrices
from sklearn.utils._param_validation import _Callables
from sklearn.utils._param_validation import _InstancesOf
from sklearn.utils._param_validation import _NoneConstraint
from sklearn.utils._param_validation import make_constraint


@pytest.mark.parametrize("interval_type", [Integral, Real])
def test_interval_range(interval_type):
    """Check the range of values depending on closed."""
    interval = Interval(interval_type, -2, 2, closed="left")
    assert -2 in interval and 2 not in interval

    interval = Interval(interval_type, -2, 2, closed="right")
    assert -2 not in interval and 2 in interval

    interval = Interval(interval_type, -2, 2, closed="both")
    assert -2 in interval and 2 in interval

    interval = Interval(interval_type, -2, 2, closed="neither")
    assert -2 not in interval and 2 not in interval


def test_interval_inf_in_bounds():
    """Check that inf is included if a bound is closed and set to None.

    Only valid for real intervals.
    """
    interval = Interval(Real, 0, None, closed="right")
    assert np.inf in interval

    interval = Interval(Real, None, 0, closed="left")
    assert -np.inf in interval


@pytest.mark.parametrize("params, error, match", [
    (
        {"type": int, "left": 0, "right": 1},
        ValueError,
        r"type must be numbers\.Integral or numbers\.Real"
    ),
    (
        {"type": Integral, "left": 1., "right": 2},
        TypeError,
        r"Expecting left to be an int for an interval over the integers"
    ),
    (
        {"type": Integral, "left": 1, "right": 2.},
        TypeError,
        r"Expecting right to be an int for an interval over the integers"
    ),
    (
        {"type": Integral, "left": None, "right": 0, "closed": "left"},
        ValueError,
        r"left can't be None when closed == left"
    ),
    (
        {"type": Integral, "left": 0, "right": None, "closed": "right"},
        ValueError,
        r"right can't be None when closed == right"
    ),
])
def test_interval_errors(params, error, match):
    """Check that informative errors are raised for invalid combination of parameters"""
    with pytest.raises(error, match=match):
        Interval(**params)


def test_stroptions():
    """Sanity check for the StrOptions constraint"""
    options = StrOptions({"a", "b", "c"}, deprecated={"c"})
    assert options.is_satisfied_by("a")
    assert options.is_satisfied_by("c")
    assert not options.is_satisfied_by("d")

    assert "(deprecated)" in repr(options)


@pytest.mark.parametrize("constraint", [
    Interval(Real, None, 0),
    Interval(Real, 0, None),
    StrOptions({"a", "b", "c"})
])
def test_generate_invalid_param_val(constraint):
    """Check that the value generated does not satisfy the constraint"""
    bad_value = constraint.generate_invalid_param_val()
    assert not constraint.is_satisfied_by(bad_value)


@pytest.mark.parametrize("constraint_declaration, value",[
    (Interval(Real, 0, 1), 0.42),
    (Interval(Integral, 0, None), 42),
    (StrOptions({"a", "b", "c"}), "b"),
    (callable, lambda x: x + 1),
    (None, None),
    ("array-like", [[1, 2], [3, 4]]),
    ("array-like", np.array([[1, 2], [3, 4]])),
    ("sparse matrix", csr_matrix([[1, 2], [3, 4]])),
    (np.random.RandomState, np.random.RandomState(0)),
    (int, 1),
])
def test_is_satisified_by(constraint_declaration, value):
    """Sanity check for the is_satisfied_by method"""
    constraint = make_constraint(constraint_declaration)
    assert constraint.is_satisfied_by(value)


@pytest.mark.parametrize("constraint_declaration, expected_constraint_class", [
    (Interval(Real, 0, 1), Interval),
    (StrOptions({"option1", "option2"}), StrOptions),
    ("array-like", _ArrayLikes),
    ("sparse matrix", _SparseMatrices),
    (None, _NoneConstraint),
    (callable, _Callables),
    (int, _InstancesOf),
])
def test_make_constraint(constraint_declaration, expected_constraint_class):
    """Check that make_constraint dispaches to the appropriate constraint class"""
    constraint = make_constraint(constraint_declaration)
    assert constraint.__class__ is expected_constraint_class


def test_make_constraint_unknown():
    with pytest.raises(ValueError, match="Unknown constraint"):
        make_constraint("not a valid constraint")
