from abc import ABC
from abc import abstractmethod
from collections import UserString
import functools
from inspect import signature
import numbers
import operator

import numpy as np
from scipy.sparse import issparse

from .validation import _is_arraylike_not_scalar


def _validate_params(parameter_constraints, params):
    """Validate types and values of given parameters.

    Parameters
    ----------
    param_specs : dict
        A dictionary where the keys are the names of the parameters to validate. Each
        entry is a tuple consisting of:
        - The parameter value to validate
        - A list of tuples describing the valid types and values consisting of:
          - A type. It can be: any type, callable or "array-like".
          - The valid values for this type. It can be a list, a set or an Interval.
            If not provided, it means any value of this type is valid.
    """
    for param_name, constraints in parameter_constraints.items():
        _validate_param(param_name, params[param_name], constraints)


def _validate_param(param_name, param_val, constraints):
    """Check if a parameter has an appropriate type and value.

    Raises a TypeError if the type of the parameter is none of the valid types or a
    ValueError if the parameter has a valid type but its value is not in the accepted
    values for this type.

    Parameters
    ----------
    param_name : str
        The name of the parameter.

    param_val : object
        The value of the parameter.

    constraints : list
        The list of valid types and accepted values for each type.
    """
    constraints = map(make_constraint, constraints)

    for constraint in constraints:
        if constraint.is_satisfied_by(param_val):
            return
    else:
        if len(constraints) == 1:
            constraints_str = f"{constraints[0]}"
        else:
            constraints_str = f"{', '.join(constraints[:-1])} or {constraints[-1]}"

        raise ValueError(
            f"{param_name} must be {constraints_str}. Got {param_val} instead.")


def make_constraint(constraint):
    """"""
    if isinstance(constraint, str) and constraint == "array-like":
        return _ArrayLikes()
    if isinstance(constraint, str) and constraint == "sparse matrix":
        return _SparseMatrices()
    if constraint is callable:
        return _Callables()
    if constraint is None:
        return _NoneConstraint()
    if isinstance(constraint, (Interval, StrOptions)):
        return constraint
    if isinstance(constraint, type):
        return _InstancesOf(constraint)
    raise ValueError("Unknown constraint type ???")


def validate_params(parameter_constraints):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            param_names = [p.name for p in signature(func).parameters.values()]

            params = {}
            for i, val in enumerate(args):
                name = param_names[i]
                params[name] = val
            for name, val in kwargs.items():
                params[name] = val

            _validate_params(parameter_constraints, params)
            return func(*args, **kwargs)
        return wrapper
    return decorator


class Constraint(ABC):

    @abstractmethod
    def is_satisfied_by(self, val):
        pass

    @abstractmethod
    def __repr__(self):
        pass

class Interval(Constraint):
    """Class representing an typed interval constraint.

    Parameters
    ----------
    type : {numbers.Integral, numbers.Real}
        The set of numbers in which to set the interval.

    left : float or int or None
        The left bound of the interval. None means left bound is -∞.

    right : float, int or None
        The right bound of the interval. None means right bound is +∞.

    closed : {"left", "right", "both", "neither"}, default="left"
        Whether the interval is open or closed. Possible choices are:

        - `"left"`: the interval is closed on the left and open on the rigth.
          It is equivalent to the interval `[ left, right )`.
        - `"right"`: the interval is closed on the right and open on the left.
          It is equivalent to the interval `( left, right ]`.
        - `"both"`: the interval is closed.
          It is equivalent to the interval `[ left, right ]`.
        - `"neither"`: the interval is open.
          It is equivalent to the interval `( left, right )`.

    Notes
    -----
    Setting a bound to `None` and setting the interval closed is valid. For instance,
    strictly speaking, `Interval(0, None, closed="both")` corresponds to
    `[0, +∞) U {+∞}`.
    """

    @_validate_params(
        {
            "type": [numbers.Integral, numbers.Real],
            "left": [None, numbers.Real],
            "right": [None, numbers.Real],
            "closed": [Options({"left", "right", "both", "neither"})],
        }
    )
    def __init__(self, type, left, right, *, closed="left"):
        self.type = type
        self.left = left
        self.right = right
        self.closed = closed

        if self.type is numbers.Integral and self.left is not None and not isinstance(self.left, numbers.Integral):
            raise TypeError("Expecting left to be an int for an interval over the integers")
        if self.type is numbers.Integral and self.left is not None and not isinstance(self.left, numbers.Integral):
            raise TypeError("Expecting left to be an int for an interval over the integers")

    def __contains__(self, val):
        left_cmp = operator.lt if self.closed in ("left", "both") else operator.le
        right_cmp = operator.gt if self.closed in ("right", "both") else operator.ge

        if self.left is not None and left_cmp(val, self.left):
            return False
        if self.right is not None and right_cmp(val, self.right):
            return False
        return True

    def is_satisfied_by(self, val):
        if not isinstance(val, self.type):
            return False

        return val in self

    def __repr__(self):
        type_str = "an int" if self.type is numbers.Integral else "a float"
        left_bracket = "[" if self.closed in ("left", "both") else "("
        left_bound = "-inf" if self.left is None else self.left
        right_bound = "inf" if self.right is None else self.right
        right_bracket = "]" if self.closed in ("right", "both") else ")"
        return (
            f"{type_str} in the range "
            f"{left_bracket}{left_bound}, {right_bound}{right_bracket}"
        )


class StrOptions(Constraint):
    @validate_params(
        {
            "options": [set],
            "deprecated": [None, set],
        }
    )
    def __init__(self, options, deprecated=None):
        self.options = options
        self.deprecated = deprecated or {}

    def is_satisfied_by(self, val):
        return val in self.options
    
    def _mark_if_deprecated(self, option):
        option_str = f"{o!r}"
        if option in self.deprecated:
            option_str = f"{option_str} (deprecated)"
        return option_str
    
    def __repr__(self):
        f"a str among {{{[self._mark__if_deprecated(o) for o in self.options]}}}"


class _InstancesOf(Constraint):
    def __init__(self, type):
        self.type = type

    def _type_name(t):
        """Convert type into humman readable string."""
        module = t.__module__
        qualname = t.__qualname__
        if module == "builtins":
            return qualname
        elif t == numbers.Real:
            return "float"
        elif t == numbers.Integral:
            return "int"
        return f"{module}.{qualname}"
    
    def is_satisfied_by(self, val):
        return isinstance(val, self.type)

    def __repr__(self):
        return f"an instance of {self._type_name(self.type)!r}"


class _ArrayLikes(Constraint):
    def is_satisfied_by(self, val):
        return _is_arraylike_not_scalar(val)

    def __repr__(self):
        return "an array-like"


class _SparseMatrices(Constraint):
    def is_satisfied_by(self, val):
        return issparse(val)

    def __repr__(self):
        return "a sparse matrix"


class _Callables(Constraint):
    def is_satisfied_by(self, val):
        return callable(val)

    def __repr__(self):
        return "a callable"


class _NoneConstraint(Constraint):
    def is_satisfied_by(self, val):
        return val is None

    def __repr__(self):
        return "None"


def generate_invalid_param_val(target_type, valid_vals):
    """Generate an instance of a certain type that is not in some given container.

    Parameters
    ----------
    target_type : str, int or float
        The type of the value to generate.

    valid_vals : container
        The ensemble of values not to generate.

    Returns
    -------
    val : str, int or float
        The generated object.
    """
    if target_type is str:
        return "this_is_obviously_an_invalid_val"

    if isinstance(valid_vals, Interval):
        if valid_vals.left is not None:
            return valid_vals.left - 1
        else:
            return valid_vals.right + 1


def get_random_state_param_spec(val=None, estimator=True):
    """Appropriate valid types and values for the random_state parameter validation.

    Parameters
    ----------
    val : object, default=None
        The value of the random_state parameter. Unused for parameter validation of
        an estimator.

    estimator: bool, default=True
        Whether this is used to validate parameters of an estimator or not.
    """
    specs = [
        (numbers.Integral, Interval(0, 2**32 - 1, closed="both")),
        (np.random.RandomState,),
        (type(None),),
    ]
    return {"random_state": specs if estimator else (val, specs)}
