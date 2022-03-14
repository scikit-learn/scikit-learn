from abc import ABC
from abc import abstractmethod
import functools
from inspect import signature
from numbers import Integral
from numbers import Real
import operator

import numpy as np
from scipy.sparse import issparse

from .validation import _is_arraylike_not_scalar


def validate_parameter_constraints(parameter_constraints, params):
    """Validate types and values of given parameters.

    Parameters
    ----------
    parameter_constraints : dict
        A dictionary `param_name: list of constraints`. A parameter is valid if it
        satisfies one of the constraints from the list. Constraints can be:
        - an Interval object, representing a continuous or discrete range of numbers
        - the string "array-like"
        - the string "sparse matrix"
        - callable
        - None, meaning that None is a valid value for the parameter
        - any type, meaning that any instance of this type is valid
        - a StrOptions object, representing a set of strings

    params : dict
        A dictionary `param_name: param_value`. The parameters to validate against the
        constraints.
    """
    for param_name, constraints in parameter_constraints.items():

        param_val = params[param_name]
        constraints = [make_constraint(constraint) for constraint in constraints]

        for constraint in constraints:
            if constraint.is_satisfied_by(param_val):
                # this constraint is satisfied, no need to check further.
                break
        else:
            # No constraint is satisfied, raise with an informative message.
            if len(constraints) == 1:
                constraints_str = f"{constraints[0]}"
            else:
                constraints_str = (
                    f"{', '.join([repr(c) for c in constraints[:-1]])} or"
                    f" {constraints[-1]}"
                )

            raise ValueError(
                f"{param_name} must be {constraints_str}. Got {param_val!r} instead."
            )


def make_constraint(constraint):
    """Convert the constraint into the appropriate Constraint object.

    Parameters
    ----------
    constraint : object
        The constraint to convert.

    Returns
    -------
    constraint : instance of _Constraint
        The converted constraint.
    """
    if isinstance(constraint, str) and constraint == "array-like":
        return _ArrayLikes()
    if isinstance(constraint, str) and constraint == "sparse matrix":
        return _SparseMatrices()
    if isinstance(constraint, str) and constraint == "random_state":
        return _RandomStates()
    if constraint is callable:
        return _Callables()
    if constraint is None:
        return _NoneConstraint()
    if isinstance(constraint, type):
        return _InstancesOf(constraint)
    if isinstance(constraint, (Interval, StrOptions)):
        return constraint
    raise ValueError(f"Unknown constraint type: {constraint}")


def validate_params(parameter_constraints):
    """Decorator to validate types and values of functions and methods.

    Parameters
    ----------
    parameter_constraints : dict
        A dictionary `param_name: list of constraints`. See the docstring of
        `validate_parameter_constraints` for a description of the accepted constraints.

    Returns
    -------
    decorated_function : function or method
        The decorated function.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            sig_params = [
                p
                for p in signature(func).parameters.values()
                if p.kind not in (p.VAR_KEYWORD, p.VAR_POSITIONAL)
            ]

            # combine the signature and the actual args and kwargs to build a map
            # param_name: param_value
            params = {}

            # First, the positional arguments for which we need to parse the signature
            # to recover their names. The signature objects holds the parameters in the
            # same order as in the function definition
            for val, param in zip(args, sig_params):
                # ignore self for methods
                if param.name == "self":
                    continue

                params[param.name] = val

            # Then, the keyword arguments which already hold both the name and the value
            for name, val in kwargs.items():
                params[name] = val

            # Finally, the parameters with have a default that are unset
            for param in [p for p in sig_params if p.default is not p.empty]:
                if param.name not in params:
                    params[param.name] = param.default

            validate_parameter_constraints(parameter_constraints, params)
            return func(*args, **kwargs)

        return wrapper

    return decorator


class _Constraint(ABC):
    """Base class for the constraint objects."""

    @abstractmethod
    def is_satisfied_by(self, val):
        """Whether or not a value satisfies the constraint.

        Parameters
        ----------
        val : object
            The value to check.

        Returns
        -------
        is_satisfied : bool
            Whether or not the constraint is satisfied by this value.
        """

    @abstractmethod
    def __repr__(self):
        """A human readable representational string of the constraint."""


class _InstancesOf(_Constraint):
    """Constraint representing instances of a given type.

    Parameters
    ----------
    type : type
        The valid type.
    """

    def __init__(self, type):
        self.type = type

    def _type_name(self, t):
        """Convert type into human readable string."""
        module = t.__module__
        qualname = t.__qualname__
        if module == "builtins":
            return qualname
        elif t == Real:
            return "float"
        elif t == Integral:
            return "int"
        return f"{module}.{qualname}"

    def is_satisfied_by(self, val):
        return isinstance(val, self.type)

    def __repr__(self):
        return f"an instance of {self._type_name(self.type)!r}"


class _NoneConstraint(_Constraint):
    """Constraint representing the None singleton."""

    def is_satisfied_by(self, val):
        return val is None

    def __repr__(self):
        return "None"


class StrOptions(_Constraint):
    """Constraint representing a set of strings.

    Parameters
    ----------
    options : set of str
        The set of valid strings.

    deprecated : set of str or None, default=None
        A subset of the `options` to mark as deprecated in the repr of the constraint.
    """

    @validate_params({"options": [set], "deprecated": [set, None]})
    def __init__(self, options, deprecated=None):
        self.options = options
        self.deprecated = deprecated or {}

    def generate_invalid_param_val(self):
        """Return a value that does not satisfy the constraint."""
        return "this_is_obviously_an_invalid_val"

    def is_satisfied_by(self, val):
        return isinstance(val, str) and val in self.options

    def _mark_if_deprecated(self, option):
        """Add a deprecated mark to an option if needed."""
        option_str = f"{option!r}"
        if option in self.deprecated:
            option_str = f"{option_str} (deprecated)"
        return option_str

    def __repr__(self):
        options_str = (
            f"{', '.join([self._mark_if_deprecated(o) for o in self.options])}"
        )
        return f"a str among {{{options_str}}}"


class Interval(_Constraint):
    """Constraint representing an typed interval.

    Parameters
    ----------
    type : {numbers.Integral, numbers.Real}
        The set of numbers in which to set the interval.

    left : float or int or None
        The left bound of the interval. None means left bound is -∞.

    right : float, int or None
        The right bound of the interval. None means right bound is +∞.

    closed : {"left", "right", "both", "neither"}
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

    def _check_params(self):
        if self.type is Integral:
            suffix = "for an interval over the integers."
            if self.left is not None and not isinstance(self.left, Integral):
                raise TypeError(f"Expecting left to be an int {suffix}")
            if self.right is not None and not isinstance(self.right, Integral):
                raise TypeError(f"Expecting right to be an int {suffix}")
            if self.left is None and self.closed in ("left", "both"):
                raise ValueError(
                    f"left can't be None when closed == {self.closed} {suffix}"
                )
            if self.right is None and self.closed in ("right", "both"):
                raise ValueError(
                    f"right can't be None when closed == {self.closed} {suffix}"
                )

        if self.right is not None and self.left is not None and self.right <= self.left:
            raise ValueError(
                f"right can't be less than left. Got left={self.left} and "
                f"right={self.right}"
            )

    @validate_params(
        {
            "type": [type],
            "left": [Integral, Real, None],
            "right": [Integral, Real, None],
            "closed": [StrOptions({"left", "right", "both", "neither"})],
        }
    )
    def __init__(self, type, left, right, closed):
        self.type = type
        self.left = left
        self.right = right
        self.closed = closed

        self._check_params()

    def __contains__(self, val):
        left_cmp = operator.lt if self.closed in ("left", "both") else operator.le
        right_cmp = operator.gt if self.closed in ("right", "both") else operator.ge

        if self.left is not None and left_cmp(val, self.left):
            return False
        if self.right is not None and right_cmp(val, self.right):
            return False
        return True

    def generate_invalid_param_val(self):
        """Return a value that does not satisfy the constraint."""
        if self.left is None and self.right is None:
            raise AttributeError

        if self.left is not None:
            return self.left - 1
        else:
            return self.right + 1

    def is_satisfied_by(self, val):
        if not isinstance(val, self.type):
            return False

        return val in self

    def __repr__(self):
        type_str = "an int" if self.type is Integral else "a float"
        left_bracket = "[" if self.closed in ("left", "both") else "("
        left_bound = "-inf" if self.left is None else self.left
        right_bound = "inf" if self.right is None else self.right
        right_bracket = "]" if self.closed in ("right", "both") else ")"
        return (
            f"{type_str} in the range "
            f"{left_bracket}{left_bound}, {right_bound}{right_bracket}"
        )


class _ArrayLikes(_Constraint):
    """Constraint representing array-likes"""

    def is_satisfied_by(self, val):
        return _is_arraylike_not_scalar(val)

    def __repr__(self):
        return "an array-like"


class _SparseMatrices(_Constraint):
    """Constraint representing sparse matrices."""

    def is_satisfied_by(self, val):
        return issparse(val)

    def __repr__(self):
        return "a sparse matrix"


class _Callables(_Constraint):
    """Constraint representing callables."""

    def is_satisfied_by(self, val):
        return callable(val)

    def __repr__(self):
        return "a callable"


class _RandomStates(_Constraint):
    """Constraint representing random states.

    Convenience class for
    [Interval(Integral, 0, 2**32 - 1, closed="both"), np.random.RandomState, None]
    """

    def __init__(self):
        self._constraints = [
            Interval(Integral, 0, 2**32 - 1, closed="both"),
            _InstancesOf(np.random.RandomState),
            _NoneConstraint(),
        ]

    def is_satisfied_by(self, val):
        return any(c.is_satisfied_by(val) for c in self._constraints)

    def __repr__(self):
        return (
            f"{', '.join([repr(c) for c in self._constraints[:-1]])} or"
            f" {self._constraints[-1]}"
        )
