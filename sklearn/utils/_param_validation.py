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


def _validate_params(parameter_constraints, params):
    """Validate types and values of given parameters.

    Parameters
    ----------
    parameter_constraints : dict
        A dictionary `param_name: list of constraints`.

    params : dict
        A dictionary `param_name: param_value`. The parameters to validate against the
        constraints.
    """
    for param_name, constraints in parameter_constraints.items():
        _validate_param(param_name, params[param_name], constraints)


def _validate_param(param_name, param_val, constraints):
    """Check if a parameter satisfies a given constraint.

    Raises a ValueError if the constraint is not satified.

    Parameters
    ----------
    param_name : str
        The name of the parameter.

    param_val : object
        The value of the parameter.

    constraints : list
        To be valid, the parameter must satisfy one of these constraints.
    """
    constraints = [make_constraint(constraint) for constraint in constraints]

    for constraint in constraints:
        if constraint.is_satisfied_by(param_val):
            # this constraint is satisfied, nothing to do.
            return
    else:
        # No constraint is satisfied, raise with an informative message.
        if len(constraints) == 1:
            constraints_str = f"{constraints[0]}"
        else:
            constraints_str = (
                f"{', '.join([repr(c) for c in constraints[:-1]])} or {constraints[-1]}"
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
    if constraint is callable:
        return _Callables()
    if constraint is None:
        return _NoneConstraint()
    if isinstance(constraint, type):
        return _InstancesOf(constraint)
    if isinstance(constraint, (Interval, StrOptions, TypeOptions)):
        return constraint
    raise ValueError(f"Unknown constraint type: {constraint}")


def validate_params(parameter_constraints):
    """Decorator to validate types and values of functions and methods.

    Parameters
    ----------
    parameter_constraints : dict
        A dictionary `param_name: list of constraints`.

    Returns
    -------
    decorated_function : function or method
        The decorated function.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            param_names = [
                p.name
                for p in signature(func).parameters.values()
                if p.kind not in (p.VAR_KEYWORD, p.VAR_POSITIONAL)
            ]

            # combine the signature and the actual args and kwargs to build a map
            # param_name: param_value
            params = {}
            for i, val in enumerate(args):
                # the signature objects holds the parameters in the same order as in the
                # function definition
                name = param_names[i]

                if name == "self":
                    # ignore self for methods
                    continue

                params[name] = val
            for name, val in kwargs.items():
                params[name] = val

            _validate_params(parameter_constraints, params)
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


class StrOptions(_Constraint):
    """Constraint representing a set of strings.
    
    Parameters
    ----------
    options : set of str
        The set of valid strings.

    deprecated : set of str or None, default=None
        A subset of the `options` to mark as deprecated in the repr of the constraint.
    """
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


class TypeOptions(_Constraint):
    """Constraint representing a set of types.
    
    Parameters
    ----------
    options : set of instances of type
        The set of valid types.
    """
    def __init__(self, options):
        self.options = options

    def generate_invalid_param_val(self):
        """Return a value that does not satisfy the constraint."""
        return type("BadType", (), {})

    def is_satisfied_by(self, val):
        return isinstance(val, type) and val in self.options

    def __repr__(self):
        return f"one of {{{self.options}}}"


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

    def __init__(self, type, left, right, *, closed="left"):
        self.type = type
        self.left = left
        self.right = right
        self.closed = closed

        if (
            self.type is Integral
            and self.left is not None
            and not isinstance(self.left, Integral)
        ):
            raise TypeError(
                "Expecting left to be an int for an interval over the integers"
            )
        if (
            self.type is Integral
            and self.left is not None
            and not isinstance(self.left, Integral)
        ):
            raise TypeError(
                "Expecting left to be an int for an interval over the integers"
            )

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
        """Convert type into humman readable string."""
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


class _NoneConstraint(_Constraint):
    """Constraint representing the None singleton."""
    def is_satisfied_by(self, val):
        return val is None

    def __repr__(self):
        return "None"


def get_random_state_param_constraints():
    """Appropriate constraints for the validation of the random_state parameter.

    Return
    ------
    constraints : list
        The constraints for random_state.
    """
    constraints = [
        Interval(Integral, 0, 2**32 - 1, closed="both"),
        np.random.RandomState,
        None,
    ]
    return {"random_state": constraints}
