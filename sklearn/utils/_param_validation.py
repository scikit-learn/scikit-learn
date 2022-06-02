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


def validate_parameter_constraints(parameter_constraints, params, caller_name):
    """Validate types and values of given parameters.

    Parameters
    ----------
    parameter_constraints : dict
        A dictionary `param_name: list of constraints`. A parameter is valid if it
        satisfies one of the constraints from the list. Constraints can be:
        - an Interval object, representing a continuous or discrete range of numbers
        - the string "array-like"
        - the string "sparse matrix"
        - the string "random state"
        - callable
        - None, meaning that None is a valid value for the parameter
        - any type, meaning that any instance of this type is valid
        - a StrOptions object, representing a set of strings

    params : dict
        A dictionary `param_name: param_value`. The parameters to validate against the
        constraints.

    caller_name : str
        The name of the estimator or function or method that called this function.
    """
    if params.keys() != parameter_constraints.keys():
        raise ValueError(
            f"The parameter constraints {list(parameter_constraints.keys())} do not "
            f"match the parameters to validate {list(params.keys())}."
        )

    for param_name, param_val in params.items():
        constraints = parameter_constraints[param_name]
        constraints = [make_constraint(constraint) for constraint in constraints]

        for constraint in constraints:
            if constraint.is_satisfied_by(param_val):
                # this constraint is satisfied, no need to check further.
                break
        else:
            # No constraint is satisfied, raise with an informative message.

            # Ignore constraints that only contains internal options that we don't want
            # to expose in the error message
            constraints = [constraint for constraint in constraints if str(constraint)]

            if len(constraints) == 1:
                constraints_str = f"{constraints[0]}"
            else:
                constraints_str = (
                    f"{', '.join([str(c) for c in constraints[:-1]])} or"
                    f" {constraints[-1]}"
                )

            raise ValueError(
                f"The {param_name!r} parameter of {caller_name} must be"
                f" {constraints_str}. Got {param_val!r} instead."
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

            func_sig = signature(func)

            # Map *args/**kwargs to the function signature
            params = func_sig.bind(*args, **kwargs)
            params.apply_defaults()

            # ignore self/cls and positional/keyword markers
            to_ignore = [
                p.name
                for p in func_sig.parameters.values()
                if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD)
            ]
            to_ignore += ["self", "cls"]
            params = {k: v for k, v in params.arguments.items() if k not in to_ignore}

            validate_parameter_constraints(
                parameter_constraints, params, caller_name=func.__qualname__
            )
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
    def __str__(self):
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

    def __str__(self):
        return f"an instance of {self._type_name(self.type)!r}"


class _NoneConstraint(_Constraint):
    """Constraint representing the None singleton."""

    def is_satisfied_by(self, val):
        return val is None

    def __str__(self):
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

    @validate_params(
        {"options": [set], "deprecated": [set, None], "internal": [set, None]}
    )
    def __init__(self, options, deprecated=None, internal=None):
        self.options = options
        self.deprecated = deprecated or set()
        self.internal = internal or set()

        if self.deprecated - self.options:
            raise ValueError("The deprecated options must be a subset of the options.")

        if self.internal - self.options:
            raise ValueError("The internal options must be a subset of the options.")

        if self.deprecated & self.internal:
            raise ValueError(
                "The deprecated and internal parameters should not overlap."
            )

    def is_satisfied_by(self, val):
        return isinstance(val, str) and val in self.options

    def _mark_if_deprecated(self, option):
        """Add a deprecated mark to an option if needed."""
        option_str = f"{option!r}"
        if option in self.deprecated:
            option_str = f"{option_str} (deprecated)"
        return option_str

    def __str__(self):
        visible_options = [o for o in self.options if o not in self.internal]

        if not visible_options:
            return ""

        options_str = (
            f"{', '.join([self._mark_if_deprecated(o) for o in visible_options])}"
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

        - `"left"`: the interval is closed on the left and open on the right.
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
    strictly speaking, `Interval(Real, 0, None, closed="both")` corresponds to
    `[0, +∞) U {+∞}`.
    """

    @validate_params(
        {
            "type": [type],
            "left": [Integral, Real, None],
            "right": [Integral, Real, None],
            "closed": [StrOptions({"left", "right", "both", "neither"})],
        }
    )
    def __init__(self, type, left, right, *, closed):
        self.type = type
        self.left = left
        self.right = right
        self.closed = closed

        self._check_params()

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

    def __str__(self):
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

    def __str__(self):
        return "an array-like"


class _SparseMatrices(_Constraint):
    """Constraint representing sparse matrices."""

    def is_satisfied_by(self, val):
        return issparse(val)

    def __str__(self):
        return "a sparse matrix"


class _Callables(_Constraint):
    """Constraint representing callables."""

    def is_satisfied_by(self, val):
        return callable(val)

    def __str__(self):
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

    def __str__(self):
        return (
            f"{', '.join([repr(c) for c in self._constraints[:-1]])} or"
            f" {self._constraints[-1]}"
        )


def generate_invalid_param_val(constraint):
    """Return a value that does not satisfy the constraint.

    This is only useful for testing purpose.

    Parameters
    ----------
    constraint : Constraint
        The constraint to generate a value for.

    Returns
    -------
    val : object
        A value that does not satisfy the constraint.
    """
    if isinstance(constraint, StrOptions):
        return f"not {' or '.join(constraint.options)}"
    elif isinstance(constraint, Interval):
        interval = constraint
        if interval.left is None and interval.right is None:
            raise NotImplementedError

        if interval.left is not None:
            return interval.left - 1
        else:
            return interval.right + 1
    else:
        raise NotImplementedError
