from collections import UserString
import numbers
import operator

import numpy as np
from scipy.sparse import issparse

from .validation import _is_arraylike_not_scalar


class Interval:
    """Class representing an interval.

    Parameters
    ----------
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

    def __init__(self, left, right, *, closed="left"):
        validate_params(
            {
                "left": (left, [(type(None),), (numbers.Real,)]),
                "right": (right, [(type(None),), (numbers.Real,)]),
                "closed": (closed, [(str, {"left", "right", "both", "neither"})]),
            }
        )

        self.left = left
        self.right = right
        self.closed = closed

    def __contains__(self, val):
        left_cmp = operator.lt if self.closed in ("left", "both") else operator.le
        right_cmp = operator.gt if self.closed in ("right", "both") else operator.ge

        if self.left is not None and left_cmp(val, self.left):
            return False
        if self.right is not None and right_cmp(val, self.right):
            return False
        return True

    def __repr__(self):
        left_bracket = "[" if self.closed in ("left", "both") else "("
        left_bound = "-inf" if self.left is None else self.left
        right_bound = "inf" if self.right is None else self.right
        right_bracket = "]" if self.closed in ("right", "both") else ")"
        return f"{left_bracket}{left_bound}, {right_bound}{right_bracket}"


def _has_type(obj, target_type):
    """Extension of isinstance to account for callable and array-like.

    Parameters
    ----------
    obj : object
        The object to test.

    target_type : type, callable or "array-like"
        The target type to check against.

    Returns
    -------
    has_type : bool
        Whether the object has the target type or not.
    """
    if isinstance(target_type, str) and target_type == "array-like":
        return _is_arraylike_not_scalar(obj)
    if isinstance(target_type, str) and target_type == "sparse matrix":
        return issparse(obj)
    if target_type is callable:
        return callable(obj)
    return isinstance(obj, target_type)


def validate_params(param_specs):
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
    for param_name, (param_val, expected_type_and_vals) in param_specs.items():
        validate_param(param_name, param_val, expected_type_and_vals)


def validate_param(param_name, param_val, expected_type_and_vals):
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

    expected_type_and_vals : list of tuples
        The list of valid types and accepted values for each type.
    """
    has_a_valid_type = False

    for target_type, *target_vals in expected_type_and_vals:
        if _has_type(param_val, target_type):
            has_a_valid_type = True

            if not target_vals:
                continue

            if param_val not in target_vals[0]:
                raise ValueError(
                    f"{param_name} of type {_type_name(target_type)} must be "
                    f"{_format_valid_vals_str(target_vals[0])}. "
                    f"Got {param_val} instead."
                )

    if not has_a_valid_type:
        valid_types = [t for t, *_ in expected_type_and_vals]
        raise TypeError(
            f"{param_name} must be {_format_valid_types_str(valid_types)}. "
            f"Got {type(param_val).__qualname__} instead."
        )


def _format_valid_vals_str(vals):
    """Convert a collection of values into human readable string."""
    if isinstance(vals, Interval):
        valid_vals_str = f"in the range {vals}"
    elif len(vals) == 1:
        valid_vals_str = f"{vals[0]!r}"
    else:
        valid_vals_str = f"one of {{{', '.join([repr(v) for v in vals])}}}"

    return valid_vals_str


def _format_valid_types_str(types):
    """Convert a list of types into human readable string."""
    valid_types_str = []
    if "array-like" in types:
        types.remove("array-like")
        valid_types_str.append("an array-like")
    if "sparse matrix" in types:
        types.remove("sparse matrix")
        valid_types_str.append("a sparse matrix")
    if callable in types:
        types.remove(callable)
        valid_types_str.append("a callable")

    if len(types) == 1:
        valid_types_str.append(f"an instance of {_type_name(types[0])}")
    else:
        valid_types_str.append(
            f"an instance of {{{', '.join(_type_name(t) for t in types)}}}"
        )

    if len(valid_types_str) == 1:
        valid_types_str = valid_types_str[0]
    else:
        valid_types_str = f"{', '.join(valid_types_str[:-1])} or {valid_types_str[-1]}"

    return valid_types_str


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


class DeprecatedParamStr(UserString):
    """Mark a valid parameter value of type str as deprecated."""

    def __repr__(self):
        return f"{self.data!r} (deprecated)"


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
