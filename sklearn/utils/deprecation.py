# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import functools
import warnings
from inspect import signature

__all__ = ["deprecated"]


class deprecated:
    """Decorator to mark a function or class as deprecated.

    Issue a warning when the function is called/the class is instantiated and
    adds a warning to the docstring.

    The optional extra argument will be appended to the deprecation message
    and the docstring. Note: to use this with the default value for extra, put
    in an empty of parentheses:

    Examples
    --------
    >>> from sklearn.utils import deprecated
    >>> deprecated()
    <sklearn.utils.deprecation.deprecated object at ...>
    >>> @deprecated()
    ... def some_function(): pass

    Parameters
    ----------
    extra : str, default=''
          To be added to the deprecation messages.
    """

    # Adapted from https://wiki.python.org/moin/PythonDecoratorLibrary,
    # but with many changes.

    def __init__(self, extra=""):
        self.extra = extra

    def __call__(self, obj):
        """Call method

        Parameters
        ----------
        obj : object
        """
        if isinstance(obj, type):
            return self._decorate_class(obj)
        elif isinstance(obj, property):
            # Note that this is only triggered properly if the `deprecated`
            # decorator is placed before the `property` decorator, like so:
            #
            # @deprecated(msg)
            # @property
            # def deprecated_attribute_(self):
            #     ...
            return self._decorate_property(obj)
        else:
            return self._decorate_fun(obj)

    def _decorate_class(self, cls):
        msg = "Class %s is deprecated" % cls.__name__
        if self.extra:
            msg += "; %s" % self.extra

        new = cls.__new__
        sig = signature(cls)

        def wrapped(cls, *args, **kwargs):
            warnings.warn(msg, category=FutureWarning)
            if new is object.__new__:
                return object.__new__(cls)

            return new(cls, *args, **kwargs)

        cls.__new__ = wrapped

        wrapped.__name__ = "__new__"
        wrapped.deprecated_original = new
        # Restore the original signature, see PEP 362.
        cls.__signature__ = sig

        return cls

    def _decorate_fun(self, fun):
        """Decorate function fun"""

        msg = "Function %s is deprecated" % fun.__name__
        if self.extra:
            msg += "; %s" % self.extra

        @functools.wraps(fun)
        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=FutureWarning)
            return fun(*args, **kwargs)

        # Add a reference to the wrapped function so that we can introspect
        # on function arguments in Python 2 (already works in Python 3)
        wrapped.__wrapped__ = fun

        return wrapped

    def _decorate_property(self, prop):
        msg = self.extra

        @property
        @functools.wraps(prop.fget)
        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=FutureWarning)
            return prop.fget(*args, **kwargs)

        return wrapped


def _is_deprecated(func):
    """Helper to check if func is wrapped by our deprecated decorator"""
    closures = getattr(func, "__closure__", [])
    if closures is None:
        closures = []
    is_deprecated = "deprecated" in "".join(
        [c.cell_contents for c in closures if isinstance(c.cell_contents, str)]
    )
    return is_deprecated


# TODO(1.8): remove force_all_finite and change the default value of ensure_all_finite
# to True (remove None without deprecation).
def _deprecate_force_all_finite(force_all_finite, ensure_all_finite):
    """Helper to deprecate force_all_finite in favor of ensure_all_finite."""
    if force_all_finite != "deprecated":
        warnings.warn(
            "'force_all_finite' was renamed to 'ensure_all_finite' in 1.6 and will be "
            "removed in 1.8.",
            FutureWarning,
        )

        if ensure_all_finite is not None:
            raise ValueError(
                "'force_all_finite' and 'ensure_all_finite' cannot be used together. "
                "Pass `ensure_all_finite` only."
            )

        return force_all_finite

    if ensure_all_finite is None:
        return True

    return ensure_all_finite
