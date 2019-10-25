import warnings
import functools
import sys

__all__ = ["deprecated"]


class deprecated:
    """Decorator to mark a function or class as deprecated.

    Issue a warning when the function is called/the class is instantiated and
    adds a warning to the docstring.

    The optional extra argument will be appended to the deprecation message
    and the docstring. Note: to use this with the default value for extra, put
    in an empty of parentheses:

    >>> from sklearn.utils import deprecated
    >>> deprecated()
    <sklearn.utils.deprecation.deprecated object at ...>

    >>> @deprecated()
    ... def some_function(): pass

    Parameters
    ----------
    extra : string
          to be added to the deprecation messages
    """

    # Adapted from https://wiki.python.org/moin/PythonDecoratorLibrary,
    # but with many changes.

    def __init__(self, extra=''):
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
            # Note that this is only triggered properly if the `property`
            # decorator comes before the `deprecated` decorator, like so:
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

        # FIXME: we should probably reset __new__ for full generality
        init = cls.__init__

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return init(*args, **kwargs)
        cls.__init__ = wrapped

        wrapped.__name__ = '__init__'
        wrapped.__doc__ = self._update_doc(init.__doc__)
        wrapped.deprecated_original = init

        return cls

    def _decorate_fun(self, fun):
        """Decorate function fun"""

        msg = "Function %s is deprecated" % fun.__name__
        if self.extra:
            msg += "; %s" % self.extra

        @functools.wraps(fun)
        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return fun(*args, **kwargs)

        wrapped.__doc__ = self._update_doc(wrapped.__doc__)
        # Add a reference to the wrapped function so that we can introspect
        # on function arguments in Python 2 (already works in Python 3)
        wrapped.__wrapped__ = fun

        return wrapped

    def _decorate_property(self, prop):
        msg = self.extra

        @property
        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return prop.fget(*args, **kwargs)

        return wrapped

    def _update_doc(self, olddoc):
        newdoc = "DEPRECATED"
        if self.extra:
            newdoc = "%s: %s" % (newdoc, self.extra)
        if olddoc:
            newdoc = "%s\n\n%s" % (newdoc, olddoc)
        return newdoc


def _is_deprecated(func):
    """Helper to check if func is wraped by our deprecated decorator"""
    closures = getattr(func, '__closure__', [])
    if closures is None:
        closures = []
    is_deprecated = ('deprecated' in ''.join([c.cell_contents
                                              for c in closures
                     if isinstance(c.cell_contents, str)]))
    return is_deprecated


def _raise_dep_warning_if_not_pytest(deprecated_path, correct_path):

    # Raise a deprecation warning with standardized deprecation message.
    # Useful because we are now deprecating # anything that isn't explicitly
    # in an __init__ file.

    # TODO: remove in 0.24 since this shouldn't be needed anymore.

    message = (
        "The {deprecated_path} module is  deprecated in version "
        "0.22 and will be removed in version 0.24. "
        "The corresponding classes / functions "
        "should instead be imported from {correct_path}. "
        "Anything that cannot be imported from {correct_path} is now "
        "part of the private API."
    ).format(deprecated_path=deprecated_path, correct_path=correct_path)

    warnings.warn(message, DeprecationWarning)
