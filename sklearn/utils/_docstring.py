"""Docstring utilities."""


class Substitution:
    """Decorate a function's or a class' docstring to perform string
    substitution on it.

    Examples
    --------
    >>> my_docstring = "my docstring to inject"
    >>> from sklearn.utils._docstring import Substitution
    >>> @Substitution(injected_docstring=my_docstring)
    ... def func(x):
    ...     "{injected_docstring}"
    ...     pass
    """

    def __init__(self, *args, **kwargs):
        if args and kwargs:
            raise AssertionError("Only positional or keyword args are allowed")

        self.params = args or kwargs

    def __call__(self, obj):
        obj.__doc__ = obj.__doc__.format(**self.params)
        return obj
