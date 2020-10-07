"""Docstring utilities."""


def inject_docstring(**kwargs):
    """Class or function decorator allowing to perform string substitution.

    Parameters
    ----------
    **kwargs
        The key corresponds to the value present in the docstring which should
        be substituted while the value is the value injected.

    Examples
    --------
    >>> my_docstring = "my docstring to inject"
    >>> from sklearn.utils._docstring import inject_docstring
    >>> @inject_docstring(key_to_replace=my_docstring)
    ... def func(x):
    ...     "{key_to_replace}"
    ...     pass
    >>> @inject_docstring(key_to_replace=my_docstring)
    ... class Klazz:
    ...     "{key_to_replace}"
    ...     pass
    """
    def decorate(obj):
        obj.__doc__ = obj.__doc__.format(**kwargs)
        return obj
    return decorate
