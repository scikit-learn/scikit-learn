"""Sphinx deprecation classes and utilities."""

from __future__ import annotations

import warnings


class RemovedInSphinx90Warning(DeprecationWarning):
    pass


class RemovedInSphinx10Warning(PendingDeprecationWarning):
    pass


RemovedInNextVersionWarning = RemovedInSphinx90Warning


def _deprecation_warning(
    module: str,
    attribute: str,
    canonical_name: str = '',
    *,
    remove: tuple[int, int],
    raises: bool = False,
) -> None:
    """Helper function for module-level deprecations using ``__getattr__``.

    :param module: The module containing a deprecated object.
    :param attribute: The name of the deprecated object.
    :param canonical_name: Optional fully-qualified name for its replacement.
    :param remove: Target version for removal.
    :param raises: Indicate whether to raise an exception instead of a warning.

    When *raises* is ``True``, an :exc:`AttributeError` is raised instead
    of emitting a warning so that it is easy to locate deprecated objects
    in tests that could suppress deprecation warnings.

    Usage::

       # deprecated name -> (object to return, canonical path or empty string, removal version)
       _DEPRECATED_OBJECTS = {
           'deprecated_name': (
               object_to_return,
               'fully_qualified_replacement_name',
               (9, 0),
           ),
       }


       def __getattr__(name: str) -> Any:
           if name not in _DEPRECATED_OBJECTS:
               msg = f'module {__name__!r} has no attribute {name!r}'
               raise AttributeError(msg)

           from sphinx.deprecation import _deprecation_warning

           deprecated_object, canonical_name, remove = _DEPRECATED_OBJECTS[name]
           _deprecation_warning(__name__, name, canonical_name, remove=remove)
           return deprecated_object
    """
    if remove == (9, 0):
        warning_class: type[Warning] = RemovedInSphinx90Warning
    elif remove == (10, 0):
        warning_class = RemovedInSphinx10Warning
    else:
        msg = f'removal version {remove!r} is invalid!'
        raise RuntimeError(msg)

    qualname = f'{module}.{attribute}'
    if canonical_name:
        message = (
            f'The alias {qualname!r} is deprecated, use {canonical_name!r} instead.'
        )
    else:
        message = f'{qualname!r} is deprecated.'

    if raises:
        raise AttributeError(message)

    message = f'{message} Check CHANGES for Sphinx API modifications.'
    warnings.warn(message, warning_class, stacklevel=3)
