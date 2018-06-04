from __future__ import absolute_import, division, print_function


class FrozenInstanceError(AttributeError):
    """
    A frozen/immutable instance has been attempted to be modified.

    It mirrors the behavior of ``namedtuples`` by using the same error message
    and subclassing :exc:`AttributeError`.

    .. versionadded:: 16.1.0
    """
    msg = "can't set attribute"
    args = [msg]


class AttrsAttributeNotFoundError(ValueError):
    """
    An ``attrs`` function couldn't find an attribute that the user asked for.

    .. versionadded:: 16.2.0
    """


class NotAnAttrsClassError(ValueError):
    """
    A non-``attrs`` class has been passed into an ``attrs`` function.

    .. versionadded:: 16.2.0
    """


class DefaultAlreadySetError(RuntimeError):
    """
    A default has been set using ``attr.ib()`` and is attempted to be reset
    using the decorator.

    .. versionadded:: 17.1.0
    """


class UnannotatedAttributeError(RuntimeError):
    """
    A class with ``auto_attribs=True`` has an ``attr.ib()`` without a type
    annotation.

    .. versionadded:: 17.3.0
    """
