"""
Commonly useful validators.
"""

from __future__ import absolute_import, division, print_function

from ._make import _AndValidator, and_, attrib, attrs


__all__ = [
    "and_",
    "in_",
    "instance_of",
    "optional",
    "provides",
]


@attrs(repr=False, slots=True, hash=True)
class _InstanceOfValidator(object):
    type = attrib()

    def __call__(self, inst, attr, value):
        """
        We use a callable class to be able to change the ``__repr__``.
        """
        if not isinstance(value, self.type):
            raise TypeError(
                "'{name}' must be {type!r} (got {value!r} that is a "
                "{actual!r})."
                .format(name=attr.name, type=self.type,
                        actual=value.__class__, value=value),
                attr, self.type, value,
            )

    def __repr__(self):
        return (
            "<instance_of validator for type {type!r}>"
            .format(type=self.type)
        )


def instance_of(type):
    """
    A validator that raises a :exc:`TypeError` if the initializer is called
    with a wrong type for this particular attribute (checks are performed using
    :func:`isinstance` therefore it's also valid to pass a tuple of types).

    :param type: The type to check for.
    :type type: type or tuple of types

    :raises TypeError: With a human readable error message, the attribute
        (of type :class:`attr.Attribute`), the expected type, and the value it
        got.
    """
    return _InstanceOfValidator(type)


@attrs(repr=False, slots=True, hash=True)
class _ProvidesValidator(object):
    interface = attrib()

    def __call__(self, inst, attr, value):
        """
        We use a callable class to be able to change the ``__repr__``.
        """
        if not self.interface.providedBy(value):
            raise TypeError(
                "'{name}' must provide {interface!r} which {value!r} "
                "doesn't."
                .format(name=attr.name, interface=self.interface, value=value),
                attr, self.interface, value,
            )

    def __repr__(self):
        return (
            "<provides validator for interface {interface!r}>"
            .format(interface=self.interface)
        )


def provides(interface):
    """
    A validator that raises a :exc:`TypeError` if the initializer is called
    with an object that does not provide the requested *interface* (checks are
    performed using ``interface.providedBy(value)`` (see `zope.interface
    <https://zopeinterface.readthedocs.io/en/latest/>`_).

    :param zope.interface.Interface interface: The interface to check for.

    :raises TypeError: With a human readable error message, the attribute
        (of type :class:`attr.Attribute`), the expected interface, and the
        value it got.
    """
    return _ProvidesValidator(interface)


@attrs(repr=False, slots=True, hash=True)
class _OptionalValidator(object):
    validator = attrib()

    def __call__(self, inst, attr, value):
        if value is None:
            return

        self.validator(inst, attr, value)

    def __repr__(self):
        return (
            "<optional validator for {what} or None>"
            .format(what=repr(self.validator))
        )


def optional(validator):
    """
    A validator that makes an attribute optional.  An optional attribute is one
    which can be set to ``None`` in addition to satisfying the requirements of
    the sub-validator.

    :param validator: A validator (or a list of validators) that is used for
        non-``None`` values.
    :type validator: callable or :class:`list` of callables.

    .. versionadded:: 15.1.0
    .. versionchanged:: 17.1.0 *validator* can be a list of validators.
    """
    if isinstance(validator, list):
        return _OptionalValidator(_AndValidator(validator))
    return _OptionalValidator(validator)


@attrs(repr=False, slots=True, hash=True)
class _InValidator(object):
    options = attrib()

    def __call__(self, inst, attr, value):
        if value not in self.options:
            raise ValueError(
                "'{name}' must be in {options!r} (got {value!r})"
                .format(name=attr.name, options=self.options, value=value)
            )

    def __repr__(self):
        return (
            "<in_ validator with options {options!r}>"
            .format(options=self.options)
        )


def in_(options):
    """
    A validator that raises a :exc:`ValueError` if the initializer is called
    with a value that does not belong in the options provided.  The check is
    performed using ``value in options``.

    :param options: Allowed options.
    :type options: list, tuple, :class:`enum.Enum`, ...

    :raises ValueError: With a human readable error message, the attribute (of
       type :class:`attr.Attribute`), the expected options, and the value it
       got.

    .. versionadded:: 17.1.0
    """
    return _InValidator(options)
