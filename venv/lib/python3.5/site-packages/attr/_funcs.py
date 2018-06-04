from __future__ import absolute_import, division, print_function

import copy

from ._compat import iteritems
from ._make import NOTHING, _obj_setattr, fields
from .exceptions import AttrsAttributeNotFoundError


def asdict(inst, recurse=True, filter=None, dict_factory=dict,
           retain_collection_types=False):
    """
    Return the ``attrs`` attribute values of *inst* as a dict.

    Optionally recurse into other ``attrs``-decorated classes.

    :param inst: Instance of an ``attrs``-decorated class.
    :param bool recurse: Recurse into classes that are also
        ``attrs``-decorated.
    :param callable filter: A callable whose return code determines whether an
        attribute or element is included (``True``) or dropped (``False``).  Is
        called with the :class:`attr.Attribute` as the first argument and the
        value as the second argument.
    :param callable dict_factory: A callable to produce dictionaries from.  For
        example, to produce ordered dictionaries instead of normal Python
        dictionaries, pass in ``collections.OrderedDict``.
    :param bool retain_collection_types: Do not convert to ``list`` when
        encountering an attribute whose type is ``tuple`` or ``set``.  Only
        meaningful if ``recurse`` is ``True``.

    :rtype: return type of *dict_factory*

    :raise attr.exceptions.NotAnAttrsClassError: If *cls* is not an ``attrs``
        class.

    ..  versionadded:: 16.0.0 *dict_factory*
    ..  versionadded:: 16.1.0 *retain_collection_types*
    """
    attrs = fields(inst.__class__)
    rv = dict_factory()
    for a in attrs:
        v = getattr(inst, a.name)
        if filter is not None and not filter(a, v):
            continue
        if recurse is True:
            if has(v.__class__):
                rv[a.name] = asdict(v, recurse=True, filter=filter,
                                    dict_factory=dict_factory)
            elif isinstance(v, (tuple, list, set)):
                cf = v.__class__ if retain_collection_types is True else list
                rv[a.name] = cf([
                    asdict(i, recurse=True, filter=filter,
                           dict_factory=dict_factory)
                    if has(i.__class__) else i
                    for i in v
                ])
            elif isinstance(v, dict):
                df = dict_factory
                rv[a.name] = df((
                    asdict(kk, dict_factory=df) if has(kk.__class__) else kk,
                    asdict(vv, dict_factory=df) if has(vv.__class__) else vv)
                    for kk, vv in iteritems(v))
            else:
                rv[a.name] = v
        else:
            rv[a.name] = v
    return rv


def astuple(inst, recurse=True, filter=None, tuple_factory=tuple,
            retain_collection_types=False):
    """
    Return the ``attrs`` attribute values of *inst* as a tuple.

    Optionally recurse into other ``attrs``-decorated classes.

    :param inst: Instance of an ``attrs``-decorated class.
    :param bool recurse: Recurse into classes that are also
        ``attrs``-decorated.
    :param callable filter: A callable whose return code determines whether an
        attribute or element is included (``True``) or dropped (``False``).  Is
        called with the :class:`attr.Attribute` as the first argument and the
        value as the second argument.
    :param callable tuple_factory: A callable to produce tuples from.  For
        example, to produce lists instead of tuples.
    :param bool retain_collection_types: Do not convert to ``list``
        or ``dict`` when encountering an attribute which type is
        ``tuple``, ``dict`` or ``set``.  Only meaningful if ``recurse`` is
        ``True``.

    :rtype: return type of *tuple_factory*

    :raise attr.exceptions.NotAnAttrsClassError: If *cls* is not an ``attrs``
        class.

    ..  versionadded:: 16.2.0
    """
    attrs = fields(inst.__class__)
    rv = []
    retain = retain_collection_types  # Very long. :/
    for a in attrs:
        v = getattr(inst, a.name)
        if filter is not None and not filter(a, v):
            continue
        if recurse is True:
            if has(v.__class__):
                rv.append(astuple(v, recurse=True, filter=filter,
                                  tuple_factory=tuple_factory,
                                  retain_collection_types=retain))
            elif isinstance(v, (tuple, list, set)):
                cf = v.__class__ if retain is True else list
                rv.append(cf([
                    astuple(j, recurse=True, filter=filter,
                            tuple_factory=tuple_factory,
                            retain_collection_types=retain)
                    if has(j.__class__) else j
                    for j in v
                ]))
            elif isinstance(v, dict):
                df = v.__class__ if retain is True else dict
                rv.append(df(
                        (
                            astuple(
                                kk,
                                tuple_factory=tuple_factory,
                                retain_collection_types=retain
                            ) if has(kk.__class__) else kk,
                            astuple(
                                vv,
                                tuple_factory=tuple_factory,
                                retain_collection_types=retain
                            ) if has(vv.__class__) else vv
                        )
                        for kk, vv in iteritems(v)))
            else:
                rv.append(v)
        else:
            rv.append(v)
    return rv if tuple_factory is list else tuple_factory(rv)


def has(cls):
    """
    Check whether *cls* is a class with ``attrs`` attributes.

    :param type cls: Class to introspect.
    :raise TypeError: If *cls* is not a class.

    :rtype: :class:`bool`
    """
    return getattr(cls, "__attrs_attrs__", None) is not None


def assoc(inst, **changes):
    """
    Copy *inst* and apply *changes*.

    :param inst: Instance of a class with ``attrs`` attributes.
    :param changes: Keyword changes in the new copy.

    :return: A copy of inst with *changes* incorporated.

    :raise attr.exceptions.AttrsAttributeNotFoundError: If *attr_name* couldn't
        be found on *cls*.
    :raise attr.exceptions.NotAnAttrsClassError: If *cls* is not an ``attrs``
        class.

    ..  deprecated:: 17.1.0
        Use :func:`evolve` instead.
    """
    import warnings
    warnings.warn("assoc is deprecated and will be removed after 2018/01.",
                  DeprecationWarning, stacklevel=2)
    new = copy.copy(inst)
    attrs = fields(inst.__class__)
    for k, v in iteritems(changes):
        a = getattr(attrs, k, NOTHING)
        if a is NOTHING:
            raise AttrsAttributeNotFoundError(
                "{k} is not an attrs attribute on {cl}."
                .format(k=k, cl=new.__class__)
            )
        _obj_setattr(new, k, v)
    return new


def evolve(inst, **changes):
    """
    Create a new instance, based on *inst* with *changes* applied.

    :param inst: Instance of a class with ``attrs`` attributes.
    :param changes: Keyword changes in the new copy.

    :return: A copy of inst with *changes* incorporated.

    :raise TypeError: If *attr_name* couldn't be found in the class
        ``__init__``.
    :raise attr.exceptions.NotAnAttrsClassError: If *cls* is not an ``attrs``
        class.

    ..  versionadded:: 17.1.0
    """
    cls = inst.__class__
    attrs = fields(cls)
    for a in attrs:
        if not a.init:
            continue
        attr_name = a.name  # To deal with private attributes.
        init_name = attr_name if attr_name[0] != "_" else attr_name[1:]
        if init_name not in changes:
            changes[init_name] = getattr(inst, attr_name)
    return cls(**changes)
