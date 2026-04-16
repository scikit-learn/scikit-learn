# SPDX-License-Identifier: MIT

"""
This module contains type definitions for the tables used in the
``pyproject.toml``.  You should either import this at type-check time only, or
make sure ``typing_extensions`` is available for Python 3.10 and below.

Documentation notice: the fields with hyphens are not shown due to a sphinx-autodoc bug.
"""

from __future__ import annotations

import re
import sys
import typing
from typing import (
    Any,
    Dict,
    List,
    Literal,
    TypedDict,
    Union,
)

import packaging.requirements
import packaging.specifiers
import packaging.version

from ._dispatch import get_name, is_typed_dict, keydispatch
from .errors import ConfigurationError, ConfigurationTypeError, SimpleErrorCollector

if sys.version_info < (3, 11):
    if typing.TYPE_CHECKING:
        from typing_extensions import Required
    else:
        try:
            from typing_extensions import Required
        except ModuleNotFoundError:
            V = typing.TypeVar("V")

            class Required:
                def __class_getitem__(cls, item: V) -> V:
                    return item
else:
    from typing import Required


__all__ = [
    "BuildSystemTable",
    "ContactTable",
    "Dynamic",
    "IncludeGroupTable",
    "LicenseTable",
    "ProjectTable",
    "PyProjectTable",
    "ReadmeTable",
    "to_project_table",
]


VALID_ENTRY_POINT = re.compile(r"^\w+(\.\w+)*$")


def __dir__() -> list[str]:
    return __all__


class ContactTable(TypedDict, total=False):
    """
    Can have either name or email.
    """

    name: str
    email: str


class LicenseTable(TypedDict, total=False):
    """
    Can have either text or file. Legacy.
    """

    text: str
    file: str


ReadmeTable = TypedDict(
    "ReadmeTable", {"file": str, "text": str, "content-type": str}, total=False
)

Dynamic = Literal[
    "authors",
    "classifiers",
    "dependencies",
    "description",
    "entry-points",
    "gui-scripts",
    "import-names",
    "import-namespaces",
    "keywords",
    "license",
    "license-files",
    "maintainers",
    "optional-dependencies",
    "readme",
    "requires-python",
    "scripts",
    "urls",
    "version",
]

ProjectTable = TypedDict(
    "ProjectTable",
    {
        "name": Required[str],
        "version": str,
        "description": str,
        "license": Union[LicenseTable, str],
        "license-files": List[str],
        "readme": Union[str, ReadmeTable],
        "requires-python": str,
        "dependencies": List[str],
        "optional-dependencies": Dict[str, List[str]],
        "entry-points": Dict[str, Dict[str, str]],
        "authors": List[ContactTable],
        "maintainers": List[ContactTable],
        "urls": Dict[str, str],
        "classifiers": List[str],
        "keywords": List[str],
        "scripts": Dict[str, str],
        "gui-scripts": Dict[str, str],
        "import-names": List[str],
        "import-namespaces": List[str],
        "dynamic": List[Dynamic],
    },
    total=False,
)

BuildSystemTable = TypedDict(
    "BuildSystemTable",
    {
        "build-backend": str,
        "requires": List[str],
        "backend-path": List[str],
    },
    total=False,
)

# total=False here because this could be
# extended in the future
IncludeGroupTable = TypedDict(
    "IncludeGroupTable",
    {"include-group": str},
    total=False,
)

PyProjectTable = TypedDict(
    "PyProjectTable",
    {
        "build-system": BuildSystemTable,
        "project": ProjectTable,
        "tool": Dict[str, Any],
        "dependency-groups": Dict[str, List[Union[str, IncludeGroupTable]]],
    },
    total=False,
)

T = typing.TypeVar("T")


def join(prefix: str, name: str) -> str:
    """
    Join a prefix and a name.
    """
    if not prefix:
        return name
    return f"{prefix}.{name!r}" if "." in name else f"{prefix}.{name}"


@keydispatch
def validate_via_prefix(
    prefix: str, data: object, error_collector: SimpleErrorCollector
) -> None:
    """
    Validate a TypedDict at runtime.
    """


@validate_via_prefix.register("project")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, dict):
        return

    if "name" not in data:
        msg = f'Field "{prefix}.name" is required if "{prefix}" is present'
        error_collector.error(ConfigurationError(msg, key=f"{prefix}.name"))


@validate_via_prefix.register(r"project\.(authors|maintainers)\[\d+\]")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, dict):
        return

    if "name" not in data and "email" not in data:
        msg = f'Field "{prefix}" must have at least one of "name" or "email" keys'
        error_collector.error(ConfigurationError(msg, key=prefix))

    extra_keys = set(data.keys()) - {"name", "email"}
    if extra_keys:
        extra_keys_list = ", ".join(f'"{k}"' for k in sorted(extra_keys))
        msg = f'Field "{prefix}" contains unexpected keys: {extra_keys_list}'
        error_collector.error(ConfigurationError(msg, key=prefix))


@validate_via_prefix.register(r"project\.license")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, dict):
        return

    if len({"text", "file"} & set(data.keys())) != 1:
        msg = f'Field "{prefix}" must have exactly one of "text" or "file" keys'
        error_collector.error(ConfigurationError(msg, key=prefix))

    extra_keys = set(data.keys()) - {"text", "file"}
    if extra_keys:
        extra_keys_list = ", ".join(f'"{k}"' for k in sorted(extra_keys))
        msg = f'Field "{prefix}" contains unexpected keys: {extra_keys_list}'
        error_collector.error(ConfigurationError(msg, key=prefix))


@validate_via_prefix.register(r"project\.readme")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, dict):
        return
    extra_keys = set(data.keys()) - {"file", "text", "content-type"}
    if extra_keys:
        extra_keys_list = ", ".join(f'"{k}"' for k in sorted(extra_keys))
        msg = f'Field "{prefix}" contains unexpected keys: {extra_keys_list}'
        error_collector.error(ConfigurationError(msg, key=prefix))
    if len({"file", "text"} & set(data.keys())) != 1:
        msg = f'Field "{prefix}" must have exactly one of "file" or "text" keys'
        error_collector.error(ConfigurationError(msg, key=prefix))
    if "content-type" not in data:
        msg = f'Field "{prefix}" is missing required key "content-type"'
        error_collector.error(ConfigurationError(msg, key=prefix))


@validate_via_prefix.register(r"project\.version")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, str):
        return
    try:
        packaging.version.Version(data)
    except packaging.version.InvalidVersion:
        msg = f'Field "{prefix}" is an invalid PEP 440 version string (got {data!r})'
        error_collector.error(ConfigurationError(msg, key=prefix))


@validate_via_prefix.register(
    r"project\.(dependencies|optional-dependencies\.[^\.]+)\[\d+\]"
)
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, str):
        return
    try:
        packaging.requirements.Requirement(data)
    except packaging.requirements.InvalidRequirement as exc:
        msg = f'Field "{prefix}" is an invalid PEP 508 requirement string {data!r} ({exc!r})'
        error_collector.error(ConfigurationError(msg, key=prefix))


@validate_via_prefix.register(r"project\.entry-points")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, dict):
        return
    for name in data:
        if not VALID_ENTRY_POINT.fullmatch(name):
            msg = (
                f'Field "{prefix}" has an invalid key, expecting a key containing'
                f" only alphanumeric, underscore, or dot characters (got {name!r})"
            )
            error_collector.error(ConfigurationTypeError(msg, key=prefix))


@validate_via_prefix.register(r"project\.requires-python")
def _(prefix: str, data: object, error_collector: SimpleErrorCollector) -> None:
    if not isinstance(data, str):
        return
    try:
        packaging.specifiers.SpecifierSet(data)
    except packaging.specifiers.InvalidSpecifier:
        msg = f'Field "{prefix}" is an invalid Python version specifier string (got {data!r})'
        error_collector.error(ConfigurationError(msg, key=prefix))


def _cast_typed_dict(
    cls: type[Any], data: object, prefix: str, error_collector: SimpleErrorCollector
) -> None:
    """
    Runtime cast for TypedDicts.
    """
    if not isinstance(data, dict):
        msg = f'Field "{prefix}" has an invalid type, expecting {get_name(cls)} (got {get_name(type(data))})'
        raise ConfigurationTypeError(msg, key=prefix)

    hints = typing.get_type_hints(cls)
    for key, typ in hints.items():
        if key in data:
            with error_collector.collect(ConfigurationError):
                _cast(
                    typ,
                    data[key],
                    join(prefix, key),
                    error_collector,
                )
        # Required keys could be enforced here on 3.11+ eventually
        # instead of in the validators


def _cast_literal(
    args: tuple[type[Any], ...],
    data: object,
    prefix: str,
    error_collector: SimpleErrorCollector,
) -> None:
    if data not in args:
        arg_names = ", ".join(repr(a) for a in args)
        msg = f'Field "{prefix}" expected one of {arg_names} (got {data!r})'
        error_collector.error(ConfigurationTypeError(msg, key=prefix))


def _cast_list(
    args: tuple[type[Any], ...],
    data: object,
    prefix: str,
    error_collector: SimpleErrorCollector,
) -> None:
    (item_type,) = args
    if not isinstance(data, list):
        msg = f'Field "{prefix}" has an invalid type, expecting list[{get_name(item_type)}] (got {get_name(type(data))})'
        raise ConfigurationTypeError(msg, key=prefix)
    for index, item in enumerate(data):
        with error_collector.collect(ConfigurationError):
            _cast(item_type, item, f"{prefix}[{index}]", error_collector)


def _cast_dict(
    args: tuple[type[Any], ...],
    data: object,
    prefix: str,
    error_collector: SimpleErrorCollector,
) -> None:
    _, value_type = args
    if not isinstance(data, dict):
        msg = f'Field "{prefix}" has an invalid type, expecting dict[str, {get_name(value_type)}] (got {get_name(type(data))})'
        raise ConfigurationTypeError(msg, key=prefix)
    for key, value in data.items():
        with error_collector.collect(ConfigurationError):
            _cast(value_type, value, join(prefix, key), error_collector)


def _cast_union(
    args: tuple[type[Any], ...],
    data: object,
    prefix: str,
    error_collector: SimpleErrorCollector,
) -> None:
    """
    Runtime cast for Union types.

    Checks parent type only (does not check for TypedDict contents), which gives
    better errors. Currently implements dicts and strs only, as that's all
    that's needed.
    """
    for arg in args:
        if arg is str and isinstance(data, str):
            return
        if (arg is dict or is_typed_dict(arg)) and isinstance(data, dict):
            _cast(arg, data, prefix, error_collector)
            return
    # No variant matched, raise a union type error
    arg_names = " | ".join(get_name(a) for a in args)
    msg = f'Field "{prefix}" does not match any of: {arg_names} (got {get_name(type(data))})'
    raise ConfigurationTypeError(msg, key=prefix)


def _cast(
    type_hint: type[Any],
    data: object,
    prefix: str,
    error_collector: SimpleErrorCollector,
) -> None:
    """
    Runtime cast for types.

    Just enough to cover the dicts above (not general or public). Calls validators as well.
    This may raise ConfigurationError even when the error collector is collecting; this is used
    to short-circuit further validation when the type is wrong.
    """
    origin = typing.get_origin(type_hint)
    # Special case Required, needed on 3.10 and older
    if origin is Required:
        (type_hint,) = typing.get_args(type_hint)
        origin = typing.get_origin(type_hint)
    args = typing.get_args(type_hint)

    # Any accepts everything, so no validation
    if type_hint is Any:  # type: ignore[comparison-overlap]
        validate_via_prefix(prefix, data, error_collector)
        return

    # TypedDict
    if is_typed_dict(type_hint):
        validate_via_prefix(prefix, data, error_collector)
        _cast_typed_dict(type_hint, data, prefix, error_collector)
        return

    if origin is typing.Literal:
        validate_via_prefix(prefix, data, error_collector)
        _cast_literal(args, data, prefix, error_collector)
    elif origin is list:
        validate_via_prefix(prefix, data, error_collector)
        _cast_list(args, data, prefix, error_collector)
    elif origin is dict:
        validate_via_prefix(prefix, data, error_collector)
        _cast_dict(args, data, prefix, error_collector)
    elif origin is typing.Union:
        # A union does not run the validator, the selected branch will
        _cast_union(args, data, prefix, error_collector)
    elif isinstance(data, origin or type_hint):
        validate_via_prefix(prefix, data, error_collector)
    else:
        msg = f'Field "{prefix}" has an invalid type, expecting {get_name(type_hint)} (got {get_name(type(data))})'
        error_collector.error(ConfigurationTypeError(msg, key=prefix))


def to_project_table(
    data: dict[str, Any], /, *, collect_errors: bool
) -> PyProjectTable:
    """
    Convert a dict to a PyProjectTable, validating types at runtime.

    Note that only the types that are affected by a TypedDict are validated;
    extra keys are ignored. If it throws, it will be an ExceptionGroup with all
    the ConfigurationTypeError found.
    """
    error_collector = SimpleErrorCollector(collect_errors=collect_errors)
    with error_collector.collect(ConfigurationError):
        _cast(PyProjectTable, data, "", error_collector)
    error_collector.finalize('Errors in "pyproject.toml"')

    return typing.cast("PyProjectTable", data)
