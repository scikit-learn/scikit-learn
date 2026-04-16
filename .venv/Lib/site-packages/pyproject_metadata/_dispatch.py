# SPDX-License-Identifier: MIT

"""
This module contains dispatchers for use in project_table.py.
"""

from __future__ import annotations

import functools
import re
import sys
import typing
from typing import Callable, Generic

if sys.version_info < (3, 10):
    if typing.TYPE_CHECKING:
        from typing_extensions import Concatenate, ParamSpec
    else:
        ParamSpec = typing.TypeVar
        Concatenate = object
else:
    from typing import Concatenate, ParamSpec


__all__ = [
    "get_name",
    "is_typed_dict",
    "keydispatch",
]


def __dir__() -> list[str]:
    return __all__


P = ParamSpec("P")
R = typing.TypeVar("R")
T = typing.TypeVar("T")


class KeyDispatcher(Generic[P, R]):
    def __init__(self, func: Callable[Concatenate[str, P], R]) -> None:
        self.default: Callable[Concatenate[str, P], R] = func
        self.registry: dict[str, Callable[Concatenate[str, P], R]] = {}
        functools.update_wrapper(self, func)

    def register(
        self, value: str
    ) -> Callable[[Callable[Concatenate[str, P], R]], Callable[Concatenate[str, P], R]]:
        """Register a function for an exact value."""

        def decorator(
            func: Callable[Concatenate[str, P], R],
        ) -> Callable[Concatenate[str, P], R]:
            self.registry[value] = func
            return func

        return decorator

    def dispatch(self, value: str) -> Callable[Concatenate[str, P], R]:
        """Return the registered implementation or the default."""
        results = [key for key in self.registry if re.fullmatch(key, value, re.ASCII)]
        result = results[0] if results else ""
        return self.registry.get(result, self.default)

    def __call__(self, value: str, *args: P.args, **kwargs: P.kwargs) -> R:
        impl = self.dispatch(value)
        return impl(value, *args, **kwargs)


def keydispatch(func: Callable[Concatenate[str, P], R]) -> KeyDispatcher[P, R]:
    """Decorate a function into a KeyDispatcher."""
    return KeyDispatcher(func)


def is_typed_dict(type_hint: object) -> bool:
    if sys.version_info >= (3, 10):
        return typing.is_typeddict(type_hint)
    return hasattr(type_hint, "__annotations__") and hasattr(type_hint, "__total__")


def get_name(type_hint: type[object]) -> str:
    """
    Get the name of a type hint as a readable modern Python type.
    """
    if origin := typing.get_origin(type_hint):
        if args := typing.get_args(type_hint):
            arg_names = ", ".join(get_name(a) for a in args)
            return f"{origin.__name__}[{arg_names}]"
        return origin.__name__  # type: ignore[no-any-return]
    return type_hint.__name__
