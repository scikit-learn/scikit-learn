from __future__ import annotations

import inspect
from collections.abc import Sequence
from functools import wraps
from typing import TYPE_CHECKING, Callable, TypeVar

from polars._utils.various import issue_warning

if TYPE_CHECKING:
    import sys
    from collections.abc import Mapping

    from polars._typing import Ambiguous

    if sys.version_info >= (3, 10):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    P = ParamSpec("P")
    T = TypeVar("T")


USE_EARLIEST_TO_AMBIGUOUS: Mapping[bool, Ambiguous] = {
    True: "earliest",
    False: "latest",
}


def issue_deprecation_warning(message: str, *, version: str) -> None:
    """
    Issue a deprecation warning.

    Parameters
    ----------
    message
        The message associated with the warning.
    version
        The Polars version number in which the warning is first issued.
        This argument is used to help developers determine when to remove the
        deprecated functionality.
    """
    issue_warning(message, DeprecationWarning)


def deprecate_function(
    message: str, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark a function as deprecated."""

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            issue_deprecation_warning(
                f"`{function.__qualname__}` is deprecated. {message}",
                version=version,
            )
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate


def deprecate_renamed_function(
    new_name: str, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark a function as deprecated due to being renamed (or moved)."""
    return deprecate_function(f"It has been renamed to `{new_name}`.", version=version)


def deprecate_renamed_parameter(
    old_name: str, new_name: str, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator to mark a function argument as deprecated due to being renamed.

    Use as follows::

        @deprecate_renamed_parameter("old_name", "new_name", version="0.20.4")
        def myfunc(new_name): ...
    """

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            _rename_keyword_argument(
                old_name, new_name, kwargs, function.__qualname__, version
            )
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate


def _rename_keyword_argument(
    old_name: str,
    new_name: str,
    kwargs: dict[str, object],
    func_name: str,
    version: str,
) -> None:
    """Rename a keyword argument of a function."""
    if old_name in kwargs:
        if new_name in kwargs:
            msg = (
                f"`{func_name!r}` received both `{old_name!r}` and `{new_name!r}` as arguments;"
                f" `{old_name!r}` is deprecated, use `{new_name!r}` instead"
            )
            raise TypeError(msg)
        issue_deprecation_warning(
            f"The argument `{old_name}` for `{func_name}` is deprecated."
            f" It has been renamed to `{new_name}`.",
            version=version,
        )
        kwargs[new_name] = kwargs.pop(old_name)


def deprecate_nonkeyword_arguments(
    allowed_args: list[str] | None = None, message: str | None = None, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator to deprecate the use of non-keyword arguments of a function.

    Parameters
    ----------
    allowed_args
        The names of some first arguments of the decorated function that are allowed to
        be given as positional arguments. Should include "self" when decorating class
        methods. If set to None (default), equal to all arguments that do not have a
        default value.
    message
        Optionally overwrite the default warning message.
    version
        The Polars version number in which the warning is first issued.
        This argument is used to help developers determine when to remove the
        deprecated functionality.
    """

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        old_sig = inspect.signature(function)

        if allowed_args is not None:
            allow_args = allowed_args
        else:
            allow_args = [
                p.name
                for p in old_sig.parameters.values()
                if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
                and p.default is p.empty
            ]

        new_params = [
            p.replace(kind=p.KEYWORD_ONLY)
            if (
                p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
                and p.name not in allow_args
            )
            else p
            for p in old_sig.parameters.values()
        ]
        new_params.sort(key=lambda p: p.kind)

        new_sig = old_sig.replace(parameters=new_params)

        num_allowed_args = len(allow_args)
        if message is None:
            msg_format = (
                f"All arguments of {function.__qualname__}{{except_args}} will be keyword-only in the next breaking release."
                " Use keyword arguments to silence this warning."
            )
            msg = msg_format.format(except_args=_format_argument_list(allow_args))
        else:
            msg = message

        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            if len(args) > num_allowed_args:
                issue_deprecation_warning(msg, version=version)
            return function(*args, **kwargs)

        wrapper.__signature__ = new_sig  # type: ignore[attr-defined]
        return wrapper

    return decorate


def _format_argument_list(allowed_args: list[str]) -> str:
    """Format allowed arguments list for use in the warning message of `deprecate_nonkeyword_arguments`."""  # noqa: W505
    if "self" in allowed_args:
        allowed_args.remove("self")
    if not allowed_args:
        return ""
    elif len(allowed_args) == 1:
        return f" except for {allowed_args[0]!r}"
    else:
        last = allowed_args[-1]
        args = ", ".join([f"{x!r}" for x in allowed_args[:-1]])
        return f" except for {args} and {last!r}"


def deprecate_parameter_as_multi_positional(
    old_name: str, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator to mark a function argument as deprecated due to being made multi-positional.

    Use as follows::

        @deprecate_parameter_as_positional("columns", version="0.20.4")
        def myfunc(*columns): ...
    """  # noqa: W505

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            try:
                arg_value = kwargs.pop(old_name)
            except KeyError:
                return function(*args, **kwargs)

            issue_deprecation_warning(
                f"Passing `{old_name}` as a keyword argument is deprecated."
                " Pass it as a positional argument instead.",
                version=version,
            )

            if not isinstance(arg_value, Sequence) or isinstance(arg_value, str):
                arg_value = (arg_value,)
            elif not isinstance(arg_value, tuple):
                arg_value = tuple(arg_value)

            args = args + arg_value  # type: ignore[assignment]
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate
