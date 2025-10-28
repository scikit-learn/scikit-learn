from __future__ import annotations

import ast
import inspect
import sys
from collections import defaultdict
from collections.abc import Sequence
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, TypeVar, get_args

from polars._typing import DeprecationType

if sys.version_info >= (3, 13):
    from warnings import deprecated
else:
    try:
        from typing_extensions import deprecated
    except ImportError:

        def deprecated(  # type: ignore[no-redef]
            message: str,
        ) -> Callable[[Callable[P, T]], Callable[P, T]]:
            return _deprecate_function(message)


from polars._utils.various import issue_warning

if TYPE_CHECKING:
    from collections.abc import Mapping

    if sys.version_info >= (3, 10):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec
    from polars._typing import Ambiguous

    P = ParamSpec("P")
    T = TypeVar("T")

USE_EARLIEST_TO_AMBIGUOUS: Mapping[bool, Ambiguous] = {
    True: "earliest",
    False: "latest",
}


def issue_deprecation_warning(message: str, *, version: str = "") -> None:
    """
    Issue a deprecation warning.

    Parameters
    ----------
    message
        The message associated with the warning.
    version
        The version in which deprecation occurred
        (if the version number was not already included in `message`).
    """
    if version:
        message = f"{message.strip()}\n(Deprecated in version {version})"
    issue_warning(message, DeprecationWarning)


def _deprecate_function(message: str) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark a function as deprecated."""

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            issue_deprecation_warning(message)
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        wrapper.__deprecated__ = message  # type: ignore[attr-defined]
        return wrapper

    return decorate


def deprecate_streaming_parameter() -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark `streaming` argument as deprecated due to being renamed."""

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            if "streaming" in kwargs:
                issue_deprecation_warning(
                    "the `streaming` parameter was deprecated in 1.25.0; use `engine` instead."
                )
                if kwargs["streaming"]:
                    kwargs["engine"] = "streaming"
                elif "engine" not in kwargs:
                    kwargs["engine"] = "in-memory"

                del kwargs["streaming"]

            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate


def deprecate_renamed_parameter(
    old_name: str, new_name: str, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator to mark a function parameter as deprecated due to being renamed.

    Use as follows:

        @deprecate_renamed_parameter("old_name", new_name="new_name")
        def myfunc(new_name): ...

    Ensure that you also update the function docstring with a note about the
    deprecation, specifically adding a `.. versionchanged:: 0.0.0` directive
    that states which parameter was renamed to which new name and in which
    version the rename happened.
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
            is_deprecated = (
                f"was deprecated in version {version}" if version else "is deprecated"
            )
            msg = (
                f"`{func_name!r}` received both `{old_name!r}` and `{new_name!r}` as arguments;"
                f" `{old_name!r}` {is_deprecated}, use `{new_name!r}` instead"
            )
            raise TypeError(msg)

        in_version = f" in version {version}" if version else ""
        issue_deprecation_warning(
            f"the argument `{old_name}` for `{func_name}` is deprecated. "
            f"It was renamed to `{new_name}`{in_version}."
        )
        kwargs[new_name] = kwargs.pop(old_name)


def deprecate_nonkeyword_arguments(
    allowed_args: list[str] | None = None, message: str | None = None, *, version: str
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator for deprecating the use of non-keyword arguments in a function.

    Use as follows:

        @deprecate_nonkeyword_arguments(allowed_args=["self", "val"], version="1.0.0")
        def myfunc(self, val: int = 0, other: int: = 0): ...

    Ensure that you also update the function docstring with a note about the
    deprecation, specifically adding a `.. versionchanged:: 0.0.0` directive
    that states that we now expect keyword args and in which version this
    update happened.

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
                f"all arguments of {function.__qualname__}{{except_args}} will be keyword-only in the next breaking release."
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
    old_name: str,
) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator to mark a function argument as deprecated due to being made multi-positional.

    Use as follows:

        @deprecate_parameter_as_multi_positional("columns")
        def myfunc(*columns): ...

    Ensure that you also update the function docstring with a note about the
    deprecation, specifically adding a `.. versionchanged:: 0.0.0` directive
    that states that we now expect positional args and in which version this
    update happened.
    """  # noqa: W505

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            try:
                arg_value = kwargs.pop(old_name)
            except KeyError:
                return function(*args, **kwargs)

            issue_deprecation_warning(
                f"passing `{old_name}` as a keyword argument is deprecated."
                " Pass it as a positional argument instead."
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


def _find_deprecated_functions(
    source: str, module_path: str
) -> defaultdict[str, list[str]]:
    tree = ast.parse(source)
    object_path: list[str] = []

    def deprecated(decorator: Any) -> str:
        if isinstance(decorator, ast.Name):
            return decorator.id if "deprecate" in decorator.id else ""
        elif isinstance(decorator, ast.Call):
            return deprecated(decorator.func)
        return ""

    def qualified_name(func_name: str) -> str:
        return ".".join([module_path, *object_path, func_name])

    results = defaultdict(list)

    class FunctionVisitor(ast.NodeVisitor):
        def visit_ClassDef(self, node: Any) -> None:
            object_path.append(node.name)
            self.generic_visit(node)
            object_path.pop()

        def visit_FunctionDef(self, node: Any) -> None:
            if any((decorator_name := deprecated(d)) for d in node.decorator_list):
                key = decorator_name.removeprefix("deprecate_").replace(
                    "deprecated", "function"
                )
                results[key].append(qualified_name(node.name))
            self.generic_visit(node)

        visit_AsyncFunctionDef = visit_FunctionDef

    FunctionVisitor().visit(tree)
    return results


def identify_deprecations(*types: DeprecationType) -> dict[str, list[str]]:
    """
    Return a dict identifying functions/methods that are deprecated in some way.

    Parameters
    ----------
    *types
        The types of deprecations to identify.
        If empty, all types are returned; recognised values are:
            - "function"
            - "renamed_parameter"
            - "streaming_parameter"
            - "nonkeyword_arguments"
            - "parameter_as_multi_positional"

    Examples
    --------
    >>> from polars._utils.deprecation import identify_deprecations
    >>> identify_deprecations("streaming_parameter")  # doctest: +IGNORE_RESULT
    {'streaming_parameter': [
        'functions.lazy.collect_all',
        'functions.lazy.collect_all_async',
        'lazyframe.frame.LazyFrame.collect',
        'lazyframe.frame.LazyFrame.collect_async',
        'lazyframe.frame.LazyFrame.explain',
        'lazyframe.frame.LazyFrame.show_graph',
    ]}
    """
    valid_types = set(get_args(DeprecationType))
    for tp in types:
        if tp not in valid_types:
            msg = (
                f"unrecognised deprecation type {tp!r}.\n"
                f"Expected one (or more) of {repr(sorted(valid_types))[1:-1]}"
            )
            raise ValueError(msg)

    package_path = Path(sys.modules["polars"].__file__).parent  # type: ignore[arg-type]
    results = defaultdict(list)

    for py_file in package_path.rglob("*.py"):
        rel_path = py_file.relative_to(package_path)
        module_path = ".".join(rel_path.parts).removesuffix(".py")
        with py_file.open("r", encoding="utf-8") as src:
            for deprecation_type, func_names in _find_deprecated_functions(
                source=src.read(),
                module_path=module_path,
            ).items():
                if deprecation_type not in valid_types:
                    # note: raising here implies we have a new deprecation function
                    # that should be added to the DeprecationType type alias
                    msg = f"unrecognised deprecation type {tp!r}.\n"
                    raise ValueError(msg)

                results[deprecation_type].extend(func_names)

    return {
        dep: sorted(results[dep])
        for dep in sorted(results)
        if not types or dep in types
    }


__all__ = [
    "deprecate_nonkeyword_arguments",
    "deprecate_parameter_as_multi_positional",
    "deprecate_renamed_parameter",
    "deprecate_streaming_parameter",
    "deprecated",
    "identify_deprecations",
]
