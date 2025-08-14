# SPDX-License-Identifier: MIT

"""
This module defines exceptions and error handling utilities. It is the
recommend path to access ``ConfiguratonError``, ``ConfigurationWarning``, and
``ExceptionGroup``. For backward compatibility, ``ConfigurationError`` is
re-exported in the top-level package.
"""

from __future__ import annotations

import builtins
import contextlib
import dataclasses
import sys
import typing
import warnings

__all__ = [
    "ConfigurationError",
    "ConfigurationWarning",
    "ExceptionGroup",
]


def __dir__() -> list[str]:
    return __all__


class ConfigurationError(Exception):
    """Error in the backend metadata. Has an optional key attribute, which will be non-None
    if the error is related to a single key in the pyproject.toml file."""

    def __init__(self, msg: str, *, key: str | None = None):
        super().__init__(msg)
        self._key = key

    @property
    def key(self) -> str | None:  # pragma: no cover
        return self._key


class ConfigurationWarning(UserWarning):
    """Warnings about backend metadata."""


if sys.version_info >= (3, 11):
    ExceptionGroup = builtins.ExceptionGroup
else:

    class ExceptionGroup(Exception):
        """A minimal implementation of `ExceptionGroup` from Python 3.11.

        Users can replace this with a more complete implementation, such as from
        the exceptiongroup backport package, if better error messages and
        integration with tooling is desired and the addition of a dependency is
        acceptable.
        """

        message: str
        exceptions: list[Exception]

        def __init__(self, message: str, exceptions: list[Exception]) -> None:
            self.message = message
            self.exceptions = exceptions

        def __repr__(self) -> str:
            return f"{self.__class__.__name__}({self.message!r}, {self.exceptions!r})"


@dataclasses.dataclass
class ErrorCollector:
    """
    Collect errors and raise them as a group at the end (if collect_errors is True),
    otherwise raise them immediately.
    """

    collect_errors: bool
    errors: list[Exception] = dataclasses.field(default_factory=list)

    def config_error(
        self,
        msg: str,
        *,
        key: str | None = None,
        got: typing.Any = None,
        got_type: type[typing.Any] | None = None,
        warn: bool = False,
        **kwargs: typing.Any,
    ) -> None:
        """Raise a configuration error, or add it to the error list."""
        msg = msg.format(key=f'"{key}"', **kwargs)
        if got is not None:
            msg = f"{msg} (got {got!r})"
        if got_type is not None:
            msg = f"{msg} (got {got_type.__name__})"

        if warn:
            warnings.warn(msg, ConfigurationWarning, stacklevel=3)
        elif self.collect_errors:
            self.errors.append(ConfigurationError(msg, key=key))
        else:
            raise ConfigurationError(msg, key=key)

    def finalize(self, msg: str) -> None:
        """Raise a group exception if there are any errors."""
        if self.errors:
            raise ExceptionGroup(msg, self.errors)

    @contextlib.contextmanager
    def collect(self) -> typing.Generator[None, None, None]:
        """Support nesting; add any grouped errors to the error list."""
        if self.collect_errors:
            try:
                yield
            except ExceptionGroup as error:
                self.errors.extend(error.exceptions)
        else:
            yield
