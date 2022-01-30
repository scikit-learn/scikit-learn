import unittest.result
from typing import Callable, TypeVar, overload
from typing_extensions import ParamSpec

_P = ParamSpec("_P")
_T = TypeVar("_T")

def installHandler() -> None: ...
def registerResult(result: unittest.result.TestResult) -> None: ...
def removeResult(result: unittest.result.TestResult) -> bool: ...
@overload
def removeHandler(method: None = ...) -> None: ...
@overload
def removeHandler(method: Callable[_P, _T]) -> Callable[_P, _T]: ...  # type: ignore
