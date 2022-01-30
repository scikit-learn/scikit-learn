import sys
import unittest.case
import unittest.loader
import unittest.result
import unittest.suite
from types import ModuleType
from typing import Any, Iterable, Protocol, Type

class _TestRunner(Protocol):
    def run(self, test: unittest.suite.TestSuite | unittest.case.TestCase) -> unittest.result.TestResult: ...

# not really documented
class TestProgram:
    result: unittest.result.TestResult
    module: None | str | ModuleType
    verbosity: int
    failfast: bool | None
    catchbreak: bool | None
    buffer: bool | None
    progName: str | None
    warnings: str | None

    if sys.version_info >= (3, 7):
        testNamePatterns: list[str] | None
    def __init__(
        self,
        module: None | str | ModuleType = ...,
        defaultTest: str | Iterable[str] | None = ...,
        argv: list[str] | None = ...,
        testRunner: Type[_TestRunner] | _TestRunner | None = ...,
        testLoader: unittest.loader.TestLoader = ...,
        exit: bool = ...,
        verbosity: int = ...,
        failfast: bool | None = ...,
        catchbreak: bool | None = ...,
        buffer: bool | None = ...,
        warnings: str | None = ...,
        *,
        tb_locals: bool = ...,
    ) -> None: ...
    def usageExit(self, msg: Any = ...) -> None: ...
    def parseArgs(self, argv: list[str]) -> None: ...
    if sys.version_info >= (3, 7):
        def createTests(self, from_discovery: bool = ..., Loader: unittest.loader.TestLoader | None = ...) -> None: ...
    else:
        def createTests(self) -> None: ...
    def runTests(self) -> None: ...  # undocumented

main = TestProgram
