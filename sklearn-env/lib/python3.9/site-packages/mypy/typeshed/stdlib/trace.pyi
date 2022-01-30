import sys
import types
from _typeshed import StrPath
from typing import Any, Callable, Mapping, Optional, Sequence, Tuple, TypeVar
from typing_extensions import ParamSpec

_T = TypeVar("_T")
_P = ParamSpec("_P")
_localtrace = Callable[[types.FrameType, str, Any], Callable[..., Any]]
_fileModuleFunction = Tuple[str, Optional[str], str]

class CoverageResults:
    def __init__(
        self,
        counts: dict[tuple[str, int], int] | None = ...,
        calledfuncs: dict[_fileModuleFunction, int] | None = ...,
        infile: StrPath | None = ...,
        callers: dict[tuple[_fileModuleFunction, _fileModuleFunction], int] | None = ...,
        outfile: StrPath | None = ...,
    ) -> None: ...  # undocumented
    def update(self, other: CoverageResults) -> None: ...
    def write_results(self, show_missing: bool = ..., summary: bool = ..., coverdir: StrPath | None = ...) -> None: ...
    def write_results_file(
        self, path: StrPath, lines: Sequence[str], lnotab: Any, lines_hit: Mapping[int, int], encoding: str | None = ...
    ) -> tuple[int, int]: ...
    def is_ignored_filename(self, filename: str) -> bool: ...  # undocumented

class Trace:
    def __init__(
        self,
        count: int = ...,
        trace: int = ...,
        countfuncs: int = ...,
        countcallers: int = ...,
        ignoremods: Sequence[str] = ...,
        ignoredirs: Sequence[str] = ...,
        infile: StrPath | None = ...,
        outfile: StrPath | None = ...,
        timing: bool = ...,
    ) -> None: ...
    def run(self, cmd: str | types.CodeType) -> None: ...
    def runctx(
        self, cmd: str | types.CodeType, globals: Mapping[str, Any] | None = ..., locals: Mapping[str, Any] | None = ...
    ) -> None: ...
    if sys.version_info >= (3, 9):
        def runfunc(self, __func: Callable[_P, _T], *args: _P.args, **kw: _P.kwargs) -> _T: ...  # type: ignore
    else:
        def runfunc(self, func: Callable[_P, _T], *args: _P.args, **kw: _P.kwargs) -> _T: ...  # type: ignore
    def file_module_function_of(self, frame: types.FrameType) -> _fileModuleFunction: ...
    def globaltrace_trackcallers(self, frame: types.FrameType, why: str, arg: Any) -> None: ...
    def globaltrace_countfuncs(self, frame: types.FrameType, why: str, arg: Any) -> None: ...
    def globaltrace_lt(self, frame: types.FrameType, why: str, arg: Any) -> None: ...
    def localtrace_trace_and_count(self, frame: types.FrameType, why: str, arg: Any) -> _localtrace: ...
    def localtrace_trace(self, frame: types.FrameType, why: str, arg: Any) -> _localtrace: ...
    def localtrace_count(self, frame: types.FrameType, why: str, arg: Any) -> _localtrace: ...
    def results(self) -> CoverageResults: ...
