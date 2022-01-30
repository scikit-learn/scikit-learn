import threading
from _typeshed import StrPath, SupportsWrite
from time import struct_time
from types import FrameType, TracebackType
from typing import (
    IO,
    Any,
    Callable,
    Dict,
    Generic,
    List,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Text,
    Tuple,
    TypeVar,
    Union,
    overload,
)

_SysExcInfoType = Union[Tuple[type, BaseException, Optional[TracebackType]], Tuple[None, None, None]]
_ExcInfoType = Union[None, bool, _SysExcInfoType]
_ArgsType = Union[Tuple[Any, ...], Mapping[str, Any]]
_FilterType = Union[Filter, Callable[[LogRecord], int]]
_Level = Union[int, Text]

raiseExceptions: bool
logThreads: bool
logMultiprocessing: bool
logProcesses: bool
_srcfile: str | None

def currentframe() -> FrameType: ...

_levelNames: Dict[int | str, str | int]  # Union[int:str, str:int]

class Filterer(object):
    filters: List[Filter]
    def __init__(self) -> None: ...
    def addFilter(self, filter: _FilterType) -> None: ...
    def removeFilter(self, filter: _FilterType) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...

class Logger(Filterer):
    name: str
    level: int
    parent: Logger | PlaceHolder
    propagate: bool
    handlers: List[Handler]
    disabled: int
    def __init__(self, name: str, level: _Level = ...) -> None: ...
    def setLevel(self, level: _Level) -> None: ...
    def isEnabledFor(self, level: int) -> bool: ...
    def getEffectiveLevel(self) -> int: ...
    def getChild(self, suffix: str) -> Logger: ...
    def debug(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def info(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def warning(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    warn = warning
    def error(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def critical(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    fatal = critical
    def log(
        self, level: int, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def exception(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def _log(
        self, level: int, msg: Any, args: _ArgsType, exc_info: _ExcInfoType | None = ..., extra: Dict[str, Any] | None = ...
    ) -> None: ...  # undocumented
    def filter(self, record: LogRecord) -> bool: ...
    def addHandler(self, hdlr: Handler) -> None: ...
    def removeHandler(self, hdlr: Handler) -> None: ...
    def findCaller(self) -> Tuple[str, int, str]: ...
    def handle(self, record: LogRecord) -> None: ...
    def makeRecord(
        self,
        name: str,
        level: int,
        fn: str,
        lno: int,
        msg: Any,
        args: _ArgsType,
        exc_info: _SysExcInfoType | None,
        func: str | None = ...,
        extra: Mapping[str, Any] | None = ...,
    ) -> LogRecord: ...

CRITICAL: int
FATAL: int
ERROR: int
WARNING: int
WARN: int
INFO: int
DEBUG: int
NOTSET: int

class Handler(Filterer):
    level: int  # undocumented
    formatter: Formatter | None  # undocumented
    lock: threading.Lock | None  # undocumented
    name: str | None  # undocumented
    def __init__(self, level: _Level = ...) -> None: ...
    def createLock(self) -> None: ...
    def acquire(self) -> None: ...
    def release(self) -> None: ...
    def setLevel(self, level: _Level) -> None: ...
    def setFormatter(self, fmt: Formatter) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...
    def flush(self) -> None: ...
    def close(self) -> None: ...
    def handle(self, record: LogRecord) -> None: ...
    def handleError(self, record: LogRecord) -> None: ...
    def format(self, record: LogRecord) -> str: ...
    def emit(self, record: LogRecord) -> None: ...

class Formatter:
    converter: Callable[[float | None], struct_time]
    _fmt: str | None
    datefmt: str | None
    def __init__(self, fmt: str | None = ..., datefmt: str | None = ...) -> None: ...
    def format(self, record: LogRecord) -> str: ...
    def formatTime(self, record: LogRecord, datefmt: str | None = ...) -> str: ...
    def formatException(self, ei: _SysExcInfoType) -> str: ...

class Filter:
    def __init__(self, name: str = ...) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...

class LogRecord:
    args: _ArgsType
    asctime: str
    created: int
    exc_info: _SysExcInfoType | None
    exc_text: str | None
    filename: str
    funcName: str
    levelname: str
    levelno: int
    lineno: int
    module: str
    msecs: int
    message: str
    msg: str
    name: str
    pathname: str
    process: int
    processName: str
    relativeCreated: int
    thread: int
    threadName: str
    def __init__(
        self,
        name: str,
        level: int,
        pathname: str,
        lineno: int,
        msg: Any,
        args: _ArgsType,
        exc_info: _SysExcInfoType | None,
        func: str | None = ...,
    ) -> None: ...
    def getMessage(self) -> str: ...

_L = TypeVar("_L", Logger, LoggerAdapter[Logger], LoggerAdapter[Any])

class LoggerAdapter(Generic[_L]):
    logger: _L
    extra: Mapping[str, Any]
    def __init__(self, logger: _L, extra: Mapping[str, Any]) -> None: ...
    def process(self, msg: Any, kwargs: MutableMapping[str, Any]) -> Tuple[Any, MutableMapping[str, Any]]: ...
    def debug(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def info(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def warning(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def error(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def exception(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def critical(
        self, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def log(
        self, level: int, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
    ) -> None: ...
    def isEnabledFor(self, level: int) -> bool: ...

@overload
def getLogger() -> Logger: ...
@overload
def getLogger(name: Text | str) -> Logger: ...
def getLoggerClass() -> type: ...
def debug(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...
def info(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...
def warning(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...

warn = warning

def error(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...
def critical(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...
def exception(msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any) -> None: ...
def log(
    level: int, msg: Any, *args: Any, exc_info: _ExcInfoType = ..., extra: Dict[str, Any] | None = ..., **kwargs: Any
) -> None: ...

fatal = critical

def disable(level: int) -> None: ...
def addLevelName(level: int, levelName: str) -> None: ...
def getLevelName(level: int | str) -> Any: ...
def makeLogRecord(dict: Mapping[str, Any]) -> LogRecord: ...
@overload
def basicConfig() -> None: ...
@overload
def basicConfig(
    *,
    filename: str | None = ...,
    filemode: str = ...,
    format: str = ...,
    datefmt: str | None = ...,
    level: _Level | None = ...,
    stream: IO[str] = ...,
) -> None: ...
def shutdown(handlerList: Sequence[Any] = ...) -> None: ...  # handlerList is undocumented
def setLoggerClass(klass: type) -> None: ...
def captureWarnings(capture: bool) -> None: ...

_StreamT = TypeVar("_StreamT", bound=SupportsWrite[str])

class StreamHandler(Handler, Generic[_StreamT]):
    stream: _StreamT  # undocumented
    @overload
    def __init__(self: StreamHandler[IO[str]], stream: None = ...) -> None: ...
    @overload
    def __init__(self: StreamHandler[_StreamT], stream: _StreamT) -> None: ...

class FileHandler(StreamHandler):
    baseFilename: str  # undocumented
    mode: str  # undocumented
    encoding: str | None  # undocumented
    delay: bool  # undocumented
    def __init__(self, filename: StrPath, mode: str = ..., encoding: str | None = ..., delay: bool = ...) -> None: ...
    def _open(self) -> IO[Any]: ...

class NullHandler(Handler): ...

class PlaceHolder:
    def __init__(self, alogger: Logger) -> None: ...
    def append(self, alogger: Logger) -> None: ...

# Below aren't in module docs but still visible

class RootLogger(Logger):
    def __init__(self, level: int) -> None: ...

root: RootLogger

BASIC_FORMAT: str
