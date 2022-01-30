import sys
import threading
from _typeshed import StrPath, SupportsWrite
from collections.abc import Callable, Iterable, Mapping, MutableMapping, Sequence
from io import TextIOWrapper
from string import Template
from time import struct_time
from types import FrameType, TracebackType
from typing import Any, ClassVar, Generic, Optional, Pattern, TextIO, Tuple, Type, TypeVar, Union, overload
from typing_extensions import Literal

_SysExcInfoType = Union[Tuple[Type[BaseException], BaseException, Optional[TracebackType]], Tuple[None, None, None]]
_ExcInfoType = Union[None, bool, _SysExcInfoType, BaseException]
_ArgsType = Union[Tuple[object, ...], Mapping[str, object]]
_FilterType = Union[Filter, Callable[[LogRecord], int]]
_Level = Union[int, str]
_FormatStyle = Literal["%", "{", "$"]

raiseExceptions: bool
logThreads: bool
logMultiprocessing: bool
logProcesses: bool
_srcfile: str | None

def currentframe() -> FrameType: ...

_levelToName: dict[int, str]
_nameToLevel: dict[str, int]

class Filterer(object):
    filters: list[Filter]
    def __init__(self) -> None: ...
    def addFilter(self, filter: _FilterType) -> None: ...
    def removeFilter(self, filter: _FilterType) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...

class Manager(object):  # undocumented
    root: RootLogger
    disable: int
    emittedNoHandlerWarning: bool
    loggerDict: dict[str, Logger | PlaceHolder]
    loggerClass: Type[Logger] | None
    logRecordFactory: Callable[..., LogRecord] | None
    def __init__(self, rootnode: RootLogger) -> None: ...
    def getLogger(self, name: str) -> Logger: ...
    def setLoggerClass(self, klass: Type[Logger]) -> None: ...
    def setLogRecordFactory(self, factory: Callable[..., LogRecord]) -> None: ...

class Logger(Filterer):
    name: str  # undocumented
    level: int  # undocumented
    parent: Logger | None  # undocumented
    propagate: bool
    handlers: list[Handler]  # undocumented
    disabled: bool  # undocumented
    root: ClassVar[RootLogger]  # undocumented
    manager: Manager  # undocumented
    def __init__(self, name: str, level: _Level = ...) -> None: ...
    def setLevel(self, level: _Level) -> None: ...
    def isEnabledFor(self, level: int) -> bool: ...
    def getEffectiveLevel(self) -> int: ...
    def getChild(self, suffix: str) -> Logger: ...
    if sys.version_info >= (3, 8):
        def debug(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def info(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def warning(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def warn(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def error(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def exception(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def critical(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def log(
            self,
            level: int,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def _log(
            self,
            level: int,
            msg: object,
            args: _ArgsType,
            exc_info: _ExcInfoType | None = ...,
            extra: Mapping[str, object] | None = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
        ) -> None: ...  # undocumented
    else:
        def debug(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def info(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def warning(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def warn(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def error(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def critical(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def log(
            self,
            level: int,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def exception(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
        ) -> None: ...
        def _log(
            self,
            level: int,
            msg: object,
            args: _ArgsType,
            exc_info: _ExcInfoType | None = ...,
            extra: Mapping[str, object] | None = ...,
            stack_info: bool = ...,
        ) -> None: ...  # undocumented
    fatal = critical
    def filter(self, record: LogRecord) -> bool: ...
    def addHandler(self, hdlr: Handler) -> None: ...
    def removeHandler(self, hdlr: Handler) -> None: ...
    if sys.version_info >= (3, 8):
        def findCaller(self, stack_info: bool = ..., stacklevel: int = ...) -> tuple[str, int, str, str | None]: ...
    else:
        def findCaller(self, stack_info: bool = ...) -> tuple[str, int, str, str | None]: ...
    def handle(self, record: LogRecord) -> None: ...
    def makeRecord(
        self,
        name: str,
        level: int,
        fn: str,
        lno: int,
        msg: object,
        args: _ArgsType,
        exc_info: _SysExcInfoType | None,
        func: str | None = ...,
        extra: Mapping[str, object] | None = ...,
        sinfo: str | None = ...,
    ) -> LogRecord: ...
    def hasHandlers(self) -> bool: ...
    def callHandlers(self, record: LogRecord) -> None: ...  # undocumented

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
    def get_name(self) -> str: ...  # undocumented
    def set_name(self, name: str) -> None: ...  # undocumented
    def createLock(self) -> None: ...
    def acquire(self) -> None: ...
    def release(self) -> None: ...
    def setLevel(self, level: _Level) -> None: ...
    def setFormatter(self, fmt: Formatter | None) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...
    def flush(self) -> None: ...
    def close(self) -> None: ...
    def handle(self, record: LogRecord) -> bool: ...
    def handleError(self, record: LogRecord) -> None: ...
    def format(self, record: LogRecord) -> str: ...
    def emit(self, record: LogRecord) -> None: ...

class Formatter:
    converter: Callable[[float | None], struct_time]
    _fmt: str | None  # undocumented
    datefmt: str | None  # undocumented
    _style: PercentStyle  # undocumented
    default_time_format: str
    if sys.version_info >= (3, 9):
        default_msec_format: str | None
    else:
        default_msec_format: str

    if sys.version_info >= (3, 8):
        def __init__(
            self, fmt: str | None = ..., datefmt: str | None = ..., style: _FormatStyle = ..., validate: bool = ...
        ) -> None: ...
    else:
        def __init__(self, fmt: str | None = ..., datefmt: str | None = ..., style: _FormatStyle = ...) -> None: ...
    def format(self, record: LogRecord) -> str: ...
    def formatTime(self, record: LogRecord, datefmt: str | None = ...) -> str: ...
    def formatException(self, ei: _SysExcInfoType) -> str: ...
    def formatMessage(self, record: LogRecord) -> str: ...  # undocumented
    def formatStack(self, stack_info: str) -> str: ...
    def usesTime(self) -> bool: ...  # undocumented

class BufferingFormatter:
    linefmt: Formatter
    def __init__(self, linefmt: Formatter | None = ...) -> None: ...
    def formatHeader(self, records: Sequence[LogRecord]) -> str: ...
    def formatFooter(self, records: Sequence[LogRecord]) -> str: ...
    def format(self, records: Sequence[LogRecord]) -> str: ...

class Filter:
    name: str  # undocumented
    nlen: int  # undocumented
    def __init__(self, name: str = ...) -> None: ...
    def filter(self, record: LogRecord) -> bool: ...

class LogRecord:
    # args can be set to None by logging.handlers.QueueHandler
    # (see https://bugs.python.org/issue44473)
    args: _ArgsType | None
    asctime: str
    created: float
    exc_info: _SysExcInfoType | None
    exc_text: str | None
    filename: str
    funcName: str
    levelname: str
    levelno: int
    lineno: int
    module: str
    msecs: float
    message: str
    msg: str
    name: str
    pathname: str
    process: int | None
    processName: str | None
    relativeCreated: float
    stack_info: str | None
    thread: int | None
    threadName: str | None
    def __init__(
        self,
        name: str,
        level: int,
        pathname: str,
        lineno: int,
        msg: object,
        args: _ArgsType | None,
        exc_info: _SysExcInfoType | None,
        func: str | None = ...,
        sinfo: str | None = ...,
    ) -> None: ...
    def getMessage(self) -> str: ...

_L = TypeVar("_L", Logger, LoggerAdapter[Logger], LoggerAdapter[Any])

class LoggerAdapter(Generic[_L]):
    logger: _L
    manager: Manager  # undocumented
    if sys.version_info >= (3, 10):
        extra: Mapping[str, object] | None
        def __init__(self, logger: _L, extra: Mapping[str, object] | None) -> None: ...
    else:
        extra: Mapping[str, object]
        def __init__(self, logger: _L, extra: Mapping[str, object]) -> None: ...
    def process(self, msg: Any, kwargs: MutableMapping[str, Any]) -> tuple[Any, MutableMapping[str, Any]]: ...
    if sys.version_info >= (3, 8):
        def debug(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def info(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def warning(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def warn(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def error(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def exception(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def critical(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def log(
            self,
            level: int,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            stacklevel: int = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
    else:
        def debug(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def info(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def warning(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def warn(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def error(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def exception(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def critical(
            self,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
        def log(
            self,
            level: int,
            msg: object,
            *args: object,
            exc_info: _ExcInfoType = ...,
            stack_info: bool = ...,
            extra: Mapping[str, object] | None = ...,
            **kwargs: object,
        ) -> None: ...
    def isEnabledFor(self, level: int) -> bool: ...
    def getEffectiveLevel(self) -> int: ...
    def setLevel(self, level: _Level) -> None: ...
    def hasHandlers(self) -> bool: ...
    def _log(
        self,
        level: int,
        msg: object,
        args: _ArgsType,
        exc_info: _ExcInfoType | None = ...,
        extra: Mapping[str, object] | None = ...,
        stack_info: bool = ...,
    ) -> None: ...  # undocumented
    @property
    def name(self) -> str: ...  # undocumented

def getLogger(name: str | None = ...) -> Logger: ...
def getLoggerClass() -> Type[Logger]: ...
def getLogRecordFactory() -> Callable[..., LogRecord]: ...

if sys.version_info >= (3, 8):
    def debug(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def info(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def warning(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def warn(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def error(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def critical(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def exception(
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...
    def log(
        level: int,
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        stacklevel: int = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...

else:
    def debug(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def info(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def warning(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def warn(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def error(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def critical(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def exception(
        msg: object, *args: object, exc_info: _ExcInfoType = ..., stack_info: bool = ..., extra: Mapping[str, object] | None = ...
    ) -> None: ...
    def log(
        level: int,
        msg: object,
        *args: object,
        exc_info: _ExcInfoType = ...,
        stack_info: bool = ...,
        extra: Mapping[str, object] | None = ...,
    ) -> None: ...

fatal = critical

if sys.version_info >= (3, 7):
    def disable(level: int = ...) -> None: ...

else:
    def disable(level: int) -> None: ...

def addLevelName(level: int, levelName: str) -> None: ...
def getLevelName(level: _Level) -> Any: ...
def makeLogRecord(dict: Mapping[str, object]) -> LogRecord: ...

if sys.version_info >= (3, 9):
    def basicConfig(
        *,
        filename: StrPath | None = ...,
        filemode: str = ...,
        format: str = ...,
        datefmt: str | None = ...,
        style: _FormatStyle = ...,
        level: _Level | None = ...,
        stream: SupportsWrite[str] | None = ...,
        handlers: Iterable[Handler] | None = ...,
        force: bool | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
    ) -> None: ...

elif sys.version_info >= (3, 8):
    def basicConfig(
        *,
        filename: StrPath | None = ...,
        filemode: str = ...,
        format: str = ...,
        datefmt: str | None = ...,
        style: _FormatStyle = ...,
        level: _Level | None = ...,
        stream: SupportsWrite[str] | None = ...,
        handlers: Iterable[Handler] | None = ...,
        force: bool = ...,
    ) -> None: ...

else:
    def basicConfig(
        *,
        filename: StrPath | None = ...,
        filemode: str = ...,
        format: str = ...,
        datefmt: str | None = ...,
        style: _FormatStyle = ...,
        level: _Level | None = ...,
        stream: SupportsWrite[str] | None = ...,
        handlers: Iterable[Handler] | None = ...,
    ) -> None: ...

def shutdown(handlerList: Sequence[Any] = ...) -> None: ...  # handlerList is undocumented
def setLoggerClass(klass: Type[Logger]) -> None: ...
def captureWarnings(capture: bool) -> None: ...
def setLogRecordFactory(factory: Callable[..., LogRecord]) -> None: ...

lastResort: StreamHandler[Any] | None

_StreamT = TypeVar("_StreamT", bound=SupportsWrite[str])

class StreamHandler(Handler, Generic[_StreamT]):
    stream: _StreamT  # undocumented
    terminator: str
    @overload
    def __init__(self: StreamHandler[TextIO], stream: None = ...) -> None: ...
    @overload
    def __init__(self: StreamHandler[_StreamT], stream: _StreamT) -> None: ...
    if sys.version_info >= (3, 7):
        def setStream(self, stream: _StreamT) -> _StreamT | None: ...

class FileHandler(StreamHandler[TextIOWrapper]):
    baseFilename: str  # undocumented
    mode: str  # undocumented
    encoding: str | None  # undocumented
    delay: bool  # undocumented
    if sys.version_info >= (3, 9):
        errors: str | None  # undocumented
        def __init__(
            self, filename: StrPath, mode: str = ..., encoding: str | None = ..., delay: bool = ..., errors: str | None = ...
        ) -> None: ...
    else:
        def __init__(self, filename: StrPath, mode: str = ..., encoding: str | None = ..., delay: bool = ...) -> None: ...
    def _open(self) -> TextIOWrapper: ...  # undocumented

class NullHandler(Handler): ...

class PlaceHolder:  # undocumented
    loggerMap: dict[Logger, None]
    def __init__(self, alogger: Logger) -> None: ...
    def append(self, alogger: Logger) -> None: ...

# Below aren't in module docs but still visible

class RootLogger(Logger):
    def __init__(self, level: int) -> None: ...

root: RootLogger

class PercentStyle(object):  # undocumented
    default_format: str
    asctime_format: str
    asctime_search: str
    if sys.version_info >= (3, 8):
        validation_pattern: Pattern[str]
    _fmt: str
    def __init__(self, fmt: str) -> None: ...
    def usesTime(self) -> bool: ...
    if sys.version_info >= (3, 8):
        def validate(self) -> None: ...
    def format(self, record: Any) -> str: ...

class StrFormatStyle(PercentStyle):  # undocumented
    fmt_spec = Any
    field_spec = Any

class StringTemplateStyle(PercentStyle):  # undocumented
    _tpl: Template

_STYLES: dict[str, tuple[PercentStyle, str]]

BASIC_FORMAT: str
