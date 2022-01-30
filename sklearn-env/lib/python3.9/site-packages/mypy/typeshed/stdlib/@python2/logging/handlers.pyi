from _typeshed import StrPath
from logging import FileHandler, Handler, LogRecord
from socket import SocketKind, SocketType
from typing import Any, ClassVar, Dict, List, Tuple

DEFAULT_TCP_LOGGING_PORT: int
DEFAULT_UDP_LOGGING_PORT: int
DEFAULT_HTTP_LOGGING_PORT: int
DEFAULT_SOAP_LOGGING_PORT: int
SYSLOG_UDP_PORT: int
SYSLOG_TCP_PORT: int

class WatchedFileHandler(FileHandler):
    dev: int
    ino: int
    def __init__(self, filename: StrPath, mode: str = ..., encoding: str | None = ..., delay: bool = ...) -> None: ...
    def _statstream(self) -> None: ...

class RotatingFileHandler(Handler):
    def __init__(
        self,
        filename: str,
        mode: str = ...,
        maxBytes: int = ...,
        backupCount: int = ...,
        encoding: str | None = ...,
        delay: bool = ...,
    ) -> None: ...
    def doRollover(self) -> None: ...

class TimedRotatingFileHandler(Handler):
    def __init__(
        self,
        filename: str,
        when: str = ...,
        interval: int = ...,
        backupCount: int = ...,
        encoding: str | None = ...,
        delay: bool = ...,
        utc: bool = ...,
    ) -> None: ...
    def doRollover(self) -> None: ...

class SocketHandler(Handler):
    retryStart: float
    retryFactor: float
    retryMax: float
    def __init__(self, host: str, port: int) -> None: ...
    def makeSocket(self, timeout: float = ...) -> SocketType: ...  # timeout is undocumented
    def makePickle(self, record: LogRecord) -> bytes: ...
    def send(self, s: bytes) -> None: ...
    def createSocket(self) -> None: ...

class DatagramHandler(SocketHandler):
    def makeSocket(self) -> SocketType: ...  # type: ignore

class SysLogHandler(Handler):
    LOG_EMERG: int
    LOG_ALERT: int
    LOG_CRIT: int
    LOG_ERR: int
    LOG_WARNING: int
    LOG_NOTICE: int
    LOG_INFO: int
    LOG_DEBUG: int

    LOG_KERN: int
    LOG_USER: int
    LOG_MAIL: int
    LOG_DAEMON: int
    LOG_AUTH: int
    LOG_SYSLOG: int
    LOG_LPR: int
    LOG_NEWS: int
    LOG_UUCP: int
    LOG_CRON: int
    LOG_AUTHPRIV: int
    LOG_FTP: int
    LOG_LOCAL0: int
    LOG_LOCAL1: int
    LOG_LOCAL2: int
    LOG_LOCAL3: int
    LOG_LOCAL4: int
    LOG_LOCAL5: int
    LOG_LOCAL6: int
    LOG_LOCAL7: int
    address: Tuple[str, int] | str  # undocumented
    unixsocket: bool  # undocumented
    socktype: SocketKind  # undocumented
    facility: int  # undocumented
    priority_names: ClassVar[Dict[str, int]]  # undocumented
    facility_names: ClassVar[Dict[str, int]]  # undocumented
    priority_map: ClassVar[Dict[str, str]]  # undocumented
    def __init__(self, address: Tuple[str, int] | str = ..., facility: int = ..., socktype: SocketKind | None = ...) -> None: ...
    def encodePriority(self, facility: int | str, priority: int | str) -> int: ...
    def mapPriority(self, levelName: str) -> str: ...

class NTEventLogHandler(Handler):
    def __init__(self, appname: str, dllname: str | None = ..., logtype: str = ...) -> None: ...
    def getEventCategory(self, record: LogRecord) -> int: ...
    # TODO correct return value?
    def getEventType(self, record: LogRecord) -> int: ...
    def getMessageID(self, record: LogRecord) -> int: ...

class SMTPHandler(Handler):
    # TODO `secure` can also be an empty tuple
    def __init__(
        self,
        mailhost: str | Tuple[str, int],
        fromaddr: str,
        toaddrs: List[str],
        subject: str,
        credentials: Tuple[str, str] | None = ...,
        secure: Tuple[str] | Tuple[str, str] | None = ...,
    ) -> None: ...
    def getSubject(self, record: LogRecord) -> str: ...

class BufferingHandler(Handler):
    buffer: List[LogRecord]
    def __init__(self, capacity: int) -> None: ...
    def shouldFlush(self, record: LogRecord) -> bool: ...

class MemoryHandler(BufferingHandler):
    def __init__(self, capacity: int, flushLevel: int = ..., target: Handler | None = ...) -> None: ...
    def setTarget(self, target: Handler) -> None: ...

class HTTPHandler(Handler):
    def __init__(self, host: str, url: str, method: str = ...) -> None: ...
    def mapLogRecord(self, record: LogRecord) -> Dict[str, Any]: ...
