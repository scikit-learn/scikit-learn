import email.message
import io
import ssl
import sys
import types
from _typeshed import Self, WriteableBuffer
from socket import socket
from typing import IO, Any, BinaryIO, Callable, Iterable, Iterator, Mapping, Protocol, Type, TypeVar, Union, overload

_DataType = Union[bytes, IO[Any], Iterable[bytes], str]
_T = TypeVar("_T")

HTTP_PORT: int
HTTPS_PORT: int

CONTINUE: int
SWITCHING_PROTOCOLS: int
PROCESSING: int

OK: int
CREATED: int
ACCEPTED: int
NON_AUTHORITATIVE_INFORMATION: int
NO_CONTENT: int
RESET_CONTENT: int
PARTIAL_CONTENT: int
MULTI_STATUS: int
IM_USED: int

MULTIPLE_CHOICES: int
MOVED_PERMANENTLY: int
FOUND: int
SEE_OTHER: int
NOT_MODIFIED: int
USE_PROXY: int
TEMPORARY_REDIRECT: int

BAD_REQUEST: int
UNAUTHORIZED: int
PAYMENT_REQUIRED: int
FORBIDDEN: int
NOT_FOUND: int
METHOD_NOT_ALLOWED: int
NOT_ACCEPTABLE: int
PROXY_AUTHENTICATION_REQUIRED: int
REQUEST_TIMEOUT: int
CONFLICT: int
GONE: int
LENGTH_REQUIRED: int
PRECONDITION_FAILED: int
REQUEST_ENTITY_TOO_LARGE: int
REQUEST_URI_TOO_LONG: int
UNSUPPORTED_MEDIA_TYPE: int
REQUESTED_RANGE_NOT_SATISFIABLE: int
EXPECTATION_FAILED: int
UNPROCESSABLE_ENTITY: int
LOCKED: int
FAILED_DEPENDENCY: int
UPGRADE_REQUIRED: int
PRECONDITION_REQUIRED: int
TOO_MANY_REQUESTS: int
REQUEST_HEADER_FIELDS_TOO_LARGE: int

INTERNAL_SERVER_ERROR: int
NOT_IMPLEMENTED: int
BAD_GATEWAY: int
SERVICE_UNAVAILABLE: int
GATEWAY_TIMEOUT: int
HTTP_VERSION_NOT_SUPPORTED: int
INSUFFICIENT_STORAGE: int
NOT_EXTENDED: int
NETWORK_AUTHENTICATION_REQUIRED: int

responses: dict[int, str]

class HTTPMessage(email.message.Message):
    def getallmatchingheaders(self, name: str) -> list[str]: ...  # undocumented

def parse_headers(fp: io.BufferedIOBase, _class: Callable[[], email.message.Message] = ...) -> HTTPMessage: ...

class HTTPResponse(io.BufferedIOBase, BinaryIO):
    msg: HTTPMessage
    headers: HTTPMessage
    version: int
    debuglevel: int
    fp: io.BufferedReader
    closed: bool
    status: int
    reason: str
    chunked: bool
    chunk_left: int | None
    length: int | None
    will_close: bool
    def __init__(self, sock: socket, debuglevel: int = ..., method: str | None = ..., url: str | None = ...) -> None: ...
    def peek(self, n: int = ...) -> bytes: ...
    def read(self, amt: int | None = ...) -> bytes: ...
    def read1(self, n: int = ...) -> bytes: ...
    def readinto(self, b: WriteableBuffer) -> int: ...
    def readline(self, limit: int = ...) -> bytes: ...  # type: ignore
    @overload
    def getheader(self, name: str) -> str | None: ...
    @overload
    def getheader(self, name: str, default: _T) -> str | _T: ...
    def getheaders(self) -> list[tuple[str, str]]: ...
    def fileno(self) -> int: ...
    def isclosed(self) -> bool: ...
    def __iter__(self) -> Iterator[bytes]: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: types.TracebackType | None
    ) -> bool | None: ...
    def info(self) -> email.message.Message: ...
    def geturl(self) -> str: ...
    def getcode(self) -> int: ...
    def begin(self) -> None: ...

# This is an API stub only for the class below, not a class itself.
# urllib.request uses it for a parameter.
class _HTTPConnectionProtocol(Protocol):
    if sys.version_info >= (3, 7):
        def __call__(
            self,
            host: str,
            port: int | None = ...,
            timeout: float = ...,
            source_address: tuple[str, int] | None = ...,
            blocksize: int = ...,
        ) -> HTTPConnection: ...
    else:
        def __call__(
            self, host: str, port: int | None = ..., timeout: float = ..., source_address: tuple[str, int] | None = ...
        ) -> HTTPConnection: ...

class HTTPConnection:
    auto_open: int  # undocumented
    debuglevel: int
    default_port: int  # undocumented
    response_class: Type[HTTPResponse]  # undocumented
    timeout: float | None
    host: str
    port: int
    sock: Any
    if sys.version_info >= (3, 7):
        def __init__(
            self,
            host: str,
            port: int | None = ...,
            timeout: float | None = ...,
            source_address: tuple[str, int] | None = ...,
            blocksize: int = ...,
        ) -> None: ...
    else:
        def __init__(
            self, host: str, port: int | None = ..., timeout: float | None = ..., source_address: tuple[str, int] | None = ...
        ) -> None: ...
    def request(
        self, method: str, url: str, body: _DataType | None = ..., headers: Mapping[str, str] = ..., *, encode_chunked: bool = ...
    ) -> None: ...
    def getresponse(self) -> HTTPResponse: ...
    def set_debuglevel(self, level: int) -> None: ...
    def set_tunnel(self, host: str, port: int | None = ..., headers: Mapping[str, str] | None = ...) -> None: ...
    def connect(self) -> None: ...
    def close(self) -> None: ...
    def putrequest(self, method: str, url: str, skip_host: bool = ..., skip_accept_encoding: bool = ...) -> None: ...
    def putheader(self, header: str, *argument: str) -> None: ...
    def endheaders(self, message_body: _DataType | None = ..., *, encode_chunked: bool = ...) -> None: ...
    def send(self, data: _DataType) -> None: ...

class HTTPSConnection(HTTPConnection):
    if sys.version_info >= (3, 7):
        def __init__(
            self,
            host: str,
            port: int | None = ...,
            key_file: str | None = ...,
            cert_file: str | None = ...,
            timeout: float | None = ...,
            source_address: tuple[str, int] | None = ...,
            *,
            context: ssl.SSLContext | None = ...,
            check_hostname: bool | None = ...,
            blocksize: int = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            host: str,
            port: int | None = ...,
            key_file: str | None = ...,
            cert_file: str | None = ...,
            timeout: float | None = ...,
            source_address: tuple[str, int] | None = ...,
            *,
            context: ssl.SSLContext | None = ...,
            check_hostname: bool | None = ...,
        ) -> None: ...

class HTTPException(Exception): ...

error = HTTPException

class NotConnected(HTTPException): ...
class InvalidURL(HTTPException): ...

class UnknownProtocol(HTTPException):
    def __init__(self, version: str) -> None: ...

class UnknownTransferEncoding(HTTPException): ...
class UnimplementedFileMode(HTTPException): ...

class IncompleteRead(HTTPException):
    def __init__(self, partial: bytes, expected: int | None = ...) -> None: ...
    partial: bytes
    expected: int | None

class ImproperConnectionState(HTTPException): ...
class CannotSendRequest(ImproperConnectionState): ...
class CannotSendHeader(ImproperConnectionState): ...
class ResponseNotReady(ImproperConnectionState): ...

class BadStatusLine(HTTPException):
    def __init__(self, line: str) -> None: ...

class LineTooLong(HTTPException):
    def __init__(self, line_type: str) -> None: ...

class RemoteDisconnected(ConnectionResetError, BadStatusLine): ...
