import ssl
import sys
from _typeshed import StrOrBytesPath
from email.message import Message
from http.client import HTTPMessage, HTTPResponse, _HTTPConnectionProtocol
from http.cookiejar import CookieJar
from typing import IO, Any, Callable, ClassVar, Mapping, NoReturn, Pattern, Sequence, Tuple, TypeVar, overload
from urllib.error import HTTPError
from urllib.response import addclosehook, addinfourl

_T = TypeVar("_T")
_UrlopenRet = Any

def urlopen(
    url: str | Request,
    data: bytes | None = ...,
    timeout: float | None = ...,
    *,
    cafile: str | None = ...,
    capath: str | None = ...,
    cadefault: bool = ...,
    context: ssl.SSLContext | None = ...,
) -> _UrlopenRet: ...
def install_opener(opener: OpenerDirector) -> None: ...
def build_opener(*handlers: BaseHandler | Callable[[], BaseHandler]) -> OpenerDirector: ...

if sys.platform == "win32":
    from nturl2path import pathname2url as pathname2url, url2pathname as url2pathname
else:
    def url2pathname(pathname: str) -> str: ...
    def pathname2url(pathname: str) -> str: ...

def getproxies() -> dict[str, str]: ...
def parse_http_list(s: str) -> list[str]: ...
def parse_keqv_list(l: list[str]) -> dict[str, str]: ...

if sys.platform == "win32" or sys.platform == "darwin":
    def proxy_bypass(host: str) -> Any: ...  # undocumented

else:
    def proxy_bypass(host: str, proxies: Mapping[str, str] | None = ...) -> Any: ...  # undocumented

class Request:
    @property
    def full_url(self) -> str: ...
    @full_url.setter
    def full_url(self, value: str) -> None: ...
    @full_url.deleter
    def full_url(self) -> None: ...
    type: str
    host: str
    origin_req_host: str
    selector: str
    data: bytes | None
    headers: dict[str, str]
    unredirected_hdrs: dict[str, str]
    unverifiable: bool
    method: str | None
    timeout: float | None  # Undocumented, only set after __init__() by OpenerDirector.open()
    def __init__(
        self,
        url: str,
        data: bytes | None = ...,
        headers: dict[str, str] = ...,
        origin_req_host: str | None = ...,
        unverifiable: bool = ...,
        method: str | None = ...,
    ) -> None: ...
    def get_method(self) -> str: ...
    def add_header(self, key: str, val: str) -> None: ...
    def add_unredirected_header(self, key: str, val: str) -> None: ...
    def has_header(self, header_name: str) -> bool: ...
    def remove_header(self, header_name: str) -> None: ...
    def get_full_url(self) -> str: ...
    def set_proxy(self, host: str, type: str) -> None: ...
    @overload
    def get_header(self, header_name: str) -> str | None: ...
    @overload
    def get_header(self, header_name: str, default: _T) -> str | _T: ...
    def header_items(self) -> list[tuple[str, str]]: ...
    def has_proxy(self) -> bool: ...

class OpenerDirector:
    addheaders: list[tuple[str, str]]
    def add_handler(self, handler: BaseHandler) -> None: ...
    def open(self, fullurl: str | Request, data: bytes | None = ..., timeout: float | None = ...) -> _UrlopenRet: ...
    def error(self, proto: str, *args: Any) -> _UrlopenRet: ...
    def close(self) -> None: ...

class BaseHandler:
    handler_order: ClassVar[int]
    parent: OpenerDirector
    def add_parent(self, parent: OpenerDirector) -> None: ...
    def close(self) -> None: ...

class HTTPDefaultErrorHandler(BaseHandler):
    def http_error_default(
        self, req: Request, fp: IO[bytes], code: int, msg: str, hdrs: HTTPMessage
    ) -> HTTPError: ...  # undocumented

class HTTPRedirectHandler(BaseHandler):
    max_redirections: ClassVar[int]  # undocumented
    max_repeats: ClassVar[int]  # undocumented
    inf_msg: ClassVar[str]  # undocumented
    def redirect_request(
        self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage, newurl: str
    ) -> Request | None: ...
    def http_error_301(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...
    def http_error_302(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...
    def http_error_303(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...
    def http_error_307(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...

class HTTPCookieProcessor(BaseHandler):
    cookiejar: CookieJar
    def __init__(self, cookiejar: CookieJar | None = ...) -> None: ...
    def http_request(self, request: Request) -> Request: ...  # undocumented
    def http_response(self, request: Request, response: HTTPResponse) -> HTTPResponse: ...  # undocumented
    def https_request(self, request: Request) -> Request: ...  # undocumented
    def https_response(self, request: Request, response: HTTPResponse) -> HTTPResponse: ...  # undocumented

class ProxyHandler(BaseHandler):
    def __init__(self, proxies: dict[str, str] | None = ...) -> None: ...
    def proxy_open(self, req: Request, proxy: str, type: str) -> _UrlopenRet | None: ...  # undocumented
    # TODO add a method for every (common) proxy protocol

class HTTPPasswordMgr:
    def add_password(self, realm: str, uri: str | Sequence[str], user: str, passwd: str) -> None: ...
    def find_user_password(self, realm: str, authuri: str) -> tuple[str | None, str | None]: ...
    def is_suburi(self, base: str, test: str) -> bool: ...  # undocumented
    def reduce_uri(self, uri: str, default_port: bool = ...) -> str: ...  # undocumented

class HTTPPasswordMgrWithDefaultRealm(HTTPPasswordMgr):
    def add_password(self, realm: str | None, uri: str | Sequence[str], user: str, passwd: str) -> None: ...
    def find_user_password(self, realm: str | None, authuri: str) -> tuple[str | None, str | None]: ...

class HTTPPasswordMgrWithPriorAuth(HTTPPasswordMgrWithDefaultRealm):
    def add_password(
        self, realm: str | None, uri: str | Sequence[str], user: str, passwd: str, is_authenticated: bool = ...
    ) -> None: ...
    def update_authenticated(self, uri: str | Sequence[str], is_authenticated: bool = ...) -> None: ...
    def is_authenticated(self, authuri: str) -> bool: ...

class AbstractBasicAuthHandler:
    rx: ClassVar[Pattern[str]]  # undocumented
    passwd: HTTPPasswordMgr
    add_password: Callable[[str, str | Sequence[str], str, str], None]
    def __init__(self, password_mgr: HTTPPasswordMgr | None = ...) -> None: ...
    def http_error_auth_reqed(self, authreq: str, host: str, req: Request, headers: HTTPMessage) -> None: ...
    def http_request(self, req: Request) -> Request: ...  # undocumented
    def http_response(self, req: Request, response: HTTPResponse) -> HTTPResponse: ...  # undocumented
    def https_request(self, req: Request) -> Request: ...  # undocumented
    def https_response(self, req: Request, response: HTTPResponse) -> HTTPResponse: ...  # undocumented
    def retry_http_basic_auth(self, host: str, req: Request, realm: str) -> _UrlopenRet | None: ...  # undocumented

class HTTPBasicAuthHandler(AbstractBasicAuthHandler, BaseHandler):
    auth_header: ClassVar[str]  # undocumented
    def http_error_401(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...

class ProxyBasicAuthHandler(AbstractBasicAuthHandler, BaseHandler):
    auth_header: ClassVar[str]
    def http_error_407(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...

class AbstractDigestAuthHandler:
    def __init__(self, passwd: HTTPPasswordMgr | None = ...) -> None: ...
    def reset_retry_count(self) -> None: ...
    def http_error_auth_reqed(self, auth_header: str, host: str, req: Request, headers: HTTPMessage) -> None: ...
    def retry_http_digest_auth(self, req: Request, auth: str) -> _UrlopenRet | None: ...
    def get_cnonce(self, nonce: str) -> str: ...
    def get_authorization(self, req: Request, chal: Mapping[str, str]) -> str: ...
    def get_algorithm_impls(self, algorithm: str) -> tuple[Callable[[str], str], Callable[[str, str], str]]: ...
    def get_entity_digest(self, data: bytes | None, chal: Mapping[str, str]) -> str | None: ...

class HTTPDigestAuthHandler(BaseHandler, AbstractDigestAuthHandler):
    auth_header: ClassVar[str]  # undocumented
    def http_error_401(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...

class ProxyDigestAuthHandler(BaseHandler, AbstractDigestAuthHandler):
    auth_header: ClassVar[str]  # undocumented
    def http_error_407(self, req: Request, fp: IO[bytes], code: int, msg: str, headers: HTTPMessage) -> _UrlopenRet | None: ...

class AbstractHTTPHandler(BaseHandler):  # undocumented
    def __init__(self, debuglevel: int = ...) -> None: ...
    def set_http_debuglevel(self, level: int) -> None: ...
    def do_request_(self, request: Request) -> Request: ...
    def do_open(self, http_class: _HTTPConnectionProtocol, req: Request, **http_conn_args: Any) -> HTTPResponse: ...

class HTTPHandler(AbstractHTTPHandler):
    def http_open(self, req: Request) -> HTTPResponse: ...
    def http_request(self, request: Request) -> Request: ...  # undocumented

class HTTPSHandler(AbstractHTTPHandler):
    def __init__(
        self, debuglevel: int = ..., context: ssl.SSLContext | None = ..., check_hostname: bool | None = ...
    ) -> None: ...
    def https_open(self, req: Request) -> HTTPResponse: ...
    def https_request(self, request: Request) -> Request: ...  # undocumented

class FileHandler(BaseHandler):
    names: ClassVar[Tuple[str, ...] | None]  # undocumented
    def file_open(self, req: Request) -> addinfourl: ...
    def get_names(self) -> Tuple[str, ...]: ...  # undocumented
    def open_local_file(self, req: Request) -> addinfourl: ...  # undocumented

class DataHandler(BaseHandler):
    def data_open(self, req: Request) -> addinfourl: ...

class ftpwrapper:  # undocumented
    def __init__(
        self, user: str, passwd: str, host: str, port: int, dirs: str, timeout: float | None = ..., persistent: bool = ...
    ) -> None: ...
    def close(self) -> None: ...
    def endtransfer(self) -> None: ...
    def file_close(self) -> None: ...
    def init(self) -> None: ...
    def real_close(self) -> None: ...
    def retrfile(self, file: str, type: str) -> tuple[addclosehook, int]: ...

class FTPHandler(BaseHandler):
    def ftp_open(self, req: Request) -> addinfourl: ...
    def connect_ftp(
        self, user: str, passwd: str, host: str, port: int, dirs: str, timeout: float
    ) -> ftpwrapper: ...  # undocumented

class CacheFTPHandler(FTPHandler):
    def setTimeout(self, t: float) -> None: ...
    def setMaxConns(self, m: int) -> None: ...
    def check_cache(self) -> None: ...  # undocumented
    def clear_cache(self) -> None: ...  # undocumented
    def connect_ftp(
        self, user: str, passwd: str, host: str, port: int, dirs: str, timeout: float
    ) -> ftpwrapper: ...  # undocumented

class UnknownHandler(BaseHandler):
    def unknown_open(self, req: Request) -> NoReturn: ...

class HTTPErrorProcessor(BaseHandler):
    def http_response(self, request: Request, response: HTTPResponse) -> _UrlopenRet: ...
    def https_response(self, request: Request, response: HTTPResponse) -> _UrlopenRet: ...

def urlretrieve(
    url: str,
    filename: StrOrBytesPath | None = ...,
    reporthook: Callable[[int, int, int], None] | None = ...,
    data: bytes | None = ...,
) -> tuple[str, HTTPMessage]: ...
def urlcleanup() -> None: ...

class URLopener:
    version: ClassVar[str]
    def __init__(self, proxies: dict[str, str] | None = ..., **x509: str) -> None: ...
    def open(self, fullurl: str, data: bytes | None = ...) -> _UrlopenRet: ...
    def open_unknown(self, fullurl: str, data: bytes | None = ...) -> _UrlopenRet: ...
    def retrieve(
        self,
        url: str,
        filename: str | None = ...,
        reporthook: Callable[[int, int, int], None] | None = ...,
        data: bytes | None = ...,
    ) -> tuple[str, Message | None]: ...
    def addheader(self, *args: tuple[str, str]) -> None: ...  # undocumented
    def cleanup(self) -> None: ...  # undocumented
    def close(self) -> None: ...  # undocumented
    def http_error(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None = ...
    ) -> _UrlopenRet: ...  # undocumented
    def http_error_default(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage
    ) -> _UrlopenRet: ...  # undocumented
    def open_data(self, url: str, data: bytes | None = ...) -> addinfourl: ...  # undocumented
    def open_file(self, url: str) -> addinfourl: ...  # undocumented
    def open_ftp(self, url: str) -> addinfourl: ...  # undocumented
    def open_http(self, url: str, data: bytes | None = ...) -> _UrlopenRet: ...  # undocumented
    def open_https(self, url: str, data: bytes | None = ...) -> _UrlopenRet: ...  # undocumented
    def open_local_file(self, url: str) -> addinfourl: ...  # undocumented
    def open_unknown_proxy(self, proxy: str, fullurl: str, data: bytes | None = ...) -> None: ...  # undocumented

class FancyURLopener(URLopener):
    def prompt_user_passwd(self, host: str, realm: str) -> tuple[str, str]: ...
    def get_user_passwd(self, host: str, realm: str, clear_cache: int = ...) -> tuple[str, str]: ...  # undocumented
    def http_error_301(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None = ...
    ) -> _UrlopenRet | addinfourl | None: ...  # undocumented
    def http_error_302(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None = ...
    ) -> _UrlopenRet | addinfourl | None: ...  # undocumented
    def http_error_303(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None = ...
    ) -> _UrlopenRet | addinfourl | None: ...  # undocumented
    def http_error_307(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None = ...
    ) -> _UrlopenRet | addinfourl | None: ...  # undocumented
    def http_error_401(
        self,
        url: str,
        fp: IO[bytes],
        errcode: int,
        errmsg: str,
        headers: HTTPMessage,
        data: bytes | None = ...,
        retry: bool = ...,
    ) -> _UrlopenRet | None: ...  # undocumented
    def http_error_407(
        self,
        url: str,
        fp: IO[bytes],
        errcode: int,
        errmsg: str,
        headers: HTTPMessage,
        data: bytes | None = ...,
        retry: bool = ...,
    ) -> _UrlopenRet | None: ...  # undocumented
    def http_error_default(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage
    ) -> addinfourl: ...  # undocumented
    def redirect_internal(
        self, url: str, fp: IO[bytes], errcode: int, errmsg: str, headers: HTTPMessage, data: bytes | None
    ) -> _UrlopenRet | None: ...  # undocumented
    def retry_http_basic_auth(self, url: str, realm: str, data: bytes | None = ...) -> _UrlopenRet | None: ...  # undocumented
    def retry_https_basic_auth(self, url: str, realm: str, data: bytes | None = ...) -> _UrlopenRet | None: ...  # undocumented
    def retry_proxy_http_basic_auth(
        self, url: str, realm: str, data: bytes | None = ...
    ) -> _UrlopenRet | None: ...  # undocumented
    def retry_proxy_https_basic_auth(
        self, url: str, realm: str, data: bytes | None = ...
    ) -> _UrlopenRet | None: ...  # undocumented
