import enum
import socket
import sys
from _typeshed import ReadableBuffer, Self, StrOrBytesPath, WriteableBuffer
from typing import Any, Callable, ClassVar, Dict, Iterable, List, NamedTuple, Optional, Set, Tuple, Type, Union, overload
from typing_extensions import Literal, TypedDict

_PCTRTT = Tuple[Tuple[str, str], ...]
_PCTRTTT = Tuple[_PCTRTT, ...]
_PeerCertRetDictType = Dict[str, Union[str, _PCTRTTT, _PCTRTT]]
_PeerCertRetType = Union[_PeerCertRetDictType, bytes, None]
_EnumRetType = List[Tuple[bytes, str, Union[Set[str], bool]]]
_PasswordType = Union[Callable[[], Union[str, bytes]], str, bytes]

_SrvnmeCbType = Callable[[Union[SSLSocket, SSLObject], Optional[str], SSLSocket], Optional[int]]

class _Cipher(TypedDict):
    aead: bool
    alg_bits: int
    auth: str
    description: str
    digest: str | None
    id: int
    kea: str
    name: str
    protocol: str
    strength_bits: int
    symmetric: str

class SSLError(OSError):
    library: str
    reason: str

class SSLZeroReturnError(SSLError): ...
class SSLWantReadError(SSLError): ...
class SSLWantWriteError(SSLError): ...
class SSLSyscallError(SSLError): ...
class SSLEOFError(SSLError): ...

if sys.version_info >= (3, 7):
    class SSLCertVerificationError(SSLError, ValueError):
        verify_code: int
        verify_message: str
    CertificateError = SSLCertVerificationError
else:
    class CertificateError(ValueError): ...

def wrap_socket(
    sock: socket.socket,
    keyfile: StrOrBytesPath | None = ...,
    certfile: StrOrBytesPath | None = ...,
    server_side: bool = ...,
    cert_reqs: int = ...,
    ssl_version: int = ...,
    ca_certs: str | None = ...,
    do_handshake_on_connect: bool = ...,
    suppress_ragged_eofs: bool = ...,
    ciphers: str | None = ...,
) -> SSLSocket: ...
def create_default_context(
    purpose: Purpose = ...,
    *,
    cafile: StrOrBytesPath | None = ...,
    capath: StrOrBytesPath | None = ...,
    cadata: str | bytes | None = ...,
) -> SSLContext: ...

if sys.version_info >= (3, 7):
    def _create_unverified_context(
        protocol: int = ...,
        *,
        cert_reqs: int = ...,
        check_hostname: bool = ...,
        purpose: Purpose = ...,
        certfile: StrOrBytesPath | None = ...,
        keyfile: StrOrBytesPath | None = ...,
        cafile: StrOrBytesPath | None = ...,
        capath: StrOrBytesPath | None = ...,
        cadata: str | bytes | None = ...,
    ) -> SSLContext: ...

else:
    def _create_unverified_context(
        protocol: int = ...,
        *,
        cert_reqs: int | None = ...,
        check_hostname: bool = ...,
        purpose: Purpose = ...,
        certfile: StrOrBytesPath | None = ...,
        keyfile: StrOrBytesPath | None = ...,
        cafile: StrOrBytesPath | None = ...,
        capath: StrOrBytesPath | None = ...,
        cadata: str | bytes | None = ...,
    ) -> SSLContext: ...

_create_default_https_context: Callable[..., SSLContext]

def RAND_bytes(__num: int) -> bytes: ...
def RAND_pseudo_bytes(__num: int) -> tuple[bytes, bool]: ...
def RAND_status() -> bool: ...
def RAND_egd(path: str) -> None: ...
def RAND_add(__s: bytes, __entropy: float) -> None: ...
def match_hostname(cert: _PeerCertRetType, hostname: str) -> None: ...
def cert_time_to_seconds(cert_time: str) -> int: ...
def get_server_certificate(addr: tuple[str, int], ssl_version: int = ..., ca_certs: str | None = ...) -> str: ...
def DER_cert_to_PEM_cert(der_cert_bytes: bytes) -> str: ...
def PEM_cert_to_DER_cert(pem_cert_string: str) -> bytes: ...

class DefaultVerifyPaths(NamedTuple):
    cafile: str
    capath: str
    openssl_cafile_env: str
    openssl_cafile: str
    openssl_capath_env: str
    openssl_capath: str

def get_default_verify_paths() -> DefaultVerifyPaths: ...

if sys.platform == "win32":
    def enum_certificates(store_name: str) -> _EnumRetType: ...
    def enum_crls(store_name: str) -> _EnumRetType: ...

class VerifyMode(enum.IntEnum):
    CERT_NONE: int
    CERT_OPTIONAL: int
    CERT_REQUIRED: int

CERT_NONE: VerifyMode
CERT_OPTIONAL: VerifyMode
CERT_REQUIRED: VerifyMode

class VerifyFlags(enum.IntFlag):
    VERIFY_DEFAULT: int
    VERIFY_CRL_CHECK_LEAF: int
    VERIFY_CRL_CHECK_CHAIN: int
    VERIFY_X509_STRICT: int
    VERIFY_X509_TRUSTED_FIRST: int

VERIFY_DEFAULT: VerifyFlags
VERIFY_CRL_CHECK_LEAF: VerifyFlags
VERIFY_CRL_CHECK_CHAIN: VerifyFlags
VERIFY_X509_STRICT: VerifyFlags
VERIFY_X509_TRUSTED_FIRST: VerifyFlags

class _SSLMethod(enum.IntEnum):
    PROTOCOL_SSLv23: int
    PROTOCOL_SSLv2: int
    PROTOCOL_SSLv3: int
    PROTOCOL_TLSv1: int
    PROTOCOL_TLSv1_1: int
    PROTOCOL_TLSv1_2: int
    PROTOCOL_TLS: int
    PROTOCOL_TLS_CLIENT: int
    PROTOCOL_TLS_SERVER: int

PROTOCOL_SSLv23: _SSLMethod
PROTOCOL_SSLv2: _SSLMethod
PROTOCOL_SSLv3: _SSLMethod
PROTOCOL_TLSv1: _SSLMethod
PROTOCOL_TLSv1_1: _SSLMethod
PROTOCOL_TLSv1_2: _SSLMethod
PROTOCOL_TLS: _SSLMethod
PROTOCOL_TLS_CLIENT: _SSLMethod
PROTOCOL_TLS_SERVER: _SSLMethod

class Options(enum.IntFlag):
    OP_ALL: int
    OP_NO_SSLv2: int
    OP_NO_SSLv3: int
    OP_NO_TLSv1: int
    OP_NO_TLSv1_1: int
    OP_NO_TLSv1_2: int
    OP_NO_TLSv1_3: int
    OP_CIPHER_SERVER_PREFERENCE: int
    OP_SINGLE_DH_USE: int
    OP_SINGLE_ECDH_USE: int
    OP_NO_COMPRESSION: int
    OP_NO_TICKET: int
    if sys.version_info >= (3, 7):
        OP_NO_RENEGOTIATION: int
    if sys.version_info >= (3, 8):
        OP_ENABLE_MIDDLEBOX_COMPAT: int

OP_ALL: Options
OP_NO_SSLv2: Options
OP_NO_SSLv3: Options
OP_NO_TLSv1: Options
OP_NO_TLSv1_1: Options
OP_NO_TLSv1_2: Options
OP_NO_TLSv1_3: Options
OP_CIPHER_SERVER_PREFERENCE: Options
OP_SINGLE_DH_USE: Options
OP_SINGLE_ECDH_USE: Options
OP_NO_COMPRESSION: Options
OP_NO_TICKET: Options
if sys.version_info >= (3, 7):
    OP_NO_RENEGOTIATION: Options
if sys.version_info >= (3, 8):
    OP_ENABLE_MIDDLEBOX_COMPAT: Options

if sys.version_info >= (3, 7):
    HAS_NEVER_CHECK_COMMON_NAME: bool
    HAS_SSLv2: bool
    HAS_SSLv3: bool
    HAS_TLSv1: bool
    HAS_TLSv1_1: bool
    HAS_TLSv1_2: bool
    HAS_TLSv1_3: bool
HAS_ALPN: bool
HAS_ECDH: bool
HAS_SNI: bool
HAS_NPN: bool
CHANNEL_BINDING_TYPES: list[str]

OPENSSL_VERSION: str
OPENSSL_VERSION_INFO: tuple[int, int, int, int, int]
OPENSSL_VERSION_NUMBER: int

class AlertDescription(enum.IntEnum):
    ALERT_DESCRIPTION_ACCESS_DENIED: int
    ALERT_DESCRIPTION_BAD_CERTIFICATE: int
    ALERT_DESCRIPTION_BAD_CERTIFICATE_HASH_VALUE: int
    ALERT_DESCRIPTION_BAD_CERTIFICATE_STATUS_RESPONSE: int
    ALERT_DESCRIPTION_BAD_RECORD_MAC: int
    ALERT_DESCRIPTION_CERTIFICATE_EXPIRED: int
    ALERT_DESCRIPTION_CERTIFICATE_REVOKED: int
    ALERT_DESCRIPTION_CERTIFICATE_UNKNOWN: int
    ALERT_DESCRIPTION_CERTIFICATE_UNOBTAINABLE: int
    ALERT_DESCRIPTION_CLOSE_NOTIFY: int
    ALERT_DESCRIPTION_DECODE_ERROR: int
    ALERT_DESCRIPTION_DECOMPRESSION_FAILURE: int
    ALERT_DESCRIPTION_DECRYPT_ERROR: int
    ALERT_DESCRIPTION_HANDSHAKE_FAILURE: int
    ALERT_DESCRIPTION_ILLEGAL_PARAMETER: int
    ALERT_DESCRIPTION_INSUFFICIENT_SECURITY: int
    ALERT_DESCRIPTION_INTERNAL_ERROR: int
    ALERT_DESCRIPTION_NO_RENEGOTIATION: int
    ALERT_DESCRIPTION_PROTOCOL_VERSION: int
    ALERT_DESCRIPTION_RECORD_OVERFLOW: int
    ALERT_DESCRIPTION_UNEXPECTED_MESSAGE: int
    ALERT_DESCRIPTION_UNKNOWN_CA: int
    ALERT_DESCRIPTION_UNKNOWN_PSK_IDENTITY: int
    ALERT_DESCRIPTION_UNRECOGNIZED_NAME: int
    ALERT_DESCRIPTION_UNSUPPORTED_CERTIFICATE: int
    ALERT_DESCRIPTION_UNSUPPORTED_EXTENSION: int
    ALERT_DESCRIPTION_USER_CANCELLED: int

ALERT_DESCRIPTION_HANDSHAKE_FAILURE: AlertDescription
ALERT_DESCRIPTION_INTERNAL_ERROR: AlertDescription
ALERT_DESCRIPTION_ACCESS_DENIED: AlertDescription
ALERT_DESCRIPTION_BAD_CERTIFICATE: AlertDescription
ALERT_DESCRIPTION_BAD_CERTIFICATE_HASH_VALUE: AlertDescription
ALERT_DESCRIPTION_BAD_CERTIFICATE_STATUS_RESPONSE: AlertDescription
ALERT_DESCRIPTION_BAD_RECORD_MAC: AlertDescription
ALERT_DESCRIPTION_CERTIFICATE_EXPIRED: AlertDescription
ALERT_DESCRIPTION_CERTIFICATE_REVOKED: AlertDescription
ALERT_DESCRIPTION_CERTIFICATE_UNKNOWN: AlertDescription
ALERT_DESCRIPTION_CERTIFICATE_UNOBTAINABLE: AlertDescription
ALERT_DESCRIPTION_CLOSE_NOTIFY: AlertDescription
ALERT_DESCRIPTION_DECODE_ERROR: AlertDescription
ALERT_DESCRIPTION_DECOMPRESSION_FAILURE: AlertDescription
ALERT_DESCRIPTION_DECRYPT_ERROR: AlertDescription
ALERT_DESCRIPTION_ILLEGAL_PARAMETER: AlertDescription
ALERT_DESCRIPTION_INSUFFICIENT_SECURITY: AlertDescription
ALERT_DESCRIPTION_NO_RENEGOTIATION: AlertDescription
ALERT_DESCRIPTION_PROTOCOL_VERSION: AlertDescription
ALERT_DESCRIPTION_RECORD_OVERFLOW: AlertDescription
ALERT_DESCRIPTION_UNEXPECTED_MESSAGE: AlertDescription
ALERT_DESCRIPTION_UNKNOWN_CA: AlertDescription
ALERT_DESCRIPTION_UNKNOWN_PSK_IDENTITY: AlertDescription
ALERT_DESCRIPTION_UNRECOGNIZED_NAME: AlertDescription
ALERT_DESCRIPTION_UNSUPPORTED_CERTIFICATE: AlertDescription
ALERT_DESCRIPTION_UNSUPPORTED_EXTENSION: AlertDescription
ALERT_DESCRIPTION_USER_CANCELLED: AlertDescription

class _ASN1Object(NamedTuple):
    nid: int
    shortname: str
    longname: str
    oid: str
    @classmethod
    def fromnid(cls: Type[Self], nid: int) -> Self: ...
    @classmethod
    def fromname(cls: Type[Self], name: str) -> Self: ...

class Purpose(_ASN1Object, enum.Enum):
    SERVER_AUTH: _ASN1Object
    CLIENT_AUTH: _ASN1Object

class SSLSocket(socket.socket):
    context: SSLContext
    server_side: bool
    server_hostname: str | None
    session: SSLSession | None
    session_reused: bool | None
    if sys.version_info < (3, 7):
        def __init__(
            self,
            sock: socket.socket | None = ...,
            keyfile: str | None = ...,
            certfile: str | None = ...,
            server_side: bool = ...,
            cert_reqs: int = ...,
            ssl_version: int = ...,
            ca_certs: str | None = ...,
            do_handshake_on_connect: bool = ...,
            family: int = ...,
            type: int = ...,
            proto: int = ...,
            fileno: int | None = ...,
            suppress_ragged_eofs: bool = ...,
            npn_protocols: Iterable[str] | None = ...,
            ciphers: str | None = ...,
            server_hostname: str | None = ...,
            _context: SSLContext | None = ...,
            _session: Any | None = ...,
        ) -> None: ...
    else:
        def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def connect(self, addr: socket._Address | bytes) -> None: ...
    def connect_ex(self, addr: socket._Address | bytes) -> int: ...
    def recv(self, buflen: int = ..., flags: int = ...) -> bytes: ...
    def recv_into(self, buffer: WriteableBuffer, nbytes: int | None = ..., flags: int = ...) -> int: ...
    def recvfrom(self, buflen: int = ..., flags: int = ...) -> tuple[bytes, socket._RetAddress]: ...
    def recvfrom_into(
        self, buffer: WriteableBuffer, nbytes: int | None = ..., flags: int = ...
    ) -> tuple[int, socket._RetAddress]: ...
    def send(self, data: ReadableBuffer, flags: int = ...) -> int: ...
    def sendall(self, data: ReadableBuffer, flags: int = ...) -> None: ...
    @overload
    def sendto(self, data: ReadableBuffer, flags_or_addr: socket._Address) -> int: ...
    @overload
    def sendto(self, data: ReadableBuffer, flags_or_addr: int | socket._Address, addr: socket._Address | None = ...) -> int: ...
    def shutdown(self, how: int) -> None: ...
    def read(self, len: int = ..., buffer: bytearray | None = ...) -> bytes: ...
    def write(self, data: bytes) -> int: ...
    def do_handshake(self, block: bool = ...) -> None: ...  # block is undocumented
    @overload
    def getpeercert(self, binary_form: Literal[False] = ...) -> _PeerCertRetDictType | None: ...
    @overload
    def getpeercert(self, binary_form: Literal[True]) -> bytes | None: ...
    @overload
    def getpeercert(self, binary_form: bool) -> _PeerCertRetType: ...
    def cipher(self) -> tuple[str, str, int] | None: ...
    def shared_ciphers(self) -> list[tuple[str, str, int]] | None: ...
    def compression(self) -> str | None: ...
    def get_channel_binding(self, cb_type: str = ...) -> bytes | None: ...
    def selected_alpn_protocol(self) -> str | None: ...
    def selected_npn_protocol(self) -> str | None: ...
    def accept(self) -> tuple[SSLSocket, socket._RetAddress]: ...
    def unwrap(self) -> socket.socket: ...
    def version(self) -> str | None: ...
    def pending(self) -> int: ...
    if sys.version_info >= (3, 8):
        def verify_client_post_handshake(self) -> None: ...

if sys.version_info >= (3, 7):
    class TLSVersion(enum.IntEnum):
        MINIMUM_SUPPORTED: int
        MAXIMUM_SUPPORTED: int
        SSLv3: int
        TLSv1: int
        TLSv1_1: int
        TLSv1_2: int
        TLSv1_3: int

class SSLContext:
    check_hostname: bool
    options: Options
    verify_flags: VerifyFlags
    verify_mode: VerifyMode
    @property
    def protocol(self) -> _SSLMethod: ...
    if sys.version_info >= (3, 7):
        hostname_checks_common_name: bool
        maximum_version: TLSVersion
        minimum_version: TLSVersion
        sni_callback: Callable[[SSLObject, str, SSLContext], None | int] | None
        sslobject_class: ClassVar[Type[SSLObject]]
        sslsocket_class: ClassVar[Type[SSLSocket]]
    if sys.version_info >= (3, 8):
        keylog_filename: str
        post_handshake_auth: bool
    def __new__(cls, protocol: int = ..., *args: Any, **kwargs: Any) -> SSLContext: ...
    def __init__(self, protocol: int = ...) -> None: ...
    def cert_store_stats(self) -> dict[str, int]: ...
    def load_cert_chain(
        self, certfile: StrOrBytesPath, keyfile: StrOrBytesPath | None = ..., password: _PasswordType | None = ...
    ) -> None: ...
    def load_default_certs(self, purpose: Purpose = ...) -> None: ...
    def load_verify_locations(
        self, cafile: StrOrBytesPath | None = ..., capath: StrOrBytesPath | None = ..., cadata: str | bytes | None = ...
    ) -> None: ...
    def get_ca_certs(self, binary_form: bool = ...) -> list[_PeerCertRetDictType] | list[bytes]: ...
    def get_ciphers(self) -> list[_Cipher]: ...
    def set_default_verify_paths(self) -> None: ...
    def set_ciphers(self, __cipherlist: str) -> None: ...
    def set_alpn_protocols(self, alpn_protocols: Iterable[str]) -> None: ...
    def set_npn_protocols(self, npn_protocols: Iterable[str]) -> None: ...
    if sys.version_info >= (3, 7):
        def set_servername_callback(self, server_name_callback: _SrvnmeCbType | None) -> None: ...
    else:
        def set_servername_callback(self, __method: _SrvnmeCbType | None) -> None: ...
    def load_dh_params(self, __path: str) -> None: ...
    def set_ecdh_curve(self, __name: str) -> None: ...
    def wrap_socket(
        self,
        sock: socket.socket,
        server_side: bool = ...,
        do_handshake_on_connect: bool = ...,
        suppress_ragged_eofs: bool = ...,
        server_hostname: str | None = ...,
        session: SSLSession | None = ...,
    ) -> SSLSocket: ...
    def wrap_bio(
        self,
        incoming: MemoryBIO,
        outgoing: MemoryBIO,
        server_side: bool = ...,
        server_hostname: str | None = ...,
        session: SSLSession | None = ...,
    ) -> SSLObject: ...
    def session_stats(self) -> dict[str, int]: ...

class SSLObject:
    context: SSLContext
    server_side: bool
    server_hostname: str | None
    session: SSLSession | None
    session_reused: bool
    if sys.version_info >= (3, 7):
        def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    else:
        def __init__(self, sslobj: Any, owner: SSLSocket | SSLObject | None = ..., session: Any | None = ...) -> None: ...
    def read(self, len: int = ..., buffer: bytearray | None = ...) -> bytes: ...
    def write(self, data: bytes) -> int: ...
    @overload
    def getpeercert(self, binary_form: Literal[False] = ...) -> _PeerCertRetDictType | None: ...
    @overload
    def getpeercert(self, binary_form: Literal[True]) -> bytes | None: ...
    @overload
    def getpeercert(self, binary_form: bool) -> _PeerCertRetType: ...
    def selected_alpn_protocol(self) -> str | None: ...
    def selected_npn_protocol(self) -> str | None: ...
    def cipher(self) -> tuple[str, str, int] | None: ...
    def shared_ciphers(self) -> list[tuple[str, str, int]] | None: ...
    def compression(self) -> str | None: ...
    def pending(self) -> int: ...
    def do_handshake(self) -> None: ...
    def unwrap(self) -> None: ...
    def version(self) -> str | None: ...
    def get_channel_binding(self, cb_type: str = ...) -> bytes | None: ...
    if sys.version_info >= (3, 8):
        def verify_client_post_handshake(self) -> None: ...

class MemoryBIO:
    pending: int
    eof: bool
    def read(self, __size: int = ...) -> bytes: ...
    def write(self, __buf: bytes) -> int: ...
    def write_eof(self) -> None: ...

class SSLSession:
    id: bytes
    time: int
    timeout: int
    ticket_lifetime_hint: int
    has_ticket: bool

class SSLErrorNumber(enum.IntEnum):
    SSL_ERROR_EOF: int
    SSL_ERROR_INVALID_ERROR_CODE: int
    SSL_ERROR_SSL: int
    SSL_ERROR_SYSCALL: int
    SSL_ERROR_WANT_CONNECT: int
    SSL_ERROR_WANT_READ: int
    SSL_ERROR_WANT_WRITE: int
    SSL_ERROR_WANT_X509_LOOKUP: int
    SSL_ERROR_ZERO_RETURN: int

SSL_ERROR_EOF: SSLErrorNumber  # undocumented
SSL_ERROR_INVALID_ERROR_CODE: SSLErrorNumber  # undocumented
SSL_ERROR_SSL: SSLErrorNumber  # undocumented
SSL_ERROR_SYSCALL: SSLErrorNumber  # undocumented
SSL_ERROR_WANT_CONNECT: SSLErrorNumber  # undocumented
SSL_ERROR_WANT_READ: SSLErrorNumber  # undocumented
SSL_ERROR_WANT_WRITE: SSLErrorNumber  # undocumented
SSL_ERROR_WANT_X509_LOOKUP: SSLErrorNumber  # undocumented
SSL_ERROR_ZERO_RETURN: SSLErrorNumber  # undocumented

def get_protocol_name(protocol_code: int) -> str: ...

if sys.version_info < (3, 9):
    AF_INET: int
PEM_FOOTER: str
PEM_HEADER: str
SOCK_STREAM: int
SOL_SOCKET: int
SO_TYPE: int
