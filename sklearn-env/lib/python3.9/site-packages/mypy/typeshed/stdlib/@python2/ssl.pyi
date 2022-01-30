import socket
import sys
from _typeshed import StrPath
from typing import Any, Callable, ClassVar, Dict, Iterable, List, NamedTuple, Optional, Set, Text, Tuple, Union, overload
from typing_extensions import Literal

_PCTRTT = Tuple[Tuple[str, str], ...]
_PCTRTTT = Tuple[_PCTRTT, ...]
_PeerCertRetDictType = Dict[str, Union[str, _PCTRTTT, _PCTRTT]]
_PeerCertRetType = Union[_PeerCertRetDictType, bytes, None]
_EnumRetType = List[Tuple[bytes, str, Union[Set[str], bool]]]
_PasswordType = Union[Callable[[], Union[str, bytes]], str, bytes]

_SC1ArgT = SSLSocket
_SrvnmeCbType = Callable[[_SC1ArgT, Optional[str], SSLSocket], Optional[int]]

class SSLError(OSError):
    library: str
    reason: str

class SSLZeroReturnError(SSLError): ...
class SSLWantReadError(SSLError): ...
class SSLWantWriteError(SSLError): ...
class SSLSyscallError(SSLError): ...
class SSLEOFError(SSLError): ...
class CertificateError(ValueError): ...

def wrap_socket(
    sock: socket.socket,
    keyfile: str | None = ...,
    certfile: str | None = ...,
    server_side: bool = ...,
    cert_reqs: int = ...,
    ssl_version: int = ...,
    ca_certs: str | None = ...,
    do_handshake_on_connect: bool = ...,
    suppress_ragged_eofs: bool = ...,
    ciphers: str | None = ...,
) -> SSLSocket: ...
def create_default_context(
    purpose: Any = ..., *, cafile: str | None = ..., capath: str | None = ..., cadata: Text | bytes | None = ...
) -> SSLContext: ...
def _create_unverified_context(
    protocol: int = ...,
    *,
    cert_reqs: int | None = ...,
    check_hostname: bool = ...,
    purpose: Any = ...,
    certfile: str | None = ...,
    keyfile: str | None = ...,
    cafile: str | None = ...,
    capath: str | None = ...,
    cadata: Text | bytes | None = ...,
) -> SSLContext: ...

_create_default_https_context: Callable[..., SSLContext]

def RAND_status() -> bool: ...
def RAND_egd(path: str) -> None: ...
def RAND_add(__s: bytes, __entropy: float) -> None: ...
def match_hostname(cert: _PeerCertRetType, hostname: str) -> None: ...
def cert_time_to_seconds(cert_time: str) -> int: ...
def get_server_certificate(addr: Tuple[str, int], ssl_version: int = ..., ca_certs: str | None = ...) -> str: ...
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

CERT_NONE: int
CERT_OPTIONAL: int
CERT_REQUIRED: int

VERIFY_DEFAULT: int
VERIFY_CRL_CHECK_LEAF: int
VERIFY_CRL_CHECK_CHAIN: int
VERIFY_X509_STRICT: int
VERIFY_X509_TRUSTED_FIRST: int

PROTOCOL_SSLv23: int
PROTOCOL_SSLv2: int
PROTOCOL_SSLv3: int
PROTOCOL_TLSv1: int
PROTOCOL_TLSv1_1: int
PROTOCOL_TLSv1_2: int
PROTOCOL_TLS: int
OP_ALL: int
OP_NO_SSLv2: int
OP_NO_SSLv3: int
OP_NO_TLSv1: int
OP_NO_TLSv1_1: int
OP_NO_TLSv1_2: int
OP_CIPHER_SERVER_PREFERENCE: int
OP_SINGLE_DH_USE: int
OP_SINGLE_ECDH_USE: int
OP_NO_COMPRESSION: int

HAS_ALPN: bool
HAS_ECDH: bool
HAS_SNI: bool
HAS_NPN: bool
CHANNEL_BINDING_TYPES: List[str]

OPENSSL_VERSION: str
OPENSSL_VERSION_INFO: Tuple[int, int, int, int, int]
OPENSSL_VERSION_NUMBER: int

ALERT_DESCRIPTION_HANDSHAKE_FAILURE: int
ALERT_DESCRIPTION_INTERNAL_ERROR: int
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
ALERT_DESCRIPTION_ILLEGAL_PARAMETER: int
ALERT_DESCRIPTION_INSUFFICIENT_SECURITY: int
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

class _ASN1Object(NamedTuple):
    nid: int
    shortname: str
    longname: str
    oid: str

class Purpose(_ASN1Object):
    SERVER_AUTH: ClassVar[Purpose]
    CLIENT_AUTH: ClassVar[Purpose]

class SSLSocket(socket.socket):
    context: SSLContext
    server_side: bool
    server_hostname: str | None
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
    def connect(self, addr: socket._Address | bytes) -> None: ...
    def connect_ex(self, addr: socket._Address | bytes) -> int: ...
    def recv(self, buflen: int = ..., flags: int = ...) -> bytes: ...
    def recv_into(self, buffer: socket._WriteBuffer, nbytes: int | None = ..., flags: int = ...) -> int: ...
    def recvfrom(self, buflen: int = ..., flags: int = ...) -> tuple[bytes, socket._RetAddress]: ...
    def recvfrom_into(
        self, buffer: socket._WriteBuffer, nbytes: int | None = ..., flags: int = ...
    ) -> tuple[int, socket._RetAddress]: ...
    @overload
    def sendto(self, data: bytes, flags_or_addr: socket._Address) -> int: ...
    @overload
    def sendto(self, data: bytes, flags_or_addr: int | socket._Address, addr: socket._Address | None = ...) -> int: ...
    def read(self, len: int = ..., buffer: bytearray | None = ...) -> bytes: ...
    def write(self, data: bytes) -> int: ...
    def do_handshake(self, block: bool = ...) -> None: ...  # block is undocumented
    @overload
    def getpeercert(self, binary_form: Literal[False] = ...) -> _PeerCertRetDictType | None: ...
    @overload
    def getpeercert(self, binary_form: Literal[True]) -> bytes | None: ...
    @overload
    def getpeercert(self, binary_form: bool) -> _PeerCertRetType: ...
    def cipher(self) -> Tuple[str, str, int] | None: ...
    def compression(self) -> str | None: ...
    def get_channel_binding(self, cb_type: str = ...) -> bytes | None: ...
    def selected_alpn_protocol(self) -> str | None: ...
    def selected_npn_protocol(self) -> str | None: ...
    def accept(self) -> Tuple[SSLSocket, socket._RetAddress]: ...
    def unwrap(self) -> socket.socket: ...
    def version(self) -> str | None: ...
    def pending(self) -> int: ...

class SSLContext:
    check_hostname: bool
    options: int
    def __new__(cls, protocol: int, *args: Any, **kwargs: Any) -> SSLContext: ...
    @property
    def protocol(self) -> int: ...
    verify_flags: int
    verify_mode: int
    def __init__(self, protocol: int) -> None: ...
    def cert_store_stats(self) -> Dict[str, int]: ...
    def load_cert_chain(self, certfile: StrPath, keyfile: StrPath | None = ..., password: _PasswordType | None = ...) -> None: ...
    def load_default_certs(self, purpose: Purpose = ...) -> None: ...
    def load_verify_locations(
        self, cafile: StrPath | None = ..., capath: StrPath | None = ..., cadata: Text | bytes | None = ...
    ) -> None: ...
    def get_ca_certs(self, binary_form: bool = ...) -> List[_PeerCertRetDictType] | List[bytes]: ...
    def set_default_verify_paths(self) -> None: ...
    def set_ciphers(self, __cipherlist: str) -> None: ...
    def set_alpn_protocols(self, alpn_protocols: Iterable[str]) -> None: ...
    def set_npn_protocols(self, npn_protocols: Iterable[str]) -> None: ...
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
    ) -> SSLSocket: ...
    def session_stats(self) -> Dict[str, int]: ...

# TODO below documented in cpython but not in docs.python.org
# taken from python 3.4
SSL_ERROR_EOF: int
SSL_ERROR_INVALID_ERROR_CODE: int
SSL_ERROR_SSL: int
SSL_ERROR_SYSCALL: int
SSL_ERROR_WANT_CONNECT: int
SSL_ERROR_WANT_READ: int
SSL_ERROR_WANT_WRITE: int
SSL_ERROR_WANT_X509_LOOKUP: int
SSL_ERROR_ZERO_RETURN: int

def get_protocol_name(protocol_code: int) -> str: ...

AF_INET: int
PEM_FOOTER: str
PEM_HEADER: str
SOCK_STREAM: int
SOL_SOCKET: int
SO_TYPE: int
