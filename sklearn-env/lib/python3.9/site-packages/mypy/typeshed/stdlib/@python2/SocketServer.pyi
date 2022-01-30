import sys
from socket import SocketType
from typing import Any, BinaryIO, Callable, ClassVar, List, Text, Tuple, Union

class BaseServer:
    address_family: int
    RequestHandlerClass: Callable[..., BaseRequestHandler]
    server_address: Tuple[str, int]
    socket: SocketType
    allow_reuse_address: bool
    request_queue_size: int
    socket_type: int
    timeout: float | None
    def __init__(self, server_address: Any, RequestHandlerClass: Callable[..., BaseRequestHandler]) -> None: ...
    def fileno(self) -> int: ...
    def handle_request(self) -> None: ...
    def serve_forever(self, poll_interval: float = ...) -> None: ...
    def shutdown(self) -> None: ...
    def server_close(self) -> None: ...
    def finish_request(self, request: bytes, client_address: Tuple[str, int]) -> None: ...
    def get_request(self) -> Tuple[SocketType, Tuple[str, int]]: ...
    def handle_error(self, request: bytes, client_address: Tuple[str, int]) -> None: ...
    def handle_timeout(self) -> None: ...
    def process_request(self, request: bytes, client_address: Tuple[str, int]) -> None: ...
    def server_activate(self) -> None: ...
    def server_bind(self) -> None: ...
    def verify_request(self, request: bytes, client_address: Tuple[str, int]) -> bool: ...

class TCPServer(BaseServer):
    def __init__(
        self,
        server_address: Tuple[str, int],
        RequestHandlerClass: Callable[..., BaseRequestHandler],
        bind_and_activate: bool = ...,
    ) -> None: ...

class UDPServer(BaseServer):
    def __init__(
        self,
        server_address: Tuple[str, int],
        RequestHandlerClass: Callable[..., BaseRequestHandler],
        bind_and_activate: bool = ...,
    ) -> None: ...

if sys.platform != "win32":
    class UnixStreamServer(BaseServer):
        def __init__(
            self,
            server_address: Text | bytes,
            RequestHandlerClass: Callable[..., BaseRequestHandler],
            bind_and_activate: bool = ...,
        ) -> None: ...
    class UnixDatagramServer(BaseServer):
        def __init__(
            self,
            server_address: Text | bytes,
            RequestHandlerClass: Callable[..., BaseRequestHandler],
            bind_and_activate: bool = ...,
        ) -> None: ...

if sys.platform != "win32":
    class ForkingMixIn:
        timeout: float | None  # undocumented
        active_children: List[int] | None  # undocumented
        max_children: int  # undocumented
        def collect_children(self) -> None: ...  # undocumented
        def handle_timeout(self) -> None: ...  # undocumented
        def process_request(self, request: bytes, client_address: Tuple[str, int]) -> None: ...

class ThreadingMixIn:
    daemon_threads: bool
    def process_request_thread(self, request: bytes, client_address: Tuple[str, int]) -> None: ...  # undocumented
    def process_request(self, request: bytes, client_address: Tuple[str, int]) -> None: ...

if sys.platform != "win32":
    class ForkingTCPServer(ForkingMixIn, TCPServer): ...
    class ForkingUDPServer(ForkingMixIn, UDPServer): ...

class ThreadingTCPServer(ThreadingMixIn, TCPServer): ...
class ThreadingUDPServer(ThreadingMixIn, UDPServer): ...

if sys.platform != "win32":
    class ThreadingUnixStreamServer(ThreadingMixIn, UnixStreamServer): ...
    class ThreadingUnixDatagramServer(ThreadingMixIn, UnixDatagramServer): ...

class BaseRequestHandler:
    # Those are technically of types, respectively:
    # * Union[SocketType, Tuple[bytes, SocketType]]
    # * Union[Tuple[str, int], str]
    # But there are some concerns that having unions here would cause
    # too much inconvenience to people using it (see
    # https://github.com/python/typeshed/pull/384#issuecomment-234649696)
    request: Any
    client_address: Any
    server: BaseServer
    def __init__(self, request: Any, client_address: Any, server: BaseServer) -> None: ...
    def setup(self) -> None: ...
    def handle(self) -> None: ...
    def finish(self) -> None: ...

class StreamRequestHandler(BaseRequestHandler):
    rbufsize: ClassVar[int]  # undocumented
    wbufsize: ClassVar[int]  # undocumented
    timeout: ClassVar[float | None]  # undocumented
    disable_nagle_algorithm: ClassVar[bool]  # undocumented
    connection: SocketType  # undocumented
    rfile: BinaryIO
    wfile: BinaryIO

class DatagramRequestHandler(BaseRequestHandler):
    packet: SocketType  # undocumented
    socket: SocketType  # undocumented
    rfile: BinaryIO
    wfile: BinaryIO
