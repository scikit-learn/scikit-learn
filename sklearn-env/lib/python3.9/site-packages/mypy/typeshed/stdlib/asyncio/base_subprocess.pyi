import subprocess
from collections import deque
from typing import IO, Any, Callable, Optional, Sequence, Tuple, Union

from . import events, futures, protocols, transports

_File = Optional[Union[int, IO[Any]]]

class BaseSubprocessTransport(transports.SubprocessTransport):

    _closed: bool  # undocumented
    _protocol: protocols.SubprocessProtocol  # undocumented
    _loop: events.AbstractEventLoop  # undocumented
    _proc: subprocess.Popen[Any] | None  # undocumented
    _pid: int | None  # undocumented
    _returncode: int | None  # undocumented
    _exit_waiters: list[futures.Future[Any]]  # undocumented
    _pending_calls: deque[tuple[Callable[..., Any], Tuple[Any, ...]]]  # undocumented
    _pipes: dict[int, _File]  # undocumented
    _finished: bool  # undocumented
    def __init__(
        self,
        loop: events.AbstractEventLoop,
        protocol: protocols.SubprocessProtocol,
        args: str | bytes | Sequence[str | bytes],
        shell: bool,
        stdin: _File,
        stdout: _File,
        stderr: _File,
        bufsize: int,
        waiter: futures.Future[Any] | None = ...,
        extra: Any | None = ...,
        **kwargs: Any,
    ) -> None: ...
    def _start(
        self,
        args: str | bytes | Sequence[str | bytes],
        shell: bool,
        stdin: _File,
        stdout: _File,
        stderr: _File,
        bufsize: int,
        **kwargs: Any,
    ) -> None: ...  # undocumented
    def set_protocol(self, protocol: protocols.BaseProtocol) -> None: ...
    def get_protocol(self) -> protocols.BaseProtocol: ...
    def is_closing(self) -> bool: ...
    def close(self) -> None: ...
    def get_pid(self) -> int | None: ...  # type: ignore
    def get_returncode(self) -> int | None: ...
    def get_pipe_transport(self, fd: int) -> _File: ...  # type: ignore
    def _check_proc(self) -> None: ...  # undocumented
    def send_signal(self, signal: int) -> None: ...  # type: ignore
    def terminate(self) -> None: ...
    def kill(self) -> None: ...
    async def _connect_pipes(self, waiter: futures.Future[Any] | None) -> None: ...  # undocumented
    def _call(self, cb: Callable[..., Any], *data: Any) -> None: ...  # undocumented
    def _pipe_connection_lost(self, fd: int, exc: BaseException | None) -> None: ...  # undocumented
    def _pipe_data_received(self, fd: int, data: bytes) -> None: ...  # undocumented
    def _process_exited(self, returncode: int) -> None: ...  # undocumented
    async def _wait(self) -> int: ...  # undocumented
    def _try_finish(self) -> None: ...  # undocumented
    def _call_connection_lost(self, exc: BaseException | None) -> None: ...  # undocumented

class WriteSubprocessPipeProto(protocols.BaseProtocol):  # undocumented
    def __init__(self, proc: BaseSubprocessTransport, fd: int) -> None: ...
    def connection_made(self, transport: transports.BaseTransport) -> None: ...
    def connection_lost(self, exc: BaseException | None) -> None: ...
    def pause_writing(self) -> None: ...
    def resume_writing(self) -> None: ...

class ReadSubprocessPipeProto(WriteSubprocessPipeProto, protocols.Protocol):  # undocumented
    def data_received(self, data: bytes) -> None: ...
