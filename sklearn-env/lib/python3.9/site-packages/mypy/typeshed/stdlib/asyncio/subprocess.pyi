import subprocess
import sys
from _typeshed import StrOrBytesPath
from asyncio import events, protocols, streams, transports
from typing import IO, Any, Callable, Union
from typing_extensions import Literal

if sys.version_info >= (3, 8):
    _ExecArg = StrOrBytesPath
else:
    _ExecArg = Union[str, bytes]

PIPE: int
STDOUT: int
DEVNULL: int

class SubprocessStreamProtocol(streams.FlowControlMixin, protocols.SubprocessProtocol):
    stdin: streams.StreamWriter | None
    stdout: streams.StreamReader | None
    stderr: streams.StreamReader | None
    def __init__(self, limit: int, loop: events.AbstractEventLoop) -> None: ...
    def connection_made(self, transport: transports.BaseTransport) -> None: ...
    def pipe_data_received(self, fd: int, data: bytes | str) -> None: ...
    def pipe_connection_lost(self, fd: int, exc: Exception | None) -> None: ...
    def process_exited(self) -> None: ...

class Process:
    stdin: streams.StreamWriter | None
    stdout: streams.StreamReader | None
    stderr: streams.StreamReader | None
    pid: int
    def __init__(
        self, transport: transports.BaseTransport, protocol: protocols.BaseProtocol, loop: events.AbstractEventLoop
    ) -> None: ...
    @property
    def returncode(self) -> int | None: ...
    async def wait(self) -> int: ...
    def send_signal(self, signal: int) -> None: ...
    def terminate(self) -> None: ...
    def kill(self) -> None: ...
    async def communicate(self, input: bytes | None = ...) -> tuple[bytes, bytes]: ...

if sys.version_info >= (3, 10):
    async def create_subprocess_shell(
        cmd: str | bytes,
        stdin: int | IO[Any] | None = ...,
        stdout: int | IO[Any] | None = ...,
        stderr: int | IO[Any] | None = ...,
        limit: int = ...,
        *,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        text: Literal[False, None] = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        executable: StrOrBytesPath | None = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: subprocess._ENV | None = ...,
        startupinfo: Any | None = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
    async def create_subprocess_exec(
        program: _ExecArg,
        *args: _ExecArg,
        stdin: int | IO[Any] | None = ...,
        stdout: int | IO[Any] | None = ...,
        stderr: int | IO[Any] | None = ...,
        limit: int = ...,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        text: bool | None = ...,
        executable: StrOrBytesPath | None = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: subprocess._ENV | None = ...,
        startupinfo: Any | None = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...

else:
    async def create_subprocess_shell(
        cmd: str | bytes,
        stdin: int | IO[Any] | None = ...,
        stdout: int | IO[Any] | None = ...,
        stderr: int | IO[Any] | None = ...,
        loop: events.AbstractEventLoop | None = ...,
        limit: int = ...,
        *,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        text: Literal[False, None] = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        executable: StrOrBytesPath | None = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: subprocess._ENV | None = ...,
        startupinfo: Any | None = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
    async def create_subprocess_exec(
        program: _ExecArg,
        *args: _ExecArg,
        stdin: int | IO[Any] | None = ...,
        stdout: int | IO[Any] | None = ...,
        stderr: int | IO[Any] | None = ...,
        loop: events.AbstractEventLoop | None = ...,
        limit: int = ...,
        # These parameters are forced to these values by BaseEventLoop.subprocess_shell
        universal_newlines: Literal[False] = ...,
        shell: Literal[True] = ...,
        bufsize: Literal[0] = ...,
        encoding: None = ...,
        errors: None = ...,
        # These parameters are taken by subprocess.Popen, which this ultimately delegates to
        text: bool | None = ...,
        executable: StrOrBytesPath | None = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: subprocess._ENV | None = ...,
        startupinfo: Any | None = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
    ) -> Process: ...
