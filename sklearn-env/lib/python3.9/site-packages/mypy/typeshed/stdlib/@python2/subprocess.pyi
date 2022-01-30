from typing import IO, Any, Callable, Generic, Mapping, Optional, Sequence, Text, Tuple, TypeVar, Union

_FILE = Union[None, int, IO[Any]]
_TXT = Union[bytes, Text]
_CMD = Union[_TXT, Sequence[_TXT]]
_ENV = Union[Mapping[bytes, _TXT], Mapping[Text, _TXT]]

# Same args as Popen.__init__
def call(
    args: _CMD,
    bufsize: int = ...,
    executable: _TXT = ...,
    stdin: _FILE = ...,
    stdout: _FILE = ...,
    stderr: _FILE = ...,
    preexec_fn: Callable[[], Any] = ...,
    close_fds: bool = ...,
    shell: bool = ...,
    cwd: _TXT | None = ...,
    env: _ENV | None = ...,
    universal_newlines: bool = ...,
    startupinfo: Any = ...,
    creationflags: int = ...,
) -> int: ...
def check_call(
    args: _CMD,
    bufsize: int = ...,
    executable: _TXT = ...,
    stdin: _FILE = ...,
    stdout: _FILE = ...,
    stderr: _FILE = ...,
    preexec_fn: Callable[[], Any] = ...,
    close_fds: bool = ...,
    shell: bool = ...,
    cwd: _TXT | None = ...,
    env: _ENV | None = ...,
    universal_newlines: bool = ...,
    startupinfo: Any = ...,
    creationflags: int = ...,
) -> int: ...

# Same args as Popen.__init__ except for stdout
def check_output(
    args: _CMD,
    bufsize: int = ...,
    executable: _TXT = ...,
    stdin: _FILE = ...,
    stderr: _FILE = ...,
    preexec_fn: Callable[[], Any] = ...,
    close_fds: bool = ...,
    shell: bool = ...,
    cwd: _TXT | None = ...,
    env: _ENV | None = ...,
    universal_newlines: bool = ...,
    startupinfo: Any = ...,
    creationflags: int = ...,
) -> bytes: ...

PIPE: int
STDOUT: int

class CalledProcessError(Exception):
    returncode: int
    # morally: _CMD
    cmd: Any
    # morally: Optional[bytes]
    output: bytes
    def __init__(self, returncode: int, cmd: _CMD, output: bytes | None = ...) -> None: ...

# We use a dummy type variable used to make Popen generic like it is in python 3
_T = TypeVar("_T", bound=bytes)

class Popen(Generic[_T]):
    stdin: IO[bytes] | None
    stdout: IO[bytes] | None
    stderr: IO[bytes] | None
    pid: int
    returncode: int
    def __new__(
        cls,
        args: _CMD,
        bufsize: int = ...,
        executable: _TXT | None = ...,
        stdin: _FILE | None = ...,
        stdout: _FILE | None = ...,
        stderr: _FILE | None = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: _TXT | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any | None = ...,
        creationflags: int = ...,
    ) -> Popen[bytes]: ...
    def poll(self) -> int | None: ...
    def wait(self) -> int: ...
    # morally: -> Tuple[Optional[bytes], Optional[bytes]]
    def communicate(self, input: _TXT | None = ...) -> Tuple[bytes, bytes]: ...
    def send_signal(self, signal: int) -> None: ...
    def terminate(self) -> None: ...
    def kill(self) -> None: ...

def list2cmdline(seq: Sequence[str]) -> str: ...  # undocumented

# Windows-only: STARTUPINFO etc.

STD_INPUT_HANDLE: Any
STD_OUTPUT_HANDLE: Any
STD_ERROR_HANDLE: Any
SW_HIDE: Any
STARTF_USESTDHANDLES: Any
STARTF_USESHOWWINDOW: Any
CREATE_NEW_CONSOLE: Any
CREATE_NEW_PROCESS_GROUP: Any
