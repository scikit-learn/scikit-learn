import sys
from _typeshed import Self, StrOrBytesPath
from types import TracebackType
from typing import IO, Any, AnyStr, Callable, Generic, Iterable, Mapping, Sequence, Type, TypeVar, Union, overload
from typing_extensions import Literal

if sys.version_info >= (3, 9):
    from types import GenericAlias

# We prefer to annotate inputs to methods (eg subprocess.check_call) with these
# union types.
# For outputs we use laborious literal based overloads to try to determine
# which specific return types to use, and prefer to fall back to Any when
# this does not work, so the caller does not have to use an assertion to confirm
# which type.
#
# For example:
#
# try:
#    x = subprocess.check_output(["ls", "-l"])
#    reveal_type(x)  # bytes, based on the overloads
# except TimeoutError as e:
#    reveal_type(e.cmd)  # Any, but morally is _CMD
_FILE = Union[None, int, IO[Any]]
_TXT = Union[bytes, str]
if sys.version_info >= (3, 8):
    _CMD = Union[StrOrBytesPath, Sequence[StrOrBytesPath]]
else:
    # Python 3.6 doesn't support _CMD being a single PathLike.
    # See: https://bugs.python.org/issue31961
    _CMD = Union[_TXT, Sequence[StrOrBytesPath]]
if sys.platform == "win32":
    _ENV = Mapping[str, str]
else:
    _ENV = Union[Mapping[bytes, StrOrBytesPath], Mapping[str, StrOrBytesPath]]

_T = TypeVar("_T")

class CompletedProcess(Generic[_T]):
    # morally: _CMD
    args: Any
    returncode: int
    # These can both be None, but requiring checks for None would be tedious
    # and writing all the overloads would be horrific.
    stdout: _T
    stderr: _T
    def __init__(self, args: _CMD, returncode: int, stdout: _T | None = ..., stderr: _T | None = ...) -> None: ...
    def check_returncode(self) -> None: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

if sys.version_info >= (3, 7):
    # Nearly the same args as for 3.6, except for capture_output and text
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        capture_output: bool = ...,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        input: str | None = ...,
        text: Literal[True],
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        capture_output: bool = ...,
        check: bool = ...,
        encoding: str,
        errors: str | None = ...,
        input: str | None = ...,
        text: bool | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        capture_output: bool = ...,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str,
        input: str | None = ...,
        text: bool | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        *,
        universal_newlines: Literal[True],
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        # where the *real* keyword only args start
        capture_output: bool = ...,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        input: str | None = ...,
        text: bool | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: Literal[False] = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        capture_output: bool = ...,
        check: bool = ...,
        encoding: None = ...,
        errors: None = ...,
        input: bytes | None = ...,
        text: Literal[None, False] = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[bytes]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        capture_output: bool = ...,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        input: _TXT | None = ...,
        text: bool | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[Any]: ...

else:
    # Nearly same args as Popen.__init__ except for timeout, input, and check
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        check: bool = ...,
        encoding: str,
        errors: str | None = ...,
        input: str | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str,
        input: str | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        *,
        universal_newlines: Literal[True],
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        # where the *real* keyword only args start
        check: bool = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        input: str | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[str]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: Literal[False] = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        check: bool = ...,
        encoding: None = ...,
        errors: None = ...,
        input: bytes | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[bytes]: ...
    @overload
    def run(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stdout: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        check: bool = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        input: _TXT | None = ...,
        timeout: float | None = ...,
    ) -> CompletedProcess[Any]: ...

# Same args as Popen.__init__
def call(
    args: _CMD,
    bufsize: int = ...,
    executable: StrOrBytesPath | None = ...,
    stdin: _FILE = ...,
    stdout: _FILE = ...,
    stderr: _FILE = ...,
    preexec_fn: Callable[[], Any] | None = ...,
    close_fds: bool = ...,
    shell: bool = ...,
    cwd: StrOrBytesPath | None = ...,
    env: _ENV | None = ...,
    universal_newlines: bool = ...,
    startupinfo: Any = ...,
    creationflags: int = ...,
    restore_signals: bool = ...,
    start_new_session: bool = ...,
    pass_fds: Any = ...,
    *,
    timeout: float | None = ...,
) -> int: ...

# Same args as Popen.__init__
def check_call(
    args: _CMD,
    bufsize: int = ...,
    executable: StrOrBytesPath = ...,
    stdin: _FILE = ...,
    stdout: _FILE = ...,
    stderr: _FILE = ...,
    preexec_fn: Callable[[], Any] | None = ...,
    close_fds: bool = ...,
    shell: bool = ...,
    cwd: StrOrBytesPath | None = ...,
    env: _ENV | None = ...,
    universal_newlines: bool = ...,
    startupinfo: Any = ...,
    creationflags: int = ...,
    restore_signals: bool = ...,
    start_new_session: bool = ...,
    pass_fds: Any = ...,
    timeout: float | None = ...,
) -> int: ...

if sys.version_info >= (3, 7):
    # 3.7 added text
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        text: Literal[True],
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str,
        errors: str | None = ...,
        text: bool | None = ...,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str,
        text: bool | None = ...,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        *,
        universal_newlines: Literal[True],
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        # where the real keyword only ones start
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        text: bool | None = ...,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: Literal[False] = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: None = ...,
        errors: None = ...,
        text: Literal[None, False] = ...,
    ) -> bytes: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
        text: bool | None = ...,
    ) -> Any: ...  # morally: -> _TXT

else:
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str,
        errors: str | None = ...,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        universal_newlines: Literal[True],
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
    ) -> str: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: Literal[False] = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: None = ...,
        errors: None = ...,
    ) -> bytes: ...
    @overload
    def check_output(
        args: _CMD,
        bufsize: int = ...,
        executable: StrOrBytesPath | None = ...,
        stdin: _FILE = ...,
        stderr: _FILE = ...,
        preexec_fn: Callable[[], Any] | None = ...,
        close_fds: bool = ...,
        shell: bool = ...,
        cwd: StrOrBytesPath | None = ...,
        env: _ENV | None = ...,
        universal_newlines: bool = ...,
        startupinfo: Any = ...,
        creationflags: int = ...,
        restore_signals: bool = ...,
        start_new_session: bool = ...,
        pass_fds: Any = ...,
        *,
        timeout: float | None = ...,
        input: _TXT | None = ...,
        encoding: str | None = ...,
        errors: str | None = ...,
    ) -> Any: ...  # morally: -> _TXT

PIPE: int
STDOUT: int
DEVNULL: int

class SubprocessError(Exception): ...

class TimeoutExpired(SubprocessError):
    def __init__(self, cmd: _CMD, timeout: float, output: _TXT | None = ..., stderr: _TXT | None = ...) -> None: ...
    # morally: _CMD
    cmd: Any
    timeout: float
    # morally: _TXT | None
    output: Any
    stdout: Any
    stderr: Any

class CalledProcessError(SubprocessError):
    returncode: int
    # morally: _CMD
    cmd: Any
    # morally: _TXT | None
    output: Any

    # morally: _TXT | None
    stdout: Any
    stderr: Any
    def __init__(self, returncode: int, cmd: _CMD, output: _TXT | None = ..., stderr: _TXT | None = ...) -> None: ...

class Popen(Generic[AnyStr]):
    args: _CMD
    stdin: IO[AnyStr] | None
    stdout: IO[AnyStr] | None
    stderr: IO[AnyStr] | None
    pid: int
    returncode: int
    universal_newlines: bool

    # Technically it is wrong that Popen provides __new__ instead of __init__
    # but this shouldn't come up hopefully?

    if sys.version_info >= (3, 7):
        # text is added in 3.7
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            text: bool | None = ...,
            encoding: str,
            errors: str | None = ...,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            text: bool | None = ...,
            encoding: str | None = ...,
            errors: str,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            *,
            universal_newlines: Literal[True],
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            # where the *real* keyword only args start
            text: bool | None = ...,
            encoding: str | None = ...,
            errors: str | None = ...,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            text: Literal[True],
            encoding: str | None = ...,
            errors: str | None = ...,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: Literal[False] = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            text: Literal[None, False] = ...,
            encoding: None = ...,
            errors: None = ...,
        ) -> Popen[bytes]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            text: bool | None = ...,
            encoding: str | None = ...,
            errors: str | None = ...,
        ) -> Popen[Any]: ...
    else:
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            encoding: str,
            errors: str | None = ...,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            encoding: str | None = ...,
            errors: str,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            *,
            universal_newlines: Literal[True],
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            # where the *real* keyword only args start
            encoding: str | None = ...,
            errors: str | None = ...,
        ) -> Popen[str]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: Literal[False] = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            encoding: None = ...,
            errors: None = ...,
        ) -> Popen[bytes]: ...
        @overload
        def __new__(
            cls,
            args: _CMD,
            bufsize: int = ...,
            executable: StrOrBytesPath | None = ...,
            stdin: _FILE | None = ...,
            stdout: _FILE | None = ...,
            stderr: _FILE | None = ...,
            preexec_fn: Callable[[], Any] | None = ...,
            close_fds: bool = ...,
            shell: bool = ...,
            cwd: StrOrBytesPath | None = ...,
            env: _ENV | None = ...,
            universal_newlines: bool = ...,
            startupinfo: Any | None = ...,
            creationflags: int = ...,
            restore_signals: bool = ...,
            start_new_session: bool = ...,
            pass_fds: Any = ...,
            *,
            encoding: str | None = ...,
            errors: str | None = ...,
        ) -> Popen[Any]: ...
    def poll(self) -> int | None: ...
    if sys.version_info >= (3, 7):
        def wait(self, timeout: float | None = ...) -> int: ...
    else:
        def wait(self, timeout: float | None = ..., endtime: float | None = ...) -> int: ...
    # Return str/bytes
    def communicate(
        self,
        input: AnyStr | None = ...,
        timeout: float | None = ...,
        # morally this should be optional
    ) -> tuple[AnyStr, AnyStr]: ...
    def send_signal(self, sig: int) -> None: ...
    def terminate(self) -> None: ...
    def kill(self) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, type: Type[BaseException] | None, value: BaseException | None, traceback: TracebackType | None
    ) -> None: ...
    if sys.version_info >= (3, 9):
        def __class_getitem__(cls, item: Any) -> GenericAlias: ...

# The result really is always a str.
def getstatusoutput(cmd: _TXT) -> tuple[int, str]: ...
def getoutput(cmd: _TXT) -> str: ...

if sys.version_info >= (3, 8):
    def list2cmdline(seq: Iterable[StrOrBytesPath]) -> str: ...  # undocumented

else:
    def list2cmdline(seq: Iterable[str]) -> str: ...  # undocumented

if sys.platform == "win32":
    class STARTUPINFO:
        if sys.version_info >= (3, 7):
            def __init__(
                self,
                *,
                dwFlags: int = ...,
                hStdInput: Any | None = ...,
                hStdOutput: Any | None = ...,
                hStdError: Any | None = ...,
                wShowWindow: int = ...,
                lpAttributeList: Mapping[str, Any] | None = ...,
            ) -> None: ...
        dwFlags: int
        hStdInput: Any | None
        hStdOutput: Any | None
        hStdError: Any | None
        wShowWindow: int
        if sys.version_info >= (3, 7):
            lpAttributeList: Mapping[str, Any]
    STD_INPUT_HANDLE: Any
    STD_OUTPUT_HANDLE: Any
    STD_ERROR_HANDLE: Any
    SW_HIDE: int
    STARTF_USESTDHANDLES: int
    STARTF_USESHOWWINDOW: int
    CREATE_NEW_CONSOLE: int
    CREATE_NEW_PROCESS_GROUP: int
    if sys.version_info >= (3, 7):
        ABOVE_NORMAL_PRIORITY_CLASS: int
        BELOW_NORMAL_PRIORITY_CLASS: int
        HIGH_PRIORITY_CLASS: int
        IDLE_PRIORITY_CLASS: int
        NORMAL_PRIORITY_CLASS: int
        REALTIME_PRIORITY_CLASS: int
        CREATE_NO_WINDOW: int
        DETACHED_PROCESS: int
        CREATE_DEFAULT_ERROR_MODE: int
        CREATE_BREAKAWAY_FROM_JOB: int
