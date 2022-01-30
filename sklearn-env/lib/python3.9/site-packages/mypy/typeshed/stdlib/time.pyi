import sys
from types import SimpleNamespace
from typing import Any, NamedTuple, Tuple
from typing_extensions import final

_TimeTuple = Tuple[int, int, int, int, int, int, int, int, int]

altzone: int
daylight: int
timezone: int
tzname: tuple[str, str]

if sys.version_info >= (3, 7):
    if sys.platform == "linux":
        CLOCK_BOOTTIME: int
    if sys.platform != "linux" and sys.platform != "win32" and sys.platform != "darwin":
        CLOCK_PROF: int  # FreeBSD, NetBSD, OpenBSD
        CLOCK_UPTIME: int  # FreeBSD, OpenBSD

if sys.platform != "win32":
    CLOCK_MONOTONIC: int
    CLOCK_MONOTONIC_RAW: int
    CLOCK_PROCESS_CPUTIME_ID: int
    CLOCK_REALTIME: int
    CLOCK_THREAD_CPUTIME_ID: int
    if sys.platform != "linux" and sys.platform != "darwin":
        CLOCK_HIGHRES: int  # Solaris only

if sys.version_info >= (3, 8) and sys.platform == "darwin":
    CLOCK_UPTIME_RAW: int

if sys.version_info >= (3, 9) and sys.platform == "linux":
    CLOCK_TAI: int

class _struct_time(NamedTuple):
    tm_year: int
    tm_mon: int
    tm_mday: int
    tm_hour: int
    tm_min: int
    tm_sec: int
    tm_wday: int
    tm_yday: int
    tm_isdst: int
    @property
    def n_fields(self) -> int: ...
    @property
    def n_sequence_fields(self) -> int: ...
    @property
    def n_unnamed_fields(self) -> int: ...

@final
class struct_time(_struct_time):
    def __init__(
        self,
        o: tuple[int, int, int, int, int, int, int, int, int]
        | tuple[int, int, int, int, int, int, int, int, int, str]
        | tuple[int, int, int, int, int, int, int, int, int, str, int],
        _arg: Any = ...,
    ) -> None: ...
    def __new__(
        cls,
        o: tuple[int, int, int, int, int, int, int, int, int]
        | tuple[int, int, int, int, int, int, int, int, int, str]
        | tuple[int, int, int, int, int, int, int, int, int, str, int],
        _arg: Any = ...,
    ) -> struct_time: ...
    @property
    def tm_zone(self) -> str: ...
    @property
    def tm_gmtoff(self) -> int: ...

def asctime(t: _TimeTuple | struct_time = ...) -> str: ...

if sys.version_info < (3, 8):
    def clock() -> float: ...

def ctime(secs: float | None = ...) -> str: ...
def gmtime(secs: float | None = ...) -> struct_time: ...
def localtime(secs: float | None = ...) -> struct_time: ...
def mktime(t: _TimeTuple | struct_time) -> float: ...
def sleep(secs: float) -> None: ...
def strftime(format: str, t: _TimeTuple | struct_time = ...) -> str: ...
def strptime(string: str, format: str = ...) -> struct_time: ...
def time() -> float: ...

if sys.platform != "win32":
    def tzset() -> None: ...  # Unix only

def get_clock_info(name: str) -> SimpleNamespace: ...
def monotonic() -> float: ...
def perf_counter() -> float: ...
def process_time() -> float: ...

if sys.platform != "win32":
    def clock_getres(clk_id: int) -> float: ...  # Unix only
    def clock_gettime(clk_id: int) -> float: ...  # Unix only
    def clock_settime(clk_id: int, time: float) -> None: ...  # Unix only

if sys.version_info >= (3, 7):
    if sys.platform != "win32":
        def clock_gettime_ns(clock_id: int) -> int: ...
        def clock_settime_ns(clock_id: int, time: int) -> int: ...
    def monotonic_ns() -> int: ...
    def perf_counter_ns() -> int: ...
    def process_time_ns() -> int: ...
    def time_ns() -> int: ...
    def thread_time() -> float: ...
    def thread_time_ns() -> int: ...
