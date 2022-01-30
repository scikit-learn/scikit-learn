from typing import IO, Any, Callable, Iterator, Optional, Tuple, Type, Union

HIGHEST_PROTOCOL: int
bytes_types: Tuple[Type[Any], ...]  # undocumented

def dump(obj: Any, file: IO[bytes], protocol: int | None = ...) -> None: ...
def dumps(obj: Any, protocol: int | None = ...) -> bytes: ...
def load(file: IO[bytes]) -> Any: ...
def loads(string: bytes) -> Any: ...

class PickleError(Exception): ...
class PicklingError(PickleError): ...
class UnpicklingError(PickleError): ...

_reducedtype = Union[
    str,
    Tuple[Callable[..., Any], Tuple[Any, ...]],
    Tuple[Callable[..., Any], Tuple[Any, ...], Any],
    Tuple[Callable[..., Any], Tuple[Any, ...], Any, Optional[Iterator[Any]]],
    Tuple[Callable[..., Any], Tuple[Any, ...], Any, Optional[Iterator[Any]], Optional[Iterator[Any]]],
]

class Pickler:
    fast: bool
    def __init__(self, file: IO[bytes], protocol: int | None = ...) -> None: ...
    def dump(self, __obj: Any) -> None: ...
    def clear_memo(self) -> None: ...
    def persistent_id(self, obj: Any) -> Any: ...

class Unpickler:
    def __init__(self, file: IO[bytes]) -> None: ...
    def load(self) -> Any: ...
    def find_class(self, __module_name: str, __global_name: str) -> Any: ...

MARK: bytes
STOP: bytes
POP: bytes
POP_MARK: bytes
DUP: bytes
FLOAT: bytes
INT: bytes
BININT: bytes
BININT1: bytes
LONG: bytes
BININT2: bytes
NONE: bytes
PERSID: bytes
BINPERSID: bytes
REDUCE: bytes
STRING: bytes
BINSTRING: bytes
SHORT_BINSTRING: bytes
UNICODE: bytes
BINUNICODE: bytes
APPEND: bytes
BUILD: bytes
GLOBAL: bytes
DICT: bytes
EMPTY_DICT: bytes
APPENDS: bytes
GET: bytes
BINGET: bytes
INST: bytes
LONG_BINGET: bytes
LIST: bytes
EMPTY_LIST: bytes
OBJ: bytes
PUT: bytes
BINPUT: bytes
LONG_BINPUT: bytes
SETITEM: bytes
TUPLE: bytes
EMPTY_TUPLE: bytes
SETITEMS: bytes
BINFLOAT: bytes

TRUE: bytes
FALSE: bytes

# protocol 2
PROTO: bytes
NEWOBJ: bytes
EXT1: bytes
EXT2: bytes
EXT4: bytes
TUPLE1: bytes
TUPLE2: bytes
TUPLE3: bytes
NEWTRUE: bytes
NEWFALSE: bytes
LONG1: bytes
LONG4: bytes

def encode_long(x: int) -> bytes: ...  # undocumented
def decode_long(data: bytes) -> int: ...  # undocumented
