from typing import IO, Any, Callable, Iterator, List, MutableMapping, Text, Tuple, Type

_Reader = Callable[[IO[bytes]], Any]

UP_TO_NEWLINE: int
TAKEN_FROM_ARGUMENT1: int
TAKEN_FROM_ARGUMENT4: int

class ArgumentDescriptor(object):
    name: str
    n: int
    reader: _Reader
    doc: str
    def __init__(self, name: str, n: int, reader: _Reader, doc: str) -> None: ...

def read_uint1(f: IO[bytes]) -> int: ...

uint1: ArgumentDescriptor

def read_uint2(f: IO[bytes]) -> int: ...

uint2: ArgumentDescriptor

def read_int4(f: IO[bytes]) -> int: ...

int4: ArgumentDescriptor

def read_stringnl(f: IO[bytes], decode: bool = ..., stripquotes: bool = ...) -> bytes | Text: ...

stringnl: ArgumentDescriptor

def read_stringnl_noescape(f: IO[bytes]) -> str: ...

stringnl_noescape: ArgumentDescriptor

def read_stringnl_noescape_pair(f: IO[bytes]) -> Text: ...

stringnl_noescape_pair: ArgumentDescriptor

def read_string1(f: IO[bytes]) -> str: ...

string1: ArgumentDescriptor

def read_string4(f: IO[bytes]) -> str: ...

string4: ArgumentDescriptor

def read_unicodestringnl(f: IO[bytes]) -> Text: ...

unicodestringnl: ArgumentDescriptor

def read_unicodestring4(f: IO[bytes]) -> Text: ...

unicodestring4: ArgumentDescriptor

def read_decimalnl_short(f: IO[bytes]) -> int: ...
def read_decimalnl_long(f: IO[bytes]) -> int: ...

decimalnl_short: ArgumentDescriptor
decimalnl_long: ArgumentDescriptor

def read_floatnl(f: IO[bytes]) -> float: ...

floatnl: ArgumentDescriptor

def read_float8(f: IO[bytes]) -> float: ...

float8: ArgumentDescriptor

def read_long1(f: IO[bytes]) -> int: ...

long1: ArgumentDescriptor

def read_long4(f: IO[bytes]) -> int: ...

long4: ArgumentDescriptor

class StackObject(object):
    name: str
    obtype: Type[Any] | Tuple[Type[Any], ...]
    doc: str
    def __init__(self, name: str, obtype: Type[Any] | Tuple[Type[Any], ...], doc: str) -> None: ...

pyint: StackObject
pylong: StackObject
pyinteger_or_bool: StackObject
pybool: StackObject
pyfloat: StackObject
pystring: StackObject
pyunicode: StackObject
pynone: StackObject
pytuple: StackObject
pylist: StackObject
pydict: StackObject
anyobject: StackObject
markobject: StackObject
stackslice: StackObject

class OpcodeInfo(object):
    name: str
    code: str
    arg: ArgumentDescriptor | None
    stack_before: List[StackObject]
    stack_after: List[StackObject]
    proto: int
    doc: str
    def __init__(
        self,
        name: str,
        code: str,
        arg: ArgumentDescriptor | None,
        stack_before: List[StackObject],
        stack_after: List[StackObject],
        proto: int,
        doc: str,
    ) -> None: ...

opcodes: List[OpcodeInfo]

def genops(pickle: bytes | IO[bytes]) -> Iterator[Tuple[OpcodeInfo, Any | None, int | None]]: ...
def optimize(p: bytes | IO[bytes]) -> bytes: ...
def dis(
    pickle: bytes | IO[bytes], out: IO[str] | None = ..., memo: MutableMapping[int, Any] | None = ..., indentlevel: int = ...
) -> None: ...
