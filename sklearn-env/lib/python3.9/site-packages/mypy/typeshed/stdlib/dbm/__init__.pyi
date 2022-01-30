from _typeshed import Self
from types import TracebackType
from typing import Iterator, MutableMapping, Tuple, Type, Union
from typing_extensions import Literal

_KeyType = Union[str, bytes]
_ValueType = Union[str, bytes]
_TFlags = Literal[
    "r",
    "w",
    "c",
    "n",
    "rf",
    "wf",
    "cf",
    "nf",
    "rs",
    "ws",
    "cs",
    "ns",
    "ru",
    "wu",
    "cu",
    "nu",
    "rfs",
    "wfs",
    "cfs",
    "nfs",
    "rfu",
    "wfu",
    "cfu",
    "nfu",
    "rsf",
    "wsf",
    "csf",
    "nsf",
    "rsu",
    "wsu",
    "csu",
    "nsu",
    "ruf",
    "wuf",
    "cuf",
    "nuf",
    "rus",
    "wus",
    "cus",
    "nus",
    "rfsu",
    "wfsu",
    "cfsu",
    "nfsu",
    "rfus",
    "wfus",
    "cfus",
    "nfus",
    "rsfu",
    "wsfu",
    "csfu",
    "nsfu",
    "rsuf",
    "wsuf",
    "csuf",
    "nsuf",
    "rufs",
    "wufs",
    "cufs",
    "nufs",
    "rusf",
    "wusf",
    "cusf",
    "nusf",
]

class _Database(MutableMapping[_KeyType, bytes]):
    def close(self) -> None: ...
    def __getitem__(self, key: _KeyType) -> bytes: ...
    def __setitem__(self, key: _KeyType, value: _ValueType) -> None: ...
    def __delitem__(self, key: _KeyType) -> None: ...
    def __iter__(self) -> Iterator[bytes]: ...
    def __len__(self) -> int: ...
    def __del__(self) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(
        self, exc_type: Type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...

class _error(Exception): ...

error = Tuple[Type[_error], Type[OSError]]

def whichdb(filename: str) -> str: ...
def open(file: str, flag: _TFlags = ..., mode: int = ...) -> _Database: ...
