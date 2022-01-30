from _csv import (
    QUOTE_ALL as QUOTE_ALL,
    QUOTE_MINIMAL as QUOTE_MINIMAL,
    QUOTE_NONE as QUOTE_NONE,
    QUOTE_NONNUMERIC as QUOTE_NONNUMERIC,
    Dialect as Dialect,
    Error as Error,
    _DialectLike,
    _reader,
    _writer,
    field_size_limit as field_size_limit,
    get_dialect as get_dialect,
    list_dialects as list_dialects,
    reader as reader,
    register_dialect as register_dialect,
    unregister_dialect as unregister_dialect,
    writer as writer,
)
from typing import (
    Any,
    Dict as _DictReadMapping,
    Generic,
    Iterable,
    Iterator,
    List,
    Mapping,
    Sequence,
    Text,
    Type,
    TypeVar,
    overload,
)

_T = TypeVar("_T")

class excel(Dialect):
    delimiter: str
    quotechar: str
    doublequote: bool
    skipinitialspace: bool
    lineterminator: str
    quoting: int

class excel_tab(excel):
    delimiter: str

class DictReader(Generic[_T], Iterator[_DictReadMapping[_T, str]]):
    fieldnames: Sequence[_T] | None
    restkey: str | None
    restval: str | None
    reader: _reader
    dialect: _DialectLike
    line_num: int
    @overload
    def __init__(
        self,
        f: Iterable[Text],
        fieldnames: Sequence[_T],
        restkey: str | None = ...,
        restval: str | None = ...,
        dialect: _DialectLike = ...,
        *args: Any,
        **kwds: Any,
    ) -> None: ...
    @overload
    def __init__(
        self: DictReader[str],
        f: Iterable[Text],
        fieldnames: Sequence[str] | None = ...,
        restkey: str | None = ...,
        restval: str | None = ...,
        dialect: _DialectLike = ...,
        *args: Any,
        **kwds: Any,
    ) -> None: ...
    def __iter__(self) -> DictReader[_T]: ...
    def next(self) -> _DictReadMapping[_T, str]: ...

class DictWriter(Generic[_T]):
    fieldnames: Sequence[_T]
    restval: Any | None
    extrasaction: str
    writer: _writer
    def __init__(
        self,
        f: Any,
        fieldnames: Sequence[_T],
        restval: Any | None = ...,
        extrasaction: str = ...,
        dialect: _DialectLike = ...,
        *args: Any,
        **kwds: Any,
    ) -> None: ...
    def writeheader(self) -> None: ...
    def writerow(self, rowdict: Mapping[_T, Any]) -> Any: ...
    def writerows(self, rowdicts: Iterable[Mapping[_T, Any]]) -> None: ...

class Sniffer(object):
    preferred: List[str]
    def __init__(self) -> None: ...
    def sniff(self, sample: str, delimiters: str | None = ...) -> Type[Dialect]: ...
    def has_header(self, sample: str) -> bool: ...
