import pyexpat.errors as errors
import pyexpat.model as model
from _typeshed import SupportsRead
from typing import Any, Callable, Optional, Tuple
from typing_extensions import final

EXPAT_VERSION: str  # undocumented
version_info: tuple[int, int, int]  # undocumented
native_encoding: str  # undocumented
features: list[tuple[str, int]]  # undocumented

class ExpatError(Exception):
    code: int
    lineno: int
    offset: int

error = ExpatError

XML_PARAM_ENTITY_PARSING_NEVER: int
XML_PARAM_ENTITY_PARSING_UNLESS_STANDALONE: int
XML_PARAM_ENTITY_PARSING_ALWAYS: int

_Model = Tuple[int, int, Optional[str], Tuple[Any, ...]]

@final
class XMLParserType(object):
    def Parse(self, __data: str | bytes, __isfinal: bool = ...) -> int: ...
    def ParseFile(self, __file: SupportsRead[bytes]) -> int: ...
    def SetBase(self, __base: str) -> None: ...
    def GetBase(self) -> str | None: ...
    def GetInputContext(self) -> bytes | None: ...
    def ExternalEntityParserCreate(self, __context: str | None, __encoding: str = ...) -> XMLParserType: ...
    def SetParamEntityParsing(self, __flag: int) -> int: ...
    def UseForeignDTD(self, __flag: bool = ...) -> None: ...
    buffer_size: int
    buffer_text: bool
    buffer_used: int
    namespace_prefixes: bool  # undocumented
    ordered_attributes: bool
    specified_attributes: bool
    ErrorByteIndex: int
    ErrorCode: int
    ErrorColumnNumber: int
    ErrorLineNumber: int
    CurrentByteIndex: int
    CurrentColumnNumber: int
    CurrentLineNumber: int
    XmlDeclHandler: Callable[[str, str | None, int], Any] | None
    StartDoctypeDeclHandler: Callable[[str, str | None, str | None, bool], Any] | None
    EndDoctypeDeclHandler: Callable[[], Any] | None
    ElementDeclHandler: Callable[[str, _Model], Any] | None
    AttlistDeclHandler: Callable[[str, str, str, str | None, bool], Any] | None
    StartElementHandler: Callable[[str, dict[str, str]], Any] | Callable[[str, list[str]], Any] | Callable[
        [str, dict[str, str], list[str]], Any
    ] | None
    EndElementHandler: Callable[[str], Any] | None
    ProcessingInstructionHandler: Callable[[str, str], Any] | None
    CharacterDataHandler: Callable[[str], Any] | None
    UnparsedEntityDeclHandler: Callable[[str, str | None, str, str | None, str], Any] | None
    EntityDeclHandler: Callable[[str, bool, str | None, str | None, str, str | None, str | None], Any] | None
    NotationDeclHandler: Callable[[str, str | None, str, str | None], Any] | None
    StartNamespaceDeclHandler: Callable[[str, str], Any] | None
    EndNamespaceDeclHandler: Callable[[str], Any] | None
    CommentHandler: Callable[[str], Any] | None
    StartCdataSectionHandler: Callable[[], Any] | None
    EndCdataSectionHandler: Callable[[], Any] | None
    DefaultHandler: Callable[[str], Any] | None
    DefaultHandlerExpand: Callable[[str], Any] | None
    NotStandaloneHandler: Callable[[], int] | None
    ExternalEntityRefHandler: Callable[[str, str | None, str | None, str | None], int] | None

def ErrorString(__code: int) -> str: ...

# intern is undocumented
def ParserCreate(
    encoding: str | None = ..., namespace_separator: str | None = ..., intern: dict[str, Any] | None = ...
) -> XMLParserType: ...
