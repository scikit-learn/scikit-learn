from sre_constants import (
    SRE_FLAG_ASCII as SRE_FLAG_ASCII,
    SRE_FLAG_DEBUG as SRE_FLAG_DEBUG,
    SRE_FLAG_DOTALL as SRE_FLAG_DOTALL,
    SRE_FLAG_IGNORECASE as SRE_FLAG_IGNORECASE,
    SRE_FLAG_LOCALE as SRE_FLAG_LOCALE,
    SRE_FLAG_MULTILINE as SRE_FLAG_MULTILINE,
    SRE_FLAG_TEMPLATE as SRE_FLAG_TEMPLATE,
    SRE_FLAG_UNICODE as SRE_FLAG_UNICODE,
    SRE_FLAG_VERBOSE as SRE_FLAG_VERBOSE,
    SRE_INFO_CHARSET as SRE_INFO_CHARSET,
    SRE_INFO_LITERAL as SRE_INFO_LITERAL,
    SRE_INFO_PREFIX as SRE_INFO_PREFIX,
    _NamedIntConstant,
)
from sre_parse import SubPattern
from typing import Any, Pattern

MAXCODE: int

def dis(code: list[_NamedIntConstant]) -> None: ...
def isstring(obj: Any) -> bool: ...
def compile(p: str | bytes | SubPattern, flags: int = ...) -> Pattern[Any]: ...
