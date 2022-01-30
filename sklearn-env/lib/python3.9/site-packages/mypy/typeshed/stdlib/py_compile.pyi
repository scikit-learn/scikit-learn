import sys
from typing import AnyStr, Type

class PyCompileError(Exception):
    exc_type_name: str
    exc_value: BaseException
    file: str
    msg: str
    def __init__(self, exc_type: Type[BaseException], exc_value: BaseException, file: str, msg: str = ...) -> None: ...

if sys.version_info >= (3, 7):
    import enum
    class PycInvalidationMode(enum.Enum):
        TIMESTAMP: int
        CHECKED_HASH: int
        UNCHECKED_HASH: int
    def _get_default_invalidation_mode() -> PycInvalidationMode: ...

if sys.version_info >= (3, 8):
    def compile(
        file: AnyStr,
        cfile: AnyStr | None = ...,
        dfile: AnyStr | None = ...,
        doraise: bool = ...,
        optimize: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
        quiet: int = ...,
    ) -> AnyStr | None: ...

elif sys.version_info >= (3, 7):
    def compile(
        file: AnyStr,
        cfile: AnyStr | None = ...,
        dfile: AnyStr | None = ...,
        doraise: bool = ...,
        optimize: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
    ) -> AnyStr | None: ...

else:
    def compile(
        file: AnyStr, cfile: AnyStr | None = ..., dfile: AnyStr | None = ..., doraise: bool = ..., optimize: int = ...
    ) -> AnyStr | None: ...

def main(args: list[str] | None = ...) -> int: ...
