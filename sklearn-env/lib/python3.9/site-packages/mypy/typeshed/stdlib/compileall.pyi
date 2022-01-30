import sys
from _typeshed import StrPath
from typing import Any, Protocol

if sys.version_info >= (3, 7):
    from py_compile import PycInvalidationMode

class _SupportsSearch(Protocol):
    def search(self, string: str) -> Any: ...

if sys.version_info >= (3, 9):
    def compile_dir(
        dir: StrPath,
        maxlevels: int | None = ...,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
        *,
        stripdir: str | None = ...,  # TODO: change to StrPath | None once https://bugs.python.org/issue40447 is resolved
        prependdir: StrPath | None = ...,
        limit_sl_dest: StrPath | None = ...,
        hardlink_dupes: bool = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
        *,
        stripdir: str | None = ...,  # TODO: change to StrPath | None once https://bugs.python.org/issue40447 is resolved
        prependdir: StrPath | None = ...,
        limit_sl_dest: StrPath | None = ...,
        hardlink_dupes: bool = ...,
    ) -> int: ...

elif sys.version_info >= (3, 7):
    def compile_dir(
        dir: StrPath,
        maxlevels: int = ...,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
    ) -> int: ...

else:
    def compile_dir(
        dir: StrPath,
        maxlevels: int = ...,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: StrPath | None = ...,
        force: bool = ...,
        rx: _SupportsSearch | None = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
    ) -> int: ...

if sys.version_info >= (3, 7):
    def compile_path(
        skip_curdir: bool = ...,
        maxlevels: int = ...,
        force: bool = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: PycInvalidationMode | None = ...,
    ) -> int: ...

else:
    def compile_path(
        skip_curdir: bool = ...,
        maxlevels: int = ...,
        force: bool = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
    ) -> int: ...
