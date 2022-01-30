from collections.abc import Sequence
from typing import Any, Pattern

from ..ccompiler import CCompiler
from ..cmd import Command

LANG_EXT: dict[str, str]

class config(Command):
    description: str
    # Tuple is full name, short name, description
    user_options: Sequence[tuple[str, str | None, str]]
    compiler: str | CCompiler
    cc: str | None
    include_dirs: Sequence[str] | None
    libraries: Sequence[str] | None
    library_dirs: Sequence[str] | None
    noisy: int
    dump_source: int
    temp_files: Sequence[str]
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
    def try_cpp(
        self,
        body: str | None = ...,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        lang: str = ...,
    ) -> bool: ...
    def search_cpp(
        self,
        pattern: Pattern[str] | str,
        body: str | None = ...,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        lang: str = ...,
    ) -> bool: ...
    def try_compile(
        self, body: str, headers: Sequence[str] | None = ..., include_dirs: Sequence[str] | None = ..., lang: str = ...
    ) -> bool: ...
    def try_link(
        self,
        body: str,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        libraries: Sequence[str] | None = ...,
        library_dirs: Sequence[str] | None = ...,
        lang: str = ...,
    ) -> bool: ...
    def try_run(
        self,
        body: str,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        libraries: Sequence[str] | None = ...,
        library_dirs: Sequence[str] | None = ...,
        lang: str = ...,
    ) -> bool: ...
    def check_func(
        self,
        func: str,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        libraries: Sequence[str] | None = ...,
        library_dirs: Sequence[str] | None = ...,
        decl: int = ...,
        call: int = ...,
    ) -> bool: ...
    def check_lib(
        self,
        library: str,
        library_dirs: Sequence[str] | None = ...,
        headers: Sequence[str] | None = ...,
        include_dirs: Sequence[str] | None = ...,
        other_libraries: list[str] = ...,
    ) -> bool: ...
    def check_header(
        self, header: str, include_dirs: Sequence[str] | None = ..., library_dirs: Sequence[str] | None = ..., lang: str = ...
    ) -> bool: ...

def dump_file(filename: str, head: Any | None = ...) -> None: ...
