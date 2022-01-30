import sys
from types import CodeType
from typing import IO, Any, Container, Iterable, Iterator, Sequence, Tuple

LOAD_CONST: int  # undocumented
IMPORT_NAME: int  # undocumented
STORE_NAME: int  # undocumented
STORE_GLOBAL: int  # undocumented
STORE_OPS: tuple[int, int]  # undocumented
EXTENDED_ARG: int  # undocumented

packagePathMap: dict[str, list[str]]  # undocumented

def AddPackagePath(packagename: str, path: str) -> None: ...

replacePackageMap: dict[str, str]  # undocumented

def ReplacePackage(oldname: str, newname: str) -> None: ...

class Module:  # undocumented
    def __init__(self, name: str, file: str | None = ..., path: str | None = ...) -> None: ...
    def __repr__(self) -> str: ...

class ModuleFinder:

    modules: dict[str, Module]
    path: list[str]  # undocumented
    badmodules: dict[str, dict[str, int]]  # undocumented
    debug: int  # undocumented
    indent: int  # undocumented
    excludes: Container[str]  # undocumented
    replace_paths: Sequence[tuple[str, str]]  # undocumented

    if sys.version_info >= (3, 8):
        def __init__(
            self,
            path: list[str] | None = ...,
            debug: int = ...,
            excludes: Container[str] | None = ...,
            replace_paths: Sequence[tuple[str, str]] | None = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            path: list[str] | None = ...,
            debug: int = ...,
            excludes: Container[str] = ...,
            replace_paths: Sequence[tuple[str, str]] = ...,
        ) -> None: ...
    def msg(self, level: int, str: str, *args: Any) -> None: ...  # undocumented
    def msgin(self, *args: Any) -> None: ...  # undocumented
    def msgout(self, *args: Any) -> None: ...  # undocumented
    def run_script(self, pathname: str) -> None: ...
    def load_file(self, pathname: str) -> None: ...  # undocumented
    def import_hook(
        self, name: str, caller: Module | None = ..., fromlist: list[str] | None = ..., level: int = ...
    ) -> Module | None: ...  # undocumented
    def determine_parent(self, caller: Module | None, level: int = ...) -> Module | None: ...  # undocumented
    def find_head_package(self, parent: Module, name: str) -> tuple[Module, str]: ...  # undocumented
    def load_tail(self, q: Module, tail: str) -> Module: ...  # undocumented
    def ensure_fromlist(self, m: Module, fromlist: Iterable[str], recursive: int = ...) -> None: ...  # undocumented
    def find_all_submodules(self, m: Module) -> Iterable[str]: ...  # undocumented
    def import_module(self, partname: str, fqname: str, parent: Module) -> Module | None: ...  # undocumented
    def load_module(self, fqname: str, fp: IO[str], pathname: str, file_info: tuple[str, str, str]) -> Module: ...  # undocumented
    def scan_opcodes(self, co: CodeType) -> Iterator[tuple[str, Tuple[Any, ...]]]: ...  # undocumented
    def scan_code(self, co: CodeType, m: Module) -> None: ...  # undocumented
    def load_package(self, fqname: str, pathname: str) -> Module: ...  # undocumented
    def add_module(self, fqname: str) -> Module: ...  # undocumented
    def find_module(
        self, name: str, path: str | None, parent: Module | None = ...
    ) -> tuple[IO[Any] | None, str | None, tuple[str, str, int]]: ...  # undocumented
    def report(self) -> None: ...
    def any_missing(self) -> list[str]: ...  # undocumented
    def any_missing_maybe(self) -> tuple[list[str], list[str]]: ...  # undocumented
    def replace_paths_in_code(self, co: CodeType) -> CodeType: ...  # undocumented
