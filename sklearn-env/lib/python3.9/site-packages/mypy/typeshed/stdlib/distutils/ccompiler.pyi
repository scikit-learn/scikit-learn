from typing import Any, Callable, Optional, Tuple, Union

_Macro = Union[Tuple[str], Tuple[str, Optional[str]]]

def gen_lib_options(
    compiler: CCompiler, library_dirs: list[str], runtime_library_dirs: list[str], libraries: list[str]
) -> list[str]: ...
def gen_preprocess_options(macros: list[_Macro], include_dirs: list[str]) -> list[str]: ...
def get_default_compiler(osname: str | None = ..., platform: str | None = ...) -> str: ...
def new_compiler(
    plat: str | None = ..., compiler: str | None = ..., verbose: int = ..., dry_run: int = ..., force: int = ...
) -> CCompiler: ...
def show_compilers() -> None: ...

class CCompiler:
    dry_run: bool
    force: bool
    verbose: bool
    output_dir: str | None
    macros: list[_Macro]
    include_dirs: list[str]
    libraries: list[str]
    library_dirs: list[str]
    runtime_library_dirs: list[str]
    objects: list[str]
    def __init__(self, verbose: int = ..., dry_run: int = ..., force: int = ...) -> None: ...
    def add_include_dir(self, dir: str) -> None: ...
    def set_include_dirs(self, dirs: list[str]) -> None: ...
    def add_library(self, libname: str) -> None: ...
    def set_libraries(self, libnames: list[str]) -> None: ...
    def add_library_dir(self, dir: str) -> None: ...
    def set_library_dirs(self, dirs: list[str]) -> None: ...
    def add_runtime_library_dir(self, dir: str) -> None: ...
    def set_runtime_library_dirs(self, dirs: list[str]) -> None: ...
    def define_macro(self, name: str, value: str | None = ...) -> None: ...
    def undefine_macro(self, name: str) -> None: ...
    def add_link_object(self, object: str) -> None: ...
    def set_link_objects(self, objects: list[str]) -> None: ...
    def detect_language(self, sources: str | list[str]) -> str | None: ...
    def find_library_file(self, dirs: list[str], lib: str, debug: bool = ...) -> str | None: ...
    def has_function(
        self,
        funcname: str,
        includes: list[str] | None = ...,
        include_dirs: list[str] | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
    ) -> bool: ...
    def library_dir_option(self, dir: str) -> str: ...
    def library_option(self, lib: str) -> str: ...
    def runtime_library_dir_option(self, dir: str) -> str: ...
    def set_executables(self, **args: str) -> None: ...
    def compile(
        self,
        sources: list[str],
        output_dir: str | None = ...,
        macros: _Macro | None = ...,
        include_dirs: list[str] | None = ...,
        debug: bool = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        depends: list[str] | None = ...,
    ) -> list[str]: ...
    def create_static_lib(
        self,
        objects: list[str],
        output_libname: str,
        output_dir: str | None = ...,
        debug: bool = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link(
        self,
        target_desc: str,
        objects: list[str],
        output_filename: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_executable(
        self,
        objects: list[str],
        output_progname: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        debug: bool = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_shared_lib(
        self,
        objects: list[str],
        output_libname: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def link_shared_object(
        self,
        objects: list[str],
        output_filename: str,
        output_dir: str | None = ...,
        libraries: list[str] | None = ...,
        library_dirs: list[str] | None = ...,
        runtime_library_dirs: list[str] | None = ...,
        export_symbols: list[str] | None = ...,
        debug: bool = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
        build_temp: str | None = ...,
        target_lang: str | None = ...,
    ) -> None: ...
    def preprocess(
        self,
        source: str,
        output_file: str | None = ...,
        macros: list[_Macro] | None = ...,
        include_dirs: list[str] | None = ...,
        extra_preargs: list[str] | None = ...,
        extra_postargs: list[str] | None = ...,
    ) -> None: ...
    def executable_filename(self, basename: str, strip_dir: int = ..., output_dir: str = ...) -> str: ...
    def library_filename(self, libname: str, lib_type: str = ..., strip_dir: int = ..., output_dir: str = ...) -> str: ...
    def object_filenames(self, source_filenames: list[str], strip_dir: int = ..., output_dir: str = ...) -> list[str]: ...
    def shared_object_filename(self, basename: str, strip_dir: int = ..., output_dir: str = ...) -> str: ...
    def execute(self, func: Callable[..., None], args: Tuple[Any, ...], msg: str | None = ..., level: int = ...) -> None: ...
    def spawn(self, cmd: list[str]) -> None: ...
    def mkpath(self, name: str, mode: int = ...) -> None: ...
    def move_file(self, src: str, dst: str) -> str: ...
    def announce(self, msg: str, level: int = ...) -> None: ...
    def warn(self, msg: str) -> None: ...
    def debug_print(self, msg: str) -> None: ...
