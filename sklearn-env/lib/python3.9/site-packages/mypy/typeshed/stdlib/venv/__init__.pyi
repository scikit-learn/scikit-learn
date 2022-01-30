import sys
from _typeshed import StrOrBytesPath
from types import SimpleNamespace
from typing import Sequence

class EnvBuilder:
    system_site_packages: bool
    clear: bool
    symlinks: bool
    upgrade: bool
    with_pip: bool
    prompt: str | None

    if sys.version_info >= (3, 9):
        def __init__(
            self,
            system_site_packages: bool = ...,
            clear: bool = ...,
            symlinks: bool = ...,
            upgrade: bool = ...,
            with_pip: bool = ...,
            prompt: str | None = ...,
            upgrade_deps: bool = ...,
        ) -> None: ...
    else:
        def __init__(
            self,
            system_site_packages: bool = ...,
            clear: bool = ...,
            symlinks: bool = ...,
            upgrade: bool = ...,
            with_pip: bool = ...,
            prompt: str | None = ...,
        ) -> None: ...
    def create(self, env_dir: StrOrBytesPath) -> None: ...
    def clear_directory(self, path: StrOrBytesPath) -> None: ...  # undocumented
    def ensure_directories(self, env_dir: StrOrBytesPath) -> SimpleNamespace: ...
    def create_configuration(self, context: SimpleNamespace) -> None: ...
    def symlink_or_copy(
        self, src: StrOrBytesPath, dst: StrOrBytesPath, relative_symlinks_ok: bool = ...
    ) -> None: ...  # undocumented
    def setup_python(self, context: SimpleNamespace) -> None: ...
    def _setup_pip(self, context: SimpleNamespace) -> None: ...  # undocumented
    def setup_scripts(self, context: SimpleNamespace) -> None: ...
    def post_setup(self, context: SimpleNamespace) -> None: ...
    def replace_variables(self, text: str, context: SimpleNamespace) -> str: ...  # undocumented
    def install_scripts(self, context: SimpleNamespace, path: str) -> None: ...
    if sys.version_info >= (3, 9):
        def upgrade_dependencies(self, context: SimpleNamespace) -> None: ...

if sys.version_info >= (3, 9):
    def create(
        env_dir: StrOrBytesPath,
        system_site_packages: bool = ...,
        clear: bool = ...,
        symlinks: bool = ...,
        with_pip: bool = ...,
        prompt: str | None = ...,
        upgrade_deps: bool = ...,
    ) -> None: ...

else:
    def create(
        env_dir: StrOrBytesPath,
        system_site_packages: bool = ...,
        clear: bool = ...,
        symlinks: bool = ...,
        with_pip: bool = ...,
        prompt: str | None = ...,
    ) -> None: ...

def main(args: Sequence[str] | None = ...) -> None: ...
