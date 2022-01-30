from distutils.cmd import Command
from typing import Text

class install(Command):
    user: bool
    prefix: Text | None
    home: Text | None
    root: Text | None
    install_lib: Text | None
    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    def run(self) -> None: ...
