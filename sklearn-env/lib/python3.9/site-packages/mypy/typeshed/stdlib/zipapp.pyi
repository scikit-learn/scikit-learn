import sys
from pathlib import Path
from typing import BinaryIO, Callable, Union

_Path = Union[str, Path, BinaryIO]

class ZipAppError(ValueError): ...

if sys.version_info >= (3, 7):
    def create_archive(
        source: _Path,
        target: _Path | None = ...,
        interpreter: str | None = ...,
        main: str | None = ...,
        filter: Callable[[Path], bool] | None = ...,
        compressed: bool = ...,
    ) -> None: ...

else:
    def create_archive(
        source: _Path, target: _Path | None = ..., interpreter: str | None = ..., main: str | None = ...
    ) -> None: ...

def get_interpreter(archive: _Path) -> str: ...
