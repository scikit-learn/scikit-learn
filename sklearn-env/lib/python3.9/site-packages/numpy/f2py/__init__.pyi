import os
import subprocess
from typing import Literal as L, Any, List, Iterable, Dict, overload, TypedDict

from numpy._pytesttester import PytestTester

class _F2PyDictBase(TypedDict):
    csrc: List[str]
    h: List[str]

class _F2PyDict(_F2PyDictBase, total=False):
    fsrc: List[str]
    ltx: List[str]

__all__: List[str]
__path__: List[str]
test: PytestTester

def run_main(comline_list: Iterable[str]) -> Dict[str, _F2PyDict]: ...

@overload
def compile(  # type: ignore[misc]
    source: str | bytes,
    modulename: str = ...,
    extra_args: str | List[str] = ...,
    verbose: bool = ...,
    source_fn: None | str | bytes | os.PathLike[Any] = ...,
    extension: L[".f", ".f90"] = ...,
    full_output: L[False] = ...,
) -> int: ...
@overload
def compile(
    source: str | bytes,
    modulename: str = ...,
    extra_args: str | List[str] = ...,
    verbose: bool = ...,
    source_fn: None | str | bytes | os.PathLike[Any] = ...,
    extension: L[".f", ".f90"] = ...,
    full_output: L[True] = ...,
) -> subprocess.CompletedProcess[bytes]: ...

def get_include() -> str: ...
