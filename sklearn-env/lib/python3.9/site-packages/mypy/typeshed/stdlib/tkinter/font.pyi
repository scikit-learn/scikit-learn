import _tkinter
import sys
import tkinter
from typing import Any, List, Tuple, Union, overload
from typing_extensions import Literal, TypedDict

NORMAL: Literal["normal"]
ROMAN: Literal["roman"]
BOLD: Literal["bold"]
ITALIC: Literal["italic"]

_FontDescription = Union[
    # "Helvetica 12"
    str,
    # A font object constructed in Python
    Font,
    # ("Helvetica", 12, BOLD)
    List[Any],
    Tuple[Any, ...],
    # A font object constructed in Tcl
    _tkinter.Tcl_Obj,
]

class _FontDict(TypedDict):
    family: str
    size: int
    weight: Literal["normal", "bold"]
    slant: Literal["roman", "italic"]
    underline: bool
    overstrike: bool

class _MetricsDict(TypedDict):
    ascent: int
    descent: int
    linespace: int
    fixed: bool

class Font:
    name: str
    delete_font: bool
    def __init__(
        self,
        # In tkinter, 'root' refers to tkinter.Tk by convention, but the code
        # actually works with any tkinter widget so we use tkinter.Misc.
        root: tkinter.Misc | None = ...,
        font: _FontDescription | None = ...,
        name: str | None = ...,
        exists: bool = ...,
        *,
        family: str = ...,
        size: int = ...,
        weight: Literal["normal", "bold"] = ...,
        slant: Literal["roman", "italic"] = ...,
        underline: bool = ...,
        overstrike: bool = ...,
    ) -> None: ...
    def __setitem__(self, key: str, value: Any) -> None: ...
    @overload
    def cget(self, option: Literal["family"]) -> str: ...
    @overload
    def cget(self, option: Literal["size"]) -> int: ...
    @overload
    def cget(self, option: Literal["weight"]) -> Literal["normal", "bold"]: ...
    @overload
    def cget(self, option: Literal["slant"]) -> Literal["roman", "italic"]: ...
    @overload
    def cget(self, option: Literal["underline", "overstrike"]) -> bool: ...
    @overload
    def cget(self, option: str) -> Any: ...
    __getitem__ = cget
    @overload
    def actual(self, option: Literal["family"], displayof: tkinter.Misc | None = ...) -> str: ...
    @overload
    def actual(self, option: Literal["size"], displayof: tkinter.Misc | None = ...) -> int: ...
    @overload
    def actual(self, option: Literal["weight"], displayof: tkinter.Misc | None = ...) -> Literal["normal", "bold"]: ...
    @overload
    def actual(self, option: Literal["slant"], displayof: tkinter.Misc | None = ...) -> Literal["roman", "italic"]: ...
    @overload
    def actual(self, option: Literal["underline", "overstrike"], displayof: tkinter.Misc | None = ...) -> bool: ...
    @overload
    def actual(self, option: None, displayof: tkinter.Misc | None = ...) -> _FontDict: ...
    @overload
    def actual(self, *, displayof: tkinter.Misc | None = ...) -> _FontDict: ...
    def config(
        self,
        *,
        family: str = ...,
        size: int = ...,
        weight: Literal["normal", "bold"] = ...,
        slant: Literal["roman", "italic"] = ...,
        underline: bool = ...,
        overstrike: bool = ...,
    ) -> _FontDict | None: ...
    configure = config
    def copy(self) -> Font: ...
    @overload
    def metrics(self, __option: Literal["ascent", "descent", "linespace"], *, displayof: tkinter.Misc | None = ...) -> int: ...
    @overload
    def metrics(self, __option: Literal["fixed"], *, displayof: tkinter.Misc | None = ...) -> bool: ...
    @overload
    def metrics(self, *, displayof: tkinter.Misc | None = ...) -> _MetricsDict: ...
    def measure(self, text: str, displayof: tkinter.Misc | None = ...) -> int: ...

def families(root: tkinter.Misc | None = ..., displayof: tkinter.Misc | None = ...) -> Tuple[str, ...]: ...
def names(root: tkinter.Misc | None = ...) -> Tuple[str, ...]: ...

if sys.version_info >= (3, 10):
    def nametofont(name: str, root: tkinter.Misc | None = ...) -> Font: ...

else:
    def nametofont(name: str) -> Font: ...
