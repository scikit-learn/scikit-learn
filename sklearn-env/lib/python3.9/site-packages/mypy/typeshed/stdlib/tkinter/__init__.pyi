import _tkinter
import sys
from _typeshed import StrOrBytesPath
from enum import Enum
from tkinter.constants import *  # comment this out to find undefined identifier names with flake8
from tkinter.font import _FontDescription
from types import TracebackType
from typing import Any, Callable, Generic, List, Mapping, Optional, Protocol, Sequence, Tuple, Type, TypeVar, Union, overload
from typing_extensions import Literal, TypedDict

# Using anything from tkinter.font in this file means that 'import tkinter'
# seems to also load tkinter.font. That's not how it actually works, but
# unfortunately not much can be done about it. https://github.com/python/typeshed/pull/4346

TclError = _tkinter.TclError
wantobjects: int
TkVersion: float
TclVersion: float
READABLE = _tkinter.READABLE
WRITABLE = _tkinter.WRITABLE
EXCEPTION = _tkinter.EXCEPTION

# Quick guide for figuring out which widget class to choose:
#   - Misc: any widget (don't use BaseWidget because Tk doesn't inherit from BaseWidget)
#   - Widget: anything that is meant to be put into another widget with e.g. pack or grid
#
# Instructions for figuring out the correct type of each widget option:
#  - See discussion on #4363.
#
#  - Find the option from the manual page of the widget. Usually the manual
#    page of a non-ttk widget has the same name as the tkinter class, in the
#    3tk section:
#
#        $ sudo apt install tk-doc
#        $ man 3tk label
#
#    Ttk manual pages tend to have ttk_ prefixed names:
#
#        $ man 3tk ttk_label
#
#    Non-GUI things like the .after() method are often in the 3tcl section:
#
#        $ sudo apt install tcl-doc
#        $ man 3tcl after
#
#    If you don't have man or apt, you can read these manual pages online:
#
#        https://www.tcl.tk/doc/
#
#    Every option has '-' in front of its name in the manual page (and in Tcl).
#    For example, there's an option named '-text' in the label manual page.
#
#  - Tkinter has some options documented in docstrings, but don't rely on them.
#    They aren't updated when a new version of Tk comes out, so the latest Tk
#    manual pages (see above) are much more likely to actually contain all
#    possible options.
#
#    Also, reading tkinter's source code typically won't help much because it
#    uses a lot of **kwargs and duck typing. Typically every argument goes into
#    self.tk.call, which is _tkinter.TkappType.call, and the return value is
#    whatever that returns. The type of that depends on how the Tcl interpreter
#    represents the return value of the executed Tcl code.
#
#  - If you think that int is an appropriate type for something, then you may
#    actually want _ScreenUnits instead.
#
#  - If you think that Callable[something] is an appropriate type for
#    something, then you may actually want Callable[something] | str,
#    because it's often possible to specify a string of Tcl code.
#
#  - If you think the correct type is Iterable[Foo] or Sequence[Foo], it is
#    probably list[Foo] | tuple[Foo, ...], disallowing other sequences such
#    as deques:
#
#        >>> tkinter.Label(font=('Helvetica', 12, collections.deque(['bold'])))
#        Traceback (most recent call last):
#          ...
#        _tkinter.TclError: unknown font style "deque(['bold'])"
#
#  - Some options can be set only in __init__, but all options are available
#    when getting their values with configure's return value or cget.
#
#  - Asks other tkinter users if you haven't worked much with tkinter.

# Some widgets have an option named -compound that accepts different values
# than the _Compound defined here. Many other options have similar things.
_Anchor = Literal["nw", "n", "ne", "w", "center", "e", "sw", "s", "se"]  # manual page: Tk_GetAnchor
_Bitmap = str  # manual page: Tk_GetBitmap
_ButtonCommand = Union[str, Callable[[], Any]]  # return value is returned from Button.invoke()
_CanvasItemId = int
_Color = str  # typically '#rrggbb', '#rgb' or color names.
_Compound = Literal["top", "left", "center", "right", "bottom", "none"]  # -compound in manual page named 'options'
_Cursor = Union[str, Tuple[str], Tuple[str, str], Tuple[str, str, str], Tuple[str, str, str, str]]  # manual page: Tk_GetCursor
_EntryValidateCommand = Union[
    Callable[[], bool], str, List[str], Tuple[str, ...]
]  # example when it's sequence:  entry['invalidcommand'] = [entry.register(print), '%P']
_GridIndex = Union[int, str, Literal["all"]]
_ImageSpec = Union[_Image, str]  # str can be from e.g. tkinter.image_names()
_Padding = Union[
    _ScreenUnits,
    Tuple[_ScreenUnits],
    Tuple[_ScreenUnits, _ScreenUnits],
    Tuple[_ScreenUnits, _ScreenUnits, _ScreenUnits],
    Tuple[_ScreenUnits, _ScreenUnits, _ScreenUnits, _ScreenUnits],
]
_Relief = Literal["raised", "sunken", "flat", "ridge", "solid", "groove"]  # manual page: Tk_GetRelief
_ScreenUnits = Union[str, float]  # manual page: Tk_GetPixels
_XYScrollCommand = Union[str, Callable[[float, float], Any]]  # -xscrollcommand and -yscrollcommand in 'options' manual page
_TakeFocusValue = Union[int, Literal[""], Callable[[str], Optional[bool]]]  # -takefocus in manual page named 'options'

class EventType(str, Enum):
    Activate: str
    ButtonPress: str
    ButtonRelease: str
    Circulate: str
    CirculateRequest: str
    ClientMessage: str
    Colormap: str
    Configure: str
    ConfigureRequest: str
    Create: str
    Deactivate: str
    Destroy: str
    Enter: str
    Expose: str
    FocusIn: str
    FocusOut: str
    GraphicsExpose: str
    Gravity: str
    KeyPress: str
    KeyRelease: str
    Keymap: str
    Leave: str
    Map: str
    MapRequest: str
    Mapping: str
    Motion: str
    MouseWheel: str
    NoExpose: str
    Property: str
    Reparent: str
    ResizeRequest: str
    Selection: str
    SelectionClear: str
    SelectionRequest: str
    Unmap: str
    VirtualEvent: str
    Visibility: str

_W = TypeVar("_W", bound="Misc")
# Events considered covariant because you should never assign to event.widget.
_W_co = TypeVar("_W_co", covariant=True, bound="Misc")

class Event(Generic[_W_co]):
    serial: int
    num: int
    focus: bool
    height: int
    width: int
    keycode: int
    state: int | str
    time: int
    x: int
    y: int
    x_root: int
    y_root: int
    char: str
    send_event: bool
    keysym: str
    keysym_num: int
    type: EventType
    widget: _W_co
    delta: int

def NoDefaultRoot(): ...

_TraceMode = Literal["array", "read", "write", "unset"]

class Variable:
    def __init__(self, master: Misc | None = ..., value: Any | None = ..., name: str | None = ...) -> None: ...
    def set(self, value: Any) -> None: ...
    initialize = set
    def get(self) -> Any: ...
    def trace_add(self, mode: _TraceMode, callback: Callable[[str, str, str], Any]) -> str: ...
    def trace_remove(self, mode: _TraceMode, cbname: str) -> None: ...
    def trace_info(self) -> list[tuple[Tuple[_TraceMode, ...], str]]: ...
    def trace_variable(self, mode, callback): ...  # deprecated
    def trace_vdelete(self, mode, cbname): ...  # deprecated
    def trace_vinfo(self): ...  # deprecated
    trace = trace_variable  # deprecated

class StringVar(Variable):
    def __init__(self, master: Misc | None = ..., value: str | None = ..., name: str | None = ...) -> None: ...
    def set(self, value: str) -> None: ...
    initialize = set
    def get(self) -> str: ...

class IntVar(Variable):
    def __init__(self, master: Misc | None = ..., value: int | None = ..., name: str | None = ...) -> None: ...
    def set(self, value: int) -> None: ...
    initialize = set
    def get(self) -> int: ...

class DoubleVar(Variable):
    def __init__(self, master: Misc | None = ..., value: float | None = ..., name: str | None = ...) -> None: ...
    def set(self, value: float) -> None: ...
    initialize = set
    def get(self) -> float: ...

class BooleanVar(Variable):
    def __init__(self, master: Misc | None = ..., value: bool | None = ..., name: str | None = ...) -> None: ...
    def set(self, value: bool) -> None: ...
    initialize = set
    def get(self) -> bool: ...

def mainloop(n: int = ...) -> None: ...

getint: Any
getdouble: Any

def getboolean(s): ...

class _GridIndexInfo(TypedDict, total=False):
    minsize: _ScreenUnits
    pad: _ScreenUnits
    uniform: str | None
    weight: int

class Misc:
    master: Misc | None
    tk: _tkinter.TkappType
    children: dict[str, Widget]
    def destroy(self) -> None: ...
    def deletecommand(self, name: str) -> None: ...
    def tk_strictMotif(self, boolean: Any | None = ...): ...
    def tk_bisque(self): ...
    def tk_setPalette(self, *args, **kw): ...
    def wait_variable(self, name: str | Variable = ...) -> None: ...
    waitvar = wait_variable
    def wait_window(self, window: Misc | None = ...) -> None: ...
    def wait_visibility(self, window: Misc | None = ...) -> None: ...
    def setvar(self, name: str = ..., value: str = ...): ...
    def getvar(self, name: str = ...): ...
    def getint(self, s): ...
    def getdouble(self, s): ...
    def getboolean(self, s): ...
    def focus_set(self) -> None: ...
    focus = focus_set
    def focus_force(self) -> None: ...
    def focus_get(self) -> Misc | None: ...
    def focus_displayof(self) -> Misc | None: ...
    def focus_lastfor(self) -> Misc | None: ...
    def tk_focusFollowsMouse(self) -> None: ...
    def tk_focusNext(self) -> Misc | None: ...
    def tk_focusPrev(self) -> Misc | None: ...
    @overload
    def after(self, ms: int, func: None = ...) -> None: ...
    @overload
    def after(self, ms: int | Literal["idle"], func: Callable[..., Any], *args: Any) -> str: ...
    # after_idle is essentially partialmethod(after, "idle")
    def after_idle(self, func: Callable[..., Any], *args: Any) -> str: ...
    def after_cancel(self, id: str) -> None: ...
    def bell(self, displayof: Literal[0] | Misc | None = ...): ...
    def clipboard_get(self, *, displayof: Misc = ..., type: str = ...) -> str: ...
    def clipboard_clear(self, *, displayof: Misc = ...) -> None: ...
    def clipboard_append(self, string: str, *, displayof: Misc = ..., format: str = ..., type: str = ...): ...
    def grab_current(self): ...
    def grab_release(self): ...
    def grab_set(self) -> None: ...
    def grab_set_global(self) -> None: ...
    def grab_status(self): ...
    def option_add(self, pattern, value, priority: Any | None = ...): ...
    def option_clear(self): ...
    def option_get(self, name, className): ...
    def option_readfile(self, fileName, priority: Any | None = ...): ...
    def selection_clear(self, **kw): ...
    def selection_get(self, **kw): ...
    def selection_handle(self, command, **kw): ...
    def selection_own(self, **kw): ...
    def selection_own_get(self, **kw): ...
    def send(self, interp, cmd, *args): ...
    def lower(self, belowThis: Any | None = ...): ...
    def tkraise(self, aboveThis: Any | None = ...): ...
    lift = tkraise
    def winfo_atom(self, name: str, displayof: Literal[0] | Misc | None = ...): ...
    def winfo_atomname(self, id: int, displayof: Literal[0] | Misc | None = ...): ...
    def winfo_cells(self) -> int: ...
    def winfo_children(self) -> list[Widget]: ...  # Widget because it can't be Toplevel or Tk
    def winfo_class(self) -> str: ...
    def winfo_colormapfull(self) -> bool: ...
    def winfo_containing(self, rootX: int, rootY: int, displayof: Literal[0] | Misc | None = ...) -> Misc | None: ...
    def winfo_depth(self) -> int: ...
    def winfo_exists(self) -> bool: ...
    def winfo_fpixels(self, number: _ScreenUnits) -> float: ...
    def winfo_geometry(self) -> str: ...
    def winfo_height(self) -> int: ...
    def winfo_id(self) -> int: ...
    def winfo_interps(self, displayof: Literal[0] | Misc | None = ...) -> Tuple[str, ...]: ...
    def winfo_ismapped(self) -> bool: ...
    def winfo_manager(self) -> str: ...
    def winfo_name(self) -> str: ...
    def winfo_parent(self) -> str: ...  # return value needs nametowidget()
    def winfo_pathname(self, id: int, displayof: Literal[0] | Misc | None = ...): ...
    def winfo_pixels(self, number: _ScreenUnits) -> int: ...
    def winfo_pointerx(self) -> int: ...
    def winfo_pointerxy(self) -> tuple[int, int]: ...
    def winfo_pointery(self) -> int: ...
    def winfo_reqheight(self) -> int: ...
    def winfo_reqwidth(self) -> int: ...
    def winfo_rgb(self, color: _Color) -> tuple[int, int, int]: ...
    def winfo_rootx(self) -> int: ...
    def winfo_rooty(self) -> int: ...
    def winfo_screen(self) -> str: ...
    def winfo_screencells(self) -> int: ...
    def winfo_screendepth(self) -> int: ...
    def winfo_screenheight(self) -> int: ...
    def winfo_screenmmheight(self) -> int: ...
    def winfo_screenmmwidth(self) -> int: ...
    def winfo_screenvisual(self) -> str: ...
    def winfo_screenwidth(self) -> int: ...
    def winfo_server(self) -> str: ...
    def winfo_toplevel(self) -> Tk | Toplevel: ...
    def winfo_viewable(self) -> bool: ...
    def winfo_visual(self) -> str: ...
    def winfo_visualid(self) -> str: ...
    def winfo_visualsavailable(self, includeids: int = ...) -> list[tuple[str, int]]: ...
    def winfo_vrootheight(self) -> int: ...
    def winfo_vrootwidth(self) -> int: ...
    def winfo_vrootx(self) -> int: ...
    def winfo_vrooty(self) -> int: ...
    def winfo_width(self) -> int: ...
    def winfo_x(self) -> int: ...
    def winfo_y(self) -> int: ...
    def update(self) -> None: ...
    def update_idletasks(self) -> None: ...
    def bindtags(self, tagList: Any | None = ...): ...
    # bind with isinstance(func, str) doesn't return anything, but all other
    # binds do. The default value of func is not str.
    @overload
    def bind(
        self,
        sequence: str | None = ...,
        func: Callable[[Event[Misc]], Any] | None = ...,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def bind(self, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    @overload
    def bind(self, *, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    # There's no way to know what type of widget bind_all and bind_class
    # callbacks will get, so those are Misc.
    @overload
    def bind_all(
        self,
        sequence: str | None = ...,
        func: Callable[[Event[Misc]], Any] | None = ...,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def bind_all(self, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    @overload
    def bind_all(self, *, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    @overload
    def bind_class(
        self,
        className: str,
        sequence: str | None = ...,
        func: Callable[[Event[Misc]], Any] | None = ...,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def bind_class(self, className: str, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    @overload
    def bind_class(self, className: str, *, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    def unbind(self, sequence: str, funcid: str | None = ...) -> None: ...
    def unbind_all(self, sequence: str) -> None: ...
    def unbind_class(self, className: str, sequence: str) -> None: ...
    def mainloop(self, n: int = ...) -> None: ...
    def quit(self): ...
    def nametowidget(self, name: str | Misc | _tkinter.Tcl_Obj) -> Any: ...
    def register(
        self, func: Callable[..., Any], subst: Callable[..., Sequence[Any]] | None = ..., needcleanup: int = ...
    ) -> str: ...
    def keys(self) -> list[str]: ...
    @overload
    def pack_propagate(self, flag: bool) -> bool | None: ...
    @overload
    def pack_propagate(self) -> None: ...
    propagate = pack_propagate
    def grid_anchor(self, anchor: _Anchor | None = ...) -> None: ...
    anchor = grid_anchor
    @overload
    def grid_bbox(
        self, column: None = ..., row: None = ..., col2: None = ..., row2: None = ...
    ) -> tuple[int, int, int, int] | None: ...
    @overload
    def grid_bbox(self, column: int, row: int, col2: None = ..., row2: None = ...) -> tuple[int, int, int, int] | None: ...
    @overload
    def grid_bbox(self, column: int, row: int, col2: int, row2: int) -> tuple[int, int, int, int] | None: ...
    bbox = grid_bbox
    def grid_columnconfigure(
        self,
        index: _GridIndex,
        cnf: _GridIndexInfo = ...,
        *,
        minsize: _ScreenUnits = ...,
        pad: _ScreenUnits = ...,
        uniform: str = ...,
        weight: int = ...,
    ) -> _GridIndexInfo | Any: ...  # can be None but annoying to check
    def grid_rowconfigure(
        self,
        index: _GridIndex,
        cnf: _GridIndexInfo = ...,
        *,
        minsize: _ScreenUnits = ...,
        pad: _ScreenUnits = ...,
        uniform: str = ...,
        weight: int = ...,
    ) -> _GridIndexInfo | Any: ...  # can be None but annoying to check
    columnconfigure = grid_columnconfigure
    rowconfigure = grid_rowconfigure
    def grid_location(self, x: _ScreenUnits, y: _ScreenUnits) -> tuple[int, int]: ...
    @overload
    def grid_propagate(self, flag: bool) -> None: ...
    @overload
    def grid_propagate(self) -> bool: ...
    def grid_size(self) -> tuple[int, int]: ...
    size = grid_size
    # Widget because Toplevel or Tk is never a slave
    def pack_slaves(self) -> list[Widget]: ...
    def grid_slaves(self, row: int | None = ..., column: int | None = ...) -> list[Widget]: ...
    def place_slaves(self) -> list[Widget]: ...
    slaves = pack_slaves
    def event_add(self, virtual: str, *sequences: str) -> None: ...
    def event_delete(self, virtual: str, *sequences: str) -> None: ...
    def event_generate(
        self,
        sequence: str,
        *,
        above: Misc | int = ...,
        borderwidth: _ScreenUnits = ...,
        button: int = ...,
        count: int = ...,
        data: Any = ...,  # anything with usable str() value
        delta: int = ...,
        detail: str = ...,
        focus: bool = ...,
        height: _ScreenUnits = ...,
        keycode: int = ...,
        keysym: str = ...,
        mode: str = ...,
        override: bool = ...,
        place: Literal["PlaceOnTop", "PlaceOnBottom"] = ...,
        root: Misc | int = ...,
        rootx: _ScreenUnits = ...,
        rooty: _ScreenUnits = ...,
        sendevent: bool = ...,
        serial: int = ...,
        state: int | str = ...,
        subwindow: Misc | int = ...,
        time: int = ...,
        warp: bool = ...,
        width: _ScreenUnits = ...,
        when: Literal["now", "tail", "head", "mark"] = ...,
        x: _ScreenUnits = ...,
        y: _ScreenUnits = ...,
    ) -> None: ...
    def event_info(self, virtual: str | None = ...) -> Tuple[str, ...]: ...
    def image_names(self) -> Tuple[str, ...]: ...
    def image_types(self) -> Tuple[str, ...]: ...
    # See #4363 and #4891
    def __setitem__(self, key: str, value: Any) -> None: ...
    def __getitem__(self, key: str) -> Any: ...
    def cget(self, key: str) -> Any: ...
    def configure(self, cnf: Any = ...) -> Any: ...
    # TODO: config is an alias of configure, but adding that here creates lots of mypy errors

class CallWrapper:
    func: Any
    subst: Any
    widget: Any
    def __init__(self, func, subst, widget): ...
    def __call__(self, *args): ...

class XView:
    @overload
    def xview(self) -> tuple[float, float]: ...
    @overload
    def xview(self, *args: Any) -> Any: ...
    def xview_moveto(self, fraction: float) -> None: ...
    @overload
    def xview_scroll(self, number: int, what: Literal["units", "pages"]) -> None: ...
    @overload
    def xview_scroll(self, number: _ScreenUnits, what: Literal["pixels"]) -> None: ...

class YView:
    @overload
    def yview(self) -> tuple[float, float]: ...
    @overload
    def yview(self, *args: Any) -> Any: ...
    def yview_moveto(self, fraction: float) -> None: ...
    @overload
    def yview_scroll(self, number: int, what: Literal["units", "pages"]) -> None: ...
    @overload
    def yview_scroll(self, number: _ScreenUnits, what: Literal["pixels"]) -> None: ...

class Wm:
    @overload
    def wm_aspect(self, minNumer: int, minDenom: int, maxNumer: int, maxDenom: int) -> None: ...
    @overload
    def wm_aspect(
        self, minNumer: None = ..., minDenom: None = ..., maxNumer: None = ..., maxDenom: None = ...
    ) -> tuple[int, int, int, int] | None: ...
    aspect = wm_aspect
    @overload
    def wm_attributes(self) -> Tuple[Any, ...]: ...
    @overload
    def wm_attributes(self, __option: str) -> Any: ...
    @overload
    def wm_attributes(self, __option: str, __value: Any, *__other_option_value_pairs: Any) -> None: ...
    attributes = wm_attributes
    def wm_client(self, name: str | None = ...) -> str: ...
    client = wm_client
    @overload
    def wm_colormapwindows(self) -> list[Misc]: ...
    @overload
    def wm_colormapwindows(self, __wlist: list[Misc] | Tuple[Misc, ...]) -> None: ...
    @overload
    def wm_colormapwindows(self, __first_wlist_item: Misc, *other_wlist_items: Misc) -> None: ...
    colormapwindows = wm_colormapwindows
    def wm_command(self, value: str | None = ...) -> str: ...
    command = wm_command
    # Some of these always return empty string, but return type is set to None to prevent accidentally using it
    def wm_deiconify(self) -> None: ...
    deiconify = wm_deiconify
    def wm_focusmodel(self, model: Any | None = ...): ...
    focusmodel = wm_focusmodel
    def wm_forget(self, window: Wm) -> None: ...
    forget = wm_forget
    def wm_frame(self): ...
    frame = wm_frame
    @overload
    def wm_geometry(self, newGeometry: None = ...) -> str: ...
    @overload
    def wm_geometry(self, newGeometry: str) -> None: ...
    geometry = wm_geometry
    def wm_grid(
        self, baseWidth: Any | None = ..., baseHeight: Any | None = ..., widthInc: Any | None = ..., heightInc: Any | None = ...
    ): ...
    grid = wm_grid
    def wm_group(self, pathName: Any | None = ...): ...
    group = wm_group
    def wm_iconbitmap(self, bitmap: Any | None = ..., default: Any | None = ...): ...
    iconbitmap = wm_iconbitmap
    def wm_iconify(self) -> None: ...
    iconify = wm_iconify
    def wm_iconmask(self, bitmap: Any | None = ...): ...
    iconmask = wm_iconmask
    def wm_iconname(self, newName: Any | None = ...): ...
    iconname = wm_iconname
    def wm_iconphoto(self, default: bool, __image1: Image, *args: Image) -> None: ...
    iconphoto = wm_iconphoto
    def wm_iconposition(self, x: Any | None = ..., y: Any | None = ...): ...
    iconposition = wm_iconposition
    def wm_iconwindow(self, pathName: Any | None = ...): ...
    iconwindow = wm_iconwindow
    def wm_manage(self, widget): ...
    manage = wm_manage
    @overload
    def wm_maxsize(self, width: None = ..., height: None = ...) -> tuple[int, int]: ...
    @overload
    def wm_maxsize(self, width: int, height: int) -> None: ...
    maxsize = wm_maxsize
    @overload
    def wm_minsize(self, width: None = ..., height: None = ...) -> tuple[int, int]: ...
    @overload
    def wm_minsize(self, width: int, height: int) -> None: ...
    minsize = wm_minsize
    @overload
    def wm_overrideredirect(self, boolean: None = ...) -> bool | None: ...  # returns True or None
    @overload
    def wm_overrideredirect(self, boolean: bool) -> None: ...
    overrideredirect = wm_overrideredirect
    def wm_positionfrom(self, who: Any | None = ...): ...
    positionfrom = wm_positionfrom
    @overload
    def wm_protocol(self, name: str, func: Callable[[], Any] | str) -> None: ...
    @overload
    def wm_protocol(self, name: str, func: None = ...) -> str: ...
    @overload
    def wm_protocol(self, name: None = ..., func: None = ...) -> Tuple[str, ...]: ...
    protocol = wm_protocol
    @overload
    def wm_resizable(self, width: None = ..., height: None = ...) -> tuple[bool, bool]: ...
    @overload
    def wm_resizable(self, width: bool, height: bool) -> None: ...
    resizable = wm_resizable
    def wm_sizefrom(self, who: Any | None = ...): ...
    sizefrom = wm_sizefrom
    @overload
    def wm_state(self, newstate: None = ...) -> str: ...
    @overload
    def wm_state(self, newstate: str) -> None: ...
    state = wm_state
    @overload
    def wm_title(self, string: None = ...) -> str: ...
    @overload
    def wm_title(self, string: str) -> None: ...
    title = wm_title
    @overload
    def wm_transient(self, master: None = ...) -> _tkinter.Tcl_Obj: ...
    @overload
    def wm_transient(self, master: Wm | _tkinter.Tcl_Obj) -> None: ...
    transient = wm_transient
    def wm_withdraw(self) -> None: ...
    withdraw = wm_withdraw

class _ExceptionReportingCallback(Protocol):
    def __call__(self, __exc: Type[BaseException], __val: BaseException, __tb: TracebackType | None) -> Any: ...

class Tk(Misc, Wm):
    master: None
    def __init__(
        # please update ttkthemes stub if you change this
        self,
        screenName: str | None = ...,
        baseName: str | None = ...,
        className: str = ...,
        useTk: bool = ...,
        sync: bool = ...,
        use: str | None = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        menu: Menu = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def loadtk(self) -> None: ...  # differs from _tkinter.TkappType.loadtk
    def destroy(self) -> None: ...
    def readprofile(self, baseName: str, className: str) -> None: ...
    report_callback_exception: _ExceptionReportingCallback
    # Tk has __getattr__ so that tk_instance.foo falls back to tk_instance.tk.foo
    # Please keep in sync with _tkinter.TkappType
    call: Callable[..., Any]
    def eval(self, __code: str) -> str: ...
    adderrorinfo: Any
    createcommand: Any
    createfilehandler: Any
    createtimerhandler: Any
    deletecommand: Any
    deletefilehandler: Any
    dooneevent: Any
    evalfile: Any
    exprboolean: Any
    exprdouble: Any
    exprlong: Any
    exprstring: Any
    getboolean: Any
    getdouble: Any
    getint: Any
    getvar: Any
    globalgetvar: Any
    globalsetvar: Any
    globalunsetvar: Any
    interpaddr: Any
    mainloop: Any
    quit: Any
    record: Any
    setvar: Any
    split: Any
    splitlist: Any
    unsetvar: Any
    wantobjects: Any
    willdispatch: Any

def Tcl(screenName: Any | None = ..., baseName: Any | None = ..., className: str = ..., useTk: bool = ...): ...

_InMiscTotal = TypedDict("_InMiscTotal", {"in": Misc})
_InMiscNonTotal = TypedDict("_InMiscNonTotal", {"in": Misc}, total=False)

class _PackInfo(_InMiscTotal):
    # 'before' and 'after' never appear in _PackInfo
    anchor: _Anchor
    expand: bool
    fill: Literal["none", "x", "y", "both"]
    side: Literal["left", "right", "top", "bottom"]
    # Paddings come out as int or tuple of int, even though any _ScreenUnits
    # can be specified in pack().
    ipadx: int
    ipady: int
    padx: int | tuple[int, int]
    pady: int | tuple[int, int]

class Pack:
    # _PackInfo is not the valid type for cnf because pad stuff accepts any
    # _ScreenUnits instead of int only. I didn't bother to create another
    # TypedDict for cnf because it appears to be a legacy thing that was
    # replaced by **kwargs.
    def pack_configure(
        self,
        cnf: Mapping[str, Any] | None = ...,
        *,
        after: Misc = ...,
        anchor: _Anchor = ...,
        before: Misc = ...,
        expand: int = ...,
        fill: Literal["none", "x", "y", "both"] = ...,
        side: Literal["left", "right", "top", "bottom"] = ...,
        ipadx: _ScreenUnits = ...,
        ipady: _ScreenUnits = ...,
        padx: _ScreenUnits | tuple[_ScreenUnits, _ScreenUnits] = ...,
        pady: _ScreenUnits | tuple[_ScreenUnits, _ScreenUnits] = ...,
        in_: Misc = ...,
        **kw: Any,  # allow keyword argument named 'in', see #4836
    ) -> None: ...
    def pack_forget(self) -> None: ...
    def pack_info(self) -> _PackInfo: ...  # errors if widget hasn't been packed
    pack = pack_configure
    forget = pack_forget
    propagate = Misc.pack_propagate
    # commented out to avoid mypy getting confused with multiple
    # inheritance and how things get overridden with different things
    # info = pack_info
    # pack_propagate = Misc.pack_propagate
    # configure = pack_configure
    # config = pack_configure
    # slaves = Misc.pack_slaves
    # pack_slaves = Misc.pack_slaves

class _PlaceInfo(_InMiscNonTotal):  # empty dict if widget hasn't been placed
    anchor: _Anchor
    bordermode: Literal["inside", "outside", "ignore"]
    width: str  # can be int()ed (even after e.g. widget.place(height='2.3c') or similar)
    height: str  # can be int()ed
    x: str  # can be int()ed
    y: str  # can be int()ed
    relheight: str  # can be float()ed if not empty string
    relwidth: str  # can be float()ed if not empty string
    relx: str  # can be float()ed if not empty string
    rely: str  # can be float()ed if not empty string

class Place:
    def place_configure(
        self,
        cnf: Mapping[str, Any] | None = ...,
        *,
        anchor: _Anchor = ...,
        bordermode: Literal["inside", "outside", "ignore"] = ...,
        width: _ScreenUnits = ...,
        height: _ScreenUnits = ...,
        x: _ScreenUnits = ...,
        y: _ScreenUnits = ...,
        # str allowed for compatibility with place_info()
        relheight: str | float = ...,
        relwidth: str | float = ...,
        relx: str | float = ...,
        rely: str | float = ...,
        in_: Misc = ...,
        **kw: Any,  # allow keyword argument named 'in', see #4836
    ) -> None: ...
    def place_forget(self) -> None: ...
    def place_info(self) -> _PlaceInfo: ...
    place = place_configure
    info = place_info
    # commented out to avoid mypy getting confused with multiple
    # inheritance and how things get overridden with different things
    # config = place_configure
    # configure = place_configure
    # forget = place_forget
    # slaves = Misc.place_slaves
    # place_slaves = Misc.place_slaves

class _GridInfo(_InMiscNonTotal):  # empty dict if widget hasn't been gridded
    column: int
    columnspan: int
    row: int
    rowspan: int
    ipadx: int
    ipady: int
    padx: int | tuple[int, int]
    pady: int | tuple[int, int]
    sticky: str  # consists of letters 'n', 's', 'w', 'e', no repeats, may be empty

class Grid:
    def grid_configure(
        self,
        cnf: Mapping[str, Any] | None = ...,
        *,
        column: int = ...,
        columnspan: int = ...,
        row: int = ...,
        rowspan: int = ...,
        ipadx: _ScreenUnits = ...,
        ipady: _ScreenUnits = ...,
        padx: _ScreenUnits | tuple[_ScreenUnits, _ScreenUnits] = ...,
        pady: _ScreenUnits | tuple[_ScreenUnits, _ScreenUnits] = ...,
        sticky: str = ...,  # consists of letters 'n', 's', 'w', 'e', may contain repeats, may be empty
        in_: Misc = ...,
        **kw: Any,  # allow keyword argument named 'in', see #4836
    ) -> None: ...
    def grid_forget(self) -> None: ...
    def grid_remove(self) -> None: ...
    def grid_info(self) -> _GridInfo: ...
    grid = grid_configure
    location = Misc.grid_location
    size = Misc.grid_size
    # commented out to avoid mypy getting confused with multiple
    # inheritance and how things get overridden with different things
    # bbox = Misc.grid_bbox
    # grid_bbox = Misc.grid_bbox
    # forget = grid_forget
    # info = grid_info
    # grid_location = Misc.grid_location
    # grid_propagate = Misc.grid_propagate
    # grid_size = Misc.grid_size
    # rowconfigure = Misc.grid_rowconfigure
    # grid_rowconfigure = Misc.grid_rowconfigure
    # grid_columnconfigure = Misc.grid_columnconfigure
    # columnconfigure = Misc.grid_columnconfigure
    # config = grid_configure
    # configure = grid_configure
    # propagate = Misc.grid_propagate
    # slaves = Misc.grid_slaves
    # grid_slaves = Misc.grid_slaves

class BaseWidget(Misc):
    master: Misc
    widgetName: Any
    def __init__(self, master, widgetName, cnf=..., kw=..., extra=...): ...
    def destroy(self) -> None: ...

# This class represents any widget except Toplevel or Tk.
class Widget(BaseWidget, Pack, Place, Grid):
    # Allow bind callbacks to take e.g. Event[Label] instead of Event[Misc].
    # Tk and Toplevel get notified for their child widgets' events, but other
    # widgets don't.
    @overload
    def bind(
        self: _W,
        sequence: str | None = ...,
        func: Callable[[Event[_W]], Any] | None = ...,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def bind(self, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    @overload
    def bind(self, *, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...

class Toplevel(BaseWidget, Wm):
    # Toplevel and Tk have the same options because they correspond to the same
    # Tcl/Tk toplevel widget. For some reason, config and configure must be
    # copy/pasted here instead of aliasing as 'config = Tk.config'.
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        class_: str = ...,
        colormap: Literal["new", ""] | Misc = ...,
        container: bool = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        menu: Menu = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        screen: str = ...,  # can't be changed after creating widget
        takefocus: _TakeFocusValue = ...,
        use: int = ...,
        visual: str | tuple[str, int] = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        menu: Menu = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Button(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,  # same as borderwidth
        bg: _Color = ...,  # same as background
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,  # same as borderwidth
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,  # same as foreground
        font: _FontDescription = ...,
        foreground: _Color = ...,
        # width and height must be int for buttons containing just text, but
        # ints are also valid _ScreenUnits
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        # We allow the textvariable to be any Variable, not necessarily
        # StringVar. This is useful for e.g. a button that displays the value
        # of an IntVar.
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def flash(self): ...
    def invoke(self) -> Any: ...

class Canvas(Widget, XView, YView):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        closeenough: float = ...,
        confine: bool = ...,
        cursor: _Cursor = ...,
        # canvas manual page has a section named COORDINATES, and the first
        # part of it describes _ScreenUnits.
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        name: str = ...,
        offset: Any = ...,  # undocumented
        relief: _Relief = ...,
        # Setting scrollregion to None doesn't reset it back to empty,
        # but setting it to () does.
        scrollregion: tuple[_ScreenUnits, _ScreenUnits, _ScreenUnits, _ScreenUnits] | tuple[()] = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        # man page says that state can be 'hidden', but it can't
        state: Literal["normal", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        width: _ScreenUnits = ...,
        xscrollcommand: _XYScrollCommand = ...,
        xscrollincrement: _ScreenUnits = ...,
        yscrollcommand: _XYScrollCommand = ...,
        yscrollincrement: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        closeenough: float = ...,
        confine: bool = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        offset: Any = ...,  # undocumented
        relief: _Relief = ...,
        scrollregion: tuple[_ScreenUnits, _ScreenUnits, _ScreenUnits, _ScreenUnits] | tuple[()] = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        state: Literal["normal", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        width: _ScreenUnits = ...,
        xscrollcommand: _XYScrollCommand = ...,
        xscrollincrement: _ScreenUnits = ...,
        yscrollcommand: _XYScrollCommand = ...,
        yscrollincrement: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def addtag(self, *args): ...  # internal method
    def addtag_above(self, newtag: str, tagOrId: str | _CanvasItemId) -> None: ...
    def addtag_all(self, newtag: str) -> None: ...
    def addtag_below(self, newtag: str, tagOrId: str | _CanvasItemId) -> None: ...
    def addtag_closest(
        self,
        newtag: str,
        x: _ScreenUnits,
        y: _ScreenUnits,
        halo: _ScreenUnits | None = ...,
        start: str | _CanvasItemId | None = ...,
    ) -> None: ...
    def addtag_enclosed(self, newtag: str, x1: _ScreenUnits, y1: _ScreenUnits, x2: _ScreenUnits, y2: _ScreenUnits) -> None: ...
    def addtag_overlapping(self, newtag: str, x1: _ScreenUnits, y1: _ScreenUnits, x2: _ScreenUnits, y2: _ScreenUnits) -> None: ...
    def addtag_withtag(self, newtag: str, tagOrId: str | _CanvasItemId) -> None: ...
    def find(self, *args): ...  # internal method
    def find_above(self, tagOrId: str | _CanvasItemId) -> Tuple[_CanvasItemId, ...]: ...
    def find_all(self) -> Tuple[_CanvasItemId, ...]: ...
    def find_below(self, tagOrId: str | _CanvasItemId) -> Tuple[_CanvasItemId, ...]: ...
    def find_closest(
        self, x: _ScreenUnits, y: _ScreenUnits, halo: _ScreenUnits | None = ..., start: str | _CanvasItemId | None = ...
    ) -> Tuple[_CanvasItemId, ...]: ...
    def find_enclosed(
        self, x1: _ScreenUnits, y1: _ScreenUnits, x2: _ScreenUnits, y2: _ScreenUnits
    ) -> Tuple[_CanvasItemId, ...]: ...
    def find_overlapping(self, x1: _ScreenUnits, y1: _ScreenUnits, x2: _ScreenUnits, y2: float) -> Tuple[_CanvasItemId, ...]: ...
    def find_withtag(self, tagOrId: str | _CanvasItemId) -> Tuple[_CanvasItemId, ...]: ...
    # Canvas.bbox() args are `str | _CanvasItemId`, but mypy rejects that
    # description because it's incompatible with Misc.bbox(), an alias for
    # Misc.grid_bbox(). Yes it is, but there's not much we can do about it.
    def bbox(self, *args: str | _CanvasItemId) -> tuple[int, int, int, int]: ...  # type: ignore
    @overload
    def tag_bind(
        self,
        tagOrId: str | int,
        sequence: str | None = ...,
        func: Callable[[Event[Canvas]], Any] | None = ...,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def tag_bind(
        self, tagOrId: str | int, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...
    ) -> None: ...
    @overload
    def tag_bind(self, tagOrId: str | int, *, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    def tag_unbind(self, tagOrId: str | int, sequence: str, funcid: str | None = ...) -> None: ...
    def canvasx(self, screenx, gridspacing: Any | None = ...): ...
    def canvasy(self, screeny, gridspacing: Any | None = ...): ...
    @overload
    def coords(self) -> list[float]: ...
    @overload
    def coords(self, __args: list[int] | list[float] | Tuple[float, ...]) -> None: ...
    @overload
    def coords(self, __x1: float, __y1: float, *args: float) -> None: ...
    # create_foo() methods accept coords as a list, a tuple, or as separate arguments.
    # Keyword arguments should be the same in each pair of overloads.
    def create_arc(self, *args, **kw) -> _CanvasItemId: ...
    def create_bitmap(self, *args, **kw) -> _CanvasItemId: ...
    def create_image(self, *args, **kw) -> _CanvasItemId: ...
    @overload
    def create_line(
        self,
        __x0: float,
        __y0: float,
        __x1: float,
        __y1: float,
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        arrow: Literal["first", "last", "both"] = ...,
        arrowshape: tuple[float, float, float] = ...,
        capstyle: Literal["round", "projecting", "butt"] = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        joinstyle: Literal["round", "bevel", "miter"] = ...,
        offset: _ScreenUnits = ...,
        smooth: bool = ...,
        splinesteps: float = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_line(
        self,
        __coords: tuple[float, float, float, float] | list[int] | list[float],
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        arrow: Literal["first", "last", "both"] = ...,
        arrowshape: tuple[float, float, float] = ...,
        capstyle: Literal["round", "projecting", "butt"] = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        joinstyle: Literal["round", "bevel", "miter"] = ...,
        offset: _ScreenUnits = ...,
        smooth: bool = ...,
        splinesteps: float = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_oval(
        self,
        __x0: float,
        __y0: float,
        __x1: float,
        __y1: float,
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_oval(
        self,
        __coords: tuple[float, float, float, float] | list[int] | list[float],
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_polygon(
        self,
        __x0: float,
        __y0: float,
        __x1: float,
        __y1: float,
        *xy_pairs: float,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        joinstyle: Literal["round", "bevel", "miter"] = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        smooth: bool = ...,
        splinesteps: float = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_polygon(
        self,
        __coords: Tuple[float, ...] | list[int] | list[float],
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        joinstyle: Literal["round", "bevel", "miter"] = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        smooth: bool = ...,
        splinesteps: float = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_rectangle(
        self,
        __x0: float,
        __y0: float,
        __x1: float,
        __y1: float,
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_rectangle(
        self,
        __coords: tuple[float, float, float, float] | list[int] | list[float],
        *,
        activedash: str | list[int] | Tuple[int, ...] = ...,
        activefill: _Color = ...,
        activeoutline: _Color = ...,
        activeoutlinestipple: _Color = ...,
        activestipple: str = ...,
        activewidth: _ScreenUnits = ...,
        dash: str | list[int] | Tuple[int, ...] = ...,
        dashoffset: _ScreenUnits = ...,
        disableddash: str | list[int] | Tuple[int, ...] = ...,
        disabledfill: _Color = ...,
        disabledoutline: _Color = ...,
        disabledoutlinestipple: _Color = ...,
        disabledstipple: _Bitmap = ...,
        disabledwidth: _ScreenUnits = ...,
        fill: _Color = ...,
        offset: _ScreenUnits = ...,
        outline: _Color = ...,
        outlineoffset: _ScreenUnits = ...,
        outlinestipple: _Bitmap = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_text(
        self,
        __x: float,
        __y: float,
        *,
        activefill: _Color = ...,
        activestipple: str = ...,
        anchor: _Anchor = ...,
        disabledfill: _Color = ...,
        disabledstipple: _Bitmap = ...,
        fill: _Color = ...,
        font: _FontDescription = ...,
        justify: Literal["left", "center", "right"] = ...,
        offset: _ScreenUnits = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        text: float | str = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_text(
        self,
        __coords: tuple[float, float] | list[int] | list[float],
        *,
        activefill: _Color = ...,
        activestipple: str = ...,
        anchor: _Anchor = ...,
        disabledfill: _Color = ...,
        disabledstipple: _Bitmap = ...,
        fill: _Color = ...,
        font: _FontDescription = ...,
        justify: Literal["left", "center", "right"] = ...,
        offset: _ScreenUnits = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        stipple: _Bitmap = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        text: float | str = ...,
        width: _ScreenUnits = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_window(
        self,
        __x: float,
        __y: float,
        *,
        anchor: _Anchor = ...,
        height: _ScreenUnits = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
        window: Widget = ...,
    ) -> _CanvasItemId: ...
    @overload
    def create_window(
        self,
        __coords: tuple[float, float] | list[int] | list[float],
        *,
        anchor: _Anchor = ...,
        height: _ScreenUnits = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
        width: _ScreenUnits = ...,
        window: Widget = ...,
    ) -> _CanvasItemId: ...
    def dchars(self, *args): ...
    def delete(self, *tagsOrCanvasIds: str | _CanvasItemId) -> None: ...
    @overload
    def dtag(self, __tag: str, __tag_to_delete: str | None = ...) -> None: ...
    @overload
    def dtag(self, __id: _CanvasItemId, __tag_to_delete: str) -> None: ...
    def focus(self, *args): ...
    def gettags(self, __tagOrId: str | _CanvasItemId) -> Tuple[str, ...]: ...
    def icursor(self, *args): ...
    def index(self, *args): ...
    def insert(self, *args): ...
    def itemcget(self, tagOrId, option): ...
    # itemconfigure kwargs depend on item type, which is not known when type checking
    def itemconfigure(
        self, tagOrId: str | _CanvasItemId, cnf: dict[str, Any] | None = ..., **kw: Any
    ) -> dict[str, tuple[str, str, str, str, str]] | None: ...
    itemconfig = itemconfigure
    def move(self, *args): ...
    if sys.version_info >= (3, 8):
        def moveto(self, tagOrId: str | _CanvasItemId, x: Literal[""] | float = ..., y: Literal[""] | float = ...) -> None: ...
    def postscript(self, cnf=..., **kw): ...
    # tkinter does:
    #    lower = tag_lower
    #    lift = tkraise = tag_raise
    #
    # But mypy doesn't like aliasing here (maybe because Misc defines the same names)
    def tag_lower(self, __first: str | _CanvasItemId, __second: str | _CanvasItemId | None = ...) -> None: ...
    def lower(self, __first: str | _CanvasItemId, __second: str | _CanvasItemId | None = ...) -> None: ...  # type: ignore
    def tag_raise(self, __first: str | _CanvasItemId, __second: str | _CanvasItemId | None = ...) -> None: ...
    def tkraise(self, __first: str | _CanvasItemId, __second: str | _CanvasItemId | None = ...) -> None: ...  # type: ignore
    def lift(self, __first: str | _CanvasItemId, __second: str | _CanvasItemId | None = ...) -> None: ...  # type: ignore
    def scale(self, *args): ...
    def scan_mark(self, x, y): ...
    def scan_dragto(self, x, y, gain: int = ...): ...
    def select_adjust(self, tagOrId, index): ...
    def select_clear(self): ...
    def select_from(self, tagOrId, index): ...
    def select_item(self): ...
    def select_to(self, tagOrId, index): ...
    def type(self, tagOrId): ...

class Checkbutton(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        offrelief: _Relief = ...,
        # The checkbutton puts a value to its variable when it's checked or
        # unchecked. We don't restrict the type of that value here, so
        # Any-typing is fine.
        #
        # I think Checkbutton shouldn't be generic, because then specifying
        # "any checkbutton regardless of what variable it uses" would be
        # difficult, and we might run into issues just like how list[float]
        # and list[int] are incompatible. Also, we would need a way to
        # specify "Checkbutton not associated with any variable", which is
        # done by setting variable to empty string (the default).
        offvalue: Any = ...,
        onvalue: Any = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        tristateimage: _ImageSpec = ...,
        tristatevalue: Any = ...,
        underline: int = ...,
        variable: Variable | Literal[""] = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        offrelief: _Relief = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        tristateimage: _ImageSpec = ...,
        tristatevalue: Any = ...,
        underline: int = ...,
        variable: Variable | Literal[""] = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def deselect(self): ...
    def flash(self): ...
    def invoke(self) -> Any: ...
    def select(self): ...
    def toggle(self): ...

_EntryIndex = Union[str, int]  # "INDICES" in manual page

class Entry(Widget, XView):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledbackground: _Color = ...,
        disabledforeground: _Color = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        invalidcommand: _EntryValidateCommand = ...,
        invcmd: _EntryValidateCommand = ...,  # same as invalidcommand
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        readonlybackground: _Color = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        show: str = ...,
        state: Literal["normal", "disabled", "readonly"] = ...,
        takefocus: _TakeFocusValue = ...,
        textvariable: Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: _EntryValidateCommand = ...,
        vcmd: _EntryValidateCommand = ...,  # same as validatecommand
        width: int = ...,
        xscrollcommand: _XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledbackground: _Color = ...,
        disabledforeground: _Color = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        invalidcommand: _EntryValidateCommand = ...,
        invcmd: _EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        readonlybackground: _Color = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        show: str = ...,
        state: Literal["normal", "disabled", "readonly"] = ...,
        takefocus: _TakeFocusValue = ...,
        textvariable: Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: _EntryValidateCommand = ...,
        vcmd: _EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: _XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def delete(self, first: _EntryIndex, last: _EntryIndex | None = ...) -> None: ...
    def get(self) -> str: ...
    def icursor(self, index: _EntryIndex) -> None: ...
    def index(self, index: _EntryIndex) -> int: ...
    def insert(self, index: _EntryIndex, string: str) -> None: ...
    def scan_mark(self, x): ...
    def scan_dragto(self, x): ...
    def selection_adjust(self, index: _EntryIndex) -> None: ...
    def selection_clear(self) -> None: ...  # type: ignore
    def selection_from(self, index: _EntryIndex) -> None: ...
    def selection_present(self) -> bool: ...
    def selection_range(self, start: _EntryIndex, end: _EntryIndex) -> None: ...
    def selection_to(self, index: _EntryIndex) -> None: ...
    select_adjust = selection_adjust
    select_clear = selection_clear
    select_from = selection_from
    select_present = selection_present
    select_range = selection_range
    select_to = selection_to

class Frame(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        class_: str = ...,
        colormap: Literal["new", ""] | Misc = ...,
        container: bool = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        visual: str | tuple[str, int] = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Label(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Listbox(Widget, XView, YView):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activestyle: Literal["dotbox", "none", "underline"] = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        exportselection: int = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: int = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        justify: Literal["left", "center", "right"] = ...,
        # There's no tkinter.ListVar, but seems like bare tkinter.Variable
        # actually works for this:
        #
        #    >>> import tkinter
        #    >>> lb = tkinter.Listbox()
        #    >>> var = lb['listvariable'] = tkinter.Variable()
        #    >>> var.set(['foo', 'bar', 'baz'])
        #    >>> lb.get(0, 'end')
        #    ('foo', 'bar', 'baz')
        listvariable: Variable = ...,
        name: str = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        # from listbox man page: "The value of the [selectmode] option may be
        # arbitrary, but the default bindings expect it to be ..."
        #
        # I have never seen anyone setting this to something else than what
        # "the default bindings expect", but let's support it anyway.
        selectmode: str = ...,
        setgrid: bool = ...,
        state: Literal["normal", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        width: int = ...,
        xscrollcommand: _XYScrollCommand = ...,
        yscrollcommand: _XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activestyle: Literal["dotbox", "none", "underline"] = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: int = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        justify: Literal["left", "center", "right"] = ...,
        listvariable: Variable = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        selectmode: str = ...,
        setgrid: bool = ...,
        state: Literal["normal", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        width: int = ...,
        xscrollcommand: _XYScrollCommand = ...,
        yscrollcommand: _XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def activate(self, index): ...
    def bbox(self, index): ...
    def curselection(self): ...
    def delete(self, first, last: Any | None = ...): ...
    def get(self, first, last: Any | None = ...): ...
    def index(self, index): ...
    def insert(self, index, *elements): ...
    def nearest(self, y): ...
    def scan_mark(self, x, y): ...
    def scan_dragto(self, x, y): ...
    def see(self, index): ...
    def selection_anchor(self, index): ...
    select_anchor: Any
    def selection_clear(self, first, last: Any | None = ...): ...  # type: ignore
    select_clear: Any
    def selection_includes(self, index): ...
    select_includes: Any
    def selection_set(self, first, last: Any | None = ...): ...
    select_set: Any
    def size(self): ...
    def itemcget(self, index, option): ...
    def itemconfigure(self, index, cnf: Any | None = ..., **kw): ...
    itemconfig: Any

_MenuIndex = Union[str, int]

class Menu(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeborderwidth: _ScreenUnits = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        name: str = ...,
        postcommand: Callable[[], Any] | str = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        takefocus: _TakeFocusValue = ...,
        tearoff: int = ...,
        # I guess tearoffcommand arguments are supposed to be widget objects,
        # but they are widget name strings. Use nametowidget() to handle the
        # arguments of tearoffcommand.
        tearoffcommand: Callable[[str, str], Any] | str = ...,
        title: str = ...,
        type: Literal["menubar", "tearoff", "normal"] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeborderwidth: _ScreenUnits = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        postcommand: Callable[[], Any] | str = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        takefocus: _TakeFocusValue = ...,
        tearoff: bool = ...,
        tearoffcommand: Callable[[str, str], Any] | str = ...,
        title: str = ...,
        type: Literal["menubar", "tearoff", "normal"] = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def tk_popup(self, x: int, y: int, entry: _MenuIndex = ...) -> None: ...
    def activate(self, index: _MenuIndex) -> None: ...
    def add(self, itemType, cnf=..., **kw): ...  # docstring says "Internal function."
    def insert(self, index, itemType, cnf=..., **kw): ...  # docstring says "Internal function."
    def add_cascade(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        label: str = ...,
        menu: Menu = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
    ) -> None: ...
    def add_checkbutton(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        label: str = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
        variable: Variable = ...,
    ) -> None: ...
    def add_command(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        label: str = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
    ) -> None: ...
    def add_radiobutton(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        label: str = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Variable = ...,
    ) -> None: ...
    def add_separator(self, cnf: dict[str, Any] | None = ..., *, background: _Color = ...) -> None: ...
    def insert_cascade(
        self,
        index: _MenuIndex,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        label: str = ...,
        menu: Menu = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
    ) -> None: ...
    def insert_checkbutton(
        self,
        index: _MenuIndex,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        label: str = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
        variable: Variable = ...,
    ) -> None: ...
    def insert_command(
        self,
        index: _MenuIndex,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        label: str = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
    ) -> None: ...
    def insert_radiobutton(
        self,
        index: _MenuIndex,
        cnf: dict[str, Any] | None = ...,
        *,
        accelerator: str = ...,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        background: _Color = ...,
        bitmap: _Bitmap = ...,
        columnbreak: int = ...,
        command: Callable[[], Any] | str = ...,
        compound: _Compound = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        hidemargin: bool = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        label: str = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Variable = ...,
    ) -> None: ...
    def insert_separator(self, index: _MenuIndex, cnf: dict[str, Any] | None = ..., *, background: _Color = ...) -> None: ...
    def delete(self, index1: _MenuIndex, index2: _MenuIndex | None = ...) -> None: ...
    def entrycget(self, index: _MenuIndex, option: str) -> Any: ...
    def entryconfigure(
        self, index: _MenuIndex, cnf: dict[str, Any] | None = ..., **kw: Any
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    entryconfig = entryconfigure
    def index(self, index: _MenuIndex) -> int | None: ...
    def invoke(self, index: _MenuIndex) -> Any: ...
    def post(self, x: int, y: int) -> None: ...
    def type(self, index: _MenuIndex) -> Literal["cascade", "checkbutton", "command", "radiobutton", "separator"]: ...
    def unpost(self) -> None: ...
    def xposition(self, index: _MenuIndex) -> int: ...
    def yposition(self, index: _MenuIndex) -> int: ...

class Menubutton(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        menu: Menu = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        menu: Menu = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        underline: int = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Message(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        anchor: _Anchor = ...,
        aspect: int = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        # there's width but no height
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        anchor: _Anchor = ...,
        aspect: int = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        justify: Literal["left", "center", "right"] = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Radiobutton(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        offrelief: _Relief = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        tristateimage: _ImageSpec = ...,
        tristatevalue: Any = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Variable | Literal[""] = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activeforeground: _Color = ...,
        anchor: _Anchor = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bitmap: _Bitmap = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: _ButtonCommand = ...,
        compound: _Compound = ...,
        cursor: _Cursor = ...,
        disabledforeground: _Color = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        image: _ImageSpec = ...,
        indicatoron: bool = ...,
        justify: Literal["left", "center", "right"] = ...,
        offrelief: _Relief = ...,
        overrelief: _Relief = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectcolor: _Color = ...,
        selectimage: _ImageSpec = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: Variable = ...,
        tristateimage: _ImageSpec = ...,
        tristatevalue: Any = ...,
        underline: int = ...,
        value: Any = ...,
        variable: Variable | Literal[""] = ...,
        width: _ScreenUnits = ...,
        wraplength: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def deselect(self): ...
    def flash(self): ...
    def invoke(self) -> Any: ...
    def select(self): ...

class Scale(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bigincrement: float = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        # don't know why the callback gets string instead of float
        command: str | Callable[[str], Any] = ...,
        cursor: _Cursor = ...,
        digits: int = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        from_: float = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        label: str = ...,
        length: _ScreenUnits = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        resolution: float = ...,
        showvalue: bool = ...,
        sliderlength: _ScreenUnits = ...,
        sliderrelief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        tickinterval: float = ...,
        to: float = ...,
        troughcolor: _Color = ...,
        variable: IntVar | DoubleVar = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        bigincrement: float = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: str | Callable[[str], Any] = ...,
        cursor: _Cursor = ...,
        digits: int = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        from_: float = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        label: str = ...,
        length: _ScreenUnits = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        resolution: float = ...,
        showvalue: bool = ...,
        sliderlength: _ScreenUnits = ...,
        sliderrelief: _Relief = ...,
        state: Literal["normal", "active", "disabled"] = ...,
        takefocus: _TakeFocusValue = ...,
        tickinterval: float = ...,
        to: float = ...,
        troughcolor: _Color = ...,
        variable: IntVar | DoubleVar = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def get(self): ...
    def set(self, value): ...
    def coords(self, value: Any | None = ...): ...
    def identify(self, x, y): ...

class Scrollbar(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activerelief: _Relief = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        # There are many ways how the command may get called. Search for
        # 'SCROLLING COMMANDS' in scrollbar man page. There doesn't seem to
        # be any way to specify an overloaded callback function, so we say
        # that it can take any args while it can't in reality.
        command: Callable[..., tuple[float, float] | None] | str = ...,
        cursor: _Cursor = ...,
        elementborderwidth: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        jump: bool = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        takefocus: _TakeFocusValue = ...,
        troughcolor: _Color = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        activerelief: _Relief = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        command: Callable[..., tuple[float, float] | None] | str = ...,
        cursor: _Cursor = ...,
        elementborderwidth: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        jump: bool = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        takefocus: _TakeFocusValue = ...,
        troughcolor: _Color = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def activate(self, index: Any | None = ...): ...
    def delta(self, deltax, deltay): ...
    def fraction(self, x, y): ...
    def identify(self, x, y): ...
    def get(self): ...
    def set(self, first, last): ...

_TextIndex = Union[_tkinter.Tcl_Obj, str, float, Misc]

class Text(Widget, XView, YView):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        autoseparators: bool = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        blockcursor: bool = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        endline: int | Literal[""] = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        # width is always int, but height is allowed to be ScreenUnits.
        # This doesn't make any sense to me, and this isn't documented.
        # The docs seem to say that both should be integers.
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        inactiveselectbackground: _Color = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertunfocussed: Literal["none", "hollow", "solid"] = ...,
        insertwidth: _ScreenUnits = ...,
        maxundo: int = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        setgrid: bool = ...,
        spacing1: _ScreenUnits = ...,
        spacing2: _ScreenUnits = ...,
        spacing3: _ScreenUnits = ...,
        startline: int | Literal[""] = ...,
        state: Literal["normal", "disabled"] = ...,
        # Literal inside Tuple doesn't actually work
        tabs: _ScreenUnits | str | Tuple[_ScreenUnits | str, ...] = ...,
        tabstyle: Literal["tabular", "wordprocessor"] = ...,
        takefocus: _TakeFocusValue = ...,
        undo: bool = ...,
        width: int = ...,
        wrap: Literal["none", "char", "word"] = ...,
        xscrollcommand: _XYScrollCommand = ...,
        yscrollcommand: _XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        autoseparators: bool = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        blockcursor: bool = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        endline: int | Literal[""] = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        inactiveselectbackground: _Color = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertunfocussed: Literal["none", "hollow", "solid"] = ...,
        insertwidth: _ScreenUnits = ...,
        maxundo: int = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        setgrid: bool = ...,
        spacing1: _ScreenUnits = ...,
        spacing2: _ScreenUnits = ...,
        spacing3: _ScreenUnits = ...,
        startline: int | Literal[""] = ...,
        state: Literal["normal", "disabled"] = ...,
        tabs: _ScreenUnits | str | Tuple[_ScreenUnits | str, ...] = ...,
        tabstyle: Literal["tabular", "wordprocessor"] = ...,
        takefocus: _TakeFocusValue = ...,
        undo: bool = ...,
        width: int = ...,
        wrap: Literal["none", "char", "word"] = ...,
        xscrollcommand: _XYScrollCommand = ...,
        yscrollcommand: _XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def bbox(self, index: _TextIndex) -> tuple[int, int, int, int] | None: ...  # type: ignore
    def compare(self, index1: _TextIndex, op: Literal["<", "<=", "==", ">=", ">", "!="], index2: _TextIndex) -> bool: ...
    def count(self, index1, index2, *args): ...  # TODO
    @overload
    def debug(self, boolean: None = ...) -> bool: ...
    @overload
    def debug(self, boolean: bool) -> None: ...
    def delete(self, index1: _TextIndex, index2: _TextIndex | None = ...) -> None: ...
    def dlineinfo(self, index: _TextIndex) -> tuple[int, int, int, int, int] | None: ...
    @overload
    def dump(
        self,
        index1: _TextIndex,
        index2: _TextIndex | None = ...,
        command: None = ...,
        *,
        all: bool = ...,
        image: bool = ...,
        mark: bool = ...,
        tag: bool = ...,
        text: bool = ...,
        window: bool = ...,
    ) -> list[tuple[str, str, str]]: ...
    @overload
    def dump(
        self,
        index1: _TextIndex,
        index2: _TextIndex | None,
        command: Callable[[str, str, str], Any] | str,
        *,
        all: bool = ...,
        image: bool = ...,
        mark: bool = ...,
        tag: bool = ...,
        text: bool = ...,
        window: bool = ...,
    ) -> None: ...
    @overload
    def dump(
        self,
        index1: _TextIndex,
        index2: _TextIndex | None = ...,
        *,
        command: Callable[[str, str, str], Any] | str,
        all: bool = ...,
        image: bool = ...,
        mark: bool = ...,
        tag: bool = ...,
        text: bool = ...,
        window: bool = ...,
    ) -> None: ...
    def edit(self, *args): ...  # docstring says "Internal method"
    @overload
    def edit_modified(self, arg: None = ...) -> bool: ...  # actually returns Literal[0, 1]
    @overload
    def edit_modified(self, arg: bool) -> None: ...  # actually returns empty string
    def edit_redo(self) -> None: ...  # actually returns empty string
    def edit_reset(self) -> None: ...  # actually returns empty string
    def edit_separator(self) -> None: ...  # actually returns empty string
    def edit_undo(self) -> None: ...  # actually returns empty string
    def get(self, index1: _TextIndex, index2: _TextIndex | None = ...) -> str: ...
    # TODO: image_* methods
    def image_cget(self, index, option): ...
    def image_configure(self, index, cnf: Any | None = ..., **kw): ...
    def image_create(self, index, cnf=..., **kw): ...
    def image_names(self): ...
    def index(self, index: _TextIndex) -> str: ...
    def insert(self, index: _TextIndex, chars: str, *args: str | list[str] | Tuple[str, ...]) -> None: ...
    @overload
    def mark_gravity(self, markName: str, direction: None = ...) -> Literal["left", "right"]: ...
    @overload
    def mark_gravity(self, markName: str, direction: Literal["left", "right"]) -> None: ...  # actually returns empty string
    def mark_names(self) -> Tuple[str, ...]: ...
    def mark_set(self, markName: str, index: _TextIndex) -> None: ...
    def mark_unset(self, *markNames: str) -> None: ...
    def mark_next(self, index: _TextIndex) -> str | None: ...
    def mark_previous(self, index: _TextIndex) -> str | None: ...
    # **kw of peer_create is same as the kwargs of Text.__init__
    def peer_create(self, newPathName: str | Text, cnf: dict[str, Any] = ..., **kw: Any) -> None: ...
    def peer_names(self) -> Tuple[_tkinter.Tcl_Obj, ...]: ...
    def replace(self, index1: _TextIndex, index2: _TextIndex, chars: str, *args: str | list[str] | Tuple[str, ...]) -> None: ...
    def scan_mark(self, x: int, y: int) -> None: ...
    def scan_dragto(self, x: int, y: int) -> None: ...
    def search(
        self,
        pattern: str,
        index: _TextIndex,
        stopindex: _TextIndex | None = ...,
        forwards: bool | None = ...,
        backwards: bool | None = ...,
        exact: bool | None = ...,
        regexp: bool | None = ...,
        nocase: bool | None = ...,
        count: Variable | None = ...,
        elide: bool | None = ...,
    ) -> str: ...  # returns empty string for not found
    def see(self, index: _TextIndex) -> None: ...
    def tag_add(self, tagName: str, index1: _TextIndex, *args: _TextIndex) -> None: ...
    # tag_bind stuff is very similar to Canvas
    @overload
    def tag_bind(
        self,
        tagName: str,
        sequence: str | None,
        func: Callable[[Event[Text]], Any] | None,
        add: Literal["", "+"] | bool | None = ...,
    ) -> str: ...
    @overload
    def tag_bind(self, tagName: str, sequence: str | None, func: str, add: Literal["", "+"] | bool | None = ...) -> None: ...
    def tag_unbind(self, tagName: str, sequence: str, funcid: str | None = ...) -> None: ...
    # allowing any string for cget instead of just Literals because there's no other way to look up tag options
    def tag_cget(self, tagName: str, option: str) -> Any: ...
    @overload
    def tag_configure(
        self,
        tagName: str,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bgstipple: _Bitmap = ...,
        borderwidth: _ScreenUnits = ...,
        border: _ScreenUnits = ...,  # alias for borderwidth
        elide: bool = ...,
        fgstipple: _Bitmap = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        justify: Literal["left", "right", "center"] = ...,
        lmargin1: _ScreenUnits = ...,
        lmargin2: _ScreenUnits = ...,
        lmargincolor: _Color = ...,
        offset: _ScreenUnits = ...,
        overstrike: bool = ...,
        overstrikefg: _Color = ...,
        relief: _Relief = ...,
        rmargin: _ScreenUnits = ...,
        rmargincolor: _Color = ...,
        selectbackground: _Color = ...,
        selectforeground: _Color = ...,
        spacing1: _ScreenUnits = ...,
        spacing2: _ScreenUnits = ...,
        spacing3: _ScreenUnits = ...,
        tabs: Any = ...,  # the exact type is kind of complicated, see manual page
        tabstyle: Literal["tabular", "wordprocessor"] = ...,
        underline: bool = ...,
        underlinefg: _Color = ...,
        wrap: Literal["none", "char", "word"] = ...,  # be careful with "none" vs None
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def tag_configure(self, tagName: str, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    tag_config = tag_configure
    def tag_delete(self, __first_tag_name: str, *tagNames: str) -> None: ...  # error if no tag names given
    def tag_lower(self, tagName: str, belowThis: str | None = ...) -> None: ...
    def tag_names(self, index: _TextIndex | None = ...) -> Tuple[str, ...]: ...
    def tag_nextrange(self, tagName: str, index1: _TextIndex, index2: _TextIndex | None = ...) -> tuple[str, str] | tuple[()]: ...
    def tag_prevrange(self, tagName: str, index1: _TextIndex, index2: _TextIndex | None = ...) -> tuple[str, str] | tuple[()]: ...
    def tag_raise(self, tagName: str, aboveThis: str | None = ...) -> None: ...
    def tag_ranges(self, tagName: str) -> Tuple[_tkinter.Tcl_Obj, ...]: ...
    # tag_remove and tag_delete are different
    def tag_remove(self, tagName: str, index1: _TextIndex, index2: _TextIndex | None = ...) -> None: ...
    # TODO: window_* methods
    def window_cget(self, index, option): ...
    def window_configure(self, index, cnf: Any | None = ..., **kw): ...
    window_config = window_configure
    def window_create(self, index, cnf=..., **kw): ...
    def window_names(self): ...
    def yview_pickplace(self, *what): ...  # deprecated

class _setit:
    def __init__(self, var, value, callback: Any | None = ...): ...
    def __call__(self, *args): ...

# manual page: tk_optionMenu
class OptionMenu(Menubutton):
    widgetName: Any
    menuname: Any
    def __init__(
        # differs from other widgets
        self,
        master: Misc | None,
        variable: StringVar,
        value: str,
        *values: str,
        # kwarg only from now on
        command: Callable[[StringVar], Any] | None = ...,
    ) -> None: ...
    # configure, config, cget are inherited from Menubutton
    # destroy and __getitem__ are overridden, signature does not change

class _Image(Protocol):
    tk: _tkinter.TkappType
    def height(self) -> int: ...
    def width(self) -> int: ...

class Image:
    name: Any
    tk: _tkinter.TkappType
    def __init__(self, imgtype, name: Any | None = ..., cnf=..., master: Misc | _tkinter.TkappType | None = ..., **kw): ...
    def __del__(self): ...
    def __setitem__(self, key, value): ...
    def __getitem__(self, key): ...
    configure: Any
    config: Any
    def height(self) -> int: ...
    def type(self): ...
    def width(self) -> int: ...

class PhotoImage(Image):
    def __init__(
        self,
        name: str | None = ...,
        cnf: dict[str, Any] = ...,
        master: Misc | _tkinter.TkappType | None = ...,
        *,
        data: str | bytes = ...,  # not same as data argument of put()
        format: str = ...,
        file: StrOrBytesPath = ...,
        gamma: float = ...,
        height: int = ...,
        palette: int | str = ...,
        width: int = ...,
    ) -> None: ...
    def configure(
        self,
        *,
        data: str | bytes = ...,
        format: str = ...,
        file: StrOrBytesPath = ...,
        gamma: float = ...,
        height: int = ...,
        palette: int | str = ...,
        width: int = ...,
    ) -> None: ...
    config = configure
    def blank(self) -> None: ...
    def cget(self, option: str) -> str: ...
    def __getitem__(self, key: str) -> str: ...  # always string: image['height'] can be '0'
    def copy(self) -> PhotoImage: ...
    def zoom(self, x: int, y: int | Literal[""] = ...) -> PhotoImage: ...
    def subsample(self, x: int, y: int | Literal[""] = ...) -> PhotoImage: ...
    def get(self, x: int, y: int) -> tuple[int, int, int]: ...
    def put(
        self,
        data: (
            str
            | list[str]
            | list[list[_Color]]
            | list[Tuple[_Color, ...]]
            | Tuple[str, ...]
            | Tuple[list[_Color], ...]
            | Tuple[Tuple[_Color, ...], ...]
        ),
        to: tuple[int, int] | None = ...,
    ) -> None: ...
    def write(self, filename: StrOrBytesPath, format: str | None = ..., from_coords: tuple[int, int] | None = ...) -> None: ...
    if sys.version_info >= (3, 8):
        def transparency_get(self, x: int, y: int) -> bool: ...
        def transparency_set(self, x: int, y: int, boolean: bool) -> None: ...

class BitmapImage(Image):
    def __init__(
        self,
        name: Any | None = ...,
        cnf: dict[str, Any] = ...,
        master: Misc | _tkinter.TkappType | None = ...,
        *,
        background: _Color = ...,
        data: str | bytes = ...,
        file: StrOrBytesPath = ...,
        foreground: _Color = ...,
        maskdata: str = ...,
        maskfile: StrOrBytesPath = ...,
    ) -> None: ...

def image_names() -> Tuple[str, ...]: ...
def image_types() -> Tuple[str, ...]: ...

class Spinbox(Widget, XView):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        buttonbackground: _Color = ...,
        buttoncursor: _Cursor = ...,
        buttondownrelief: _Relief = ...,
        buttonuprelief: _Relief = ...,
        # percent substitutions don't seem to be supported, it's similar to Entry's validation stuff
        command: Callable[[], Any] | str | list[str] | Tuple[str, ...] = ...,
        cursor: _Cursor = ...,
        disabledbackground: _Color = ...,
        disabledforeground: _Color = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        format: str = ...,
        from_: float = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        increment: float = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        invalidcommand: _EntryValidateCommand = ...,
        invcmd: _EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        readonlybackground: _Color = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        state: Literal["normal", "disabled", "readonly"] = ...,
        takefocus: _TakeFocusValue = ...,
        textvariable: Variable = ...,
        to: float = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: _EntryValidateCommand = ...,
        vcmd: _EntryValidateCommand = ...,
        values: list[str] | Tuple[str, ...] = ...,
        width: int = ...,
        wrap: bool = ...,
        xscrollcommand: _XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        activebackground: _Color = ...,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        buttonbackground: _Color = ...,
        buttoncursor: _Cursor = ...,
        buttondownrelief: _Relief = ...,
        buttonuprelief: _Relief = ...,
        command: Callable[[], Any] | str | list[str] | Tuple[str, ...] = ...,
        cursor: _Cursor = ...,
        disabledbackground: _Color = ...,
        disabledforeground: _Color = ...,
        exportselection: bool = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        format: str = ...,
        from_: float = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        increment: float = ...,
        insertbackground: _Color = ...,
        insertborderwidth: _ScreenUnits = ...,
        insertofftime: int = ...,
        insertontime: int = ...,
        insertwidth: _ScreenUnits = ...,
        invalidcommand: _EntryValidateCommand = ...,
        invcmd: _EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        readonlybackground: _Color = ...,
        relief: _Relief = ...,
        repeatdelay: int = ...,
        repeatinterval: int = ...,
        selectbackground: _Color = ...,
        selectborderwidth: _ScreenUnits = ...,
        selectforeground: _Color = ...,
        state: Literal["normal", "disabled", "readonly"] = ...,
        takefocus: _TakeFocusValue = ...,
        textvariable: Variable = ...,
        to: float = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: _EntryValidateCommand = ...,
        vcmd: _EntryValidateCommand = ...,
        values: list[str] | Tuple[str, ...] = ...,
        width: int = ...,
        wrap: bool = ...,
        xscrollcommand: _XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def bbox(self, index): ...
    def delete(self, first, last: Any | None = ...): ...
    def get(self): ...
    def icursor(self, index): ...
    def identify(self, x, y): ...
    def index(self, index): ...
    def insert(self, index, s): ...
    # spinbox.invoke("asdf") gives error mentioning .invoke("none"), but it's not documented
    def invoke(self, element: Literal["none", "buttonup", "buttondown"]) -> Literal[""]: ...
    def scan(self, *args): ...
    def scan_mark(self, x): ...
    def scan_dragto(self, x): ...
    def selection(self, *args: Any) -> Tuple[int, ...]: ...
    def selection_adjust(self, index): ...
    def selection_clear(self): ...
    def selection_element(self, element: Any | None = ...): ...
    if sys.version_info >= (3, 8):
        def selection_from(self, index: int) -> None: ...
        def selection_present(self) -> None: ...
        def selection_range(self, start: int, end: int) -> None: ...
        def selection_to(self, index: int) -> None: ...

class LabelFrame(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        class_: str = ...,
        colormap: Literal["new", ""] | Misc = ...,
        container: bool = ...,  # undocumented
        cursor: _Cursor = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        # 'ne' and 'en' are valid labelanchors, but only 'ne' is a valid _Anchor.
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: Misc = ...,
        name: str = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        visual: str | tuple[str, int] = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        fg: _Color = ...,
        font: _FontDescription = ...,
        foreground: _Color = ...,
        height: _ScreenUnits = ...,
        highlightbackground: _Color = ...,
        highlightcolor: _Color = ...,
        highlightthickness: _ScreenUnits = ...,
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: Misc = ...,
        padx: _ScreenUnits = ...,
        pady: _ScreenUnits = ...,
        relief: _Relief = ...,
        takefocus: _TakeFocusValue = ...,
        text: float | str = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class PanedWindow(Widget):
    def __init__(
        self,
        master: Misc | None = ...,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        handlepad: _ScreenUnits = ...,
        handlesize: _ScreenUnits = ...,
        height: _ScreenUnits = ...,
        name: str = ...,
        opaqueresize: bool = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        proxybackground: _Color = ...,
        proxyborderwidth: _ScreenUnits = ...,
        proxyrelief: _Relief = ...,
        relief: _Relief = ...,
        sashcursor: _Cursor = ...,
        sashpad: _ScreenUnits = ...,
        sashrelief: _Relief = ...,
        sashwidth: _ScreenUnits = ...,
        showhandle: bool = ...,
        width: _ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: _Color = ...,
        bd: _ScreenUnits = ...,
        bg: _Color = ...,
        border: _ScreenUnits = ...,
        borderwidth: _ScreenUnits = ...,
        cursor: _Cursor = ...,
        handlepad: _ScreenUnits = ...,
        handlesize: _ScreenUnits = ...,
        height: _ScreenUnits = ...,
        opaqueresize: bool = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        proxybackground: _Color = ...,
        proxyborderwidth: _ScreenUnits = ...,
        proxyrelief: _Relief = ...,
        relief: _Relief = ...,
        sashcursor: _Cursor = ...,
        sashpad: _ScreenUnits = ...,
        sashrelief: _Relief = ...,
        sashwidth: _ScreenUnits = ...,
        showhandle: bool = ...,
        width: _ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def add(self, child: Widget, **kw): ...
    def remove(self, child): ...
    forget: Any
    def identify(self, x, y): ...
    def proxy(self, *args): ...
    def proxy_coord(self): ...
    def proxy_forget(self): ...
    def proxy_place(self, x, y): ...
    def sash(self, *args): ...
    def sash_coord(self, index): ...
    def sash_mark(self, index): ...
    def sash_place(self, index, x, y): ...
    def panecget(self, child, option): ...
    def paneconfigure(self, tagOrId, cnf: Any | None = ..., **kw): ...
    paneconfig: Any
    def panes(self): ...
