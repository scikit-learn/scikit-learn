import _tkinter
import sys
import tkinter
from tkinter.font import _FontDescription
from typing import Any, Callable, Tuple, Union, overload
from typing_extensions import Literal, TypedDict

def tclobjs_to_py(adict): ...
def setup_master(master: Any | None = ...): ...

# from ttk_widget (aka ttk::widget) manual page, differs from tkinter._Compound
_TtkCompound = Literal["text", "image", tkinter._Compound]

class Style:
    master: Any
    tk: _tkinter.TkappType
    def __init__(self, master: tkinter.Misc | None = ...): ...
    def configure(self, style, query_opt: Any | None = ..., **kw): ...
    def map(self, style, query_opt: Any | None = ..., **kw): ...
    def lookup(self, style, option, state: Any | None = ..., default: Any | None = ...): ...
    def layout(self, style, layoutspec: Any | None = ...): ...
    def element_create(self, elementname, etype, *args, **kw): ...
    def element_names(self): ...
    def element_options(self, elementname): ...
    def theme_create(self, themename, parent: Any | None = ..., settings: Any | None = ...): ...
    def theme_settings(self, themename, settings): ...
    def theme_names(self) -> Tuple[str, ...]: ...
    @overload
    def theme_use(self, themename: str) -> None: ...
    @overload
    def theme_use(self, themename: None = ...) -> str: ...

class Widget(tkinter.Widget):
    def __init__(self, master: tkinter.Misc | None, widgetname, kw: Any | None = ...): ...
    def identify(self, x, y): ...
    def instate(self, statespec, callback: Any | None = ..., *args, **kw): ...
    def state(self, statespec: Any | None = ...): ...

class Button(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        default: Literal["normal", "active", "disabled"] = ...,
        image: tkinter._ImageSpec = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self) -> Any: ...

class Checkbutton(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        # Seems like variable can be empty string, but actually setting it to
        # empty string segfaults before Tcl 8.6.9. Search for ttk::checkbutton
        # here: https://sourceforge.net/projects/tcl/files/Tcl/8.6.9/tcltk-release-notes-8.6.9.txt/view
        variable: tkinter.Variable = ...,
        width: int | Literal[""] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        offvalue: Any = ...,
        onvalue: Any = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        variable: tkinter.Variable = ...,
        width: int | Literal[""] = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self) -> Any: ...

class Entry(Widget, tkinter.Entry):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        widget: str | None = ...,
        *,
        background: tkinter._Color = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Entry().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        show: str = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def config(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    def bbox(self, index): ...
    def identify(self, x, y): ...
    def validate(self): ...

class Combobox(Entry):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        background: tkinter._Color = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,  # undocumented
        foreground: tkinter._Color = ...,  # undocumented
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,  # undocumented
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        postcommand: Callable[[], Any] | str = ...,
        show: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,  # undocumented
        validatecommand: tkinter._EntryValidateCommand = ...,  # undocumented
        values: list[str] | Tuple[str, ...] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,  # undocumented
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        postcommand: Callable[[], Any] | str = ...,
        show: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        values: list[str] | Tuple[str, ...] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Combobox().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        background: tkinter._Color = ...,
        cursor: tkinter._Cursor = ...,
        exportselection: bool = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        height: int = ...,
        invalidcommand: tkinter._EntryValidateCommand = ...,
        justify: Literal["left", "center", "right"] = ...,
        postcommand: Callable[[], Any] | str = ...,
        show: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        textvariable: tkinter.Variable = ...,
        validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
        validatecommand: tkinter._EntryValidateCommand = ...,
        values: list[str] | Tuple[str, ...] = ...,
        width: int = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def config(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    def current(self, newindex: Any | None = ...): ...
    def set(self, value): ...

class Frame(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Label(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        anchor: tkinter._Anchor = ...,
        background: tkinter._Color = ...,
        border: tkinter._ScreenUnits = ...,  # alias for borderwidth
        borderwidth: tkinter._ScreenUnits = ...,  # undocumented
        class_: str = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        image: tkinter._ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
        wraplength: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        anchor: tkinter._Anchor = ...,
        background: tkinter._Color = ...,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        font: _FontDescription = ...,
        foreground: tkinter._Color = ...,
        image: tkinter._ImageSpec = ...,
        justify: Literal["left", "center", "right"] = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
        wraplength: tkinter._ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Labelframe(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,  # undocumented
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: tkinter.Misc = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,  # undocumented
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        underline: int = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        border: tkinter._ScreenUnits = ...,
        borderwidth: tkinter._ScreenUnits = ...,
        cursor: tkinter._Cursor = ...,
        height: tkinter._ScreenUnits = ...,
        labelanchor: Literal["nw", "n", "ne", "en", "e", "es", "se", "s", "sw", "ws", "w", "wn"] = ...,
        labelwidget: tkinter.Misc = ...,
        padding: tkinter._Padding = ...,
        relief: tkinter._Relief = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        underline: int = ...,
        width: tkinter._ScreenUnits = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

LabelFrame = Labelframe

class Menubutton(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        image: tkinter._ImageSpec = ...,
        menu: tkinter.Menu = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        direction: Literal["above", "below", "left", "right", "flush"] = ...,
        image: tkinter._ImageSpec = ...,
        menu: tkinter.Menu = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        width: int | Literal[""] = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Notebook(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        padding: tkinter._Padding = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def add(
        self,
        child: tkinter.Widget,
        *,
        state: Literal["normal", "disabled", "hidden"] = ...,
        sticky: str = ...,  # consists of letters 'n', 's', 'w', 'e', no repeats, may be empty
        padding: tkinter._Padding = ...,
        text: str = ...,
        image: Any = ...,  # Sequence of an image name, followed by zero or more (sequences of one or more state names followed by an image name)
        compound: tkinter._Compound = ...,
        underline: int = ...,
    ) -> None: ...
    def forget(self, tab_id): ...
    def hide(self, tab_id): ...
    def identify(self, x, y): ...
    def index(self, tab_id): ...
    def insert(self, pos, child, **kw): ...
    def select(self, tab_id: Any | None = ...): ...
    def tab(self, tab_id, option: Any | None = ..., **kw): ...
    def tabs(self): ...
    def enable_traversal(self) -> None: ...

class Panedwindow(Widget, tkinter.PanedWindow):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        # width and height for tkinter.ttk.Panedwindow are int but for tkinter.PanedWindow they are screen units
        height: int = ...,
        name: str = ...,
        orient: Literal["vertical", "horizontal"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> None: ...
    def add(self, child: tkinter.Widget, *, weight: int = ..., **kw) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Panedwindow().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        height: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        width: int = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def config(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    forget: Any
    def insert(self, pos, child, **kw): ...
    def pane(self, pane, option: Any | None = ..., **kw): ...
    def sashpos(self, index, newpos: Any | None = ...): ...

PanedWindow = Panedwindow

class Progressbar(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        length: tkinter._ScreenUnits = ...,
        maximum: float = ...,
        mode: Literal["determinate", "indeterminate"] = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        phase: int = ...,  # docs say read-only but assigning int to this works
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        value: float = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        length: tkinter._ScreenUnits = ...,
        maximum: float = ...,
        mode: Literal["determinate", "indeterminate"] = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        phase: int = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        value: float = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def start(self, interval: Any | None = ...): ...
    def step(self, amount: Any | None = ...): ...
    def stop(self): ...

class Radiobutton(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        name: str = ...,
        padding: Any = ...,  # undocumented
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        value: Any = ...,
        variable: tkinter.Variable | Literal[""] = ...,
        width: int | Literal[""] = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: tkinter._ButtonCommand = ...,
        compound: _TtkCompound = ...,
        cursor: tkinter._Cursor = ...,
        image: tkinter._ImageSpec = ...,
        padding: Any = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        text: float | str = ...,
        textvariable: tkinter.Variable = ...,
        underline: int = ...,
        value: Any = ...,
        variable: tkinter.Variable | Literal[""] = ...,
        width: int | Literal[""] = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def invoke(self) -> Any: ...

class Scale(Widget, tkinter.Scale):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        command: str | Callable[[str], Any] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,  # undocumented
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: str | Callable[[str], Any] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Scale().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: str | Callable[[str], Any] = ...,
        cursor: tkinter._Cursor = ...,
        from_: float = ...,
        length: tkinter._ScreenUnits = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        state: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        to: float = ...,
        value: float = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def config(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    def get(self, x: Any | None = ..., y: Any | None = ...): ...

class Scrollbar(Widget, tkinter.Scrollbar):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        command: Callable[..., tuple[float, float] | None] | str = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload  # type: ignore
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: Callable[..., tuple[float, float] | None] | str = ...,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    # config must be copy/pasted, otherwise ttk.Scrollbar().config is mypy error (don't know why)
    @overload  # type: ignore
    def config(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        command: Callable[..., tuple[float, float] | None] | str = ...,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def config(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...

class Separator(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        orient: Literal["horizontal", "vertical"] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

class Sizegrip(Widget):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        cursor: tkinter._Cursor = ...,
        name: str = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        cursor: tkinter._Cursor = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure

if sys.version_info >= (3, 7):
    class Spinbox(Entry):
        def __init__(
            self,
            master: tkinter.Misc | None = ...,
            *,
            background: tkinter._Color = ...,  # undocumented
            class_: str = ...,
            command: Callable[[], Any] | str | list[str] | Tuple[str, ...] = ...,
            cursor: tkinter._Cursor = ...,
            exportselection: bool = ...,  # undocumented
            font: _FontDescription = ...,  # undocumented
            foreground: tkinter._Color = ...,  # undocumented
            format: str = ...,
            from_: float = ...,
            increment: float = ...,
            invalidcommand: tkinter._EntryValidateCommand = ...,  # undocumented
            justify: Literal["left", "center", "right"] = ...,  # undocumented
            name: str = ...,
            show: Any = ...,  # undocumented
            state: str = ...,
            style: str = ...,
            takefocus: tkinter._TakeFocusValue = ...,
            textvariable: tkinter.Variable = ...,  # undocumented
            to: float = ...,
            validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
            validatecommand: tkinter._EntryValidateCommand = ...,
            values: list[str] | Tuple[str, ...] = ...,
            width: int = ...,  # undocumented
            wrap: bool = ...,
            xscrollcommand: tkinter._XYScrollCommand = ...,
        ) -> None: ...
        @overload  # type: ignore
        def configure(
            self,
            cnf: dict[str, Any] | None = ...,
            *,
            background: tkinter._Color = ...,
            command: Callable[[], Any] | str | list[str] | Tuple[str, ...] = ...,
            cursor: tkinter._Cursor = ...,
            exportselection: bool = ...,
            font: _FontDescription = ...,
            foreground: tkinter._Color = ...,
            format: str = ...,
            from_: float = ...,
            increment: float = ...,
            invalidcommand: tkinter._EntryValidateCommand = ...,
            justify: Literal["left", "center", "right"] = ...,
            show: Any = ...,
            state: str = ...,
            style: str = ...,
            takefocus: tkinter._TakeFocusValue = ...,
            textvariable: tkinter.Variable = ...,
            to: float = ...,
            validate: Literal["none", "focus", "focusin", "focusout", "key", "all"] = ...,
            validatecommand: tkinter._EntryValidateCommand = ...,
            values: list[str] | Tuple[str, ...] = ...,
            width: int = ...,
            wrap: bool = ...,
            xscrollcommand: tkinter._XYScrollCommand = ...,
        ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
        @overload
        def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
        config = configure  # type: ignore
        def set(self, value: Any) -> None: ...

class _TreeviewItemDict(TypedDict):
    text: str
    image: list[str] | Literal[""]  # no idea why it's wrapped in list
    values: list[Any] | Literal[""]
    open: bool  # actually 0 or 1
    tags: list[str] | Literal[""]

class _TreeviewTagDict(TypedDict):
    # There is also 'text' and 'anchor', but they don't seem to do anything, using them is likely a bug
    foreground: tkinter._Color
    background: tkinter._Color
    font: _FontDescription
    image: Literal[""] | str  # not wrapped in list :D

class _TreeviewHeaderDict(TypedDict):
    text: str
    image: list[str]
    anchor: tkinter._Anchor
    command: str
    state: str  # Doesn't seem to appear anywhere else than in these dicts

class _TreeviewColumnDict(TypedDict):
    width: int
    minwidth: int
    stretch: bool  # actually 0 or 1
    anchor: tkinter._Anchor
    id: str

_TreeviewColumnId = Union[int, str]  # manual page: "COLUMN IDENTIFIERS"

class Treeview(Widget, tkinter.XView, tkinter.YView):
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        *,
        class_: str = ...,
        columns: str | list[str] | Tuple[str, ...] = ...,
        cursor: tkinter._Cursor = ...,
        displaycolumns: str | list[str] | Tuple[str, ...] | list[int] | Tuple[int, ...] | Literal["#all"] = ...,
        height: int = ...,
        name: str = ...,
        padding: tkinter._Padding = ...,
        selectmode: Literal["extended", "browse", "none"] = ...,
        # list/tuple of Literal don't actually work in mypy
        #
        # 'tree headings' is same as ['tree', 'headings'], and I wouldn't be
        # surprised if someone is using it.
        show: Literal["tree", "headings", "tree headings", ""] | list[str] | Tuple[str, ...] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
        yscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> None: ...
    @overload
    def configure(
        self,
        cnf: dict[str, Any] | None = ...,
        *,
        columns: str | list[str] | Tuple[str, ...] = ...,
        cursor: tkinter._Cursor = ...,
        displaycolumns: str | list[str] | Tuple[str, ...] | list[int] | Tuple[int, ...] | Literal["#all"] = ...,
        height: int = ...,
        padding: tkinter._Padding = ...,
        selectmode: Literal["extended", "browse", "none"] = ...,
        show: Literal["tree", "headings", "tree headings", ""] | list[str] | Tuple[str, ...] = ...,
        style: str = ...,
        takefocus: tkinter._TakeFocusValue = ...,
        xscrollcommand: tkinter._XYScrollCommand = ...,
        yscrollcommand: tkinter._XYScrollCommand = ...,
    ) -> dict[str, tuple[str, str, str, Any, Any]] | None: ...
    @overload
    def configure(self, cnf: str) -> tuple[str, str, str, Any, Any]: ...
    config = configure
    def bbox(self, item, column: _TreeviewColumnId | None = ...) -> tuple[int, int, int, int] | Literal[""]: ...  # type: ignore
    def get_children(self, item: str | None = ...) -> Tuple[str, ...]: ...
    def set_children(self, item: str, *newchildren: str) -> None: ...
    @overload
    def column(self, column: _TreeviewColumnId, option: Literal["width", "minwidth"]) -> int: ...
    @overload
    def column(self, column: _TreeviewColumnId, option: Literal["stretch"]) -> bool: ...  # actually 0 or 1
    @overload
    def column(self, column: _TreeviewColumnId, option: Literal["anchor"]) -> _tkinter.Tcl_Obj: ...
    @overload
    def column(self, column: _TreeviewColumnId, option: Literal["id"]) -> str: ...
    @overload
    def column(self, column: _TreeviewColumnId, option: str) -> Any: ...
    @overload
    def column(
        self,
        column: _TreeviewColumnId,
        option: None = ...,
        *,
        width: int = ...,
        minwidth: int = ...,
        stretch: bool = ...,
        anchor: tkinter._Anchor = ...,
        # id is read-only
    ) -> _TreeviewColumnDict | None: ...
    def delete(self, *items: str) -> None: ...
    def detach(self, *items: str) -> None: ...
    def exists(self, item: str) -> bool: ...
    @overload  # type: ignore
    def focus(self, item: None = ...) -> str: ...  # can return empty string
    @overload
    def focus(self, item: str) -> Literal[""]: ...
    @overload
    def heading(self, column: _TreeviewColumnId, option: Literal["text"]) -> str: ...
    @overload
    def heading(self, column: _TreeviewColumnId, option: Literal["image"]) -> tuple[str]: ...
    @overload
    def heading(self, column: _TreeviewColumnId, option: Literal["anchor"]) -> _tkinter.Tcl_Obj: ...
    @overload
    def heading(self, column: _TreeviewColumnId, option: Literal["command"]) -> str: ...
    @overload
    def heading(self, column: _TreeviewColumnId, option: str) -> Any: ...
    @overload
    def heading(
        self,
        column: _TreeviewColumnId,
        option: None = ...,
        *,
        text: str = ...,
        image: tkinter._ImageSpec = ...,
        anchor: tkinter._Anchor = ...,
        command: str | Callable[[], Any] = ...,
    ) -> _TreeviewHeaderDict | None: ...
    def identify(self, component, x, y): ...
    def identify_row(self, y: int) -> str: ...
    def identify_column(self, x: int) -> str: ...
    def identify_region(self, x: int, y: int) -> Literal["heading", "separator", "tree", "cell", "nothing"]: ...
    def identify_element(self, x: int, y: int) -> str: ...  # don't know what possible return values are
    def index(self, item: str) -> int: ...
    def insert(
        self,
        parent: str,
        index: int | Literal["end"],
        iid: str | None = ...,
        *,
        id: str = ...,  # same as iid
        text: str = ...,
        image: tkinter._ImageSpec = ...,
        values: list[Any] | Tuple[Any, ...] = ...,
        open: bool = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
    ) -> str: ...
    @overload
    def item(self, item: str, option: Literal["text"]) -> str: ...
    @overload
    def item(self, item: str, option: Literal["image"]) -> tuple[str] | Literal[""]: ...
    @overload
    def item(self, item: str, option: Literal["values"]) -> Tuple[Any, ...] | Literal[""]: ...
    @overload
    def item(self, item: str, option: Literal["open"]) -> bool: ...  # actually 0 or 1
    @overload
    def item(self, item: str, option: Literal["tags"]) -> Tuple[str, ...] | Literal[""]: ...
    @overload
    def item(self, item: str, option: str) -> Any: ...
    @overload
    def item(
        self,
        item: str,
        option: None = ...,
        *,
        text: str = ...,
        image: tkinter._ImageSpec = ...,
        values: list[Any] | Tuple[Any, ...] | Literal[""] = ...,
        open: bool = ...,
        tags: str | list[str] | Tuple[str, ...] = ...,
    ) -> _TreeviewItemDict | None: ...
    def move(self, item: str, parent: str, index: int) -> None: ...
    reattach = move
    def next(self, item: str) -> str: ...  # returning empty string means last item
    def parent(self, item: str) -> str: ...
    def prev(self, item: str) -> str: ...  # returning empty string means first item
    def see(self, item: str) -> None: ...
    if sys.version_info >= (3, 8):
        def selection(self) -> Tuple[str, ...]: ...
    else:
        def selection(self, selop: Any | None = ..., items: Any | None = ...) -> Tuple[str, ...]: ...
    def selection_set(self, items: str | list[str] | Tuple[str, ...]) -> None: ...
    def selection_add(self, items: str | list[str] | Tuple[str, ...]) -> None: ...
    def selection_remove(self, items: str | list[str] | Tuple[str, ...]) -> None: ...
    def selection_toggle(self, items: str | list[str] | Tuple[str, ...]) -> None: ...
    @overload
    def set(self, item: str, column: None = ..., value: None = ...) -> dict[str, Any]: ...
    @overload
    def set(self, item: str, column: _TreeviewColumnId, value: None = ...) -> Any: ...
    @overload
    def set(self, item: str, column: _TreeviewColumnId, value: Any) -> Literal[""]: ...
    # There's no tag_unbind() or 'add' argument for whatever reason.
    # Also, it's 'callback' instead of 'func' here.
    @overload
    def tag_bind(
        self, tagname: str, sequence: str | None = ..., callback: Callable[[tkinter.Event[Treeview]], Any] | None = ...
    ) -> str: ...
    @overload
    def tag_bind(self, tagname: str, sequence: str | None, callback: str) -> None: ...
    @overload
    def tag_bind(self, tagname: str, *, callback: str) -> None: ...
    @overload
    def tag_configure(self, tagname: str, option: Literal["foreground", "background"]) -> tkinter._Color: ...
    @overload
    def tag_configure(self, tagname: str, option: Literal["font"]) -> _FontDescription: ...
    @overload
    def tag_configure(self, tagname: str, option: Literal["image"]) -> str: ...
    @overload
    def tag_configure(
        self,
        tagname: str,
        option: None = ...,
        *,
        # There is also 'text' and 'anchor', but they don't seem to do anything, using them is likely a bug
        foreground: tkinter._Color = ...,
        background: tkinter._Color = ...,
        font: _FontDescription = ...,
        image: tkinter._ImageSpec = ...,
    ) -> _TreeviewTagDict | Any: ...  # can be None but annoying to check
    @overload
    def tag_has(self, tagname: str, item: None = ...) -> Tuple[str, ...]: ...
    @overload
    def tag_has(self, tagname: str, item: str) -> bool: ...

class LabeledScale(Frame):
    label: Any
    scale: Any
    # TODO: don't any-type **kw. That goes to Frame.__init__.
    def __init__(
        self,
        master: tkinter.Misc | None = ...,
        variable: tkinter.IntVar | tkinter.DoubleVar | None = ...,
        from_: float = ...,
        to: float = ...,
        *,
        compound: Literal["top"] | Literal["bottom"] = ...,
        **kw: Any,
    ) -> None: ...
    # destroy is overridden, signature does not change
    value: Any

class OptionMenu(Menubutton):
    def __init__(
        self,
        master,
        variable,
        default: str | None = ...,
        *values: str,
        # rest of these are keyword-only because *args syntax used above
        style: str = ...,
        direction: Literal["above"] | Literal["below"] | Literal["left"] | Literal["right"] | Literal["flush"] = ...,
        command: Callable[[tkinter.StringVar], Any] | None = ...,
    ) -> None: ...
    # configure, config, cget, destroy are inherited from Menubutton
    # destroy and __setitem__ are overridden, signature does not change
    def set_menu(self, default: Any | None = ..., *values): ...
