from tkinter.commondialog import Dialog
from typing import Any, ClassVar

ERROR: str
INFO: str
QUESTION: str
WARNING: str
ABORTRETRYIGNORE: str
OK: str
OKCANCEL: str
RETRYCANCEL: str
YESNO: str
YESNOCANCEL: str
ABORT: str
RETRY: str
IGNORE: str
CANCEL: str
YES: str
NO: str

class Message(Dialog):
    command: ClassVar[str]

def showinfo(title: str | None = ..., message: str | None = ..., **options: Any) -> str: ...
def showwarning(title: str | None = ..., message: str | None = ..., **options: Any) -> str: ...
def showerror(title: str | None = ..., message: str | None = ..., **options: Any) -> str: ...
def askquestion(title: str | None = ..., message: str | None = ..., **options: Any) -> str: ...
def askokcancel(title: str | None = ..., message: str | None = ..., **options: Any) -> bool: ...
def askyesno(title: str | None = ..., message: str | None = ..., **options: Any) -> bool: ...
def askyesnocancel(title: str | None = ..., message: str | None = ..., **options: Any) -> bool | None: ...
def askretrycancel(title: str | None = ..., message: str | None = ..., **options: Any) -> bool: ...
