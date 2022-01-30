from tkinter.commondialog import Dialog
from typing import Any, ClassVar

class Chooser(Dialog):
    command: ClassVar[str]

def askcolor(color: str | bytes | None = ..., **options: Any) -> tuple[None, None] | tuple[tuple[float, float, float], str]: ...
