from tkinter import Frame, Misc, Scrollbar, Text
from typing import Any

# The methods from Pack, Place, and Grid are dynamically added over the parent's impls
class ScrolledText(Text):
    frame: Frame
    vbar: Scrollbar
    def __init__(self, master: Misc | None = ..., **kwargs: Any) -> None: ...
