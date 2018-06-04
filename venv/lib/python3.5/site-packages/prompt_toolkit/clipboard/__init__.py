from .base import Clipboard, ClipboardData
from .in_memory import InMemoryClipboard


# We are not importing `PyperclipClipboard` here, because it would require the
# `pyperclip` module to be present.

#from .pyperclip import PyperclipClipboard
