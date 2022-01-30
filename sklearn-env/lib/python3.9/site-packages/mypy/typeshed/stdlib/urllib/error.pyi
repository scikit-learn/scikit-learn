from email.message import Message
from typing import IO
from urllib.response import addinfourl

# Stubs for urllib.error

class URLError(IOError):
    reason: str | BaseException
    def __init__(self, reason: str | BaseException, filename: str | None = ...) -> None: ...

class HTTPError(URLError, addinfourl):
    code: int
    def __init__(self, url: str, code: int, msg: str, hdrs: Message, fp: IO[bytes] | None) -> None: ...

class ContentTooShortError(URLError):
    content: tuple[str, Message]
    def __init__(self, message: str, content: tuple[str, Message]) -> None: ...
