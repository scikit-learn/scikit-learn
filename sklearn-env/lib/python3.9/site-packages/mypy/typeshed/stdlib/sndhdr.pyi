from _typeshed import StrOrBytesPath
from typing import NamedTuple

class SndHeaders(NamedTuple):
    filetype: str
    framerate: int
    nchannels: int
    nframes: int
    sampwidth: int | str

def what(filename: StrOrBytesPath) -> SndHeaders | None: ...
def whathdr(filename: StrOrBytesPath) -> SndHeaders | None: ...
