from typing import BinaryIO, Optional, Dict, Any, Sequence, overload, Literal
from ..typing import ImageResource
import enum

EXAMPLE_IMAGES: Dict[str, str]
RETURN_BYTES = "<bytes>"
URI_BYTES = 1
URI_FILE = 2
URI_FILENAME = 3
URI_ZIPPED = 4
URI_HTTP = 5
URI_FTP = 6

class IOMode(str, enum.Enum):
    read = "r"
    write = "w"

class ImageMode(str, enum.Enum):
    single_image = "i"
    multi_image = "I"
    single_volume = "v"
    multi_volume = "V"
    any_mode = "?"

@enum.unique
class Mode(str, enum.Enum):
    read_single_image = "ri"
    read_multi_image = "rI"
    read_single_volume = "rv"
    read_multi_volume = "rV"
    read_any = "r?"
    write_single_image = "wi"
    write_multi_image = "wI"
    write_single_volume = "wv"
    write_multi_volume = "wV"
    write_any = "w?"

    @classmethod
    def _missing_(cls, value: Any) -> Mode: ...
    @property
    def io_mode(self) -> IOMode: ...
    @property
    def image_mode(self) -> ImageMode: ...

class InitializationError(Exception): ...

class Request(object):
    _uri_type: int
    raw_uri: ImageResource

    @property
    def filename(self) -> str: ...
    @property
    def extension(self) -> str: ...
    @property
    def format_hint(self) -> Optional[str]: ...
    @format_hint.setter
    def format_hint(self, format: str) -> None: ...
    @property
    def mode(self) -> Mode: ...
    @property
    def kwargs(self) -> Dict[str, Any]: ...
    @property
    def firstbytes(self) -> bytes: ...
    def __init__(
        self,
        uri: ImageResource,
        mode: str,
        *,
        extension: str = None,
        format_hint: str = None,
        **kwargs
    ) -> None: ...
    def _parse_uri(self, uri: ImageResource) -> None: ...
    def get_file(self) -> BinaryIO: ...
    def get_local_filename(self) -> str: ...
    def finish(self) -> None: ...
    def get_result(self) -> Optional[bytes]: ...
    def _read_first_bytes(self, N: int = 256) -> bytes: ...

def read_n_bytes(f: BinaryIO, N: int) -> bytes: ...

class SeekableFileObject:
    def __init__(self, f: BinaryIO) -> None: ...
    def read(self, n: int = None) -> bytes: ...
    def tell(self) -> int: ...
    def seek(self, i: int, mode: int = 0) -> int: ...
    def close(self) -> None: ...
    def isatty(self) -> bool: ...
    def seekable(self) -> bool: ...
