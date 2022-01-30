from codecs import CodecInfo
from typing import Any

class CodecRegistryError(LookupError, SystemError): ...

def normalize_encoding(encoding: str | bytes) -> str: ...
def search_function(encoding: str) -> CodecInfo | None: ...

# Needed for submodules
def __getattr__(name: str) -> Any: ...  # incomplete
