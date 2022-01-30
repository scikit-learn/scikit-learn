from typing import Any

def __getattr__(name: str) -> Any: ...  # incomplete

class DocumentLS(Any): ...  # type: ignore
class DOMImplementationLS(Any): ...  # type: ignore
