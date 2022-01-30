from _typeshed.xml import DOMImplementation
from typing import Callable, Dict, Iterable, Tuple

well_known_implementations: Dict[str, str]
registered: Dict[str, Callable[[], DOMImplementation]]

def registerDOMImplementation(name: str, factory: Callable[[], DOMImplementation]) -> None: ...
def getDOMImplementation(name: str | None = ..., features: str | Iterable[Tuple[str, str | None]] = ...) -> DOMImplementation: ...
