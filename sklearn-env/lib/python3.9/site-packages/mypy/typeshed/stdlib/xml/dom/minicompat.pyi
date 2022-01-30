from typing import Any, Iterable, List, Tuple, Type, TypeVar

_T = TypeVar("_T")

StringTypes: tuple[Type[str]]

class NodeList(List[_T]):
    length: int
    def item(self, index: int) -> _T | None: ...

class EmptyNodeList(Tuple[Any, ...]):
    length: int
    def item(self, index: int) -> None: ...
    def __add__(self, other: Iterable[_T]) -> NodeList[_T]: ...  # type: ignore
    def __radd__(self, other: Iterable[_T]) -> NodeList[_T]: ...

def defproperty(klass: Type[Any], name: str, doc: str) -> None: ...
