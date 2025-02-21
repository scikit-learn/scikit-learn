from enum import Enum

class _AutoStringNameEnum(Enum):
    def __hash__(self) -> int: ...

class JoinStyle(str, _AutoStringNameEnum):
    miter: str
    round: str
    bevel: str
    @staticmethod
    def demo() -> None: ...

class CapStyle(str, _AutoStringNameEnum):
    butt: str
    projecting: str
    round: str
    @staticmethod
    def demo() -> None: ...
