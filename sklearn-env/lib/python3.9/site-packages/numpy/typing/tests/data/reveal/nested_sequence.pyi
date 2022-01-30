from typing import Sequence, Tuple, List, Any
import numpy.typing as npt

a: Sequence[int]
b: Sequence[Sequence[int]]
c: Sequence[Sequence[Sequence[int]]]
d: Sequence[Sequence[Sequence[Sequence[int]]]]
e: Sequence[bool]
f: Tuple[int, ...]
g: List[int]
h: Sequence[Any]

def func(a: npt._NestedSequence[int]) -> None:
    ...

reveal_type(func(a))  # E: None
reveal_type(func(b))  # E: None
reveal_type(func(c))  # E: None
reveal_type(func(d))  # E: None
reveal_type(func(e))  # E: None
reveal_type(func(f))  # E: None
reveal_type(func(g))  # E: None
reveal_type(func(h))  # E: None
reveal_type(func(range(15)))  # E: None
