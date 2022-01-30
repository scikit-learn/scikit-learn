from typing import Sequence, Tuple, List
import numpy.typing as npt

a: Sequence[float]
b: List[complex]
c: Tuple[str, ...]
d: int
e: str

def func(a: npt._NestedSequence[int]) -> None:
    ...

reveal_type(func(a))  # E: incompatible type
reveal_type(func(b))  # E: incompatible type
reveal_type(func(c))  # E: incompatible type
reveal_type(func(d))  # E: incompatible type
reveal_type(func(e))  # E: incompatible type
