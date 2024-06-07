from collections.abc import Sequence
from typing import SupportsIndex, Union

_Shape = tuple[int, ...]

# Anything that can be coerced to a shape tuple
_ShapeLike = Union[SupportsIndex, Sequence[SupportsIndex]]
