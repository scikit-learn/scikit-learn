"""Static typing helpers — better implementations in the stub file."""

from collections.abc import Callable
from types import ModuleType
from typing import Any

Array = object
ArrayLike = object
ArrayNamespace = ModuleType
DType = object
Device = object
Key = object
GetIndex = object
Graph = object
NumPyObject = object
SchedulerGetCallable = object
SetIndex = object

TypeIs = Any

__all__ = [
    "Array",
    "ArrayLike",
    "ArrayNamespace",
    "DType",
    "Device",
    "GetIndex",
    "Graph",
    "Key",
    "NumPyObject",
    "SchedulerGetCallable",
    "SetIndex",
    "TypeIs",
    "override",
]


def override(func: Callable[..., Any]) -> Callable[..., Any]:  # numpydoc ignore=GL08
    return func
