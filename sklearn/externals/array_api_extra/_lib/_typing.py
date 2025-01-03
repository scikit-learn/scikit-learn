"""Static typing helpers."""

from types import ModuleType
from typing import Any

# To be changed to a Protocol later (see data-apis/array-api#589)
Array = Any  # type: ignore[no-any-explicit]
Device = Any  # type: ignore[no-any-explicit]
Index = Any  # type: ignore[no-any-explicit]

__all__ = ["Array", "Device", "Index", "ModuleType"]
