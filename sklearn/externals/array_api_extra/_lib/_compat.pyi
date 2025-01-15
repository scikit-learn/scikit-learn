"""Static type stubs for `_compat.py`."""

# https://github.com/scikit-learn/scikit-learn/pull/27910#issuecomment-2568023972
from __future__ import annotations

from types import ModuleType

from ._typing import Array, Device

# pylint: disable=missing-class-docstring,unused-argument

class ArrayModule(ModuleType):
    def device(self, x: Array, /) -> Device: ...

def array_namespace(
    *xs: Array,
    api_version: str | None = None,
    use_compat: bool | None = None,
) -> ArrayModule: ...
def device(x: Array, /) -> Device: ...
def is_jax_array(x: object, /) -> bool: ...
def is_writeable_array(x: object, /) -> bool: ...
