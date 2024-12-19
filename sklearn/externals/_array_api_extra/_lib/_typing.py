from __future__ import annotations  # https://github.com/pylint-dev/pylint/pull/9990

import typing
from types import ModuleType
from typing import Any

if typing.TYPE_CHECKING:
    from typing_extensions import override

    # To be changed to a Protocol later (see data-apis/array-api#589)
    Array = Any  # type: ignore[no-any-explicit]
    Device = Any  # type: ignore[no-any-explicit]
else:

    def no_op_decorator(f):  # pyright: ignore[reportUnreachable]
        return f

    override = no_op_decorator

__all__ = ["ModuleType", "override"]
if typing.TYPE_CHECKING:
    __all__ += ["Array", "Device"]
