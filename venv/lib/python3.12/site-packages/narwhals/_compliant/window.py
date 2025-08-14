from __future__ import annotations

from typing import TYPE_CHECKING, Generic

from narwhals._compliant.typing import NativeExprT_co

if TYPE_CHECKING:
    from collections.abc import Sequence

__all__ = ["WindowInputs"]


class WindowInputs(Generic[NativeExprT_co]):
    __slots__ = ("order_by", "partition_by")

    def __init__(
        self, partition_by: Sequence[str | NativeExprT_co], order_by: Sequence[str]
    ) -> None:
        self.partition_by = partition_by
        self.order_by = order_by
