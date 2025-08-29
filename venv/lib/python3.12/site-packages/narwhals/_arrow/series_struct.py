from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow.compute as pc

from narwhals._arrow.utils import ArrowSeriesNamespace

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries


class ArrowSeriesStructNamespace(ArrowSeriesNamespace):
    def field(self, name: str) -> ArrowSeries:
        return self.with_native(pc.struct_field(self.native, name)).alias(name)
