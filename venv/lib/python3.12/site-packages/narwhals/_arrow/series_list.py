from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import ArrowSeriesNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries


class ArrowSeriesListNamespace(ArrowSeriesNamespace):
    def len(self) -> ArrowSeries:
        return self.with_native(pc.list_value_length(self.native).cast(pa.uint32()))

    unique = not_implemented()
    contains = not_implemented()
