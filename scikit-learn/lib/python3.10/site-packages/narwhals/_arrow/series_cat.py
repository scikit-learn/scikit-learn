from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow as pa

from narwhals._arrow.utils import ArrowSeriesNamespace
from narwhals._compliant.any_namespace import CatNamespace

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries
    from narwhals._arrow.typing import Incomplete


class ArrowSeriesCatNamespace(ArrowSeriesNamespace, CatNamespace["ArrowSeries"]):
    def get_categories(self) -> ArrowSeries:
        # NOTE: Should be `list[pa.DictionaryArray]`, but `DictionaryArray` has no attributes
        chunks: Incomplete = self.native.chunks
        return self.with_native(pa.concat_arrays(x.dictionary for x in chunks).unique())
