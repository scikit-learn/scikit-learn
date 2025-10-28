from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant.any_namespace import StructNamespace
from narwhals._pandas_like.utils import PandasLikeSeriesNamespace

if TYPE_CHECKING:
    from narwhals._pandas_like.series import PandasLikeSeries


class PandasLikeSeriesStructNamespace(
    PandasLikeSeriesNamespace, StructNamespace["PandasLikeSeries"]
):
    def field(self, name: str) -> PandasLikeSeries:
        return self.with_native(self.native.struct.field(name)).alias(name)
