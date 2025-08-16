from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant.any_namespace import CatNamespace
from narwhals._pandas_like.utils import PandasLikeSeriesNamespace

if TYPE_CHECKING:
    from narwhals._pandas_like.series import PandasLikeSeries


class PandasLikeSeriesCatNamespace(
    PandasLikeSeriesNamespace, CatNamespace["PandasLikeSeries"]
):
    def get_categories(self) -> PandasLikeSeries:
        s = self.native
        return self.with_native(type(s)(s.cat.categories, name=s.name))
