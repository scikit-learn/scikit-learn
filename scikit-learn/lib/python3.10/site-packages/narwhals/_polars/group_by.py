from __future__ import annotations

from typing import TYPE_CHECKING, cast

from narwhals._utils import is_sequence_of

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from polars.dataframe.group_by import GroupBy as NativeGroupBy
    from polars.lazyframe.group_by import LazyGroupBy as NativeLazyGroupBy

    from narwhals._polars.dataframe import PolarsDataFrame, PolarsLazyFrame
    from narwhals._polars.expr import PolarsExpr


class PolarsGroupBy:
    _compliant_frame: PolarsDataFrame
    _grouped: NativeGroupBy

    @property
    def compliant(self) -> PolarsDataFrame:
        return self._compliant_frame

    def __init__(
        self,
        df: PolarsDataFrame,
        keys: Sequence[PolarsExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._keys = list(keys)
        self._compliant_frame = df.drop_nulls(keys) if drop_null_keys else df
        self._grouped = (
            self.compliant.native.group_by(keys)
            if is_sequence_of(keys, str)
            else self.compliant.native.group_by(arg.native for arg in keys)
        )

    def agg(self, *aggs: PolarsExpr) -> PolarsDataFrame:
        agg_result = self._grouped.agg(arg.native for arg in aggs)
        return self.compliant._with_native(agg_result)

    def __iter__(self) -> Iterator[tuple[tuple[str, ...], PolarsDataFrame]]:
        for key, df in self._grouped:
            yield tuple(cast("str", key)), self.compliant._with_native(df)


class PolarsLazyGroupBy:
    _compliant_frame: PolarsLazyFrame
    _grouped: NativeLazyGroupBy

    @property
    def compliant(self) -> PolarsLazyFrame:
        return self._compliant_frame

    def __init__(
        self,
        df: PolarsLazyFrame,
        keys: Sequence[PolarsExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._keys = list(keys)
        self._compliant_frame = df.drop_nulls(keys) if drop_null_keys else df
        self._grouped = (
            self.compliant.native.group_by(keys)
            if is_sequence_of(keys, str)
            else self.compliant.native.group_by(arg.native for arg in keys)
        )

    def agg(self, *aggs: PolarsExpr) -> PolarsLazyFrame:
        agg_result = self._grouped.agg(arg.native for arg in aggs)
        return self.compliant._with_native(agg_result)
