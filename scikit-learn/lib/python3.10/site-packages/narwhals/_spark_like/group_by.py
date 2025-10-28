from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._sql.group_by import SQLGroupBy

if TYPE_CHECKING:
    from collections.abc import Sequence

    from sqlframe.base.column import Column  # noqa: F401

    from narwhals._spark_like.dataframe import SparkLikeLazyFrame
    from narwhals._spark_like.expr import SparkLikeExpr


class SparkLikeLazyGroupBy(SQLGroupBy["SparkLikeLazyFrame", "SparkLikeExpr", "Column"]):
    def __init__(
        self,
        df: SparkLikeLazyFrame,
        keys: Sequence[SparkLikeExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        frame, self._keys, self._output_key_names = self._parse_keys(df, keys=keys)
        self._compliant_frame = frame.drop_nulls(self._keys) if drop_null_keys else frame

    def agg(self, *exprs: SparkLikeExpr) -> SparkLikeLazyFrame:
        result = (
            self.compliant.native.groupBy(*self._keys).agg(*agg_columns)
            if (agg_columns := tuple(self._evaluate_exprs(exprs)))
            else self.compliant.native.select(*self._keys).dropDuplicates()
        )

        return self.compliant._with_native(result).rename(
            dict(zip(self._keys, self._output_key_names))
        )
