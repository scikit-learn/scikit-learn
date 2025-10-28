from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._compliant import CompliantSelector, LazySelectorNamespace
from narwhals._spark_like.expr import SparkLikeExpr

if TYPE_CHECKING:
    from sqlframe.base.column import Column  # noqa: F401

    from narwhals._spark_like.dataframe import SparkLikeLazyFrame  # noqa: F401
    from narwhals._spark_like.expr import SparkWindowFunction


class SparkLikeSelectorNamespace(LazySelectorNamespace["SparkLikeLazyFrame", "Column"]):
    @property
    def _selector(self) -> type[SparkLikeSelector]:
        return SparkLikeSelector


class SparkLikeSelector(CompliantSelector["SparkLikeLazyFrame", "Column"], SparkLikeExpr):  # type: ignore[misc]
    _window_function: SparkWindowFunction | None = None

    def _to_expr(self) -> SparkLikeExpr:
        return SparkLikeExpr(
            self._call,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
            implementation=self._implementation,
        )
