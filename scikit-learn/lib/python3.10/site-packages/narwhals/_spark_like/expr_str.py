from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

from narwhals._spark_like.utils import strptime_to_pyspark_format
from narwhals._sql.expr_str import SQLExprStringNamespace
from narwhals._utils import _is_naive_format, not_implemented, requires

if TYPE_CHECKING:
    from sqlframe.base.column import Column

    from narwhals._spark_like.expr import SparkLikeExpr


class SparkLikeExprStringNamespace(SQLExprStringNamespace["SparkLikeExpr"]):
    def to_datetime(self, format: str | None) -> SparkLikeExpr:
        F = self.compliant._F
        if not format:
            function = F.to_timestamp
        elif _is_naive_format(format):
            function = partial(
                F.to_timestamp_ntz, format=F.lit(strptime_to_pyspark_format(format))
            )
        else:
            format = strptime_to_pyspark_format(format)
            function = partial(F.to_timestamp, format=format)
        return self.compliant._with_elementwise(
            lambda expr: function(F.replace(expr, F.lit("T"), F.lit(" ")))
        )

    def to_date(self, format: str | None) -> SparkLikeExpr:
        F = self.compliant._F
        return self.compliant._with_elementwise(
            lambda expr: F.to_date(expr, format=strptime_to_pyspark_format(format))
        )

    def to_titlecase(self) -> SparkLikeExpr:
        impl = self.compliant._implementation
        sqlframe_required_version = (3, 43, 1)
        if (
            impl.is_sqlframe()
            and (version := impl._backend_version()) < sqlframe_required_version
        ):  # pragma: no cover
            required_str = requires._unparse_version(sqlframe_required_version)
            found_str = requires._unparse_version(version)
            msg = (
                f"`str.to_titlecase` is only available in 'sqlframe>={required_str}', "
                f"found version {found_str!r}."
            )
            raise NotImplementedError(msg)

        def _to_titlecase(expr: Column) -> Column:
            F = self.compliant._F
            lower_expr = F.lower(expr)
            extract_expr = F.regexp_extract_all(
                lower_expr, regexp=F.lit(r"[a-z0-9]*[^a-z0-9]*"), idx=0
            )
            capitalized_expr = F.transform(extract_expr, f=F.initcap)
            return F.array_join(capitalized_expr, delimiter="")

        return self.compliant._with_elementwise(_to_titlecase)

    replace = not_implemented()
