from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._duckdb.utils import F, col, concat_str, lit
from narwhals._sql.expr_str import SQLExprStringNamespace
from narwhals._utils import not_implemented, requires

if TYPE_CHECKING:
    from duckdb import Expression

    from narwhals._duckdb.expr import DuckDBExpr


class DuckDBExprStringNamespace(SQLExprStringNamespace["DuckDBExpr"]):
    def to_datetime(self, format: str | None) -> DuckDBExpr:
        if format is None:
            msg = "Cannot infer format with DuckDB backend, please specify `format` explicitly."
            raise NotImplementedError(msg)

        return self.compliant._with_elementwise(
            lambda expr: F("strptime", expr, lit(format))
        )

    def to_date(self, format: str | None) -> DuckDBExpr:
        if format is not None:
            return self.to_datetime(format=format).dt.date()

        compliant_expr = self.compliant
        return compliant_expr.cast(compliant_expr._version.dtypes.Date())

    @requires.backend_version((1, 2))
    def to_titlecase(self) -> DuckDBExpr:
        from narwhals._duckdb.utils import lambda_expr

        def _to_titlecase(expr: Expression) -> Expression:
            extract_expr = F(
                "regexp_extract_all", F("lower", expr), lit(r"[a-z0-9]*[^a-z0-9]*")
            )
            elem = col("_")
            capitalize = lambda_expr(
                elem,
                concat_str(
                    F("upper", F("array_extract", elem, lit(1))),
                    F("substring", elem, lit(2)),
                ),
            )
            capitalized_expr = F("list_transform", extract_expr, capitalize)
            return F("list_aggregate", capitalized_expr, lit("string_agg"), lit(""))

        return self.compliant._with_elementwise(_to_titlecase)

    replace = not_implemented()
