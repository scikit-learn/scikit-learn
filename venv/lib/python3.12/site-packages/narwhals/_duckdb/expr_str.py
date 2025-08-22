from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._duckdb.utils import F, lit
from narwhals._sql.expr_str import SQLExprStringNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
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

    replace = not_implemented()
