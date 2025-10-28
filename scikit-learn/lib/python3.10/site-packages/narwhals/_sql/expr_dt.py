from __future__ import annotations

from typing import Any, Generic

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import DateTimeNamespace
from narwhals._sql.typing import SQLExprT


class SQLExprDateTimeNamesSpace(
    LazyExprNamespace[SQLExprT], DateTimeNamespace[SQLExprT], Generic[SQLExprT]
):
    def _function(self, name: str, *args: Any) -> SQLExprT:
        return self.compliant._function(name, *args)  # type: ignore[no-any-return]

    def year(self) -> SQLExprT:
        return self.compliant._with_elementwise(lambda expr: self._function("year", expr))

    def month(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("month", expr)
        )

    def day(self) -> SQLExprT:
        return self.compliant._with_elementwise(lambda expr: self._function("day", expr))

    def hour(self) -> SQLExprT:
        return self.compliant._with_elementwise(lambda expr: self._function("hour", expr))

    def minute(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("minute", expr)
        )

    def second(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("second", expr)
        )

    def ordinal_day(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("dayofyear", expr)
        )

    def date(self) -> SQLExprT:
        return self.compliant._with_elementwise(
            lambda expr: self._function("to_date", expr)
        )
