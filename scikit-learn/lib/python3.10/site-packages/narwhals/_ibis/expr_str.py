from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable

from ibis.expr.datatypes import Timestamp

from narwhals._sql.expr_str import SQLExprStringNamespace
from narwhals._utils import _is_naive_format, not_implemented

if TYPE_CHECKING:
    import ibis.expr.types as ir
    from typing_extensions import TypeAlias

    from narwhals._ibis.expr import IbisExpr

IntoStringValue: TypeAlias = "str | ir.StringValue"


class IbisExprStringNamespace(SQLExprStringNamespace["IbisExpr"]):
    def strip_chars(self, characters: str | None) -> IbisExpr:
        if characters is not None:
            msg = "Ibis does not support `characters` argument in `str.strip_chars`"
            raise NotImplementedError(msg)

        return self.compliant._with_callable(lambda expr: expr.strip())

    def _replace_all(
        self, pattern: IntoStringValue, value: IntoStringValue
    ) -> Callable[..., ir.StringValue]:
        def fn(expr: ir.StringColumn) -> ir.StringValue:
            return expr.re_replace(pattern, value)

        return fn

    def _replace_all_literal(
        self, pattern: IntoStringValue, value: IntoStringValue
    ) -> Callable[..., ir.StringValue]:
        def fn(expr: ir.StringColumn) -> ir.StringValue:
            return expr.replace(pattern, value)  # pyright: ignore[reportArgumentType]

        return fn

    def replace_all(
        self, pattern: str, value: str | IbisExpr, *, literal: bool
    ) -> IbisExpr:
        fn = self._replace_all_literal if literal else self._replace_all
        if isinstance(value, str):
            return self.compliant._with_callable(fn(pattern, value))
        return self.compliant._with_elementwise(
            lambda expr, value: fn(pattern, value)(expr), value=value
        )

    def _to_datetime(self, format: str) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.StringColumn) -> ir.TimestampValue:
            return expr.as_timestamp(format)

        return fn

    def _to_datetime_naive(self, format: str) -> Callable[..., ir.TimestampValue]:
        def fn(expr: ir.StringColumn) -> ir.TimestampValue:
            dtype: Any = Timestamp(timezone=None)
            return expr.as_timestamp(format).cast(dtype)

        return fn

    def to_datetime(self, format: str | None) -> IbisExpr:
        if format is None:
            msg = "Cannot infer format with Ibis backend"
            raise NotImplementedError(msg)
        fn = self._to_datetime_naive if _is_naive_format(format) else self._to_datetime
        return self.compliant._with_callable(fn(format))

    def to_date(self, format: str | None) -> IbisExpr:
        if format is None:
            msg = "Cannot infer format with Ibis backend"
            raise NotImplementedError(msg)

        def fn(expr: ir.StringColumn) -> ir.DateValue:
            return expr.as_date(format)

        return self.compliant._with_callable(fn)

    replace = not_implemented()
    to_titlecase = not_implemented()
