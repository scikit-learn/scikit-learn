from __future__ import annotations

from typing import TYPE_CHECKING

import dask.dataframe as dd

from narwhals._compliant import LazyExprNamespace
from narwhals._compliant.any_namespace import StringNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    import dask.dataframe.dask_expr as dx

    from narwhals._dask.expr import DaskExpr


class DaskExprStringNamespace(LazyExprNamespace["DaskExpr"], StringNamespace["DaskExpr"]):
    def len_chars(self) -> DaskExpr:
        return self.compliant._with_callable(lambda expr: expr.str.len(), "len")

    def replace(self, pattern: str, value: str, *, literal: bool, n: int) -> DaskExpr:
        def _replace(
            expr: dx.Series, pattern: str, value: str, *, literal: bool, n: int
        ) -> dx.Series:
            try:
                return expr.str.replace(  # pyright: ignore[reportAttributeAccessIssue]
                    pattern, value, regex=not literal, n=n
                )
            except TypeError as e:
                if not isinstance(value, str):
                    msg = "dask backed `Expr.str.replace` only supports str replacement values"
                    raise TypeError(msg) from e
                raise

        return self.compliant._with_callable(
            _replace, "replace", pattern=pattern, value=value, literal=literal, n=n
        )

    def replace_all(self, pattern: str, value: str, *, literal: bool) -> DaskExpr:
        def _replace_all(
            expr: dx.Series, pattern: str, value: str, *, literal: bool
        ) -> dx.Series:
            try:
                return expr.str.replace(  # pyright: ignore[reportAttributeAccessIssue]
                    pattern, value, regex=not literal, n=-1
                )
            except TypeError as e:
                if not isinstance(value, str):
                    msg = "dask backed `Expr.str.replace_all` only supports str replacement values."
                    raise TypeError(msg) from e
                raise

        return self.compliant._with_callable(
            _replace_all, "replace", pattern=pattern, value=value, literal=literal
        )

    def strip_chars(self, characters: str | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, characters: expr.str.strip(characters),
            "strip",
            characters=characters,
        )

    def starts_with(self, prefix: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, prefix: expr.str.startswith(prefix), "starts_with", prefix=prefix
        )

    def ends_with(self, suffix: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, suffix: expr.str.endswith(suffix), "ends_with", suffix=suffix
        )

    def contains(self, pattern: str, *, literal: bool) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, pattern, literal: expr.str.contains(
                pat=pattern, regex=not literal
            ),
            "contains",
            pattern=pattern,
            literal=literal,
        )

    def slice(self, offset: int, length: int | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, offset, length: expr.str.slice(
                start=offset, stop=offset + length if length else None
            ),
            "slice",
            offset=offset,
            length=length,
        )

    def split(self, by: str) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, by: expr.str.split(pat=by), "split", by=by
        )

    def to_datetime(self, format: str | None) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, format: dd.to_datetime(expr, format=format),
            "to_datetime",
            format=format,
        )

    def to_uppercase(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.upper(), "to_uppercase"
        )

    def to_lowercase(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.lower(), "to_lowercase"
        )

    def to_titlecase(self) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr: expr.str.title(), "to_titlecase"
        )

    def zfill(self, width: int) -> DaskExpr:
        return self.compliant._with_callable(
            lambda expr, width: expr.str.zfill(width), "zfill", width=width
        )

    to_date = not_implemented()
