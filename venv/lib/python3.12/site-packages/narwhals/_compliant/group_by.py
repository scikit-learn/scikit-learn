from __future__ import annotations

import re
from itertools import chain
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Protocol, TypeVar

from narwhals._compliant.typing import (
    CompliantDataFrameAny,
    CompliantDataFrameT_co,
    CompliantExprT_contra,
    CompliantFrameT,
    CompliantFrameT_co,
    CompliantLazyFrameAny,
    DepthTrackingExprAny,
    DepthTrackingExprT_contra,
    EagerExprT_contra,
    NarwhalsAggregation,
)
from narwhals._utils import is_sequence_of

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence

    from narwhals._compliant.expr import CompliantExpr

    _SameFrameT = TypeVar("_SameFrameT", CompliantDataFrameAny, CompliantLazyFrameAny)


__all__ = ["CompliantGroupBy", "DepthTrackingGroupBy", "EagerGroupBy"]

NativeAggregationT_co = TypeVar(
    "NativeAggregationT_co", bound="str | Callable[..., Any]", covariant=True
)


_RE_LEAF_NAME: re.Pattern[str] = re.compile(r"(\w+->)")


def _evaluate_aliases(
    frame: CompliantFrameT, exprs: Iterable[CompliantExpr[CompliantFrameT, Any]], /
) -> list[str]:
    it = (expr._evaluate_aliases(frame) for expr in exprs)
    return list(chain.from_iterable(it))


class CompliantGroupBy(Protocol[CompliantFrameT_co, CompliantExprT_contra]):
    _compliant_frame: Any

    @property
    def compliant(self) -> CompliantFrameT_co:
        return self._compliant_frame  # type: ignore[no-any-return]

    def __init__(
        self,
        compliant_frame: CompliantFrameT_co,
        keys: Sequence[CompliantExprT_contra] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None: ...

    def agg(self, *exprs: CompliantExprT_contra) -> CompliantFrameT_co: ...


class DataFrameGroupBy(
    CompliantGroupBy[CompliantDataFrameT_co, CompliantExprT_contra],
    Protocol[CompliantDataFrameT_co, CompliantExprT_contra],
):
    def __iter__(self) -> Iterator[tuple[Any, CompliantDataFrameT_co]]: ...


class ParseKeysGroupBy(
    CompliantGroupBy[CompliantFrameT_co, CompliantExprT_contra],
    Protocol[CompliantFrameT_co, CompliantExprT_contra],
):
    def _parse_keys(
        self,
        compliant_frame: _SameFrameT,
        keys: Sequence[CompliantExprT_contra] | Sequence[str],
    ) -> tuple[_SameFrameT, list[str], list[str]]:
        if is_sequence_of(keys, str):
            keys_str = list(keys)
            return compliant_frame, keys_str, keys_str.copy()
        return self._parse_expr_keys(compliant_frame, keys=keys)

    @staticmethod
    def _parse_expr_keys(
        compliant_frame: _SameFrameT, keys: Sequence[CompliantExprT_contra]
    ) -> tuple[_SameFrameT, list[str], list[str]]:
        """Parses key expressions to set up `.agg` operation with correct information.

        Since keys are expressions, it's possible to alias any such key to match
        other dataframe column names.

        In order to match polars behavior and not overwrite columns when evaluating keys:

        - We evaluate what the output key names should be, in order to remap temporary column
            names to the expected ones, and to exclude those from unnamed expressions in
            `.agg(...)` context (see https://github.com/narwhals-dev/narwhals/pull/2325#issuecomment-2800004520)
        - Create temporary names for evaluated key expressions that are guaranteed to have
            no overlap with any existing column name.
        - Add these temporary columns to the compliant dataframe.
        """
        tmp_name_length = max(len(str(c)) for c in compliant_frame.columns) + 1

        def _temporary_name(key: str) -> str:
            # 5 is the length of `__tmp`
            key_str = str(key)  # pandas allows non-string column names :sob:
            return f"_{key_str}_tmp{'_' * (tmp_name_length - len(key_str) - 5)}"

        output_names = _evaluate_aliases(compliant_frame, keys)

        safe_keys = [
            # multi-output expression cannot have duplicate names, hence it's safe to suffix
            key.name.map(_temporary_name)
            if (metadata := key._metadata) and metadata.expansion_kind.is_multi_output()
            # otherwise it's single named and we can use Expr.alias
            else key.alias(_temporary_name(new_name))
            for key, new_name in zip(keys, output_names)
        ]
        return (
            compliant_frame.with_columns(*safe_keys),
            _evaluate_aliases(compliant_frame, safe_keys),
            output_names,
        )


class DepthTrackingGroupBy(
    ParseKeysGroupBy[CompliantFrameT_co, DepthTrackingExprT_contra],
    Protocol[CompliantFrameT_co, DepthTrackingExprT_contra, NativeAggregationT_co],
):
    """`CompliantGroupBy` variant, deals with `Eager` and other backends that utilize `CompliantExpr._depth`."""

    _REMAP_AGGS: ClassVar[Mapping[NarwhalsAggregation, Any]]
    """Mapping from `narwhals` to native representation.

    Note:
    - `Dask` *may* return a `Callable` instead of a `str` referring to one.
    """

    def _ensure_all_simple(self, exprs: Sequence[DepthTrackingExprT_contra]) -> None:
        for expr in exprs:
            if not self._is_simple(expr):
                name = self.compliant._implementation.name.lower()
                msg = (
                    f"Non-trivial complex aggregation found.\n\n"
                    f"Hint: you were probably trying to apply a non-elementary aggregation with a"
                    f"{name!r} table.\n"
                    "Please rewrite your query such that group-by aggregations "
                    "are elementary. For example, instead of:\n\n"
                    "    df.group_by('a').agg(nw.col('b').round(2).mean())\n\n"
                    "use:\n\n"
                    "    df.with_columns(nw.col('b').round(2)).group_by('a').agg(nw.col('b').mean())\n\n"
                )
                raise ValueError(msg)

    @classmethod
    def _is_simple(cls, expr: DepthTrackingExprAny, /) -> bool:
        """Return `True` is we can efficiently use `expr` in a native `group_by` context."""
        return expr._is_elementary() and cls._leaf_name(expr) in cls._REMAP_AGGS

    @classmethod
    def _remap_expr_name(
        cls, name: NarwhalsAggregation | Any, /
    ) -> NativeAggregationT_co:
        """Replace `name`, with some native representation.

        Arguments:
            name: Name of a `nw.Expr` aggregation method.

        Returns:
            A native compatible representation.
        """
        return cls._REMAP_AGGS.get(name, name)

    @classmethod
    def _leaf_name(cls, expr: DepthTrackingExprAny, /) -> NarwhalsAggregation | Any:
        """Return the last function name in the chain defined by `expr`."""
        return _RE_LEAF_NAME.sub("", expr._function_name)


class EagerGroupBy(
    DepthTrackingGroupBy[
        CompliantDataFrameT_co, EagerExprT_contra, NativeAggregationT_co
    ],
    DataFrameGroupBy[CompliantDataFrameT_co, EagerExprT_contra],
    Protocol[CompliantDataFrameT_co, EagerExprT_contra, NativeAggregationT_co],
): ...
