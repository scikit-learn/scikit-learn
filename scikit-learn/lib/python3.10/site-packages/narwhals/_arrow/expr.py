from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.series import ArrowSeries
from narwhals._compliant import EagerExpr
from narwhals._expression_parsing import evaluate_output_names_and_aliases
from narwhals._utils import (
    Implementation,
    generate_temporary_column_name,
    not_implemented,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from narwhals._arrow.dataframe import ArrowDataFrame
    from narwhals._arrow.namespace import ArrowNamespace
    from narwhals._compliant.typing import AliasNames, EvalNames, EvalSeries, ScalarKwargs
    from narwhals._expression_parsing import ExprMetadata
    from narwhals._utils import Version, _LimitedContext


class ArrowExpr(EagerExpr["ArrowDataFrame", ArrowSeries]):
    _implementation: Implementation = Implementation.PYARROW

    def __init__(
        self,
        call: EvalSeries[ArrowDataFrame, ArrowSeries],
        *,
        depth: int,
        function_name: str,
        evaluate_output_names: EvalNames[ArrowDataFrame],
        alias_output_names: AliasNames | None,
        version: Version,
        scalar_kwargs: ScalarKwargs | None = None,
        implementation: Implementation | None = None,
    ) -> None:
        self._call = call
        self._depth = depth
        self._function_name = function_name
        self._depth = depth
        self._evaluate_output_names = evaluate_output_names
        self._alias_output_names = alias_output_names
        self._version = version
        self._scalar_kwargs = scalar_kwargs or {}
        self._metadata: ExprMetadata | None = None

    @classmethod
    def from_column_names(
        cls: type[Self],
        evaluate_column_names: EvalNames[ArrowDataFrame],
        /,
        *,
        context: _LimitedContext,
        function_name: str = "",
    ) -> Self:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            try:
                return [
                    ArrowSeries(
                        df.native[column_name], name=column_name, version=df._version
                    )
                    for column_name in evaluate_column_names(df)
                ]
            except KeyError as e:
                if error := df._check_columns_exist(evaluate_column_names(df)):
                    raise error from e
                raise

        return cls(
            func,
            depth=0,
            function_name=function_name,
            evaluate_output_names=evaluate_column_names,
            alias_output_names=None,
            version=context._version,
        )

    @classmethod
    def from_column_indices(cls, *column_indices: int, context: _LimitedContext) -> Self:
        def func(df: ArrowDataFrame) -> list[ArrowSeries]:
            tbl = df.native
            cols = df.columns
            return [
                ArrowSeries.from_native(tbl[i], name=cols[i], context=df)
                for i in column_indices
            ]

        return cls(
            func,
            depth=0,
            function_name="nth",
            evaluate_output_names=cls._eval_names_indices(column_indices),
            alias_output_names=None,
            version=context._version,
        )

    def __narwhals_namespace__(self) -> ArrowNamespace:
        from narwhals._arrow.namespace import ArrowNamespace

        return ArrowNamespace(version=self._version)

    def _reuse_series_extra_kwargs(
        self, *, returns_scalar: bool = False
    ) -> dict[str, Any]:
        return {"_return_py_scalar": False} if returns_scalar else {}

    def over(self, partition_by: Sequence[str], order_by: Sequence[str]) -> Self:
        meta = self._metadata
        if partition_by and meta is not None and not meta.is_scalar_like:
            msg = "Only aggregation or literal operations are supported in grouped `over` context for PyArrow."
            raise NotImplementedError(msg)

        if not partition_by:
            # e.g. `nw.col('a').cum_sum().order_by(key)`
            # which we can always easily support, as it doesn't require grouping.
            assert order_by  # noqa: S101

            def func(df: ArrowDataFrame) -> Sequence[ArrowSeries]:
                token = generate_temporary_column_name(8, df.columns)
                df = df.with_row_index(token, order_by=None).sort(
                    *order_by, descending=False, nulls_last=False
                )
                results = self(df.drop([token], strict=True))
                if meta is not None and meta.is_scalar_like:
                    # We need to broadcast the results to the original size, since
                    # `over` is a length-preserving operation.
                    size = len(df)
                    return [s._with_native(pa.repeat(s.item(), size)) for s in results]

                # TODO(marco): is there a way to do this efficiently without
                # doing 2 sorts? Here we're sorting the dataframe and then
                # again calling `sort_indices`. `ArrowSeries.scatter` would also sort.
                sorting_indices = pc.sort_indices(df.get_column(token).native)
                return [s._with_native(s.native.take(sorting_indices)) for s in results]
        else:

            def func(df: ArrowDataFrame) -> Sequence[ArrowSeries]:
                if order_by:
                    df = df.sort(*order_by, descending=False, nulls_last=False)

                output_names, aliases = evaluate_output_names_and_aliases(self, df, [])
                if overlap := set(output_names).intersection(partition_by):
                    # E.g. `df.select(nw.all().sum().over('a'))`. This is well-defined,
                    # we just don't support it yet.
                    msg = (
                        f"Column names {overlap} appear in both expression output names and in `over` keys.\n"
                        "This is not yet supported."
                    )
                    raise NotImplementedError(msg)

                tmp = df.group_by(partition_by, drop_null_keys=False).agg(self)
                tmp = df.simple_select(*partition_by).join(
                    tmp,
                    how="left",
                    left_on=partition_by,
                    right_on=partition_by,
                    suffix="_right",
                )
                return [tmp.get_column(alias) for alias in aliases]

        return self.__class__(
            func,
            depth=self._depth + 1,
            function_name=self._function_name + "->over",
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            version=self._version,
        )

    ewm_mean = not_implemented()
