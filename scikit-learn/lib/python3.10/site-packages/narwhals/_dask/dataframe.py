from __future__ import annotations

from typing import TYPE_CHECKING, Any

import dask.dataframe as dd

from narwhals._dask.utils import add_row_index, evaluate_exprs
from narwhals._expression_parsing import ExprKind
from narwhals._pandas_like.utils import native_to_narwhals_dtype, select_columns_by_name
from narwhals._typing_compat import assert_never
from narwhals._utils import (
    Implementation,
    ValidateBackendVersion,
    _remap_full_join_keys,
    check_column_names_are_unique,
    check_columns_exist,
    generate_temporary_column_name,
    not_implemented,
    parse_columns_to_drop,
    zip_strict,
)
from narwhals.typing import CompliantLazyFrame

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import dask.dataframe.dask_expr as dx
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny
    from narwhals._dask.expr import DaskExpr
    from narwhals._dask.group_by import DaskLazyGroupBy
    from narwhals._dask.namespace import DaskNamespace
    from narwhals._typing import _EagerAllowedImpl
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dataframe import LazyFrame
    from narwhals.dtypes import DType
    from narwhals.exceptions import ColumnNotFoundError
    from narwhals.typing import AsofJoinStrategy, JoinStrategy, UniqueKeepStrategy

Incomplete: TypeAlias = "Any"
"""Using `_pandas_like` utils with `_dask`.

Typing this correctly will complicate the `_pandas_like`-side.
Very low priority until `dask` adds typing.
"""


class DaskLazyFrame(
    CompliantLazyFrame["DaskExpr", "dd.DataFrame", "LazyFrame[dd.DataFrame]"],
    ValidateBackendVersion,
):
    _implementation = Implementation.DASK

    def __init__(
        self,
        native_dataframe: dd.DataFrame,
        *,
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        self._native_frame: dd.DataFrame = native_dataframe
        self._version = version
        self._cached_schema: dict[str, DType] | None = None
        self._cached_columns: list[str] | None = None
        if validate_backend_version:
            self._validate_backend_version()

    @staticmethod
    def _is_native(obj: dd.DataFrame | Any) -> TypeIs[dd.DataFrame]:
        return isinstance(obj, dd.DataFrame)

    @classmethod
    def from_native(cls, data: dd.DataFrame, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version)

    def to_narwhals(self) -> LazyFrame[dd.DataFrame]:
        return self._version.lazyframe(self, level="lazy")

    def __native_namespace__(self) -> ModuleType:
        if self._implementation is Implementation.DASK:
            return self._implementation.to_native_namespace()

        msg = f"Expected dask, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def __narwhals_namespace__(self) -> DaskNamespace:
        from narwhals._dask.namespace import DaskNamespace

        return DaskNamespace(version=self._version)

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version)

    def _with_native(self, df: Any) -> Self:
        return self.__class__(df, version=self._version)

    def _check_columns_exist(self, subset: Sequence[str]) -> ColumnNotFoundError | None:
        return check_columns_exist(subset, available=self.columns)

    def _iter_columns(self) -> Iterator[dx.Series]:
        for _col, ser in self.native.items():  # noqa: PERF102
            yield ser

    def with_columns(self, *exprs: DaskExpr) -> Self:
        new_series = evaluate_exprs(self, *exprs)
        return self._with_native(self.native.assign(**dict(new_series)))

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        result = self.native.compute(**kwargs)

        if backend is None or backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                result,
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.POLARS:
            import polars as pl  # ignore-banned-import

            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                pl.from_pandas(result),
                validate_backend_version=True,
                version=self._version,
            )

        if backend is Implementation.PYARROW:
            import pyarrow as pa  # ignore-banned-import

            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                pa.Table.from_pandas(result),
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    @property
    def columns(self) -> list[str]:
        if self._cached_columns is None:
            self._cached_columns = (
                list(self.schema)
                if self._cached_schema is not None
                else self.native.columns.tolist()
            )
        return self._cached_columns

    def filter(self, predicate: DaskExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        mask = predicate(self)[0]
        return self._with_native(self.native.loc[mask])

    def simple_select(self, *column_names: str) -> Self:
        df: Incomplete = self.native
        native = select_columns_by_name(df, list(column_names), self._implementation)
        return self._with_native(native)

    def aggregate(self, *exprs: DaskExpr) -> Self:
        new_series = evaluate_exprs(self, *exprs)
        df = dd.concat([val.rename(name) for name, val in new_series], axis=1)
        return self._with_native(df)

    def select(self, *exprs: DaskExpr) -> Self:
        new_series = evaluate_exprs(self, *exprs)
        df: Incomplete = self.native
        df = select_columns_by_name(
            df.assign(**dict(new_series)),
            [s[0] for s in new_series],
            self._implementation,
        )
        return self._with_native(df)

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        if subset is None:
            return self._with_native(self.native.dropna())
        plx = self.__narwhals_namespace__()
        mask = ~plx.any_horizontal(plx.col(*subset).is_null(), ignore_nulls=True)
        return self.filter(mask)

    @property
    def schema(self) -> dict[str, DType]:
        if self._cached_schema is None:
            native_dtypes = self.native.dtypes
            self._cached_schema = {
                col: native_to_narwhals_dtype(
                    native_dtypes[col], self._version, self._implementation
                )
                for col in self.native.columns
            }
        return self._cached_schema

    def collect_schema(self) -> dict[str, DType]:
        return self.schema

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        to_drop = parse_columns_to_drop(self, columns, strict=strict)

        return self._with_native(self.native.drop(columns=to_drop))

    def with_row_index(self, name: str, order_by: Sequence[str] | None) -> Self:
        # Implementation is based on the following StackOverflow reply:
        # https://stackoverflow.com/questions/60831518/in-dask-how-does-one-add-a-range-of-integersauto-increment-to-a-new-column/60852409#60852409
        if order_by is None:
            return self._with_native(add_row_index(self.native, name))
        plx = self.__narwhals_namespace__()
        columns = self.columns
        const_expr = plx.lit(value=1, dtype=None).alias(name).broadcast(ExprKind.LITERAL)
        row_index_expr = (
            plx.col(name).cum_sum(reverse=False).over(partition_by=[], order_by=order_by)
            - 1
        )
        return self.with_columns(const_expr).select(row_index_expr, plx.col(*columns))

    def rename(self, mapping: Mapping[str, str]) -> Self:
        return self._with_native(self.native.rename(columns=mapping))

    def head(self, n: int) -> Self:
        return self._with_native(self.native.head(n=n, compute=False, npartitions=-1))

    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        order_by: Sequence[str] | None,
    ) -> Self:
        if subset and (error := self._check_columns_exist(subset)):
            raise error
        if keep == "none":
            subset = subset or self.columns
            token = generate_temporary_column_name(
                n_bytes=8, columns=subset, prefix="count_"
            )
            ser = self.native.groupby(subset).size().rename(token)
            ser = ser[ser == 1]
            unique = ser.reset_index().drop(columns=token)
            result = self.native.merge(unique, on=subset, how="inner")
        else:
            mapped_keep = {"any": "first"}.get(keep, keep)
            if order_by:
                native = self.sort(*order_by, descending=False, nulls_last=False).native
            else:
                native = self.native
            result = native.drop_duplicates(subset=subset, keep=mapped_keep)
        return self._with_native(result)

    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        if isinstance(descending, bool):
            ascending: bool | list[bool] = not descending
        else:
            ascending = [not d for d in descending]
        position = "last" if nulls_last else "first"
        return self._with_native(
            self.native.sort_values(list(by), ascending=ascending, na_position=position)
        )

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        df = self.native
        schema = self.schema
        by = list(by)
        if isinstance(reverse, bool) and all(schema[x].is_numeric() for x in by):
            if reverse:
                return self._with_native(df.nsmallest(k, by))
            return self._with_native(df.nlargest(k, by))
        if isinstance(reverse, bool):
            reverse = [reverse] * len(by)
        return self._with_native(
            df.sort_values(by, ascending=list(reverse)).head(
                n=k, compute=False, npartitions=-1
            )
        )

    def _join_inner(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> dd.DataFrame:
        return self.native.merge(
            other.native,
            left_on=left_on,
            right_on=right_on,
            how="inner",
            suffixes=("", suffix),
        )

    def _join_left(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> dd.DataFrame:
        result_native = self.native.merge(
            other.native,
            how="left",
            left_on=left_on,
            right_on=right_on,
            suffixes=("", suffix),
        )
        extra = [
            right_key if right_key not in self.columns else f"{right_key}{suffix}"
            for left_key, right_key in zip_strict(left_on, right_on)
            if right_key != left_key
        ]
        return result_native.drop(columns=extra)

    def _join_full(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> dd.DataFrame:
        # dask does not retain keys post-join
        # we must append the suffix to each key before-hand

        right_on_mapper = _remap_full_join_keys(left_on, right_on, suffix)
        other_native = other.native.rename(columns=right_on_mapper)
        check_column_names_are_unique(other_native.columns)
        right_suffixed = list(right_on_mapper.values())
        return self.native.merge(
            other_native,
            left_on=left_on,
            right_on=right_suffixed,
            how="outer",
            suffixes=("", suffix),
        )

    def _join_cross(self, other: Self, *, suffix: str) -> dd.DataFrame:
        key_token = generate_temporary_column_name(
            n_bytes=8, columns=(*self.columns, *other.columns), prefix="cross_join_key_"
        )
        return (
            self.native.assign(**{key_token: 0})
            .merge(
                other.native.assign(**{key_token: 0}),
                how="inner",
                left_on=key_token,
                right_on=key_token,
                suffixes=("", suffix),
            )
            .drop(columns=key_token)
        )

    def _join_semi(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str]
    ) -> dd.DataFrame:
        other_native = self._join_filter_rename(
            other=other,
            columns_to_select=list(right_on),
            columns_mapping=dict(zip(right_on, left_on)),
        )
        return self.native.merge(
            other_native, how="inner", left_on=left_on, right_on=left_on
        )

    def _join_anti(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str]
    ) -> dd.DataFrame:
        indicator_token = generate_temporary_column_name(
            n_bytes=8, columns=(*self.columns, *other.columns), prefix="join_indicator_"
        )
        other_native = self._join_filter_rename(
            other=other,
            columns_to_select=list(right_on),
            columns_mapping=dict(zip(right_on, left_on)),
        )
        df = self.native.merge(
            other_native,
            how="left",
            indicator=indicator_token,  # pyright: ignore[reportArgumentType]
            left_on=left_on,
            right_on=left_on,
        )
        return df[df[indicator_token] == "left_only"].drop(columns=[indicator_token])

    def _join_filter_rename(
        self, other: Self, columns_to_select: list[str], columns_mapping: dict[str, str]
    ) -> dd.DataFrame:
        """Helper function to avoid creating extra columns and row duplication.

        Used in `"anti"` and `"semi`" join's.

        Notice that a native object is returned.
        """
        other_native: Incomplete = other.native
        # rename to avoid creating extra columns in join
        return (
            select_columns_by_name(other_native, columns_to_select, self._implementation)
            .rename(columns=columns_mapping)
            .drop_duplicates()
        )

    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        if how == "cross":
            result = self._join_cross(other=other, suffix=suffix)

        elif left_on is None or right_on is None:  # pragma: no cover
            raise ValueError(left_on, right_on)

        elif how == "inner":
            result = self._join_inner(
                other=other, left_on=left_on, right_on=right_on, suffix=suffix
            )
        elif how == "anti":
            result = self._join_anti(other=other, left_on=left_on, right_on=right_on)
        elif how == "semi":
            result = self._join_semi(other=other, left_on=left_on, right_on=right_on)
        elif how == "left":
            result = self._join_left(
                other=other, left_on=left_on, right_on=right_on, suffix=suffix
            )
        elif how == "full":
            result = self._join_full(
                other=other, left_on=left_on, right_on=right_on, suffix=suffix
            )
        else:
            assert_never(how)
        return self._with_native(result)

    def join_asof(
        self,
        other: Self,
        *,
        left_on: str,
        right_on: str,
        by_left: Sequence[str] | None,
        by_right: Sequence[str] | None,
        strategy: AsofJoinStrategy,
        suffix: str,
    ) -> Self:
        plx = self.__native_namespace__()
        return self._with_native(
            plx.merge_asof(
                self.native,
                other.native,
                left_on=left_on,
                right_on=right_on,
                left_by=by_left,
                right_by=by_right,
                direction=strategy,
                suffixes=("", suffix),
            )
        )

    def group_by(
        self, keys: Sequence[str] | Sequence[DaskExpr], *, drop_null_keys: bool
    ) -> DaskLazyGroupBy:
        from narwhals._dask.group_by import DaskLazyGroupBy

        return DaskLazyGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def tail(self, n: int) -> Self:  # pragma: no cover
        native_frame = self.native
        n_partitions = native_frame.npartitions

        if n_partitions == 1:
            return self._with_native(self.native.tail(n=n, compute=False))
        msg = (
            "`LazyFrame.tail` is not supported for Dask backend with multiple partitions."
        )
        raise NotImplementedError(msg)

    def gather_every(self, n: int, offset: int) -> Self:
        row_index_token = generate_temporary_column_name(
            n_bytes=8, columns=self.columns, prefix="row_index_"
        )
        plx = self.__narwhals_namespace__()
        return (
            self.with_row_index(row_index_token, order_by=None)
            .filter(
                (plx.col(row_index_token) >= offset)
                & ((plx.col(row_index_token) - offset) % n == 0)
            )
            .drop([row_index_token], strict=False)
        )

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        return self._with_native(
            self.native.melt(
                id_vars=index,
                value_vars=on,
                var_name=variable_name,
                value_name=value_name,
            )
        )

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        self.native.to_parquet(file)

    explode = not_implemented()
