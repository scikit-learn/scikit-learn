from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from itertools import chain, product
from typing import TYPE_CHECKING, Any, Callable, Literal, cast, overload

import numpy as np

from narwhals._compliant import EagerDataFrame
from narwhals._pandas_like.series import PANDAS_TO_NUMPY_DTYPE_MISSING, PandasLikeSeries
from narwhals._pandas_like.utils import (
    align_and_extract_native,
    get_dtype_backend,
    import_array_module,
    iter_dtype_backends,
    native_to_narwhals_dtype,
    object_native_to_narwhals_dtype,
    rename,
    select_columns_by_name,
    set_index,
)
from narwhals._typing_compat import assert_never
from narwhals._utils import (
    Implementation,
    _into_arrow_table,
    _remap_full_join_keys,
    check_column_names_are_unique,
    exclude_column_names,
    generate_temporary_column_name,
    parse_columns_to_drop,
    scale_bytes,
    zip_strict,
)
from narwhals.dependencies import is_pandas_like_dataframe
from narwhals.exceptions import InvalidOperationError, ShapeError

if TYPE_CHECKING:
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import polars as pl
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny, CompliantLazyFrameAny
    from narwhals._pandas_like.expr import PandasLikeExpr
    from narwhals._pandas_like.group_by import PandasLikeGroupBy
    from narwhals._pandas_like.namespace import PandasLikeNamespace
    from narwhals._spark_like.utils import SparkSession
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import _EagerAllowedImpl, _LazyAllowedImpl
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.typing import (
        AsofJoinStrategy,
        DTypeBackend,
        IntoSchema,
        JoinStrategy,
        PivotAgg,
        SizedMultiIndexSelector,
        SizedMultiNameSelector,
        SizeUnit,
        UniqueKeepStrategy,
        _2DArray,
        _SliceIndex,
        _SliceName,
    )

    Constructor: TypeAlias = Callable[..., pd.DataFrame]


CLASSICAL_NUMPY_DTYPES: frozenset[np.dtype[Any]] = frozenset(
    [
        np.dtype("float64"),
        np.dtype("float32"),
        np.dtype("int64"),
        np.dtype("int32"),
        np.dtype("int16"),
        np.dtype("int8"),
        np.dtype("uint64"),
        np.dtype("uint32"),
        np.dtype("uint16"),
        np.dtype("uint8"),
        np.dtype("bool"),
        np.dtype("datetime64[s]"),
        np.dtype("datetime64[ms]"),
        np.dtype("datetime64[us]"),
        np.dtype("datetime64[ns]"),
        np.dtype("timedelta64[s]"),
        np.dtype("timedelta64[ms]"),
        np.dtype("timedelta64[us]"),
        np.dtype("timedelta64[ns]"),
        np.dtype("object"),
    ]
)


class PandasLikeDataFrame(
    EagerDataFrame["PandasLikeSeries", "PandasLikeExpr", "Any", "pd.Series[Any]"]
):
    def __init__(
        self,
        native_dataframe: Any,
        *,
        implementation: Implementation,
        version: Version,
        validate_column_names: bool,
        validate_backend_version: bool = False,
    ) -> None:
        self._native_frame = native_dataframe
        self._implementation = implementation
        self._version = version
        if validate_column_names:
            check_column_names_are_unique(native_dataframe.columns)
        if validate_backend_version:
            self._validate_backend_version()

    @classmethod
    def from_arrow(cls, data: IntoArrowTable, /, *, context: _LimitedContext) -> Self:
        implementation = context._implementation
        tbl = _into_arrow_table(data, context)
        if implementation.is_pandas():
            native = tbl.to_pandas()
        elif implementation.is_modin():
            # NOTE: Function moved + deprecated (0.26.0), then old path removed (0.31.0)
            # https://github.com/modin-project/modin/pull/6806
            # https://github.com/modin-project/modin/pull/7274
            if implementation._backend_version() >= (0, 26, 0):
                from modin.pandas.io import from_arrow as mpd_from_arrow
            else:  # pragma: no cover
                from modin.pandas.utils import (
                    from_arrow as mpd_from_arrow,  # pyright: ignore[reportAttributeAccessIssue]
                )
            native = mpd_from_arrow(tbl)
        elif implementation.is_cudf():  # pragma: no cover
            native = implementation.to_native_namespace().DataFrame.from_arrow(tbl)
        else:  # pragma: no cover
            msg = "congratulations, you entered unreachable code - please report a bug"
            raise AssertionError(msg)
        return cls.from_native(native, context=context)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self:
        from narwhals.schema import Schema

        implementation = context._implementation
        ns = implementation.to_native_namespace()
        Series = cast("type[pd.Series[Any]]", ns.Series)
        DataFrame = cast("type[pd.DataFrame]", ns.DataFrame)
        aligned_data: dict[str, pd.Series[Any] | Any] = {}
        left_most: PandasLikeSeries | None = None
        for name, series in data.items():
            if isinstance(series, Series):
                compliant = PandasLikeSeries.from_native(series, context=context)
                if left_most is None:
                    left_most = compliant
                    aligned_data[name] = series
                else:
                    aligned_data[name] = align_and_extract_native(left_most, compliant)[1]
            else:
                aligned_data[name] = series
        if aligned_data or not schema:
            native = DataFrame.from_dict(aligned_data)
        else:
            native = DataFrame.from_dict({col: [] for col in schema})
        if schema:
            backend: Iterable[DTypeBackend] | None = None
            if aligned_data:
                backend = iter_dtype_backends(native.dtypes, implementation)
            native = native.astype(Schema(schema).to_pandas(backend))
        return cls.from_native(native, context=context)

    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Mapping[str, Any]],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self:
        from narwhals.schema import Schema

        implementation = context._implementation
        ns = implementation.to_native_namespace()
        DataFrame = cast("type[pd.DataFrame]", ns.DataFrame)
        if data or not schema:
            native = DataFrame.from_records(data)
        else:
            native = DataFrame.from_dict({col: [] for col in schema})
        if schema:
            backend: Iterable[DTypeBackend] | None = None
            if data:
                backend = iter_dtype_backends(native.dtypes, implementation)
            native = native.astype(Schema(schema).to_pandas(backend))
        return cls.from_native(native, context=context)

    @staticmethod
    def _is_native(obj: Any) -> TypeIs[Any]:
        return is_pandas_like_dataframe(obj)  # pragma: no cover

    @classmethod
    def from_native(cls, data: Any, /, *, context: _LimitedContext) -> Self:
        return cls(
            data,
            implementation=context._implementation,
            version=context._version,
            validate_column_names=True,
        )

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | Sequence[str] | None,
    ) -> Self:
        from narwhals.schema import Schema

        implementation = context._implementation
        DataFrame: Constructor = implementation.to_native_namespace().DataFrame
        if isinstance(schema, (Mapping, Schema)):
            it: Iterable[DTypeBackend] = (
                get_dtype_backend(native_type, implementation)
                for native_type in schema.values()
            )
            native = DataFrame(data, columns=schema.keys()).astype(
                Schema(schema).to_pandas(it)
            )
        else:
            native = DataFrame(data, columns=cls._numpy_column_names(data, schema))
        return cls.from_native(native, context=context)

    def __narwhals_dataframe__(self) -> Self:
        return self

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def __narwhals_namespace__(self) -> PandasLikeNamespace:
        from narwhals._pandas_like.namespace import PandasLikeNamespace

        return PandasLikeNamespace(self._implementation, version=self._version)

    def __native_namespace__(self) -> ModuleType:
        if self._implementation in {
            Implementation.PANDAS,
            Implementation.MODIN,
            Implementation.CUDF,
        }:
            return self._implementation.to_native_namespace()

        msg = f"Expected pandas/modin/cudf, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def __len__(self) -> int:
        return len(self.native)

    def _with_version(self, version: Version) -> Self:
        return self.__class__(
            self.native,
            implementation=self._implementation,
            version=version,
            validate_column_names=False,
        )

    def _with_native(self, df: Any, *, validate_column_names: bool = True) -> Self:
        return self.__class__(
            df,
            implementation=self._implementation,
            version=self._version,
            validate_column_names=validate_column_names,
        )

    def _extract_comparand(self, other: PandasLikeSeries) -> pd.Series[Any]:
        index = self.native.index
        if other._broadcast:
            s = other.native
            return type(s)(s.iloc[0], index=index, dtype=s.dtype, name=s.name)
        if (len_other := len(other)) != (len_idx := len(index)):
            msg = f"Expected object of length {len_idx}, got: {len_other}."
            raise ShapeError(msg)
        if other.native.index is not index:
            return set_index(other.native, index, implementation=other._implementation)
        return other.native

    @property
    def _array_funcs(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            import numpy as np

            return np
        return import_array_module(self._implementation)

    def get_column(self, name: str) -> PandasLikeSeries:
        return PandasLikeSeries.from_native(self.native[name], context=self)

    def __array__(self, dtype: Any = None, *, copy: bool | None = None) -> _2DArray:
        return self.to_numpy(dtype=dtype, copy=copy)

    def _gather(self, rows: SizedMultiIndexSelector[pd.Series[Any]]) -> Self:
        items = list(rows) if isinstance(rows, tuple) else rows
        return self._with_native(self.native.iloc[items, :])

    def _gather_slice(self, rows: _SliceIndex | range) -> Self:
        return self._with_native(
            self.native.iloc[slice(rows.start, rows.stop, rows.step), :],
            validate_column_names=False,
        )

    def _select_slice_name(self, columns: _SliceName) -> Self:
        start = (
            self.native.columns.get_loc(columns.start)
            if columns.start is not None
            else None
        )
        stop = (
            self.native.columns.get_loc(columns.stop) + 1
            if columns.stop is not None
            else None
        )
        selector = slice(start, stop, columns.step)
        return self._with_native(
            self.native.iloc[:, selector], validate_column_names=False
        )

    def _select_slice_index(self, columns: _SliceIndex | range) -> Self:
        return self._with_native(
            self.native.iloc[:, columns], validate_column_names=False
        )

    def _select_multi_index(
        self, columns: SizedMultiIndexSelector[pd.Series[Any]]
    ) -> Self:
        columns = list(columns) if isinstance(columns, tuple) else columns
        return self._with_native(
            self.native.iloc[:, columns], validate_column_names=False
        )

    def _select_multi_name(self, columns: SizedMultiNameSelector[pd.Series[Any]]) -> Self:
        return self._with_native(self.native.loc[:, columns])

    # --- properties ---
    @property
    def columns(self) -> list[str]:
        return self.native.columns.tolist()

    @overload
    def rows(self, *, named: Literal[True]) -> list[dict[str, Any]]: ...

    @overload
    def rows(self, *, named: Literal[False]) -> list[tuple[Any, ...]]: ...

    @overload
    def rows(self, *, named: bool) -> list[tuple[Any, ...]] | list[dict[str, Any]]: ...

    def rows(self, *, named: bool) -> list[tuple[Any, ...]] | list[dict[str, Any]]:
        if not named:
            # cuDF does not support itertuples. But it does support to_dict!
            if self._implementation is Implementation.CUDF:
                # Extract the row values from the named rows
                return [tuple(row.values()) for row in self.rows(named=True)]

            return list(self.native.itertuples(index=False, name=None))

        return self.native.to_dict(orient="records")

    def iter_columns(self) -> Iterator[PandasLikeSeries]:
        for _name, series in self.native.items():  # noqa: PERF102
            yield PandasLikeSeries.from_native(series, context=self)

    _iter_columns = iter_columns

    def iter_rows(
        self, *, named: bool, buffer_size: int
    ) -> Iterator[tuple[Any, ...]] | Iterator[dict[str, Any]]:
        # The param ``buffer_size`` is only here for compatibility with the Polars API
        # and has no effect on the output.
        if not named:
            yield from self.native.itertuples(index=False, name=None)
        else:
            col_names = self.native.columns
            for row in self.native.itertuples(index=False):
                yield dict(zip(col_names, row))

    @property
    def schema(self) -> dict[str, DType]:
        native_dtypes = self.native.dtypes
        return {
            col: native_to_narwhals_dtype(
                native_dtypes[col], self._version, self._implementation
            )
            if native_dtypes[col] != "object"
            else object_native_to_narwhals_dtype(
                self.native[col], self._version, self._implementation
            )
            for col in self.native.columns
        }

    def collect_schema(self) -> dict[str, DType]:
        return self.schema

    # --- reshape ---
    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(
            select_columns_by_name(self.native, list(column_names), self._implementation),
            validate_column_names=False,
        )

    def select(self, *exprs: PandasLikeExpr) -> Self:
        new_series = self._evaluate_into_exprs(*exprs)
        if not new_series:
            # return empty dataframe, like Polars does
            return self._with_native(type(self.native)(), validate_column_names=False)
        new_series = new_series[0]._align_full_broadcast(*new_series)
        namespace = self.__narwhals_namespace__()
        df = namespace._concat_horizontal([s.native for s in new_series])
        # `concat` creates a new object, so fine to modify `.columns.name` inplace.
        df.columns.name = self.native.columns.name
        return self._with_native(df, validate_column_names=True)

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        if subset is None:
            return self._with_native(
                self.native.dropna(axis=0), validate_column_names=False
            )
        plx = self.__narwhals_namespace__()
        mask = ~plx.any_horizontal(plx.col(*subset).is_null(), ignore_nulls=True)
        return self.filter(mask)

    def estimated_size(self, unit: SizeUnit) -> int | float:
        sz = self.native.memory_usage(deep=True).sum()
        return scale_bytes(sz, unit=unit)

    def with_row_index(self, name: str, order_by: Sequence[str] | None) -> Self:
        plx = self.__narwhals_namespace__()
        if order_by is None:
            size = len(self)
            data = self._array_funcs.arange(size)

            row_index = plx._expr._from_series(
                plx._series.from_iterable(
                    data, context=self, index=self.native.index, name=name
                )
            )
        else:
            rank = plx.col(order_by[0]).rank(method="ordinal", descending=False)
            row_index = (rank.over(partition_by=[], order_by=order_by) - 1).alias(name)
        return self.select(row_index, plx.all())

    def row(self, index: int) -> tuple[Any, ...]:
        return tuple(x for x in self.native.iloc[index])

    def filter(self, predicate: PandasLikeExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        mask = self._evaluate_into_exprs(predicate)[0]
        mask_native = self._extract_comparand(mask)
        return self._with_native(
            self.native.loc[mask_native], validate_column_names=False
        )

    def with_columns(self, *exprs: PandasLikeExpr) -> Self:
        columns = self._evaluate_into_exprs(*exprs)
        if not columns and len(self) == 0:
            return self
        name_columns: dict[str, PandasLikeSeries] = {s.name: s for s in columns}
        to_concat = []
        # Make sure to preserve column order
        for name in self.native.columns:
            if name in name_columns:
                series = self._extract_comparand(name_columns.pop(name))
            else:
                series = self.native[name]
            to_concat.append(series)
        to_concat.extend(self._extract_comparand(s) for s in name_columns.values())
        namespace = self.__narwhals_namespace__()
        df = namespace._concat_horizontal(to_concat)
        # `concat` creates a new object, so fine to modify `.columns.name` inplace.
        df.columns.name = self.native.columns.name
        return self._with_native(df, validate_column_names=False)

    def rename(self, mapping: Mapping[str, str]) -> Self:
        return self._with_native(
            rename(self.native, columns=mapping, implementation=self._implementation)
        )

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        to_drop = parse_columns_to_drop(self, columns, strict=strict)
        return self._with_native(
            self.native.drop(columns=to_drop), validate_column_names=False
        )

    # --- transform ---
    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        df = self.native
        if isinstance(descending, bool):
            ascending: bool | list[bool] = not descending
        else:
            ascending = [not d for d in descending]
        na_position = "last" if nulls_last else "first"
        return self._with_native(
            df.sort_values(list(by), ascending=ascending, na_position=na_position),
            validate_column_names=False,
        )

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        df = self.native
        schema = self.schema
        if isinstance(reverse, bool) and all(schema[x].is_numeric() for x in by):
            if reverse:
                return self._with_native(df.nsmallest(k, by))
            return self._with_native(df.nlargest(k, by))
        return self._with_native(
            df.sort_values(list(by), ascending=reverse).head(k),
            validate_column_names=False,
        )

    # --- convert ---
    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if backend is None:
            return PandasLikeDataFrame(
                self.native,
                implementation=self._implementation,
                version=self._version,
                validate_column_names=False,
            )

        if backend is Implementation.PANDAS:
            kwds: dict[str, Any] = {
                "implementation": Implementation.PANDAS,
                "version": self._version,
                "validate_column_names": False,
            }
            if backend is not self._implementation:
                kwds.update(validate_backend_version=True)
            return PandasLikeDataFrame(self.to_pandas(), **kwds)

        if backend is Implementation.PYARROW:
            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                native_dataframe=self.to_arrow(),
                validate_backend_version=True,
                version=self._version,
                validate_column_names=False,
            )

        if backend is Implementation.POLARS:
            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                df=self.to_polars(), validate_backend_version=True, version=self._version
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    # --- actions ---
    def group_by(
        self, keys: Sequence[str] | Sequence[PandasLikeExpr], *, drop_null_keys: bool
    ) -> PandasLikeGroupBy:
        from narwhals._pandas_like.group_by import PandasLikeGroupBy

        return PandasLikeGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def _join_inner(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> pd.DataFrame:
        return self.native.merge(
            other.native,
            left_on=left_on,
            right_on=right_on,
            how="inner",
            suffixes=("", suffix),
        )

    def _join_left(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> pd.DataFrame:
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
        # NOTE: Keep `inplace=True` to avoid making a redundant copy.
        # This may need updating, depending on https://github.com/pandas-dev/pandas/pull/51466/files
        result_native.drop(columns=extra, inplace=True)  # noqa: PD002
        return result_native

    def _join_full(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str], suffix: str
    ) -> pd.DataFrame:
        # Pandas coalesces keys in full joins unless there's no collision
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

    def _join_cross(self, other: Self, *, suffix: str) -> pd.DataFrame:
        implementation = self._implementation
        backend_version = self._backend_version
        if (implementation.is_modin() or implementation.is_cudf()) or (
            implementation.is_pandas() and backend_version < (1, 4)
        ):
            key_token = generate_temporary_column_name(
                n_bytes=8, columns=(*self.columns, *other.columns)
            )
            result_native = self.native.assign(**{key_token: 0}).merge(
                other.native.assign(**{key_token: 0}),
                how="inner",
                left_on=key_token,
                right_on=key_token,
                suffixes=("", suffix),
            )
            # NOTE: Keep `inplace=True` to avoid making a redundant copy.
            # This may need updating, depending on https://github.com/pandas-dev/pandas/pull/51466/files
            result_native.drop(columns=key_token, inplace=True)  # noqa: PD002
            return result_native
        return self.native.merge(other.native, how="cross", suffixes=("", suffix))

    def _join_semi(
        self, other: Self, *, left_on: Sequence[str], right_on: Sequence[str]
    ) -> pd.DataFrame:
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
    ) -> pd.DataFrame:
        implementation = self._implementation

        if implementation.is_cudf():
            return self.native.merge(
                other.native, how="leftanti", left_on=left_on, right_on=right_on
            )

        indicator_token = generate_temporary_column_name(
            n_bytes=8, columns=(*self.columns, *other.columns)
        )

        other_native = self._join_filter_rename(
            other=other,
            columns_to_select=list(right_on),
            columns_mapping=dict(zip(right_on, left_on)),
        )
        result_native = self.native.merge(
            other_native,
            # TODO(FBruzzesi): See https://github.com/modin-project/modin/issues/7384
            how="left" if implementation.is_pandas() else "outer",
            indicator=indicator_token,
            left_on=left_on,
            right_on=left_on,
        ).loc[lambda t: t[indicator_token] == "left_only"]
        # NOTE: Keep `inplace=True` to avoid making a redundant copy.
        # This may need updating, depending on https://github.com/pandas-dev/pandas/pull/51466/files
        result_native.drop(columns=indicator_token, inplace=True)  # noqa: PD002
        return result_native

    def _join_filter_rename(
        self, other: Self, columns_to_select: list[str], columns_mapping: dict[str, str]
    ) -> pd.DataFrame:
        """Helper function to avoid creating extra columns and row duplication.

        Used in `"anti"` and `"semi`" join's.

        Notice that a native object is returned.
        """
        implementation = self._implementation
        return rename(
            select_columns_by_name(
                other.native,
                column_names=columns_to_select,
                implementation=implementation,
            ),
            columns=columns_mapping,
            implementation=implementation,
        ).drop_duplicates()

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

    # --- partial reduction ---

    def head(self, n: int) -> Self:
        return self._with_native(self.native.head(n), validate_column_names=False)

    def tail(self, n: int) -> Self:
        return self._with_native(self.native.tail(n), validate_column_names=False)

    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        maintain_order: bool | None = None,
        order_by: Sequence[str] | None,
    ) -> Self:
        # The param `maintain_order` is only here for compatibility with the Polars API
        # and has no effect on the output.
        mapped_keep = {"none": False, "any": "first"}.get(keep, keep)
        if subset and (error := self._check_columns_exist(subset)):
            raise error
        if order_by and maintain_order:
            token = generate_temporary_column_name(8, self.columns)
            res = (
                self.with_row_index(token, order_by=None)
                .sort(*order_by, nulls_last=False, descending=False)
                .native.drop_duplicates(subset or self.columns, keep=mapped_keep)
                .sort_values(token)
            )
            res.drop(columns=token, inplace=True)  # noqa: PD002
        elif order_by:
            res = self.sort(
                *order_by, nulls_last=False, descending=False
            ).native.drop_duplicates(subset or self.columns, keep=mapped_keep)
        else:
            res = self.native.drop_duplicates(subset or self.columns, keep=mapped_keep)
        return self._with_native(res, validate_column_names=False)

    # --- lazy-only ---
    def lazy(
        self,
        backend: _LazyAllowedImpl | None = None,
        *,
        session: SparkSession | None = None,
    ) -> CompliantLazyFrameAny:
        pandas_df = self.to_pandas()
        if backend is None:
            return self
        if backend is Implementation.DUCKDB:
            import duckdb  # ignore-banned-import

            from narwhals._duckdb.dataframe import DuckDBLazyFrame

            return DuckDBLazyFrame(
                df=duckdb.table("pandas_df"),
                validate_backend_version=True,
                version=self._version,
            )
        if backend is Implementation.POLARS:
            import polars as pl  # ignore-banned-import

            from narwhals._polars.dataframe import PolarsLazyFrame

            return PolarsLazyFrame(
                df=pl.from_pandas(pandas_df).lazy(),
                validate_backend_version=True,
                version=self._version,
            )
        if backend is Implementation.DASK:
            import dask.dataframe as dd  # ignore-banned-import

            from narwhals._dask.dataframe import DaskLazyFrame

            return DaskLazyFrame(
                native_dataframe=dd.from_pandas(pandas_df),
                validate_backend_version=True,
                version=self._version,
            )
        if backend is Implementation.IBIS:
            import ibis  # ignore-banned-import

            from narwhals._ibis.dataframe import IbisLazyFrame

            return IbisLazyFrame(
                ibis.memtable(pandas_df, columns=self.columns),
                validate_backend_version=True,
                version=self._version,
            )

        if backend.is_spark_like():
            from narwhals._spark_like.dataframe import SparkLikeLazyFrame

            if session is None:
                msg = "Spark like backends require `session` to be not None."
                raise ValueError(msg)

            return SparkLikeLazyFrame(
                session.createDataFrame(pandas_df),
                version=self._version,
                implementation=backend,
                validate_backend_version=True,
            )

        raise AssertionError  # pragma: no cover

    @property
    def shape(self) -> tuple[int, int]:
        return self.native.shape

    def to_dict(self, *, as_series: bool) -> dict[str, Any]:
        if as_series:
            return {
                col: PandasLikeSeries.from_native(self.native[col], context=self)
                for col in self.columns
            }
        return self.native.to_dict(orient="list")

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _2DArray:
        native_dtypes = self.native.dtypes

        if copy is None:
            # pandas default differs from Polars, but cuDF default is True
            copy = self._implementation is Implementation.CUDF

        if native_dtypes.isin(CLASSICAL_NUMPY_DTYPES).all():
            # Fast path, no conversions necessary.
            if dtype is not None:
                return self.native.to_numpy(dtype=dtype, copy=copy)
            return self.native.to_numpy(copy=copy)

        dtype_datetime = self._version.dtypes.Datetime
        to_convert = [
            key
            for key, val in self.schema.items()
            if isinstance(val, dtype_datetime) and val.time_zone is not None
        ]
        if to_convert:
            df = self.with_columns(
                self.__narwhals_namespace__()
                .col(*to_convert)
                .dt.convert_time_zone("UTC")
                .dt.replace_time_zone(None)
            ).native
        else:
            df = self.native

        if dtype is not None:
            return df.to_numpy(dtype=dtype, copy=copy)

        # pandas return `object` dtype for nullable dtypes if dtype=None,
        # so we cast each Series to numpy and let numpy find a common dtype.
        # If there aren't any dtypes where `to_numpy()` is "broken" (i.e. it
        # returns Object) then we just call `to_numpy()` on the DataFrame.
        for col_dtype in native_dtypes:
            if str(col_dtype) in PANDAS_TO_NUMPY_DTYPE_MISSING:
                arr: Any = np.hstack(
                    [
                        self.get_column(col).to_numpy(copy=copy, dtype=None)[:, None]
                        for col in self.columns
                    ]
                )
                return arr
        return df.to_numpy(copy=copy)

    def to_pandas(self) -> pd.DataFrame:
        if self._implementation is Implementation.PANDAS:
            return self.native
        if self._implementation is Implementation.CUDF:
            return self.native.to_pandas()
        if self._implementation is Implementation.MODIN:
            return self.native._to_pandas()
        msg = f"Unknown implementation: {self._implementation}"  # pragma: no cover
        raise AssertionError(msg)

    def to_polars(self) -> pl.DataFrame:
        import polars as pl  # ignore-banned-import

        return pl.from_pandas(self.to_pandas())

    def write_parquet(self, file: str | Path | BytesIO) -> None:
        self.native.to_parquet(file)

    @overload
    def write_csv(self, file: None) -> str: ...

    @overload
    def write_csv(self, file: str | Path | BytesIO) -> None: ...

    def write_csv(self, file: str | Path | BytesIO | None) -> str | None:
        return self.native.to_csv(file, index=False)

    # --- descriptive ---
    def is_unique(self) -> PandasLikeSeries:
        return PandasLikeSeries.from_native(
            ~self.native.duplicated(keep=False), context=self
        )

    def item(self, row: int | None, column: int | str | None) -> Any:
        if row is None and column is None:
            if self.shape != (1, 1):
                msg = (
                    "can only call `.item()` if the dataframe is of shape (1, 1),"
                    " or if explicit row/col values are provided;"
                    f" frame has shape {self.shape!r}"
                )
                raise ValueError(msg)
            return self.native.iloc[0, 0]

        if row is None or column is None:
            msg = "cannot call `.item()` with only one of `row` or `column`"
            raise ValueError(msg)

        _col = self.columns.index(column) if isinstance(column, str) else column
        return self.native.iloc[row, _col]

    def clone(self) -> Self:
        return self._with_native(self.native.copy(), validate_column_names=False)

    def gather_every(self, n: int, offset: int) -> Self:
        return self._with_native(self.native.iloc[offset::n], validate_column_names=False)

    def _pivot_into_index_values(
        self,
        on: Sequence[str],
        index: Sequence[str] | None,
        values: Sequence[str] | None,
        /,
    ) -> tuple[Sequence[str], Sequence[str]]:
        index = index or (
            exclude_column_names(self, {*on, *values})
            if values
            else exclude_column_names(self, on)
        )
        values = values or exclude_column_names(self, {*on, *index})
        return index, values

    @staticmethod
    def _pivot_multi_on_name(unique_values: tuple[str, ...], /) -> str:
        LB, RB, Q = "{", "}", '"'  # noqa: N806
        body = '","'.join(unique_values)
        return f"{LB}{Q}{body}{Q}{RB}"

    @staticmethod
    def _pivot_single_on_names(
        column_names: Iterable[str], n_values: int, separator: str, /
    ) -> list[str]:
        if n_values > 1:
            return [separator.join(col).strip() for col in column_names]
        return [col[-1] for col in column_names]

    def _pivot_multi_on_names(
        self,
        column_names: Iterable[tuple[str, ...]],
        n_on: int,
        n_values: int,
        separator: str,
        /,
    ) -> Iterator[str]:
        if n_values > 1:
            for col in column_names:
                names = col[-n_on:]
                prefix = col[0]
                yield separator.join((prefix, self._pivot_multi_on_name(names)))
        else:
            for col in column_names:
                yield self._pivot_multi_on_name(col[-n_on:])

    def _pivot_remap_column_names(
        self, column_names: Iterable[Any], *, n_on: int, n_values: int, separator: str
    ) -> list[str]:
        """Reformat output column names from a native pivot operation, to match `polars`.

        Note:
            `column_names` is a `pd.MultiIndex`, but not in the stubs.
        """
        if n_on == 1:
            return self._pivot_single_on_names(column_names, n_values, separator)
        return list(self._pivot_multi_on_names(column_names, n_on, n_values, separator))

    def _pivot_table(
        self,
        on: Sequence[str],
        index: Sequence[str],
        values: Sequence[str],
        aggregate_function: Literal[
            "min", "max", "first", "last", "sum", "mean", "median"
        ],
        /,
    ) -> Any:
        kwds: dict[Any, Any] = (
            {} if self._implementation is Implementation.CUDF else {"observed": True}
        )
        return self.native.pivot_table(
            values=values,
            index=index,
            columns=on,
            aggfunc=aggregate_function,
            margins=False,
            **kwds,
        )

    def _pivot(
        self,
        on: Sequence[str],
        index: Sequence[str],
        values: Sequence[str],
        aggregate_function: PivotAgg | None,
        /,
    ) -> pd.DataFrame:
        if aggregate_function is None:
            return self.native.pivot(columns=on, index=index, values=values)
        if aggregate_function == "len":
            return (
                self.native.groupby([*on, *index], as_index=False)
                .agg(dict.fromkeys(values, "size"))
                .pivot(columns=on, index=index, values=values)
            )
        return self._pivot_table(on, index, values, aggregate_function)

    def pivot(
        self,
        on: Sequence[str],
        *,
        index: Sequence[str] | None,
        values: Sequence[str] | None,
        aggregate_function: PivotAgg | None,
        sort_columns: bool,
        separator: str,
    ) -> Self:
        implementation = self._implementation
        if implementation.is_modin():
            msg = "pivot is not supported for Modin backend due to https://github.com/modin-project/modin/issues/7409."
            raise NotImplementedError(msg)

        index, values = self._pivot_into_index_values(on, index, values)
        result = self._pivot(on, index, values, aggregate_function)

        # Select the columns in the right order
        uniques = (
            (
                self.get_column(col)
                .unique()
                .sort(descending=False, nulls_last=False)
                .to_list()
                for col in on
            )
            if sort_columns
            else (self.get_column(col).unique().to_list() for col in on)
        )
        ordered_cols = list(product(values, *chain(uniques)))
        result = result.loc[:, ordered_cols]
        columns = result.columns
        remapped = self._pivot_remap_column_names(
            columns, n_on=len(on), n_values=len(values), separator=separator
        )
        result.columns = remapped  # type: ignore[assignment]
        result.columns.names = [""]
        return self._with_native(result.reset_index())

    def to_arrow(self) -> Any:
        if self._implementation is Implementation.CUDF:
            return self.native.to_arrow(preserve_index=False)

        import pyarrow as pa  # ignore-banned-import()

        return pa.Table.from_pandas(self.native)

    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self:
        return self._with_native(
            self.native.sample(
                n=n, frac=fraction, replace=with_replacement, random_state=seed
            ),
            validate_column_names=False,
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

    def explode(self, columns: Sequence[str]) -> Self:
        dtypes = self._version.dtypes

        schema = self.collect_schema()
        for col_to_explode in columns:
            dtype = schema[col_to_explode]

            if dtype != dtypes.List:
                msg = (
                    f"`explode` operation not supported for dtype `{dtype}`, "
                    "expected List type"
                )
                raise InvalidOperationError(msg)

        if len(columns) == 1:
            return self._with_native(
                self.native.explode(columns[0]), validate_column_names=False
            )
        native_frame = self.native
        anchor_series = native_frame[columns[0]].list.len()

        if not all(
            (native_frame[col_name].list.len() == anchor_series).all()
            for col_name in columns[1:]
        ):
            msg = "exploded columns must have matching element counts"
            raise ShapeError(msg)

        original_columns = self.columns
        other_columns = [c for c in original_columns if c not in columns]

        exploded_frame = native_frame[[*other_columns, columns[0]]].explode(columns[0])
        exploded_series = [
            native_frame[col_name].explode().to_frame() for col_name in columns[1:]
        ]

        plx = self.__native_namespace__()
        return self._with_native(
            plx.concat([exploded_frame, *exploded_series], axis=1)[original_columns],
            validate_column_names=False,
        )
