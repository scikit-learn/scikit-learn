from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence, Sized
from itertools import chain
from typing import TYPE_CHECKING, Any, Literal, Protocol, TypeVar, overload

from narwhals._compliant.typing import (
    CompliantDataFrameAny,
    CompliantExprT_contra,
    CompliantLazyFrameAny,
    CompliantSeriesT,
    EagerExprT,
    EagerSeriesT,
    NativeDataFrameT,
    NativeLazyFrameT,
    NativeSeriesT,
)
from narwhals._translate import (
    ArrowConvertible,
    DictConvertible,
    FromNative,
    NumpyConvertible,
    ToNarwhals,
    ToNarwhalsT_co,
)
from narwhals._typing_compat import assert_never
from narwhals._utils import (
    ValidateBackendVersion,
    Version,
    _StoresNative,
    check_columns_exist,
    is_compliant_series,
    is_index_selector,
    is_range,
    is_sequence_like,
    is_sized_multi_index_selector,
    is_slice_index,
    is_slice_none,
)

if TYPE_CHECKING:
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import Self, TypeAlias

    from narwhals._compliant.group_by import CompliantGroupBy, DataFrameGroupBy
    from narwhals._compliant.namespace import EagerNamespace
    from narwhals._spark_like.utils import SparkSession
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import _EagerAllowedImpl, _LazyAllowedImpl
    from narwhals._utils import Implementation, _LimitedContext
    from narwhals.dataframe import DataFrame
    from narwhals.dtypes import DType
    from narwhals.exceptions import ColumnNotFoundError
    from narwhals.typing import (
        AsofJoinStrategy,
        IntoSchema,
        JoinStrategy,
        MultiColSelector,
        MultiIndexSelector,
        PivotAgg,
        SingleIndexSelector,
        SizedMultiIndexSelector,
        SizedMultiNameSelector,
        SizeUnit,
        UniqueKeepStrategy,
        _2DArray,
        _SliceIndex,
        _SliceName,
    )

    Incomplete: TypeAlias = Any

__all__ = ["CompliantDataFrame", "CompliantFrame", "CompliantLazyFrame", "EagerDataFrame"]

T = TypeVar("T")

_ToDict: TypeAlias = "dict[str, CompliantSeriesT] | dict[str, list[Any]]"  # noqa: PYI047

_NativeFrameT = TypeVar("_NativeFrameT")


class CompliantFrame(
    _StoresNative[_NativeFrameT],
    FromNative[_NativeFrameT],
    ToNarwhals[ToNarwhalsT_co],
    Protocol[CompliantExprT_contra, _NativeFrameT, ToNarwhalsT_co],
):
    """Common parts of `DataFrame`, `LazyFrame`."""

    _native_frame: _NativeFrameT
    _implementation: Implementation
    _version: Version

    def __native_namespace__(self) -> ModuleType: ...
    def __narwhals_namespace__(self) -> Any: ...
    def _with_native(self, df: _NativeFrameT) -> Self: ...
    def _with_version(self, version: Version) -> Self: ...
    @classmethod
    def from_native(cls, data: _NativeFrameT, /, *, context: _LimitedContext) -> Self: ...
    @property
    def columns(self) -> Sequence[str]: ...
    @property
    def native(self) -> _NativeFrameT:
        return self._native_frame

    @property
    def schema(self) -> Mapping[str, DType]: ...

    def collect_schema(self) -> Mapping[str, DType]: ...
    def drop(self, columns: Sequence[str], *, strict: bool) -> Self: ...
    def drop_nulls(self, subset: Sequence[str] | None) -> Self: ...
    def explode(self, columns: Sequence[str]) -> Self: ...
    def filter(self, predicate: CompliantExprT_contra | Incomplete) -> Self: ...
    def group_by(
        self,
        keys: Sequence[str] | Sequence[CompliantExprT_contra],
        *,
        drop_null_keys: bool,
    ) -> CompliantGroupBy[Self, CompliantExprT_contra]: ...
    def head(self, n: int) -> Self: ...
    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self: ...
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
    ) -> Self: ...
    def rename(self, mapping: Mapping[str, str]) -> Self: ...
    def select(self, *exprs: CompliantExprT_contra) -> Self: ...
    def simple_select(self, *column_names: str) -> Self:
        """`select` where all args are column names."""
        ...

    def sort(
        self, *by: str, descending: bool | Sequence[bool], nulls_last: bool
    ) -> Self: ...
    def tail(self, n: int) -> Self: ...
    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        order_by: Sequence[str] | None,
    ) -> Self: ...
    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self: ...
    def with_columns(self, *exprs: CompliantExprT_contra) -> Self: ...
    def with_row_index(self, name: str, order_by: Sequence[str]) -> Self: ...


class CompliantDataFrame(
    NumpyConvertible["_2DArray", "_2DArray"],
    DictConvertible["_ToDict[CompliantSeriesT]", Mapping[str, Any]],
    ArrowConvertible["pa.Table", "IntoArrowTable"],
    Sized,
    CompliantFrame[CompliantExprT_contra, NativeDataFrameT, ToNarwhalsT_co],
    Protocol[CompliantSeriesT, CompliantExprT_contra, NativeDataFrameT, ToNarwhalsT_co],
):
    def __narwhals_dataframe__(self) -> Self: ...
    @classmethod
    def from_arrow(cls, data: IntoArrowTable, /, *, context: _LimitedContext) -> Self: ...
    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self: ...
    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Mapping[str, Any]],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self: ...
    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | Sequence[str] | None,
    ) -> Self: ...
    def __array__(self, dtype: Any, *, copy: bool | None) -> _2DArray: ...
    def __getitem__(
        self,
        item: tuple[
            SingleIndexSelector | MultiIndexSelector[CompliantSeriesT],
            MultiColSelector[CompliantSeriesT],
        ],
    ) -> Self: ...

    @property
    def shape(self) -> tuple[int, int]: ...
    def clone(self) -> Self: ...
    def estimated_size(self, unit: SizeUnit) -> int | float: ...
    def gather_every(self, n: int, offset: int) -> Self: ...
    def get_column(self, name: str) -> CompliantSeriesT: ...
    def group_by(
        self,
        keys: Sequence[str] | Sequence[CompliantExprT_contra],
        *,
        drop_null_keys: bool,
    ) -> DataFrameGroupBy[Self, Any]: ...
    def item(self, row: int | None, column: int | str | None) -> Any: ...
    def iter_columns(self) -> Iterator[CompliantSeriesT]: ...
    def iter_rows(
        self, *, named: bool, buffer_size: int
    ) -> Iterator[tuple[Any, ...]] | Iterator[Mapping[str, Any]]: ...
    def is_unique(self) -> CompliantSeriesT: ...
    def lazy(
        self, backend: _LazyAllowedImpl | None, *, session: SparkSession | None
    ) -> CompliantLazyFrameAny: ...
    def pivot(
        self,
        on: Sequence[str],
        *,
        index: Sequence[str] | None,
        values: Sequence[str] | None,
        aggregate_function: PivotAgg | None,
        sort_columns: bool,
        separator: str,
    ) -> Self: ...
    def row(self, index: int) -> tuple[Any, ...]: ...
    def rows(
        self, *, named: bool
    ) -> Sequence[tuple[Any, ...]] | Sequence[Mapping[str, Any]]: ...
    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self: ...
    def to_arrow(self) -> pa.Table: ...
    def to_pandas(self) -> pd.DataFrame: ...
    def to_polars(self) -> pl.DataFrame: ...
    @overload
    def to_dict(self, *, as_series: Literal[True]) -> dict[str, CompliantSeriesT]: ...
    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...
    def to_dict(
        self, *, as_series: bool
    ) -> dict[str, CompliantSeriesT] | dict[str, list[Any]]: ...
    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        maintain_order: bool | None = None,
        order_by: Sequence[str] | None,
    ) -> Self: ...
    def with_row_index(self, name: str, order_by: Sequence[str] | None) -> Self: ...
    @overload
    def write_csv(self, file: None) -> str: ...
    @overload
    def write_csv(self, file: str | Path | BytesIO) -> None: ...
    def write_csv(self, file: str | Path | BytesIO | None) -> str | None: ...
    def write_parquet(self, file: str | Path | BytesIO) -> None: ...


class CompliantLazyFrame(
    CompliantFrame[CompliantExprT_contra, NativeLazyFrameT, ToNarwhalsT_co],
    Protocol[CompliantExprT_contra, NativeLazyFrameT, ToNarwhalsT_co],
):
    def __narwhals_lazyframe__(self) -> Self: ...
    # `LazySelectorNamespace._iter_columns` depends
    def _iter_columns(self) -> Iterator[Any]: ...
    def aggregate(self, *exprs: CompliantExprT_contra) -> Self:
        """`select` where all args are aggregations or literals.

        (so, no broadcasting is necessary).
        """
        ...

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny: ...
    def sink_parquet(self, file: str | Path | BytesIO) -> None: ...


class EagerDataFrame(
    CompliantDataFrame[
        EagerSeriesT, EagerExprT, NativeDataFrameT, "DataFrame[NativeDataFrameT]"
    ],
    CompliantLazyFrame[EagerExprT, "Incomplete", "DataFrame[NativeDataFrameT]"],
    ValidateBackendVersion,
    Protocol[EagerSeriesT, EagerExprT, NativeDataFrameT, NativeSeriesT],
):
    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    def __narwhals_namespace__(
        self,
    ) -> EagerNamespace[
        Self, EagerSeriesT, EagerExprT, NativeDataFrameT, NativeSeriesT
    ]: ...

    def to_narwhals(self) -> DataFrame[NativeDataFrameT]:
        return self._version.dataframe(self, level="full")

    def aggregate(self, *exprs: EagerExprT) -> Self:  # pyright: ignore[reportIncompatibleMethodOverride]
        # NOTE: Ignore intermittent [False Negative] (1)
        # Method "aggregate" overrides class "CompliantLazyFrame" in an incompatible manner
        # Keyword parameter "exprs" type mismatch: base parameter is type "EagerExprT@EagerDataFrame", override parameter is type "EagerExprT@EagerDataFrame"
        #  Type "EagerExprT@EagerDataFrame" is not assignable to type "EagerExprT@EagerDataFrame"
        # NOTE: Ignore intermittent [False Negative] (2)
        # Argument of type "EagerExprT@EagerDataFrame" cannot be assigned to parameter "exprs" of type "EagerExprT@EagerDataFrame" in function "select"
        #  Type "EagerExprT@EagerDataFrame" is not assignable to type "EagerExprT@EagerDataFrame"
        return self.select(*exprs)  # pyright: ignore[reportArgumentType]

    def _with_native(
        self, df: NativeDataFrameT, *, validate_column_names: bool = True
    ) -> Self: ...

    def _check_columns_exist(self, subset: Sequence[str]) -> ColumnNotFoundError | None:
        return check_columns_exist(subset, available=self.columns)

    def _evaluate_expr(self, expr: EagerExprT, /) -> EagerSeriesT:
        """Evaluate `expr` and ensure it has a **single** output."""
        result: Sequence[EagerSeriesT] = expr(self)
        assert len(result) == 1  # debug assertion  # noqa: S101
        return result[0]

    def _evaluate_into_exprs(self, *exprs: EagerExprT) -> Sequence[EagerSeriesT]:
        # NOTE: Ignore intermittent [False Negative]
        # Argument of type "EagerExprT@EagerDataFrame" cannot be assigned to parameter "expr" of type "EagerExprT@EagerDataFrame" in function "_evaluate_into_expr"
        #  Type "EagerExprT@EagerDataFrame" is not assignable to type "EagerExprT@EagerDataFrame"
        return tuple(
            chain.from_iterable(self._evaluate_into_expr(expr) for expr in exprs)  # pyright: ignore[reportArgumentType]
        )

    def _evaluate_into_expr(self, expr: EagerExprT, /) -> Sequence[EagerSeriesT]:
        """Return list of raw columns.

        For eager backends we alias operations at each step.

        As a safety precaution, here we can check that the expected result names match those
        we were expecting from the various `evaluate_output_names` / `alias_output_names` calls.

        Note that for PySpark / DuckDB, we are less free to liberally set aliases whenever we want.
        """
        aliases = expr._evaluate_aliases(self)
        result = expr(self)
        if list(aliases) != (
            result_aliases := [s.name for s in result]
        ):  # pragma: no cover
            msg = f"Safety assertion failed, expected {aliases}, got {result_aliases}"
            raise AssertionError(msg)
        return result

    def _extract_comparand(self, other: EagerSeriesT, /) -> Any:
        """Extract native Series, broadcasting to `len(self)` if necessary."""
        ...

    @staticmethod
    def _numpy_column_names(
        data: _2DArray, columns: Sequence[str] | None, /
    ) -> list[str]:
        return list(columns or (f"column_{x}" for x in range(data.shape[1])))

    def _gather(self, rows: SizedMultiIndexSelector[NativeSeriesT]) -> Self: ...
    def _gather_slice(self, rows: _SliceIndex | range) -> Self: ...
    def _select_multi_index(
        self, columns: SizedMultiIndexSelector[NativeSeriesT]
    ) -> Self: ...
    def _select_multi_name(
        self, columns: SizedMultiNameSelector[NativeSeriesT]
    ) -> Self: ...
    def _select_slice_index(self, columns: _SliceIndex | range) -> Self: ...
    def _select_slice_name(self, columns: _SliceName) -> Self: ...
    def __getitem__(  # noqa: C901, PLR0912
        self,
        item: tuple[
            SingleIndexSelector | MultiIndexSelector[EagerSeriesT],
            MultiColSelector[EagerSeriesT],
        ],
    ) -> Self:
        rows, columns = item
        compliant = self
        if not is_slice_none(columns):
            if isinstance(columns, Sized) and len(columns) == 0:
                return compliant.select()
            if is_index_selector(columns):
                if is_slice_index(columns) or is_range(columns):
                    compliant = compliant._select_slice_index(columns)
                elif is_compliant_series(columns):
                    compliant = self._select_multi_index(columns.native)
                else:
                    compliant = compliant._select_multi_index(columns)
            elif isinstance(columns, slice):
                compliant = compliant._select_slice_name(columns)
            elif is_compliant_series(columns):
                compliant = self._select_multi_name(columns.native)
            elif is_sequence_like(columns):
                compliant = self._select_multi_name(columns)
            else:
                assert_never(columns)

        if not is_slice_none(rows):
            if isinstance(rows, int):
                compliant = compliant._gather([rows])
            elif isinstance(rows, (slice, range)):
                compliant = compliant._gather_slice(rows)
            elif is_compliant_series(rows):
                compliant = compliant._gather(rows.native)
            elif is_sized_multi_index_selector(rows):
                compliant = compliant._gather(rows)
            else:
                assert_never(rows)

        return compliant

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        return self.write_parquet(file)
