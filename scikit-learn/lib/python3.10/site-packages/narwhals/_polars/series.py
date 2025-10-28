from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar, cast, overload

import polars as pl

from narwhals._polars.utils import (
    BACKEND_VERSION,
    SERIES_ACCEPTS_PD_INDEX,
    SERIES_RESPECTS_DTYPE,
    PolarsAnyNamespace,
    PolarsCatNamespace,
    PolarsDateTimeNamespace,
    PolarsListNamespace,
    PolarsStringNamespace,
    PolarsStructNamespace,
    catch_polars_exception,
    extract_args_kwargs,
    extract_native,
    narwhals_to_native_dtype,
    native_to_narwhals_dtype,
)
from narwhals._utils import Implementation, requires
from narwhals.dependencies import is_numpy_array_1d, is_pandas_index

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from types import ModuleType
    from typing import Literal, TypeVar

    import pandas as pd
    import pyarrow as pa
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._polars.dataframe import Method, PolarsDataFrame
    from narwhals._polars.namespace import PolarsNamespace
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.series import Series
    from narwhals.typing import (
        Into1DArray,
        IntoDType,
        ModeKeepStrategy,
        MultiIndexSelector,
        NonNestedLiteral,
        NumericLiteral,
        PythonLiteral,
        _1DArray,
    )

    T = TypeVar("T")
    IncludeBreakpoint: TypeAlias = Literal[False, True]

Incomplete: TypeAlias = Any

# Series methods where PolarsSeries just defers to Polars.Series directly.
INHERITED_METHODS = frozenset(
    [
        "__add__",
        "__and__",
        "__floordiv__",
        "__invert__",
        "__iter__",
        "__mod__",
        "__mul__",
        "__or__",
        "__pow__",
        "__radd__",
        "__rand__",
        "__rfloordiv__",
        "__rmod__",
        "__rmul__",
        "__ror__",
        "__rsub__",
        "__rtruediv__",
        "__sub__",
        "__truediv__",
        "abs",
        "all",
        "any",
        "arg_max",
        "arg_min",
        "arg_true",
        "ceil",
        "clip",
        "count",
        "cum_max",
        "cum_min",
        "cum_prod",
        "cum_sum",
        "diff",
        "drop_nulls",
        "exp",
        "fill_null",
        "fill_nan",
        "filter",
        "floor",
        "gather_every",
        "head",
        "is_between",
        "is_close",
        "is_duplicated",
        "is_empty",
        "is_finite",
        "is_first_distinct",
        "is_in",
        "is_last_distinct",
        "is_null",
        "is_sorted",
        "is_unique",
        "item",
        "kurtosis",
        "len",
        "log",
        "max",
        "mean",
        "min",
        "mode",
        "n_unique",
        "null_count",
        "quantile",
        "rank",
        "round",
        "sample",
        "shift",
        "skew",
        "sqrt",
        "std",
        "sum",
        "tail",
        "to_arrow",
        "to_frame",
        "to_list",
        "to_pandas",
        "unique",
        "var",
        "zip_with",
    ]
)


class PolarsSeries:
    _implementation: Implementation = Implementation.POLARS
    _native_series: pl.Series
    _version: Version

    _HIST_EMPTY_SCHEMA: ClassVar[Mapping[IncludeBreakpoint, Sequence[str]]] = {
        True: ["breakpoint", "count"],
        False: ["count"],
    }

    def __init__(self, series: pl.Series, *, version: Version) -> None:
        self._native_series = series
        self._version = version

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    def __repr__(self) -> str:  # pragma: no cover
        return "PolarsSeries"

    def __narwhals_namespace__(self) -> PolarsNamespace:
        from narwhals._polars.namespace import PolarsNamespace

        return PolarsNamespace(version=self._version)

    def __narwhals_series__(self) -> Self:
        return self

    def __native_namespace__(self) -> ModuleType:
        if self._implementation is Implementation.POLARS:
            return self._implementation.to_native_namespace()

        msg = f"Expected polars, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version)

    @classmethod
    def from_iterable(
        cls,
        data: Iterable[Any],
        *,
        context: _LimitedContext,
        name: str = "",
        dtype: IntoDType | None = None,
    ) -> Self:
        version = context._version
        dtype_pl = narwhals_to_native_dtype(dtype, version) if dtype else None
        values: Incomplete = data
        if SERIES_RESPECTS_DTYPE:
            native = pl.Series(name, values, dtype=dtype_pl)
        else:  # pragma: no cover
            if (not SERIES_ACCEPTS_PD_INDEX) and is_pandas_index(values):
                values = values.to_series()
            native = pl.Series(name, values)
            if dtype_pl:
                native = native.cast(dtype_pl)
        return cls.from_native(native, context=context)

    @staticmethod
    def _is_native(obj: pl.Series | Any) -> TypeIs[pl.Series]:
        return isinstance(obj, pl.Series)

    @classmethod
    def from_native(cls, data: pl.Series, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version)

    @classmethod
    def from_numpy(cls, data: Into1DArray, /, *, context: _LimitedContext) -> Self:
        native = pl.Series(data if is_numpy_array_1d(data) else [data])
        return cls.from_native(native, context=context)

    def to_narwhals(self) -> Series[pl.Series]:
        return self._version.series(self, level="full")

    def _with_native(self, series: pl.Series) -> Self:
        return self.__class__(series, version=self._version)

    @overload
    def _from_native_object(self, series: pl.Series) -> Self: ...

    @overload
    def _from_native_object(self, series: pl.DataFrame) -> PolarsDataFrame: ...

    @overload
    def _from_native_object(self, series: T) -> T: ...

    def _from_native_object(
        self, series: pl.Series | pl.DataFrame | T
    ) -> Self | PolarsDataFrame | T:
        if self._is_native(series):
            return self._with_native(series)
        if isinstance(series, pl.DataFrame):
            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame.from_native(series, context=self)
        # scalar
        return series

    def __getattr__(self, attr: str) -> Any:
        if attr not in INHERITED_METHODS:
            msg = f"{self.__class__.__name__} has not attribute '{attr}'."
            raise AttributeError(msg)

        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            return self._from_native_object(getattr(self.native, attr)(*pos, **kwds))

        return func

    def __len__(self) -> int:
        return len(self.native)

    def __rfloordiv__(self, other: Any) -> PolarsSeries:
        if self._backend_version < (1, 10, 0):
            name = self.name
            ns = self.__narwhals_namespace__()
            return (
                self.to_frame()
                .select((ns.col(name).__rfloordiv__(other)).alias(name))
                .get_column(name)
            )
        return self._with_native(self.native.__rfloordiv__(extract_native(other)))

    @property
    def name(self) -> str:
        return self.native.name

    @property
    def dtype(self) -> DType:
        return native_to_narwhals_dtype(self.native.dtype, self._version)

    @property
    def native(self) -> pl.Series:
        return self._native_series

    def alias(self, name: str) -> Self:
        return self._from_native_object(self.native.alias(name))

    def __getitem__(self, item: MultiIndexSelector[Self]) -> Any | Self:
        if isinstance(item, PolarsSeries):
            return self._from_native_object(self.native.__getitem__(item.native))
        return self._from_native_object(self.native.__getitem__(item))

    def cast(self, dtype: IntoDType) -> Self:
        dtype_pl = narwhals_to_native_dtype(dtype, self._version)
        return self._with_native(self.native.cast(dtype_pl))

    @requires.backend_version((1,))
    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self:
        ser = self.native
        dtype = (
            narwhals_to_native_dtype(return_dtype, self._version)
            if return_dtype
            else None
        )
        return self._with_native(ser.replace_strict(old, new, return_dtype=dtype))

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _1DArray:
        return self.__array__(dtype, copy=copy)

    def __array__(self, dtype: Any, *, copy: bool | None) -> _1DArray:
        if self._backend_version < (0, 20, 29):
            return self.native.__array__(dtype=dtype)
        return self.native.__array__(dtype=dtype, copy=copy)

    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_native(self.native.__eq__(extract_native(other)))

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_native(self.native.__ne__(extract_native(other)))

    # NOTE: These need to be anything that can't match `PolarsExpr`, due to overload order
    def __ge__(self, other: Self) -> Self:
        return self._with_native(self.native.__ge__(extract_native(other)))

    def __gt__(self, other: Self) -> Self:
        return self._with_native(self.native.__gt__(extract_native(other)))

    def __le__(self, other: Self) -> Self:
        return self._with_native(self.native.__le__(extract_native(other)))

    def __lt__(self, other: Self) -> Self:
        return self._with_native(self.native.__lt__(extract_native(other)))

    def __rpow__(self, other: PolarsSeries | Any) -> Self:
        result = self.native.__rpow__(extract_native(other))
        if self._backend_version < (1, 16, 1):
            # Explicitly set alias to work around https://github.com/pola-rs/polars/issues/20071
            result = result.alias(self.name)
        return self._with_native(result)

    def is_nan(self) -> Self:
        try:
            native_is_nan = self.native.is_nan()
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None
        if self._backend_version < (1, 18):  # pragma: no cover
            select = pl.when(self.native.is_not_null()).then(native_is_nan)
            return self._with_native(pl.select(select)[self.name])
        return self._with_native(native_is_nan)

    def median(self) -> Any:
        from narwhals.exceptions import InvalidOperationError

        if not self.dtype.is_numeric():
            msg = "`median` operation not supported for non-numeric input type."
            raise InvalidOperationError(msg)

        return self.native.median()

    def to_dummies(self, *, separator: str, drop_first: bool) -> PolarsDataFrame:
        from narwhals._polars.dataframe import PolarsDataFrame

        if self._backend_version < (0, 20, 15):
            has_nulls = self.native.is_null().any()
            result = self.native.to_dummies(separator=separator)
            output_columns = result.columns
            if drop_first:
                _ = output_columns.pop(int(has_nulls))

            result = result.select(output_columns)
        else:
            result = self.native.to_dummies(separator=separator, drop_first=drop_first)
        result = result.with_columns(pl.all().cast(pl.Int8))
        return PolarsDataFrame.from_native(result, context=self)

    def ewm_mean(
        self,
        *,
        com: float | None,
        span: float | None,
        half_life: float | None,
        alpha: float | None,
        adjust: bool,
        min_samples: int,
        ignore_nulls: bool,
    ) -> Self:
        extra_kwargs = (
            {"min_periods": min_samples}
            if self._backend_version < (1, 21, 0)
            else {"min_samples": min_samples}
        )

        native_result = self.native.ewm_mean(
            com=com,
            span=span,
            half_life=half_life,
            alpha=alpha,
            adjust=adjust,
            ignore_nulls=ignore_nulls,
            **extra_kwargs,
        )
        if self._backend_version < (1,):  # pragma: no cover
            return self._with_native(
                pl.select(
                    pl.when(~self.native.is_null()).then(native_result).otherwise(None)
                )[self.native.name]
            )

        return self._with_native(native_result)

    @requires.backend_version((1,))
    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        extra_kwargs: dict[str, Any] = (
            {"min_periods": min_samples}
            if self._backend_version < (1, 21, 0)
            else {"min_samples": min_samples}
        )
        return self._with_native(
            self.native.rolling_var(
                window_size=window_size, center=center, ddof=ddof, **extra_kwargs
            )
        )

    @requires.backend_version((1,))
    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        extra_kwargs: dict[str, Any] = (
            {"min_periods": min_samples}
            if self._backend_version < (1, 21, 0)
            else {"min_samples": min_samples}
        )
        return self._with_native(
            self.native.rolling_std(
                window_size=window_size, center=center, ddof=ddof, **extra_kwargs
            )
        )

    def rolling_sum(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        extra_kwargs: dict[str, Any] = (
            {"min_periods": min_samples}
            if self._backend_version < (1, 21, 0)
            else {"min_samples": min_samples}
        )
        return self._with_native(
            self.native.rolling_sum(
                window_size=window_size, center=center, **extra_kwargs
            )
        )

    def rolling_mean(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        extra_kwargs: dict[str, Any] = (
            {"min_periods": min_samples}
            if self._backend_version < (1, 21, 0)
            else {"min_samples": min_samples}
        )
        return self._with_native(
            self.native.rolling_mean(
                window_size=window_size, center=center, **extra_kwargs
            )
        )

    def sort(self, *, descending: bool, nulls_last: bool) -> Self:
        if self._backend_version < (0, 20, 6):
            result = self.native.sort(descending=descending)

            if nulls_last:
                is_null = result.is_null()
                result = pl.concat([result.filter(~is_null), result.filter(is_null)])
        else:
            result = self.native.sort(descending=descending, nulls_last=nulls_last)

        return self._with_native(result)

    def scatter(self, indices: int | Sequence[int], values: Any) -> Self:
        s = self.native.clone().scatter(indices, extract_native(values))
        return self._with_native(s)

    def value_counts(
        self, *, sort: bool, parallel: bool, name: str | None, normalize: bool
    ) -> PolarsDataFrame:
        from narwhals._polars.dataframe import PolarsDataFrame

        if self._backend_version < (1, 0, 0):
            value_name_ = name or ("proportion" if normalize else "count")

            result = self.native.value_counts(sort=sort, parallel=parallel).select(
                **{
                    (self.native.name): pl.col(self.native.name),
                    value_name_: pl.col("count") / pl.sum("count")
                    if normalize
                    else pl.col("count"),
                }
            )
        else:
            result = self.native.value_counts(
                sort=sort, parallel=parallel, name=name, normalize=normalize
            )
        return PolarsDataFrame.from_native(result, context=self)

    def cum_count(self, *, reverse: bool) -> Self:
        return self._with_native(self.native.cum_count(reverse=reverse))

    def __contains__(self, other: Any) -> bool:
        try:
            return self.native.__contains__(other)
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float,
        rel_tol: float,
        nans_equal: bool,
    ) -> PolarsSeries:
        if self._backend_version < (1, 32, 0):
            name = self.name
            ns = self.__narwhals_namespace__()
            other_expr = (
                ns.lit(other.native, None) if isinstance(other, PolarsSeries) else other
            )
            expr = ns.col(name).is_close(
                other_expr, abs_tol=abs_tol, rel_tol=rel_tol, nans_equal=nans_equal
            )
            return self.to_frame().select(expr).get_column(name)
        other_series = other.native if isinstance(other, PolarsSeries) else other
        result = self.native.is_close(
            other_series, abs_tol=abs_tol, rel_tol=rel_tol, nans_equal=nans_equal
        )
        return self._with_native(result)

    def mode(self, *, keep: ModeKeepStrategy) -> Self:
        result = self.native.mode()
        return self._with_native(result.head(1) if keep == "any" else result)

    def hist_from_bins(
        self, bins: list[float], *, include_breakpoint: bool
    ) -> PolarsDataFrame:
        if len(bins) <= 1:
            native = pl.DataFrame(schema=self._HIST_EMPTY_SCHEMA[include_breakpoint])
        elif self.native.is_empty():
            if include_breakpoint:
                native = (
                    pl.Series(bins[1:])
                    .to_frame("breakpoint")
                    .with_columns(count=pl.lit(0, pl.Int64))
                )
            else:
                native = pl.select(count=pl.zeros(len(bins) - 1, pl.Int64))
        else:
            return self._hist_from_data(
                bins=bins, bin_count=None, include_breakpoint=include_breakpoint
            )
        return self.__narwhals_namespace__()._dataframe.from_native(native, context=self)

    def hist_from_bin_count(
        self, bin_count: int, *, include_breakpoint: bool
    ) -> PolarsDataFrame:
        if bin_count == 0:
            native = pl.DataFrame(schema=self._HIST_EMPTY_SCHEMA[include_breakpoint])
        elif self.native.is_empty():
            if include_breakpoint:
                native = pl.select(
                    breakpoint=pl.int_range(1, bin_count + 1) / bin_count,
                    count=pl.lit(0, pl.Int64),
                )
            else:
                native = pl.select(count=pl.zeros(bin_count, pl.Int64))
        else:
            count: int | None
            if BACKEND_VERSION < (1, 15):  # pragma: no cover
                count = None
                bins = self._bins_from_bin_count(bin_count=bin_count)
            else:
                count = bin_count
                bins = None
            return self._hist_from_data(
                bins=bins,  # type: ignore[arg-type]
                bin_count=count,
                include_breakpoint=include_breakpoint,
            )
        return self.__narwhals_namespace__()._dataframe.from_native(native, context=self)

    def _bins_from_bin_count(self, bin_count: int) -> pl.Series:  # pragma: no cover
        """Prepare bins based on backend version compatibility.

        polars <1.15 does not adjust the bins when they have equivalent min/max
        polars <1.5 with bin_count=...
        returns bins that range from -inf to +inf and has bin_count + 1 bins.
          for compat: convert `bin_count=` call to `bins=`
        """
        lower = cast("float", self.native.min())
        upper = cast("float", self.native.max())

        if lower == upper:
            lower -= 0.5
            upper += 0.5

        width = (upper - lower) / bin_count
        return pl.int_range(0, bin_count + 1, eager=True) * width + lower

    def _hist_from_data(
        self, bins: list[float] | None, bin_count: int | None, *, include_breakpoint: bool
    ) -> PolarsDataFrame:
        """Calculate histogram from non-empty data and post-process the results based on the backend version."""
        from narwhals._polars.dataframe import PolarsDataFrame

        series = self.native

        # Polars inconsistently handles NaN values when computing histograms
        #   against predefined bins: https://github.com/pola-rs/polars/issues/21082
        if BACKEND_VERSION < (1, 15) or bins is not None:
            series = series.fill_nan(None)

        df = series.hist(
            bins,
            bin_count=bin_count,
            include_category=False,
            include_breakpoint=include_breakpoint,
        )

        # Apply post-processing corrections

        # Handle column naming
        if not include_breakpoint:
            col_name = df.columns[0]
            df = df.select(pl.col(col_name).alias("count"))
        elif BACKEND_VERSION < (1, 0):  # pragma: no cover
            df = df.rename({"break_point": "breakpoint"})

        if bins is not None:  # pragma: no cover
            # polars<1.6 implicitly adds -inf and inf to either end of bins
            if BACKEND_VERSION < (1, 6):
                r = pl.int_range(0, len(df))
                df = df.filter((r > 0) & (r < len(df) - 1))
            # polars<1.27 makes the lowest bin a left/right closed interval
            if BACKEND_VERSION < (1, 27):
                df = (
                    df.slice(0, 1)
                    .with_columns(pl.col("count") + ((pl.lit(series) == bins[0]).sum()))
                    .vstack(df.slice(1))
                )

        return PolarsDataFrame.from_native(df, context=self)

    def to_polars(self) -> pl.Series:
        return self.native

    def first(self) -> PythonLiteral:
        if self._backend_version < (1, 10):  # pragma: no cover
            return self.native.item(0) if len(self) else None
        return self.native.first()  # type: ignore[return-value]

    def last(self) -> PythonLiteral:
        if self._backend_version < (1, 10):  # pragma: no cover
            return self.native.item(-1) if len(self) else None
        return self.native.last()  # type: ignore[return-value]

    @property
    def dt(self) -> PolarsSeriesDateTimeNamespace:
        return PolarsSeriesDateTimeNamespace(self)

    @property
    def str(self) -> PolarsSeriesStringNamespace:
        return PolarsSeriesStringNamespace(self)

    @property
    def cat(self) -> PolarsSeriesCatNamespace:
        return PolarsSeriesCatNamespace(self)

    @property
    def struct(self) -> PolarsSeriesStructNamespace:
        return PolarsSeriesStructNamespace(self)

    __add__: Method[Self]
    __and__: Method[Self]
    __floordiv__: Method[Self]
    __invert__: Method[Self]
    __iter__: Method[Iterator[Any]]
    __mod__: Method[Self]
    __mul__: Method[Self]
    __or__: Method[Self]
    __pow__: Method[Self]
    __radd__: Method[Self]
    __rand__: Method[Self]
    __rmod__: Method[Self]
    __rmul__: Method[Self]
    __ror__: Method[Self]
    __rsub__: Method[Self]
    __rtruediv__: Method[Self]
    __sub__: Method[Self]
    __truediv__: Method[Self]
    abs: Method[Self]
    all: Method[bool]
    any: Method[bool]
    arg_max: Method[int]
    arg_min: Method[int]
    arg_true: Method[Self]
    ceil: Method[Self]
    clip: Method[Self]
    count: Method[int]
    cum_max: Method[Self]
    cum_min: Method[Self]
    cum_prod: Method[Self]
    cum_sum: Method[Self]
    diff: Method[Self]
    drop_nulls: Method[Self]
    exp: Method[Self]
    fill_null: Method[Self]
    fill_nan: Method[Self]
    filter: Method[Self]
    floor: Method[Self]
    gather_every: Method[Self]
    head: Method[Self]
    is_between: Method[Self]
    is_duplicated: Method[Self]
    is_empty: Method[bool]
    is_finite: Method[Self]
    is_first_distinct: Method[Self]
    is_in: Method[Self]
    is_last_distinct: Method[Self]
    is_null: Method[Self]
    is_sorted: Method[bool]
    is_unique: Method[Self]
    item: Method[Any]
    kurtosis: Method[float | None]
    len: Method[int]
    log: Method[Self]
    max: Method[Any]
    mean: Method[float]
    min: Method[Any]
    n_unique: Method[int]
    null_count: Method[int]
    quantile: Method[float]
    rank: Method[Self]
    round: Method[Self]
    sample: Method[Self]
    shift: Method[Self]
    skew: Method[float | None]
    sqrt: Method[Self]
    std: Method[float]
    sum: Method[float]
    tail: Method[Self]
    to_arrow: Method[pa.Array[Any]]
    to_frame: Method[PolarsDataFrame]
    to_list: Method[list[Any]]
    to_pandas: Method[pd.Series[Any]]
    unique: Method[Self]
    var: Method[float]
    zip_with: Method[Self]

    @property
    def list(self) -> PolarsSeriesListNamespace:
        return PolarsSeriesListNamespace(self)


class PolarsSeriesNamespace(PolarsAnyNamespace[PolarsSeries, pl.Series]):
    def __init__(self, series: PolarsSeries) -> None:
        self._series = series

    @property
    def compliant(self) -> PolarsSeries:
        return self._series

    @property
    def native(self) -> pl.Series:
        return self._series.native

    @property
    def name(self) -> str:
        return self.compliant.name

    def __narwhals_namespace__(self) -> PolarsNamespace:
        return self.compliant.__narwhals_namespace__()

    def to_frame(self) -> PolarsDataFrame:
        return self.compliant.to_frame()


class PolarsSeriesDateTimeNamespace(
    PolarsSeriesNamespace, PolarsDateTimeNamespace[PolarsSeries, pl.Series]
): ...


class PolarsSeriesStringNamespace(
    PolarsSeriesNamespace, PolarsStringNamespace[PolarsSeries, pl.Series]
):
    def to_titlecase(self) -> PolarsSeries:
        name = self.name
        ns = self.__narwhals_namespace__()
        return self.to_frame().select(ns.col(name).str.to_titlecase()).get_column(name)

    def zfill(self, width: int) -> PolarsSeries:
        name = self.name
        ns = self.__narwhals_namespace__()
        return self.to_frame().select(ns.col(name).str.zfill(width)).get_column(name)


class PolarsSeriesCatNamespace(
    PolarsSeriesNamespace, PolarsCatNamespace[PolarsSeries, pl.Series]
): ...


class PolarsSeriesListNamespace(
    PolarsSeriesNamespace, PolarsListNamespace[PolarsSeries, pl.Series]
):
    def len(self) -> PolarsSeries:
        name = self.name
        ns = self.__narwhals_namespace__()
        return self.to_frame().select(ns.col(name).list.len()).get_column(name)

    def contains(self, item: NonNestedLiteral) -> PolarsSeries:
        name = self.name
        ns = self.__narwhals_namespace__()
        return self.to_frame().select(ns.col(name).list.contains(item)).get_column(name)


class PolarsSeriesStructNamespace(
    PolarsSeriesNamespace, PolarsStructNamespace[PolarsSeries, pl.Series]
): ...
