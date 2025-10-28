from __future__ import annotations

import operator
import warnings
from typing import TYPE_CHECKING, Any, Callable, Literal, cast

import numpy as np

from narwhals._compliant import EagerSeries, EagerSeriesHist
from narwhals._pandas_like.series_cat import PandasLikeSeriesCatNamespace
from narwhals._pandas_like.series_dt import PandasLikeSeriesDateTimeNamespace
from narwhals._pandas_like.series_list import PandasLikeSeriesListNamespace
from narwhals._pandas_like.series_str import PandasLikeSeriesStringNamespace
from narwhals._pandas_like.series_struct import PandasLikeSeriesStructNamespace
from narwhals._pandas_like.utils import (
    align_and_extract_native,
    get_dtype_backend,
    import_array_module,
    narwhals_to_native_dtype,
    native_to_narwhals_dtype,
    object_native_to_narwhals_dtype,
    rename,
    select_columns_by_name,
    set_index,
)
from narwhals._typing_compat import assert_never
from narwhals._utils import Implementation, is_list_of, parse_version
from narwhals.dependencies import is_numpy_array_1d, is_pandas_like_series
from narwhals.exceptions import InvalidOperationError

if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable, Iterator, Mapping, Sequence
    from types import ModuleType

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._arrow.typing import ChunkedArrayAny
    from narwhals._compliant.series import HistData
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals._pandas_like.namespace import PandasLikeNamespace
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.typing import (
        ClosedInterval,
        FillNullStrategy,
        Into1DArray,
        IntoDType,
        ModeKeepStrategy,
        NonNestedLiteral,
        NumericLiteral,
        PythonLiteral,
        RankMethod,
        RollingInterpolationMethod,
        SizedMultiIndexSelector,
        TemporalLiteral,
        _1DArray,
        _SliceIndex,
    )

    PandasHistData: TypeAlias = "HistData[pd.Series[Any], list[float]]"


PANDAS_TO_NUMPY_DTYPE_NO_MISSING = {
    "Int64": "int64",
    "int64[pyarrow]": "int64",
    "Int32": "int32",
    "int32[pyarrow]": "int32",
    "Int16": "int16",
    "int16[pyarrow]": "int16",
    "Int8": "int8",
    "int8[pyarrow]": "int8",
    "UInt64": "uint64",
    "uint64[pyarrow]": "uint64",
    "UInt32": "uint32",
    "uint32[pyarrow]": "uint32",
    "UInt16": "uint16",
    "uint16[pyarrow]": "uint16",
    "UInt8": "uint8",
    "uint8[pyarrow]": "uint8",
    "Float64": "float64",
    "float64[pyarrow]": "float64",
    "Float32": "float32",
    "float32[pyarrow]": "float32",
}
PANDAS_TO_NUMPY_DTYPE_MISSING = {
    "Int64": "float64",
    "int64[pyarrow]": "float64",
    "Int32": "float64",
    "int32[pyarrow]": "float64",
    "Int16": "float64",
    "int16[pyarrow]": "float64",
    "Int8": "float64",
    "int8[pyarrow]": "float64",
    "UInt64": "float64",
    "uint64[pyarrow]": "float64",
    "UInt32": "float64",
    "uint32[pyarrow]": "float64",
    "UInt16": "float64",
    "uint16[pyarrow]": "float64",
    "UInt8": "float64",
    "uint8[pyarrow]": "float64",
    "Float64": "float64",
    "float64[pyarrow]": "float64",
    "Float32": "float32",
    "float32[pyarrow]": "float32",
}


class PandasLikeSeries(EagerSeries[Any]):
    def __init__(
        self, native_series: Any, *, implementation: Implementation, version: Version
    ) -> None:
        self._name = native_series.name
        self._native_series = native_series
        self._implementation = implementation
        self._version = version
        # Flag which indicates if, in the final step before applying an operation,
        # the single value behind the PandasLikeSeries should be extract and treated
        # as a Scalar. For example, in `nw.col('a') - nw.lit(3)`, the latter would
        # become a Series of length 1. Rather that doing a full broadcast so it matches
        # the length of the whole dataframe, we just extract the scalar.
        self._broadcast = False

    @property
    def native(self) -> Any:
        return self._native_series

    def __native_namespace__(self) -> ModuleType:
        if self._implementation.is_pandas_like():
            return self._implementation.to_native_namespace()

        msg = f"Expected pandas/modin/cudf, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def __narwhals_namespace__(self) -> PandasLikeNamespace:
        from narwhals._pandas_like.namespace import PandasLikeNamespace

        return PandasLikeNamespace(self._implementation, self._version)

    def _gather(self, rows: SizedMultiIndexSelector[pd.Series[Any]]) -> Self:
        rows = list(rows) if isinstance(rows, tuple) else rows
        return self._with_native(self.native.iloc[rows])

    def _gather_slice(self, rows: _SliceIndex | range) -> Self:
        return self._with_native(
            self.native.iloc[slice(rows.start, rows.stop, rows.step)]
        )

    def _with_version(self, version: Version) -> Self:
        return self.__class__(
            self.native, implementation=self._implementation, version=version
        )

    def _with_native(self, series: Any, *, preserve_broadcast: bool = False) -> Self:
        result = self.__class__(
            series, implementation=self._implementation, version=self._version
        )
        if preserve_broadcast:
            result._broadcast = self._broadcast
        return result

    @classmethod
    def from_iterable(
        cls,
        data: Iterable[Any],
        *,
        context: _LimitedContext,
        name: str = "",
        dtype: IntoDType | None = None,
        index: Any = None,
    ) -> Self:
        implementation = context._implementation
        version = context._version
        ns = implementation.to_native_namespace()
        kwds: dict[str, Any] = {}
        if dtype:
            kwds["dtype"] = narwhals_to_native_dtype(dtype, None, implementation, version)
        else:
            if implementation.is_pandas():
                kwds["copy"] = False
            if index is not None and len(index):
                kwds["index"] = index
        return cls.from_native(ns.Series(data, name=name, **kwds), context=context)

    @staticmethod
    def _is_native(obj: Any) -> TypeIs[Any]:
        return is_pandas_like_series(obj)  # pragma: no cover

    @classmethod
    def from_native(cls, data: Any, /, *, context: _LimitedContext) -> Self:
        return cls(data, implementation=context._implementation, version=context._version)

    @classmethod
    def from_numpy(cls, data: Into1DArray, /, *, context: _LimitedContext) -> Self:
        implementation = context._implementation
        arr = data if is_numpy_array_1d(data) else [data]
        native = implementation.to_native_namespace().Series(arr, name="")
        return cls.from_native(native, context=context)

    @classmethod
    def _align_full_broadcast(cls, *series: Self) -> Sequence[Self]:
        Series = series[0].__native_namespace__().Series
        lengths = [len(s) for s in series]
        max_length = max(lengths)
        idx = series[lengths.index(max_length)].native.index
        reindexed = []
        for s in series:
            if s._broadcast:
                native = Series(
                    s.native.iloc[0], index=idx, name=s.name, dtype=s.native.dtype
                )
                compliant = s._with_native(native)
            elif s.native.index is not idx:
                native = set_index(s.native, idx, implementation=s._implementation)
                compliant = s._with_native(native)
            else:
                compliant = s
            reindexed.append(compliant)
        return reindexed

    @property
    def name(self) -> str:
        return self._name

    @property
    def dtype(self) -> DType:
        native_dtype = self.native.dtype
        return (
            native_to_narwhals_dtype(native_dtype, self._version, self._implementation)
            if native_dtype != "object"
            else object_native_to_narwhals_dtype(
                self.native, self._version, self._implementation
            )
        )

    @property
    def _array_funcs(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            import numpy as np

            return np
        return import_array_module(self._implementation)

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
        ser = self.native
        mask_na = ser.isna()
        if self._implementation is Implementation.CUDF:
            if (min_samples == 0 and not ignore_nulls) or (not mask_na.any()):
                result = ser.ewm(
                    com=com, span=span, halflife=half_life, alpha=alpha, adjust=adjust
                ).mean()
            else:
                msg = (
                    "cuDF only supports `ewm_mean` when there are no missing values "
                    "or when both `min_period=0` and `ignore_nulls=False`"
                )
                raise NotImplementedError(msg)
        else:
            result = ser.ewm(
                com, span, half_life, alpha, min_samples, adjust, ignore_na=ignore_nulls
            ).mean()
        result[mask_na] = None
        return self._with_native(result)

    def scatter(self, indices: int | Sequence[int], values: Any) -> Self:
        if isinstance(values, self.__class__):
            values = set_index(
                values.native,
                self.native.index[indices],
                implementation=self._implementation,
            )
        s = self.native.copy(deep=True)
        s.iloc[indices] = values
        s.name = self.name
        return self._with_native(s)

    def _scatter_in_place(self, indices: Self, values: Self) -> None:
        # Scatter, modifying original Series. Use with care!
        implementation = self._implementation
        backend_version = self._backend_version
        values_native = set_index(
            values.native,
            self.native.index[indices.native],
            implementation=implementation,
        )
        if implementation is Implementation.PANDAS and parse_version(np) < (2,):
            values_native = values_native.copy()  # pragma: no cover
        min_pd_version = (1, 2)
        if implementation is Implementation.PANDAS and backend_version < min_pd_version:
            self.native.iloc[indices.native.values] = values_native  # noqa: PD011
        else:
            self.native.iloc[indices.native] = values_native

    def cast(self, dtype: IntoDType) -> Self:
        if self.dtype == dtype and self.native.dtype != "object":
            # Avoid dealing with pandas' type-system if we can. Note that it's only
            # safe to do this if we're not starting with object dtype, see tests/expr_and_series/cast_test.py::test_cast_object_pandas
            # for an example of why.
            return self._with_native(self.native, preserve_broadcast=True)
        pd_dtype = narwhals_to_native_dtype(
            dtype,
            dtype_backend=get_dtype_backend(self.native.dtype, self._implementation),
            implementation=self._implementation,
            version=self._version,
        )
        return self._with_native(self.native.astype(pd_dtype), preserve_broadcast=True)

    def item(self, index: int | None = None) -> Any:
        # cuDF doesn't have Series.item().
        if index is None:
            if len(self) != 1:
                msg = (
                    "can only call '.item()' if the Series is of length 1,"
                    f" or an explicit index is provided (Series is of length {len(self)})"
                )
                raise ValueError(msg)
            return self.native.iloc[0]
        return self.native.iloc[index]

    def to_frame(self) -> PandasLikeDataFrame:
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        return PandasLikeDataFrame(
            self.native.to_frame(),
            implementation=self._implementation,
            version=self._version,
            validate_column_names=False,
        )

    def to_list(self) -> list[Any]:
        is_cudf = self._implementation.is_cudf()
        return self.native.to_arrow().to_pylist() if is_cudf else self.native.to_list()

    def is_between(
        self, lower_bound: Any, upper_bound: Any, closed: ClosedInterval
    ) -> Self:
        ser = self.native
        _, lower_bound = align_and_extract_native(self, lower_bound)
        _, upper_bound = align_and_extract_native(self, upper_bound)
        if closed == "left":
            res = ser.ge(lower_bound) & ser.lt(upper_bound)
        elif closed == "right":
            res = ser.gt(lower_bound) & ser.le(upper_bound)
        elif closed == "none":
            res = ser.gt(lower_bound) & ser.lt(upper_bound)
        elif closed == "both":
            res = ser.ge(lower_bound) & ser.le(upper_bound)
        else:
            assert_never(closed)
        return self._with_native(res).alias(ser.name)

    def is_in(self, other: Any) -> Self:
        return self._with_native(self.native.isin(other))

    def arg_true(self) -> Self:
        ser = self.native
        size = len(ser)
        data = self._array_funcs.arange(size)
        result = ser.__class__(data, name=ser.name, index=ser.index).loc[ser]
        return self._with_native(result)

    def arg_min(self) -> int:
        return self.native.argmin()

    def arg_max(self) -> int:
        return self.native.argmax()

    # Binary comparisons

    def filter(self, predicate: Any) -> Self:
        if not is_list_of(predicate, bool):
            _, other_native = align_and_extract_native(self, predicate)
        else:
            other_native = predicate
        return self._with_native(self.native.loc[other_native]).alias(self.name)

    def first(self) -> PythonLiteral:
        return self.native.iloc[0] if len(self.native) else None

    def last(self) -> PythonLiteral:
        return self.native.iloc[-1] if len(self.native) else None

    def _with_binary(self, op: Callable[..., PandasLikeSeries], other: Any) -> Self:
        ser, other_native = align_and_extract_native(self, other)
        preserve_broadcast = self._broadcast and getattr(other, "_broadcast", True)
        return self._with_native(
            op(ser, other_native), preserve_broadcast=preserve_broadcast
        ).alias(self.name)

    def _with_binary_right(self, op: Callable[..., PandasLikeSeries], other: Any) -> Self:
        return self._with_binary(lambda x, y: op(y, x), other).alias(self.name)

    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_binary(operator.eq, other)

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_binary(operator.ne, other)

    def __ge__(self, other: Any) -> Self:
        return self._with_binary(operator.ge, other)

    def __gt__(self, other: Any) -> Self:
        return self._with_binary(operator.gt, other)

    def __le__(self, other: Any) -> Self:
        return self._with_binary(operator.le, other)

    def __lt__(self, other: Any) -> Self:
        return self._with_binary(operator.lt, other)

    def __and__(self, other: Any) -> Self:
        return self._with_binary(operator.and_, other)

    def __rand__(self, other: Any) -> Self:
        return self._with_binary_right(operator.and_, other)

    def __or__(self, other: Any) -> Self:
        return self._with_binary(operator.or_, other)

    def __ror__(self, other: Any) -> Self:
        return self._with_binary_right(operator.or_, other)

    def __add__(self, other: Any) -> Self:
        return self._with_binary(operator.add, other)

    def __radd__(self, other: Any) -> Self:
        return self._with_binary_right(operator.add, other)

    def __sub__(self, other: Any) -> Self:
        return self._with_binary(operator.sub, other)

    def __rsub__(self, other: Any) -> Self:
        return self._with_binary_right(operator.sub, other)

    def __mul__(self, other: Any) -> Self:
        return self._with_binary(operator.mul, other)

    def __rmul__(self, other: Any) -> Self:
        return self._with_binary_right(operator.mul, other)

    def __truediv__(self, other: Any) -> Self:
        return self._with_binary(operator.truediv, other)

    def __rtruediv__(self, other: Any) -> Self:
        return self._with_binary_right(operator.truediv, other)

    def __floordiv__(self, other: Any) -> Self:
        return self._with_binary(operator.floordiv, other)

    def __rfloordiv__(self, other: Any) -> Self:
        return self._with_binary_right(operator.floordiv, other)

    def __pow__(self, other: Any) -> Self:
        return self._with_binary(operator.pow, other)

    def __rpow__(self, other: Any) -> Self:
        return self._with_binary_right(operator.pow, other)

    def __mod__(self, other: Any) -> Self:
        return self._with_binary(operator.mod, other)

    def __rmod__(self, other: Any) -> Self:
        return self._with_binary_right(operator.mod, other)

    # Unary

    def __invert__(self) -> Self:
        return self._with_native(~self.native)

    # Reductions

    def any(self) -> bool:
        return self.native.any()

    def all(self) -> bool:
        return self.native.all()

    def min(self) -> Any:
        return self.native.min()

    def max(self) -> Any:
        return self.native.max()

    def sum(self) -> float:
        return self.native.sum()

    def count(self) -> int:
        return self.native.count()

    def mean(self) -> float:
        return self.native.mean()

    def median(self) -> float:
        if not self.dtype.is_numeric():
            msg = "`median` operation not supported for non-numeric input type."
            raise InvalidOperationError(msg)
        return self.native.median()

    def std(self, *, ddof: int) -> float:
        return self.native.std(ddof=ddof)

    def var(self, *, ddof: int) -> float:
        return self.native.var(ddof=ddof)

    def skew(self) -> float | None:
        ser_not_null = self.native.dropna()
        if len(ser_not_null) == 0:
            return None
        if len(ser_not_null) == 1:
            return float("nan")
        if len(ser_not_null) == 2:
            return 0.0
        m = ser_not_null - ser_not_null.mean()
        m2 = (m**2).mean()
        m3 = (m**3).mean()
        return m3 / (m2**1.5) if m2 != 0 else float("nan")

    def kurtosis(self) -> float | None:
        ser_not_null = self.native.dropna()
        if len(ser_not_null) == 0:
            return None
        if len(ser_not_null) == 1:
            return float("nan")
        m = ser_not_null - ser_not_null.mean()
        m2 = (m**2).mean()
        m4 = (m**4).mean()
        return m4 / (m2**2) - 3.0 if m2 != 0 else float("nan")

    def len(self) -> int:
        return len(self.native)

    # Transformations

    def is_null(self) -> Self:
        return self._with_native(self.native.isna(), preserve_broadcast=True)

    def is_nan(self) -> Self:
        ser = self.native
        if not self.dtype.is_numeric():
            msg = f"`.is_nan` only supported for numeric dtype and not {self.dtype}, did you mean `.is_null`?"
            raise InvalidOperationError(msg)
        # If/when pandas exposes an API which distinguishes NaN vs null, use that.
        return self._with_native(ser != ser, preserve_broadcast=True)  # noqa: PLR0124

    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self:
        ser = self.native
        kwargs = (
            {"downcast": False}
            if self._implementation is Implementation.PANDAS
            and self._backend_version < (3,)
            else {}
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "The 'downcast' keyword .*is deprecated", category=FutureWarning
            )
            if value is not None:
                _, native_value = align_and_extract_native(self, value)
                res_ser = self._with_native(
                    ser.fillna(value=native_value, **kwargs), preserve_broadcast=True
                )
            else:
                res_ser = self._with_native(
                    ser.ffill(limit=limit, **kwargs)
                    if strategy == "forward"
                    else ser.bfill(limit=limit, **kwargs),
                    preserve_broadcast=True,
                )
        return res_ser

    def fill_nan(self, value: float | None) -> Self:
        if not self.dtype.is_numeric():  # pragma: no cover
            msg = f"`.fill_nan` only supported for numeric dtype and not {self.dtype}, did you mean `.fill_null`?"
            raise InvalidOperationError(msg)
        s = self.native
        fill = s.array.dtype.na_value if value is None else value
        # If/when pandas exposes an API which distinguishes NaN vs null, use that.
        mask = s != s  # noqa: PLR0124
        # Carefully use `inplace`, as `mask` isn't provided by the user.
        mask.fillna(False, inplace=True)  # noqa: PD002
        return self._with_native(s.mask(mask, fill), preserve_broadcast=True)

    def drop_nulls(self) -> Self:
        return self._with_native(self.native.dropna())

    def n_unique(self) -> int:
        return self.native.nunique(dropna=False)

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
            )
        )

    def abs(self) -> Self:
        return self._with_native(self.native.abs())

    def cum_sum(self, *, reverse: bool) -> Self:
        result = (
            self.native.cumsum(skipna=True)
            if not reverse
            else self.native[::-1].cumsum(skipna=True)[::-1]
        )
        return self._with_native(result)

    def unique(self, *, maintain_order: bool = True) -> Self:
        """Pandas always maintains order, as per its docstring.

        > Uniques are returned in order of appearance.
        """
        return self._with_native(type(self.native)(self.native.unique(), name=self.name))

    def diff(self) -> Self:
        return self._with_native(self.native.diff())

    def shift(self, n: int) -> Self:
        return self._with_native(self.native.shift(n))

    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> PandasLikeSeries:
        tmp_name = f"{self.name}_tmp"
        dtype_backend = get_dtype_backend(self.native.dtype, self._implementation)
        dtype = (
            narwhals_to_native_dtype(
                return_dtype, dtype_backend, self._implementation, self._version
            )
            if return_dtype
            else None
        )
        namespace = self.__native_namespace__()
        other = namespace.DataFrame(
            {self.name: old, tmp_name: namespace.Series(new, dtype=dtype)}
        )
        result = self._with_native(
            self.native.to_frame().merge(other, on=self.name, how="left")[tmp_name]
        ).alias(self.name)
        if result.is_null().sum() != self.is_null().sum():
            msg = (
                "replace_strict did not replace all non-null values.\n\n"
                f"The following did not get replaced: {self.filter(~self.is_null() & result.is_null()).unique(maintain_order=False).to_list()}"
            )
            raise ValueError(msg)
        return result

    def sort(self, *, descending: bool, nulls_last: bool) -> PandasLikeSeries:
        na_position = "last" if nulls_last else "first"
        return self._with_native(
            self.native.sort_values(ascending=not descending, na_position=na_position)
        ).alias(self.name)

    def alias(self, name: str | Hashable) -> Self:
        if name != self.name:
            return self._with_native(
                rename(self.native, name, implementation=self._implementation),
                preserve_broadcast=True,
            )
        return self

    def __array__(self, dtype: Any, *, copy: bool | None) -> _1DArray:
        # pandas used to always return object dtype for nullable dtypes.
        # So, we intercept __array__ and pass to `to_numpy` ourselves to make
        # sure an appropriate numpy dtype is returned.
        return self.to_numpy(dtype=dtype, copy=copy)

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _1DArray:
        # the default is meant to be None, but pandas doesn't allow it?
        # https://numpy.org/doc/stable/reference/generated/numpy.ndarray.__array__.html
        dtypes = self._version.dtypes
        if isinstance(self.dtype, dtypes.Datetime) and self.dtype.time_zone is not None:
            s = self.dt.convert_time_zone("UTC").dt.replace_time_zone(None).native
        else:
            s = self.native

        has_missing = s.isna().any()
        kwargs: dict[Any, Any] = {"copy": copy or self._implementation.is_cudf()}
        if has_missing and str(s.dtype) in PANDAS_TO_NUMPY_DTYPE_MISSING:
            kwargs.update({"na_value": float("nan")})
            dtype = dtype or PANDAS_TO_NUMPY_DTYPE_MISSING[str(s.dtype)]
        if not has_missing and str(s.dtype) in PANDAS_TO_NUMPY_DTYPE_NO_MISSING:
            dtype = dtype or PANDAS_TO_NUMPY_DTYPE_NO_MISSING[str(s.dtype)]
        return s.to_numpy(dtype=dtype, **kwargs)

    def to_pandas(self) -> pd.Series[Any]:
        if self._implementation is Implementation.PANDAS:
            return self.native
        if self._implementation is Implementation.CUDF:  # pragma: no cover
            return self.native.to_pandas()
        if self._implementation is Implementation.MODIN:
            return self.native._to_pandas()
        msg = f"Unknown implementation: {self._implementation}"  # pragma: no cover
        raise AssertionError(msg)

    def to_polars(self) -> pl.Series:
        import polars as pl  # ignore-banned-import

        return pl.from_pandas(self.to_pandas())

    # --- descriptive ---
    def is_unique(self) -> Self:
        return self._with_native(~self.native.duplicated(keep=False)).alias(self.name)

    def null_count(self) -> int:
        return self.native.isna().sum()

    def is_first_distinct(self) -> Self:
        return self._with_native(~self.native.duplicated(keep="first")).alias(self.name)

    def is_last_distinct(self) -> Self:
        return self._with_native(~self.native.duplicated(keep="last")).alias(self.name)

    def is_sorted(self, *, descending: bool) -> bool:
        if not isinstance(descending, bool):
            msg = f"argument 'descending' should be boolean, found {type(descending)}"
            raise TypeError(msg)

        if descending:
            return self.native.is_monotonic_decreasing
        return self.native.is_monotonic_increasing

    def value_counts(
        self, *, sort: bool, parallel: bool, name: str | None, normalize: bool
    ) -> PandasLikeDataFrame:
        """Parallel is unused, exists for compatibility."""
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        index_name_ = "index" if self._name is None else self._name
        value_name_ = name or ("proportion" if normalize else "count")
        val_count = self.native.value_counts(
            dropna=False, sort=False, normalize=normalize
        ).reset_index()

        val_count.columns = [index_name_, value_name_]

        if sort:
            val_count = val_count.sort_values(value_name_, ascending=False)

        return PandasLikeDataFrame.from_native(val_count, context=self)

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> float:
        return self.native.quantile(q=quantile, interpolation=interpolation)

    def zip_with(self, mask: Any, other: Any) -> Self:
        ser = self.native
        _, mask = align_and_extract_native(self, mask)
        _, other = align_and_extract_native(self, other)
        res = ser.where(mask, other)
        return self._with_native(res)

    def head(self, n: int) -> Self:
        return self._with_native(self.native.head(n))

    def tail(self, n: int) -> Self:
        return self._with_native(self.native.tail(n))

    def round(self, decimals: int) -> Self:
        return self._with_native(self.native.round(decimals=decimals))

    def floor(self) -> Self:
        native = self.native
        native_cls = type(native)
        implementation = self._implementation
        if get_dtype_backend(native.dtype, implementation=implementation) == "pyarrow":
            import pyarrow.compute as pc

            from narwhals._arrow.utils import native_to_narwhals_dtype

            ca = native.array._pa_array
            result_arr = cast("ChunkedArrayAny", pc.floor(ca))
            nw_dtype = native_to_narwhals_dtype(result_arr.type, self._version)
            out_dtype = narwhals_to_native_dtype(
                nw_dtype, "pyarrow", self._implementation, self._version
            )
            result_native = native_cls(
                result_arr, dtype=out_dtype, index=native.index, name=native.name
            )
        else:
            array_funcs = self._array_funcs
            result_arr = array_funcs.floor(self.native)
            result_native = (
                native_cls(result_arr, index=native.index, name=native.name)
                if implementation.is_cudf()
                else result_arr
            )
        return self._with_native(result_native)

    def ceil(self) -> Self:
        native = self.native
        native_cls = type(native)
        implementation = self._implementation
        if get_dtype_backend(native.dtype, implementation=implementation) == "pyarrow":
            import pyarrow.compute as pc

            from narwhals._arrow.utils import native_to_narwhals_dtype

            ca = native.array._pa_array
            result_arr = cast("ChunkedArrayAny", pc.ceil(ca))
            nw_dtype = native_to_narwhals_dtype(result_arr.type, self._version)
            out_dtype = narwhals_to_native_dtype(
                nw_dtype, "pyarrow", self._implementation, self._version
            )
            result_native = native_cls(
                result_arr, dtype=out_dtype, index=native.index, name=native.name
            )
        else:
            array_funcs = self._array_funcs
            result_arr = array_funcs.ceil(self.native)
            result_native = (
                native_cls(result_arr, index=native.index, name=native.name)
                if implementation.is_cudf()
                else result_arr
            )
        return self._with_native(result_native)

    def to_dummies(self, *, separator: str, drop_first: bool) -> PandasLikeDataFrame:
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        plx = self.__native_namespace__()
        series = self.native
        name = str(self._name) if self._name else ""

        null_col_pl = f"{name}{separator}null"

        has_nulls = series.isna().any()
        result = plx.get_dummies(
            series,
            prefix=name,
            prefix_sep=separator,
            drop_first=drop_first,
            # Adds a null column at the end, depending on whether or not there are any.
            dummy_na=has_nulls,
            dtype="int8",
        )
        if has_nulls:
            *cols, null_col_pd = list(result.columns)
            output_order = [null_col_pd, *cols]
            result = rename(
                select_columns_by_name(result, output_order, self._implementation),
                columns={null_col_pd: null_col_pl},
                implementation=self._implementation,
            )
        return PandasLikeDataFrame.from_native(result, context=self)

    def gather_every(self, n: int, offset: int) -> Self:
        return self._with_native(self.native.iloc[offset::n])

    def clip(
        self,
        lower_bound: Self | NumericLiteral | TemporalLiteral | None,
        upper_bound: Self | NumericLiteral | TemporalLiteral | None,
    ) -> Self:
        _, lower = (
            align_and_extract_native(self, lower_bound)
            if lower_bound is not None
            else (None, None)
        )
        _, upper = (
            align_and_extract_native(self, upper_bound)
            if upper_bound is not None
            else (None, None)
        )
        impl = self._implementation
        kwargs: dict[str, Any] = {"axis": 0} if impl.is_modin() else {}
        result = self.native

        if not impl.is_pandas():
            # Workaround for both cudf and modin when clipping with a series
            #   * cudf: https://github.com/rapidsai/cudf/issues/17682
            #   * modin: https://github.com/modin-project/modin/issues/7415
            if self._is_native(lower):
                result = result.where(result >= lower, lower)
                lower = None
            if self._is_native(upper):
                result = result.where(result <= upper, upper)
                upper = None

        return self._with_native(result.clip(lower, upper, **kwargs))

    def to_arrow(self) -> pa.Array[Any]:
        if self._implementation is Implementation.CUDF:
            return self.native.to_arrow()

        import pyarrow as pa  # ignore-banned-import()

        return pa.Array.from_pandas(self.native)

    def mode(self, *, keep: ModeKeepStrategy) -> Self:
        result = self.native.mode()
        result.name = self.name
        return self._with_native(result.head(1) if keep == "any" else result)

    def cum_count(self, *, reverse: bool) -> Self:
        not_na_series = ~self.native.isna()
        result = (
            not_na_series.cumsum()
            if not reverse
            else len(self) - not_na_series.cumsum() + not_na_series - 1
        )
        return self._with_native(result)

    def cum_min(self, *, reverse: bool) -> Self:
        result = (
            self.native.cummin(skipna=True)
            if not reverse
            else self.native[::-1].cummin(skipna=True)[::-1]
        )
        return self._with_native(result)

    def cum_max(self, *, reverse: bool) -> Self:
        result = (
            self.native.cummax(skipna=True)
            if not reverse
            else self.native[::-1].cummax(skipna=True)[::-1]
        )
        return self._with_native(result)

    def cum_prod(self, *, reverse: bool) -> Self:
        result = (
            self.native.cumprod(skipna=True)
            if not reverse
            else self.native[::-1].cumprod(skipna=True)[::-1]
        )
        return self._with_native(result)

    def rolling_sum(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        result = self.native.rolling(
            window=window_size, min_periods=min_samples, center=center
        ).sum()
        return self._with_native(result)

    def rolling_mean(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        result = self.native.rolling(
            window=window_size, min_periods=min_samples, center=center
        ).mean()
        return self._with_native(result)

    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        result = self.native.rolling(
            window=window_size, min_periods=min_samples, center=center
        ).var(ddof=ddof)
        return self._with_native(result)

    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        result = self.native.rolling(
            window=window_size, min_periods=min_samples, center=center
        ).std(ddof=ddof)
        return self._with_native(result)

    def __iter__(self) -> Iterator[Any]:
        yield from self.native.__iter__()

    def __contains__(self, other: Any) -> bool:
        return self.native.isna().any() if other is None else (self.native == other).any()

    def is_finite(self) -> Self:
        s = self.native
        return self._with_native((s > float("-inf")) & (s < float("inf")))

    def rank(self, method: RankMethod, *, descending: bool) -> Self:
        pd_method = "first" if method == "ordinal" else method
        name = self.name
        if (
            self._implementation is Implementation.PANDAS
            and self._backend_version < (3,)
            and get_dtype_backend(self.native.dtype, self._implementation)
            == "numpy_nullable"
            and self.dtype.is_integer()
            and (null_mask := self.is_null()).any()
        ):
            # crazy workaround for the case of `na_option="keep"` and nullable
            # integer dtypes. This should be supported in pandas > 3.0
            # https://github.com/pandas-dev/pandas/issues/56976
            mask_name = f"{name}_is_null"
            plx = self.__narwhals_namespace__()
            df = (
                self.to_frame()
                .with_columns(plx._expr._from_series(null_mask).alias(mask_name))
                .native
            )
            return self._with_native(
                df.groupby(mask_name)
                .rank(
                    method=pd_method,
                    na_option="keep",
                    ascending=not descending,
                    pct=False,
                )
                .iloc[:, 0]
            ).alias(self.name)
        return self._with_native(
            self.native.rank(
                method=pd_method, na_option="keep", ascending=not descending, pct=False
            )
        )

    def hist_from_bins(
        self, bins: list[float], *, include_breakpoint: bool
    ) -> PandasLikeDataFrame:
        return (
            _PandasHist.from_series(self, include_breakpoint=include_breakpoint)
            .with_bins(bins)
            .to_frame()
        )

    def hist_from_bin_count(
        self, bin_count: int, *, include_breakpoint: bool
    ) -> PandasLikeDataFrame:
        return (
            _PandasHist.from_series(self, include_breakpoint=include_breakpoint)
            .with_bin_count(bin_count)
            .to_frame()
        )

    def log(self, base: float) -> Self:
        native = self.native
        native_cls = type(native)
        implementation = self._implementation

        if get_dtype_backend(native.dtype, implementation=implementation) == "pyarrow":
            import pyarrow.compute as pc

            from narwhals._arrow.utils import native_to_narwhals_dtype

            ca = native.array._pa_array
            result_arr = cast("ChunkedArrayAny", pc.logb(ca, base))
            nw_dtype = native_to_narwhals_dtype(result_arr.type, self._version)
            out_dtype = narwhals_to_native_dtype(
                nw_dtype, "pyarrow", self._implementation, self._version
            )
            result_native = native_cls(
                result_arr, dtype=out_dtype, index=native.index, name=native.name
            )
        else:
            array_funcs = self._array_funcs
            result_arr = array_funcs.log(native) / array_funcs.log(base)
            result_native = (
                native_cls(result_arr, index=native.index, name=native.name)
                if implementation.is_cudf()
                else result_arr
            )

        return self._with_native(result_native)

    def exp(self) -> Self:
        native = self.native
        native_cls = type(native)
        implementation = self._implementation

        if get_dtype_backend(native.dtype, implementation=implementation) == "pyarrow":
            import pyarrow.compute as pc

            from narwhals._arrow.utils import native_to_narwhals_dtype

            ca = native.array._pa_array
            result_arr = cast("ChunkedArrayAny", pc.exp(ca))
            nw_dtype = native_to_narwhals_dtype(result_arr.type, self._version)
            out_dtype = narwhals_to_native_dtype(
                nw_dtype, "pyarrow", self._implementation, self._version
            )
            result_native = native_cls(
                result_arr, dtype=out_dtype, index=native.index, name=native.name
            )
        else:
            result_arr = self._array_funcs.exp(native)
            result_native = (
                native_cls(result_arr, index=native.index, name=native.name)
                if implementation.is_cudf()
                else result_arr
            )

        return self._with_native(result_native)

    def sqrt(self) -> Self:
        return self._with_native(self.native.pow(0.5))

    @property
    def str(self) -> PandasLikeSeriesStringNamespace:
        return PandasLikeSeriesStringNamespace(self)

    @property
    def dt(self) -> PandasLikeSeriesDateTimeNamespace:
        return PandasLikeSeriesDateTimeNamespace(self)

    @property
    def cat(self) -> PandasLikeSeriesCatNamespace:
        return PandasLikeSeriesCatNamespace(self)

    @property
    def list(self) -> PandasLikeSeriesListNamespace:
        if not hasattr(self.native, "list"):
            msg = "Series must be of PyArrow List type to support list namespace."
            raise TypeError(msg)
        return PandasLikeSeriesListNamespace(self)

    @property
    def struct(self) -> PandasLikeSeriesStructNamespace:
        if not hasattr(self.native, "struct"):
            msg = "Series must be of PyArrow Struct type to support struct namespace."
            raise TypeError(msg)
        return PandasLikeSeriesStructNamespace(self)


class _PandasHist(EagerSeriesHist["pd.Series[Any]", "list[float]"]):
    _series: PandasLikeSeries

    def to_frame(self) -> PandasLikeDataFrame:
        from_native = self._series.__narwhals_namespace__()._dataframe.from_native
        DataFrame = self._series.__native_namespace__().DataFrame
        return from_native(DataFrame(self._data), context=self._series)

    # NOTE: *Could* be handled at narwhals-level
    def is_empty_series(self) -> bool:
        return self._series.count() < 1

    # NOTE: *Could* be handled at narwhals-level, **iff** we add `nw.repeat`, `nw.linear_space`
    # See https://github.com/narwhals-dev/narwhals/pull/2839#discussion_r2215630696
    def series_empty(self, arg: int | list[float], /) -> PandasHistData:
        count = self._zeros(arg)
        if self._breakpoint:
            return {"breakpoint": self._calculate_breakpoint(arg), "count": count}
        return {"count": count}

    def _zeros(self, arg: int | list[float], /) -> _1DArray:
        zeros = self._series._array_funcs.zeros
        return zeros(arg) if isinstance(arg, int) else zeros(len(arg) - 1)

    # NOTE: Based on `pl.Expr.cut`
    def _cut(
        self,
        breaks: list[float] | _1DArray,
        *,
        labels: Sequence[str] | None = None,
        closed: Literal["left", "right"] = "right",
    ) -> pd.Series[Any]:
        # NOTE: Polars 1.27.0 always includes the lowest bin
        cut = self._series.__native_namespace__().cut
        return cut(
            self.native,
            bins=breaks,
            right=closed == "right",
            labels=labels,
            include_lowest=True,
        )

    def _linear_space(
        self,
        start: float,
        end: float,
        num_samples: int,
        *,
        closed: Literal["both", "none"] = "both",
    ) -> _1DArray:
        return self._series._array_funcs.linspace(
            start=start, stop=end, num=num_samples, endpoint=closed == "both"
        )

    def _calculate_bins(self, bin_count: int) -> _1DArray:
        """Prepare bins for histogram calculation from bin_count."""
        lower, upper = self.native.min(), self.native.max()
        if lower == upper:
            lower -= 0.5
            upper += 0.5
        return self._linear_space(lower, upper, bin_count + 1)

    def _calculate_hist(self, bins: list[float] | _1DArray) -> PandasHistData:
        # pandas (2.2.*) .value_counts(bins=[...]) adjusts the lowest bin which should not
        #   happen since the bins were explicitly passed in.
        categories = self._cut(bins)
        # modin (0.32.0) .value_counts(...) silently drops bins with empty observations,
        #   .reindex is necessary to restore these bins.
        count = categories.value_counts(dropna=True, sort=False).reindex(
            categories.cat.categories, fill_value=0
        )
        count.reset_index(drop=True, inplace=True)  # noqa: PD002
        if self._breakpoint:
            return {"breakpoint": bins[1:], "count": count}
        return {"count": count}
