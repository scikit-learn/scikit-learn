from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Literal, cast, overload

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.series_cat import ArrowSeriesCatNamespace
from narwhals._arrow.series_dt import ArrowSeriesDateTimeNamespace
from narwhals._arrow.series_list import ArrowSeriesListNamespace
from narwhals._arrow.series_str import ArrowSeriesStringNamespace
from narwhals._arrow.series_struct import ArrowSeriesStructNamespace
from narwhals._arrow.utils import (
    cast_for_truediv,
    chunked_array,
    extract_native,
    floordiv_compat,
    is_array_or_scalar,
    lit,
    narwhals_to_native_dtype,
    native_to_narwhals_dtype,
    nulls_like,
    pad_series,
    zeros,
)
from narwhals._compliant import EagerSeries, EagerSeriesHist
from narwhals._expression_parsing import ExprKind
from narwhals._typing_compat import assert_never
from narwhals._utils import (
    Implementation,
    generate_temporary_column_name,
    is_list_of,
    not_implemented,
)
from narwhals.dependencies import is_numpy_array_1d
from narwhals.exceptions import InvalidOperationError, ShapeError

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from types import ModuleType

    import pandas as pd
    import polars as pl
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._arrow.dataframe import ArrowDataFrame
    from narwhals._arrow.namespace import ArrowNamespace
    from narwhals._arrow.typing import (  # type: ignore[attr-defined]
        ArrayAny,
        ArrayOrChunkedArray,
        ArrayOrScalar,
        ChunkedArrayAny,
        Incomplete,
        NullPlacement,
        Order,
        ScalarAny,
        TieBreaker,
        _AsPyType,
        _BasicDataType,
    )
    from narwhals._compliant.series import HistData
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
        _2DArray,
        _SliceIndex,
    )

    ArrowHistData: TypeAlias = (
        "HistData[ChunkedArrayAny, list[ScalarAny] | pa.Int64Array | list[float]]"
    )


# TODO @dangotbanned: move into `_arrow.utils`
# Lots of modules are importing inline
@overload
def maybe_extract_py_scalar(
    value: pa.Scalar[_BasicDataType[_AsPyType]],
    return_py_scalar: bool,  # noqa: FBT001
) -> _AsPyType: ...


@overload
def maybe_extract_py_scalar(
    value: pa.Scalar[pa.StructType],
    return_py_scalar: bool,  # noqa: FBT001
) -> list[dict[str, Any]]: ...


@overload
def maybe_extract_py_scalar(
    value: pa.Scalar[pa.ListType[_BasicDataType[_AsPyType]]],
    return_py_scalar: bool,  # noqa: FBT001
) -> list[_AsPyType]: ...


@overload
def maybe_extract_py_scalar(
    value: pa.Scalar[Any] | Any,
    return_py_scalar: bool,  # noqa: FBT001
) -> Any: ...


def maybe_extract_py_scalar(value: Any, return_py_scalar: bool) -> Any:  # noqa: FBT001
    if TYPE_CHECKING:
        return value.as_py()
    if return_py_scalar:
        return getattr(value, "as_py", lambda: value)()
    return value


class ArrowSeries(EagerSeries["ChunkedArrayAny"]):
    _implementation = Implementation.PYARROW

    def __init__(
        self, native_series: ChunkedArrayAny, *, name: str, version: Version
    ) -> None:
        self._name = name
        self._native_series: ChunkedArrayAny = native_series
        self._version = version
        self._broadcast = False

    @property
    def native(self) -> ChunkedArrayAny:
        return self._native_series

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, name=self._name, version=version)

    def _with_native(
        self, series: ArrayOrScalar, *, preserve_broadcast: bool = False
    ) -> Self:
        result = self.from_native(chunked_array(series), name=self.name, context=self)
        if preserve_broadcast:
            result._broadcast = self._broadcast
        return result

    def _with_binary(self, op: Callable[..., ArrayOrScalar], other: Any) -> Self:
        ser, other_native = extract_native(self, other)
        preserve_broadcast = self._broadcast and getattr(other, "_broadcast", True)
        return self._with_native(
            op(ser, other_native), preserve_broadcast=preserve_broadcast
        ).alias(self.name)

    def _with_binary_right(self, op: Callable[..., ArrayOrScalar], other: Any) -> Self:
        return self._with_binary(lambda x, y: op(y, x), other).alias(self.name)

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
        if dtype is not None:
            dtype_pa: pa.DataType | None = narwhals_to_native_dtype(dtype, version)
            if is_array_or_scalar(data):
                data = data.cast(dtype_pa)
                dtype_pa = None
            native = data if cls._is_native(data) else chunked_array([data], dtype_pa)
        else:
            native = chunked_array([data])
        return cls.from_native(native, context=context, name=name)

    def _from_scalar(self, value: Any) -> Self:
        if hasattr(value, "as_py"):
            value = value.as_py()
        return super()._from_scalar(value)

    @staticmethod
    def _is_native(obj: ChunkedArrayAny | Any) -> TypeIs[ChunkedArrayAny]:
        return isinstance(obj, pa.ChunkedArray)

    @classmethod
    def from_native(
        cls, data: ChunkedArrayAny, /, *, context: _LimitedContext, name: str = ""
    ) -> Self:
        return cls(data, version=context._version, name=name)

    @classmethod
    def from_numpy(cls, data: Into1DArray, /, *, context: _LimitedContext) -> Self:
        return cls.from_iterable(
            data if is_numpy_array_1d(data) else [data], context=context
        )

    @classmethod
    def _align_full_broadcast(cls, *series: Self) -> Sequence[Self]:
        lengths = [len(s) for s in series]
        max_length = max(lengths)
        fast_path = all(_len == max_length for _len in lengths)
        if fast_path:
            return series
        reshaped = []
        for s in series:
            if s._broadcast:
                compliant = s._with_native(pa.repeat(s.native[0], max_length))
            elif (actual_len := len(s)) != max_length:
                msg = f"Expected object of length {max_length}, got {actual_len}."
                raise ShapeError(msg)
            else:
                compliant = s
            reshaped.append(compliant)
        return reshaped

    def __narwhals_namespace__(self) -> ArrowNamespace:
        from narwhals._arrow.namespace import ArrowNamespace

        return ArrowNamespace(version=self._version)

    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_binary(pc.equal, other)

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_binary(pc.not_equal, other)

    def __ge__(self, other: Any) -> Self:
        return self._with_binary(pc.greater_equal, other)

    def __gt__(self, other: Any) -> Self:
        return self._with_binary(pc.greater, other)

    def __le__(self, other: Any) -> Self:
        return self._with_binary(pc.less_equal, other)

    def __lt__(self, other: Any) -> Self:
        return self._with_binary(pc.less, other)

    def __and__(self, other: Any) -> Self:
        return self._with_binary(pc.and_kleene, other)

    def __rand__(self, other: Any) -> Self:
        return self._with_binary_right(pc.and_kleene, other)

    def __or__(self, other: Any) -> Self:
        return self._with_binary_right(pc.or_kleene, other)

    def __ror__(self, other: Any) -> Self:
        return self._with_binary_right(pc.or_kleene, other)

    def __add__(self, other: Any) -> Self:
        return self._with_binary(pc.add, other)

    def __radd__(self, other: Any) -> Self:
        return self._with_binary_right(pc.add, other)

    def __sub__(self, other: Any) -> Self:
        return self._with_binary(pc.subtract, other)

    def __rsub__(self, other: Any) -> Self:
        return self._with_binary_right(pc.subtract, other)

    def __mul__(self, other: Any) -> Self:
        return self._with_binary(pc.multiply, other)

    def __rmul__(self, other: Any) -> Self:
        return self._with_binary_right(pc.multiply, other)

    def __pow__(self, other: Any) -> Self:
        return self._with_binary(pc.power, other)

    def __rpow__(self, other: Any) -> Self:
        return self._with_binary_right(pc.power, other)

    def __floordiv__(self, other: Any) -> Self:
        return self._with_binary(floordiv_compat, other)

    def __rfloordiv__(self, other: Any) -> Self:
        return self._with_binary_right(floordiv_compat, other)

    def __truediv__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: pc.divide(*cast_for_truediv(x, y)), other)

    def __rtruediv__(self, other: Any) -> Self:
        return self._with_binary_right(
            lambda x, y: pc.divide(*cast_for_truediv(x, y)), other
        )

    def __mod__(self, other: Any) -> Self:
        preserve_broadcast = self._broadcast and getattr(other, "_broadcast", True)
        floor_div = (self // other).native
        ser, other = extract_native(self, other)
        res = pc.subtract(ser, pc.multiply(floor_div, other))
        return self._with_native(res, preserve_broadcast=preserve_broadcast)

    def __rmod__(self, other: Any) -> Self:
        preserve_broadcast = self._broadcast and getattr(other, "_broadcast", True)
        floor_div = (other // self).native
        ser, other = extract_native(self, other)
        res = pc.subtract(other, pc.multiply(floor_div, ser))
        return self._with_native(res, preserve_broadcast=preserve_broadcast)

    def __invert__(self) -> Self:
        return self._with_native(pc.invert(self.native), preserve_broadcast=True)

    @property
    def _type(self) -> pa.DataType:
        return self.native.type

    def len(self, *, _return_py_scalar: bool = True) -> int:
        return maybe_extract_py_scalar(len(self.native), _return_py_scalar)

    def filter(self, predicate: ArrowSeries | list[bool | None]) -> Self:
        other_native: Any
        if not is_list_of(predicate, bool):
            _, other_native = extract_native(self, predicate)
        else:
            other_native = predicate
        return self._with_native(self.native.filter(other_native))

    def first(self, *, _return_py_scalar: bool = True) -> PythonLiteral:
        result = self.native[0] if len(self.native) else None
        return maybe_extract_py_scalar(result, _return_py_scalar)

    def last(self, *, _return_py_scalar: bool = True) -> PythonLiteral:
        ca = self.native
        result = ca[height - 1] if (height := len(ca)) else None
        return maybe_extract_py_scalar(result, _return_py_scalar)

    def mean(self, *, _return_py_scalar: bool = True) -> float:
        return maybe_extract_py_scalar(pc.mean(self.native), _return_py_scalar)

    def median(self, *, _return_py_scalar: bool = True) -> float:
        if not self.dtype.is_numeric():
            msg = "`median` operation not supported for non-numeric input type."
            raise InvalidOperationError(msg)

        return maybe_extract_py_scalar(
            pc.approximate_median(self.native), _return_py_scalar
        )

    def min(self, *, _return_py_scalar: bool = True) -> Any:
        return maybe_extract_py_scalar(pc.min(self.native), _return_py_scalar)

    def max(self, *, _return_py_scalar: bool = True) -> Any:
        return maybe_extract_py_scalar(pc.max(self.native), _return_py_scalar)

    def arg_min(self, *, _return_py_scalar: bool = True) -> int:
        index_min = pc.index(self.native, pc.min(self.native))
        return maybe_extract_py_scalar(index_min, _return_py_scalar)

    def arg_max(self, *, _return_py_scalar: bool = True) -> int:
        index_max = pc.index(self.native, pc.max(self.native))
        return maybe_extract_py_scalar(index_max, _return_py_scalar)

    def sum(self, *, _return_py_scalar: bool = True) -> float:
        return maybe_extract_py_scalar(
            pc.sum(self.native, min_count=0), _return_py_scalar
        )

    def drop_nulls(self) -> Self:
        return self._with_native(self.native.drop_null())

    def shift(self, n: int) -> Self:
        if n > 0:
            arrays = [nulls_like(n, self), *self.native[:-n].chunks]
        elif n < 0:
            arrays = [*self.native[-n:].chunks, nulls_like(-n, self)]
        else:
            return self._with_native(self.native)
        return self._with_native(pa.concat_arrays(arrays))

    def std(self, *, ddof: int, _return_py_scalar: bool = True) -> float:
        return maybe_extract_py_scalar(
            pc.stddev(self.native, ddof=ddof), _return_py_scalar
        )

    def var(self, *, ddof: int, _return_py_scalar: bool = True) -> float:
        return maybe_extract_py_scalar(
            pc.variance(self.native, ddof=ddof), _return_py_scalar
        )

    def skew(self, *, _return_py_scalar: bool = True) -> float | None:
        ser_not_null = self.native.drop_null()
        if len(ser_not_null) == 0:
            return None
        if len(ser_not_null) == 1:
            return float("nan")
        if len(ser_not_null) == 2:
            return 0.0
        m = pc.subtract(ser_not_null, pc.mean(ser_not_null))
        m2 = pc.mean(pc.power(m, lit(2)))
        m3 = pc.mean(pc.power(m, lit(3)))
        biased_population_skewness = pc.divide(m3, pc.power(m2, lit(1.5)))
        return maybe_extract_py_scalar(biased_population_skewness, _return_py_scalar)

    def kurtosis(self, *, _return_py_scalar: bool = True) -> float | None:
        ser_not_null = self.native.drop_null()
        if len(ser_not_null) == 0:
            return None
        if len(ser_not_null) == 1:
            return float("nan")
        m = pc.subtract(ser_not_null, pc.mean(ser_not_null))
        m2 = pc.mean(pc.power(m, lit(2)))
        m4 = pc.mean(pc.power(m, lit(4)))
        k = pc.subtract(pc.divide(m4, pc.power(m2, lit(2))), lit(3))
        return maybe_extract_py_scalar(k, _return_py_scalar)

    def count(self, *, _return_py_scalar: bool = True) -> int:
        return maybe_extract_py_scalar(pc.count(self.native), _return_py_scalar)

    def n_unique(self, *, _return_py_scalar: bool = True) -> int:
        return maybe_extract_py_scalar(
            pc.count(self.native.unique(), mode="all"), _return_py_scalar
        )

    def __native_namespace__(self) -> ModuleType:
        if self._implementation is Implementation.PYARROW:
            return self._implementation.to_native_namespace()

        msg = f"Expected pyarrow, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    @property
    def name(self) -> str:
        return self._name

    def _gather(self, rows: SizedMultiIndexSelector[ChunkedArrayAny]) -> Self:
        if len(rows) == 0:
            return self._with_native(self.native.slice(0, 0))
        if self._backend_version < (18,) and isinstance(rows, tuple):
            rows = list(rows)
        return self._with_native(self.native.take(rows))

    def _gather_slice(self, rows: _SliceIndex | range) -> Self:
        start = rows.start or 0
        stop = rows.stop if rows.stop is not None else len(self.native)
        if start < 0:
            start = len(self.native) + start
        if stop < 0:
            stop = len(self.native) + stop
        if rows.step is not None and rows.step != 1:
            msg = "Slicing with step is not supported on PyArrow tables"
            raise NotImplementedError(msg)
        return self._with_native(self.native.slice(start, stop - start))

    def scatter(self, indices: int | Sequence[int], values: Any) -> Self:
        import numpy as np  # ignore-banned-import

        values_native: ArrayAny
        if isinstance(indices, int):
            indices_native = pa.array([indices])
            values_native = pa.array([values])
        else:
            # TODO(unassigned): we may also want to let `indices` be a Series.
            # https://github.com/narwhals-dev/narwhals/issues/2155
            indices_native = pa.array(indices)
            if isinstance(values, self.__class__):
                values_native = values.native.combine_chunks()
            else:
                # NOTE: Requires fixes in https://github.com/zen-xu/pyarrow-stubs/pull/209
                pa_array: Incomplete = pa.array
                values_native = pa_array(values)

        sorting_indices = pc.sort_indices(indices_native)
        indices_native = indices_native.take(sorting_indices)
        values_native = values_native.take(sorting_indices)

        mask: _1DArray = np.zeros(self.len(), dtype=bool)
        mask[indices_native] = True
        # NOTE: Multiple issues
        # - Missing `values` type
        # - `mask` accepts a `np.ndarray`, but not mentioned in stubs
        # - Missing `replacements` type
        # - Missing return type
        pc_replace_with_mask: Incomplete = pc.replace_with_mask
        return self._with_native(pc_replace_with_mask(self.native, mask, values_native))

    def to_list(self) -> list[Any]:
        return self.native.to_pylist()

    def __array__(self, dtype: Any = None, *, copy: bool | None = None) -> _1DArray:
        return self.native.__array__(dtype=dtype, copy=copy)

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _1DArray:
        return self.native.to_numpy()

    def alias(self, name: str) -> Self:
        result = self.__class__(self.native, name=name, version=self._version)
        result._broadcast = self._broadcast
        return result

    @property
    def dtype(self) -> DType:
        return native_to_narwhals_dtype(self.native.type, self._version)

    def abs(self) -> Self:
        return self._with_native(pc.abs(self.native))

    def cum_sum(self, *, reverse: bool) -> Self:
        cum_sum = pc.cumulative_sum
        result = (
            cum_sum(self.native, skip_nulls=True)
            if not reverse
            else cum_sum(self.native[::-1], skip_nulls=True)[::-1]
        )
        return self._with_native(result)

    def round(self, decimals: int) -> Self:
        return self._with_native(
            pc.round(self.native, decimals, round_mode="half_towards_infinity")
        )

    def floor(self) -> Self:
        return self._with_native(pc.floor(self.native))

    def ceil(self) -> Self:
        return self._with_native(pc.ceil(self.native))

    def diff(self) -> Self:
        return self._with_native(pc.pairwise_diff(self.native.combine_chunks()))

    def any(self, *, _return_py_scalar: bool = True) -> bool:
        return maybe_extract_py_scalar(
            pc.any(self.native, min_count=0), _return_py_scalar
        )

    def all(self, *, _return_py_scalar: bool = True) -> bool:
        return maybe_extract_py_scalar(
            pc.all(self.native, min_count=0), _return_py_scalar
        )

    def is_between(
        self, lower_bound: Any, upper_bound: Any, closed: ClosedInterval
    ) -> Self:
        _, lower_bound = extract_native(self, lower_bound)
        _, upper_bound = extract_native(self, upper_bound)
        if closed == "left":
            ge = pc.greater_equal(self.native, lower_bound)
            lt = pc.less(self.native, upper_bound)
            res = pc.and_kleene(ge, lt)
        elif closed == "right":
            gt = pc.greater(self.native, lower_bound)
            le = pc.less_equal(self.native, upper_bound)
            res = pc.and_kleene(gt, le)
        elif closed == "none":
            gt = pc.greater(self.native, lower_bound)
            lt = pc.less(self.native, upper_bound)
            res = pc.and_kleene(gt, lt)
        elif closed == "both":
            ge = pc.greater_equal(self.native, lower_bound)
            le = pc.less_equal(self.native, upper_bound)
            res = pc.and_kleene(ge, le)
        else:
            assert_never(closed)
        return self._with_native(res)

    def is_null(self) -> Self:
        return self._with_native(self.native.is_null(), preserve_broadcast=True)

    def is_nan(self) -> Self:
        return self._with_native(pc.is_nan(self.native), preserve_broadcast=True)

    def cast(self, dtype: IntoDType) -> Self:
        data_type = narwhals_to_native_dtype(dtype, self._version)
        return self._with_native(pc.cast(self.native, data_type), preserve_broadcast=True)

    def null_count(self, *, _return_py_scalar: bool = True) -> int:
        return maybe_extract_py_scalar(self.native.null_count, _return_py_scalar)

    def head(self, n: int) -> Self:
        if n >= 0:
            return self._with_native(self.native.slice(0, n))
        num_rows = len(self)
        return self._with_native(self.native.slice(0, max(0, num_rows + n)))

    def tail(self, n: int) -> Self:
        if n >= 0:
            num_rows = len(self)
            return self._with_native(self.native.slice(max(0, num_rows - n)))
        return self._with_native(self.native.slice(abs(n)))

    def is_in(self, other: Any) -> Self:
        if self._is_native(other):
            value_set: ArrayOrChunkedArray = other
        else:
            value_set = pa.array(other)
        return self._with_native(pc.is_in(self.native, value_set=value_set))

    def arg_true(self) -> Self:
        import numpy as np  # ignore-banned-import

        res = np.flatnonzero(self.native)
        return self.from_iterable(res, name=self.name, context=self)

    def item(self, index: int | None = None) -> Any:
        if index is None:
            if len(self) != 1:
                msg = (
                    "can only call '.item()' if the Series is of length 1,"
                    f" or an explicit index is provided (Series is of length {len(self)})"
                )
                raise ValueError(msg)
            return maybe_extract_py_scalar(self.native[0], return_py_scalar=True)
        return maybe_extract_py_scalar(self.native[index], return_py_scalar=True)

    def value_counts(
        self, *, sort: bool, parallel: bool, name: str | None, normalize: bool
    ) -> ArrowDataFrame:
        """Parallel is unused, exists for compatibility."""
        from narwhals._arrow.dataframe import ArrowDataFrame

        index_name_ = "index" if self._name is None else self._name
        value_name_ = name or ("proportion" if normalize else "count")

        val_counts = pc.value_counts(self.native)
        values = val_counts.field("values")
        counts = cast("ChunkedArrayAny", val_counts.field("counts"))

        if normalize:
            arrays = [values, pc.divide(*cast_for_truediv(counts, pc.sum(counts)))]
        else:
            arrays = [values, counts]

        val_count = pa.Table.from_arrays(arrays, names=[index_name_, value_name_])

        if sort:
            val_count = val_count.sort_by([(value_name_, "descending")])

        return ArrowDataFrame(
            val_count, version=self._version, validate_column_names=True
        )

    def zip_with(self, mask: Self, other: Self) -> Self:
        cond = mask.native.combine_chunks()
        return self._with_native(pc.if_else(cond, self.native, other.native))

    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self:
        import numpy as np  # ignore-banned-import

        num_rows = len(self)
        if n is None and fraction is not None:
            n = int(num_rows * fraction)

        rng = np.random.default_rng(seed=seed)
        idx = np.arange(num_rows)
        mask = rng.choice(idx, size=n, replace=with_replacement)
        return self._with_native(self.native.take(mask))

    def fill_nan(self, value: float | None) -> Self:
        result = pc.if_else(pc.is_nan(self.native), value, self.native)
        return self._with_native(result, preserve_broadcast=True)

    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self:
        import numpy as np  # ignore-banned-import

        def fill_aux(
            arr: ChunkedArrayAny, limit: int, direction: FillNullStrategy | None
        ) -> ArrayAny:
            # this algorithm first finds the indices of the valid values to fill all the null value positions
            # then it calculates the distance of each new index and the original index
            # if the distance is equal to or less than the limit and the original value is null, it is replaced
            valid_mask = pc.is_valid(arr)
            indices = pa.array(np.arange(len(arr)), type=pa.int64())
            if direction == "forward":
                valid_index = np.maximum.accumulate(np.where(valid_mask, indices, -1))
                distance = indices - valid_index
            else:
                valid_index = np.minimum.accumulate(
                    np.where(valid_mask[::-1], indices[::-1], len(arr))
                )[::-1]
                distance = valid_index - indices
            return pc.if_else(
                pc.and_(pc.is_null(arr), pc.less_equal(distance, lit(limit))),  # pyright: ignore[reportArgumentType, reportCallIssue]
                arr.take(valid_index),
                arr,
            )

        if value is not None:
            _, native_value = extract_native(self, value)
            series: ArrayOrScalar = pc.fill_null(self.native, native_value)
        elif limit is None:
            fill_func = (
                pc.fill_null_forward if strategy == "forward" else pc.fill_null_backward
            )
            series = fill_func(self.native)
        else:
            series = fill_aux(self.native, limit, strategy)
        return self._with_native(series, preserve_broadcast=True)

    def to_frame(self) -> ArrowDataFrame:
        from narwhals._arrow.dataframe import ArrowDataFrame

        df = pa.Table.from_arrays([self.native], names=[self.name])
        return ArrowDataFrame(df, version=self._version, validate_column_names=False)

    def to_pandas(self) -> pd.Series[Any]:
        import pandas as pd  # ignore-banned-import()

        return pd.Series(self.native, name=self.name)

    def to_polars(self) -> pl.Series:
        import polars as pl  # ignore-banned-import

        return cast("pl.Series", pl.from_arrow(self.native))

    def is_unique(self) -> ArrowSeries:
        return self.to_frame().is_unique().alias(self.name)

    def is_first_distinct(self) -> Self:
        import numpy as np  # ignore-banned-import

        row_number = pa.array(np.arange(len(self)))
        col_token = generate_temporary_column_name(n_bytes=8, columns=[self.name])
        first_distinct_index = (
            pa.Table.from_arrays([self.native], names=[self.name])
            .append_column(col_token, row_number)
            .group_by(self.name)
            .aggregate([(col_token, "min")])
            .column(f"{col_token}_min")
        )

        return self._with_native(pc.is_in(row_number, first_distinct_index))

    def is_last_distinct(self) -> Self:
        import numpy as np  # ignore-banned-import

        row_number = pa.array(np.arange(len(self)))
        col_token = generate_temporary_column_name(n_bytes=8, columns=[self.name])
        last_distinct_index = (
            pa.Table.from_arrays([self.native], names=[self.name])
            .append_column(col_token, row_number)
            .group_by(self.name)
            .aggregate([(col_token, "max")])
            .column(f"{col_token}_max")
        )

        return self._with_native(pc.is_in(row_number, last_distinct_index))

    def is_sorted(self, *, descending: bool) -> bool:
        if not isinstance(descending, bool):
            msg = f"argument 'descending' should be boolean, found {type(descending)}"
            raise TypeError(msg)
        if descending:
            result = pc.all(pc.greater_equal(self.native[:-1], self.native[1:]))
        else:
            result = pc.all(pc.less_equal(self.native[:-1], self.native[1:]))
        return maybe_extract_py_scalar(result, return_py_scalar=True)

    def unique(self, *, maintain_order: bool = True) -> Self:
        # TODO(marco): `pc.unique` seems to always maintain order, is that guaranteed?
        return self._with_native(self.native.unique())

    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self:
        # https://stackoverflow.com/a/79111029/4451315
        idxs = pc.index_in(self.native, pa.array(old))
        result_native = pc.take(pa.array(new), idxs)
        if return_dtype is not None:
            result_native.cast(narwhals_to_native_dtype(return_dtype, self._version))
        result = self._with_native(result_native)
        if result.is_null().sum() != self.is_null().sum():
            msg = (
                "replace_strict did not replace all non-null values.\n\n"
                "The following did not get replaced: "
                f"{self.filter(~self.is_null() & result.is_null()).unique(maintain_order=False).to_list()}"
            )
            raise ValueError(msg)
        return result

    def sort(self, *, descending: bool, nulls_last: bool) -> Self:
        order: Order = "descending" if descending else "ascending"
        null_placement: NullPlacement = "at_end" if nulls_last else "at_start"
        sorted_indices = pc.array_sort_indices(
            self.native, order=order, null_placement=null_placement
        )
        return self._with_native(self.native.take(sorted_indices))

    def to_dummies(self, *, separator: str, drop_first: bool) -> ArrowDataFrame:
        import numpy as np  # ignore-banned-import

        from narwhals._arrow.dataframe import ArrowDataFrame

        name = self._name
        # NOTE: stub is missing attributes (https://arrow.apache.org/docs/python/generated/pyarrow.DictionaryArray.html)
        da: Incomplete = self.native.combine_chunks().dictionary_encode("encode")

        columns: _2DArray = np.zeros((len(da.dictionary), len(da)), np.int8)
        columns[da.indices, np.arange(len(da))] = 1
        null_col_pa, null_col_pl = f"{name}{separator}None", f"{name}{separator}null"
        cols = [
            {null_col_pa: null_col_pl}.get(
                f"{name}{separator}{v}", f"{name}{separator}{v}"
            )
            for v in da.dictionary
        ]

        output_order = (
            [
                null_col_pl,
                *sorted([c for c in cols if c != null_col_pl])[int(drop_first) :],
            ]
            if null_col_pl in cols
            else sorted(cols)[int(drop_first) :]
        )
        return ArrowDataFrame(
            pa.Table.from_arrays(columns, names=cols),
            version=self._version,
            validate_column_names=True,
        ).simple_select(*output_order)

    def quantile(
        self,
        quantile: float,
        interpolation: RollingInterpolationMethod,
        *,
        _return_py_scalar: bool = True,
    ) -> float:
        return maybe_extract_py_scalar(
            pc.quantile(self.native, q=quantile, interpolation=interpolation)[0],
            _return_py_scalar,
        )

    def gather_every(self, n: int, offset: int = 0) -> Self:
        return self._with_native(self.native[offset::n])

    def clip(
        self,
        lower_bound: Self | NumericLiteral | TemporalLiteral | None,
        upper_bound: Self | NumericLiteral | TemporalLiteral | None,
    ) -> Self:
        _, lower = (
            extract_native(self, lower_bound) if lower_bound is not None else (None, None)
        )
        _, upper = (
            extract_native(self, upper_bound) if upper_bound is not None else (None, None)
        )

        if lower is None:
            return self._with_native(pc.min_element_wise(self.native, upper))
        if upper is None:
            return self._with_native(pc.max_element_wise(self.native, lower))
        return self._with_native(
            pc.max_element_wise(pc.min_element_wise(self.native, upper), lower)
        )

    def to_arrow(self) -> ArrayAny:
        return self.native.combine_chunks()

    def mode(self, *, keep: ModeKeepStrategy) -> ArrowSeries:
        plx = self.__narwhals_namespace__()
        col_token = generate_temporary_column_name(n_bytes=8, columns=[self.name])
        counts = self.value_counts(
            name=col_token, normalize=False, sort=False, parallel=False
        )
        result = counts.filter(
            plx.col(col_token)
            == plx.col(col_token).max().broadcast(kind=ExprKind.AGGREGATION)
        ).get_column(self.name)
        return result.head(1) if keep == "any" else result

    def is_finite(self) -> Self:
        return self._with_native(pc.is_finite(self.native))

    def cum_count(self, *, reverse: bool) -> Self:
        dtypes = self._version.dtypes
        return (~self.is_null()).cast(dtypes.UInt32()).cum_sum(reverse=reverse)

    def cum_min(self, *, reverse: bool) -> Self:
        result = (
            pc.cumulative_min(self.native, skip_nulls=True)
            if not reverse
            else pc.cumulative_min(self.native[::-1], skip_nulls=True)[::-1]
        )
        return self._with_native(result)

    def cum_max(self, *, reverse: bool) -> Self:
        result = (
            pc.cumulative_max(self.native, skip_nulls=True)
            if not reverse
            else pc.cumulative_max(self.native[::-1], skip_nulls=True)[::-1]
        )
        return self._with_native(result)

    def cum_prod(self, *, reverse: bool) -> Self:
        result = (
            pc.cumulative_prod(self.native, skip_nulls=True)
            if not reverse
            else pc.cumulative_prod(self.native[::-1], skip_nulls=True)[::-1]
        )
        return self._with_native(result)

    def rolling_sum(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        min_samples = min_samples if min_samples is not None else window_size
        padded_series, offset = pad_series(self, window_size=window_size, center=center)

        cum_sum = padded_series.cum_sum(reverse=False).fill_null(
            value=None, strategy="forward", limit=None
        )
        rolling_sum = (
            cum_sum
            - cum_sum.shift(window_size).fill_null(value=0, strategy=None, limit=None)
            if window_size != 0
            else cum_sum
        )

        valid_count = padded_series.cum_count(reverse=False)
        count_in_window = valid_count - valid_count.shift(window_size).fill_null(
            value=0, strategy=None, limit=None
        )

        result = self._with_native(
            pc.if_else((count_in_window >= min_samples).native, rolling_sum.native, None)
        )
        return result._gather_slice(slice(offset, None))

    def rolling_mean(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        min_samples = min_samples if min_samples is not None else window_size
        padded_series, offset = pad_series(self, window_size=window_size, center=center)

        cum_sum = padded_series.cum_sum(reverse=False).fill_null(
            value=None, strategy="forward", limit=None
        )
        rolling_sum = (
            cum_sum
            - cum_sum.shift(window_size).fill_null(value=0, strategy=None, limit=None)
            if window_size != 0
            else cum_sum
        )

        valid_count = padded_series.cum_count(reverse=False)
        count_in_window = valid_count - valid_count.shift(window_size).fill_null(
            value=0, strategy=None, limit=None
        )

        result = (
            self._with_native(
                pc.if_else(
                    (count_in_window >= min_samples).native, rolling_sum.native, None
                )
            )
            / count_in_window
        )
        return result._gather_slice(slice(offset, None))

    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        min_samples = min_samples if min_samples is not None else window_size
        padded_series, offset = pad_series(self, window_size=window_size, center=center)

        cum_sum = padded_series.cum_sum(reverse=False).fill_null(
            value=None, strategy="forward", limit=None
        )
        rolling_sum = (
            cum_sum
            - cum_sum.shift(window_size).fill_null(value=0, strategy=None, limit=None)
            if window_size != 0
            else cum_sum
        )

        cum_sum_sq = (
            pow(padded_series, 2)
            .cum_sum(reverse=False)
            .fill_null(value=None, strategy="forward", limit=None)
        )
        rolling_sum_sq = (
            cum_sum_sq
            - cum_sum_sq.shift(window_size).fill_null(value=0, strategy=None, limit=None)
            if window_size != 0
            else cum_sum_sq
        )

        valid_count = padded_series.cum_count(reverse=False)
        count_in_window = valid_count - valid_count.shift(window_size).fill_null(
            value=0, strategy=None, limit=None
        )

        result = self._with_native(
            pc.if_else(
                (count_in_window >= min_samples).native,
                (rolling_sum_sq - (rolling_sum**2 / count_in_window)).native,
                None,
            )
        ) / self._with_native(pc.max_element_wise((count_in_window - ddof).native, 0))

        return result._gather_slice(slice(offset, None, None))

    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        return (
            self.rolling_var(
                window_size=window_size, min_samples=min_samples, center=center, ddof=ddof
            )
            ** 0.5
        )

    def rank(self, method: RankMethod, *, descending: bool) -> Self:
        if method == "average":
            msg = (
                "`rank` with `method='average' is not supported for pyarrow backend. "
                "The available methods are {'min', 'max', 'dense', 'ordinal'}."
            )
            raise ValueError(msg)

        sort_keys: Order = "descending" if descending else "ascending"
        tiebreaker: TieBreaker = "first" if method == "ordinal" else method

        native_series: ArrayOrChunkedArray
        if self._backend_version < (14, 0, 0):  # pragma: no cover
            native_series = self.native.combine_chunks()
        else:
            native_series = self.native

        null_mask = pc.is_null(native_series)

        rank = pc.rank(native_series, sort_keys=sort_keys, tiebreaker=tiebreaker)

        result = pc.if_else(null_mask, lit(None, rank.type), rank)
        return self._with_native(result)

    def hist_from_bins(
        self, bins: list[float], *, include_breakpoint: bool
    ) -> ArrowDataFrame:
        return (
            _ArrowHist.from_series(self, include_breakpoint=include_breakpoint)
            .with_bins(bins)
            .to_frame()
        )

    def hist_from_bin_count(
        self, bin_count: int, *, include_breakpoint: bool
    ) -> ArrowDataFrame:
        return (
            _ArrowHist.from_series(self, include_breakpoint=include_breakpoint)
            .with_bin_count(bin_count)
            .to_frame()
        )

    def __iter__(self) -> Iterator[Any]:
        for x in self.native:
            yield maybe_extract_py_scalar(x, return_py_scalar=True)

    def __contains__(self, other: Any) -> bool:
        from pyarrow import (
            ArrowInvalid,  # ignore-banned-imports
            ArrowNotImplementedError,  # ignore-banned-imports
            ArrowTypeError,  # ignore-banned-imports
        )

        try:
            other_ = lit(other) if other is not None else lit(None, type=self._type)
            return maybe_extract_py_scalar(
                pc.is_in(other_, self.native), return_py_scalar=True
            )
        except (ArrowInvalid, ArrowNotImplementedError, ArrowTypeError) as exc:
            msg = f"Unable to compare other of type {type(other)} with series of type {self.dtype}."
            raise InvalidOperationError(msg) from exc

    def log(self, base: float) -> Self:
        return self._with_native(pc.logb(self.native, lit(base)))

    def exp(self) -> Self:
        return self._with_native(pc.exp(self.native))

    def sqrt(self) -> Self:
        return self._with_native(pc.sqrt(self.native))

    @property
    def dt(self) -> ArrowSeriesDateTimeNamespace:
        return ArrowSeriesDateTimeNamespace(self)

    @property
    def cat(self) -> ArrowSeriesCatNamespace:
        return ArrowSeriesCatNamespace(self)

    @property
    def str(self) -> ArrowSeriesStringNamespace:
        return ArrowSeriesStringNamespace(self)

    @property
    def list(self) -> ArrowSeriesListNamespace:
        return ArrowSeriesListNamespace(self)

    @property
    def struct(self) -> ArrowSeriesStructNamespace:
        return ArrowSeriesStructNamespace(self)

    ewm_mean = not_implemented()


class _ArrowHist(
    EagerSeriesHist["ChunkedArrayAny", "list[ScalarAny] | pa.Int64Array | list[float]"]
):
    _series: ArrowSeries

    def to_frame(self) -> ArrowDataFrame:
        # NOTE: Constructor typing is too strict for `TypedDict`
        table: Incomplete = pa.Table.from_pydict
        from_native = self._series.__narwhals_namespace__()._dataframe.from_native
        return from_native(table(self._data), context=self._series)

    # NOTE: *Could* be handled at narwhals-level
    def is_empty_series(self) -> bool:
        # NOTE: `ChunkedArray.combine_chunks` returns the concrete array type
        # Stubs say `Array[pa.BooleanScalar]`, which is missing properties
        # https://github.com/zen-xu/pyarrow-stubs/blob/6bedee748bc74feb8513b24bf43d64b24c7fddc8/pyarrow-stubs/__lib_pxi/array.pyi#L2395-L2399
        is_null = self.native.is_null(nan_is_null=True)
        arr = cast("pa.BooleanArray", is_null.combine_chunks())
        return arr.false_count == 0

    # NOTE: *Could* be handled at narwhals-level, **iff** we add `nw.repeat`, `nw.linear_space`
    # See https://github.com/narwhals-dev/narwhals/pull/2839#discussion_r2215630696
    def series_empty(self, arg: int | list[float], /) -> ArrowHistData:
        count = self._zeros(arg)
        if self._breakpoint:
            return {"breakpoint": self._calculate_breakpoint(arg), "count": count}
        return {"count": count}

    def _zeros(self, arg: int | list[float], /) -> pa.Int64Array:
        return zeros(arg) if isinstance(arg, int) else zeros(len(arg) - 1)

    def _linear_space(
        self,
        start: float,
        end: float,
        num_samples: int,
        *,
        closed: Literal["both", "none"] = "both",
    ) -> _1DArray:
        from numpy import linspace  # ignore-banned-import

        return linspace(start=start, stop=end, num=num_samples, endpoint=closed == "both")

    def _calculate_bins(self, bin_count: int) -> _1DArray:
        """Prepare bins for histogram calculation from bin_count."""
        d = pc.min_max(self.native)
        lower, upper = d["min"].as_py(), d["max"].as_py()
        if lower == upper:
            lower -= 0.5
            upper += 0.5
        return self._linear_space(lower, upper, bin_count + 1)

    def _calculate_hist(self, bins: list[float] | _1DArray) -> ArrowHistData:
        ser = self.native
        # NOTE: `mypy` refuses to resolve `ndarray.__getitem__`
        # Previously annotated as `list[float]`, but
        # - wasn't accurate to how we implemented it
        # - `pa.scalar` overloads fail to match on `float | np.float64` (but runtime is fine)
        bins = cast("list[float]", bins)
        # Handle single bin case
        if len(bins) == 2:
            is_between_bins = pc.and_(
                pc.greater_equal(ser, lit(bins[0])), pc.less_equal(ser, lit(bins[1]))
            )
            count = pc.sum(is_between_bins.cast(pa.uint8()))
            if self._breakpoint:
                return {"breakpoint": [bins[-1]], "count": [count]}
            return {"count": [count]}

        # Handle multiple bins
        import numpy as np  # ignore-banned-import

        bin_indices = np.searchsorted(bins, ser, side="left")
        # lowest bin is inclusive
        bin_indices = pc.if_else(pc.equal(ser, lit(bins[0])), 1, bin_indices)

        # Align unique categories and counts appropriately
        obs_cats, obs_counts = np.unique(bin_indices, return_counts=True)
        obj_cats = np.arange(1, len(bins))
        counts = np.zeros_like(obj_cats)
        counts[np.isin(obj_cats, obs_cats)] = obs_counts[np.isin(obs_cats, obj_cats)]

        if self._breakpoint:
            return {"breakpoint": bins[1:], "count": counts}
        return {"count": counts}
