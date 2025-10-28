from __future__ import annotations

from typing import TYPE_CHECKING, Any, Generic, Literal, Protocol

from narwhals._compliant.any_namespace import (
    CatNamespace,
    DateTimeNamespace,
    ListNamespace,
    StringNamespace,
    StructNamespace,
)
from narwhals._compliant.column import CompliantColumn
from narwhals._compliant.typing import (
    CompliantSeriesT_co,
    EagerDataFrameAny,
    EagerSeriesT_co,
    NativeSeriesT,
    NativeSeriesT_co,
)
from narwhals._translate import FromIterable, FromNative, NumpyConvertible, ToNarwhals
from narwhals._typing_compat import TypeVar, assert_never
from narwhals._utils import (
    _StoresCompliant,
    _StoresNative,
    is_compliant_series,
    is_sized_multi_index_selector,
    unstable,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence
    from types import ModuleType

    import pandas as pd
    import polars as pl
    import pyarrow as pa
    from typing_extensions import NotRequired, Self, TypedDict

    from narwhals._compliant.dataframe import CompliantDataFrame
    from narwhals._compliant.namespace import EagerNamespace
    from narwhals._utils import Implementation, Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.series import Series
    from narwhals.typing import (
        Into1DArray,
        IntoDType,
        MultiIndexSelector,
        PythonLiteral,
        RollingInterpolationMethod,
        SizedMultiIndexSelector,
        _1DArray,
        _SliceIndex,
    )

    class HistData(TypedDict, Generic[NativeSeriesT, "_CountsT_co"]):
        breakpoint: NotRequired[list[float] | _1DArray | list[Any]]
        count: NativeSeriesT | _1DArray | _CountsT_co | list[Any]


_CountsT_co = TypeVar("_CountsT_co", bound="Iterable[Any]", covariant=True)

__all__ = [
    "CompliantSeries",
    "EagerSeries",
    "EagerSeriesCatNamespace",
    "EagerSeriesDateTimeNamespace",
    "EagerSeriesHist",
    "EagerSeriesListNamespace",
    "EagerSeriesNamespace",
    "EagerSeriesStringNamespace",
    "EagerSeriesStructNamespace",
]


class CompliantSeries(
    NumpyConvertible["_1DArray", "Into1DArray"],
    FromIterable,
    FromNative[NativeSeriesT],
    ToNarwhals["Series[NativeSeriesT]"],
    CompliantColumn,
    Protocol[NativeSeriesT],
):
    # NOTE: `narwhals`
    _implementation: Implementation

    @property
    def native(self) -> NativeSeriesT: ...
    def __narwhals_series__(self) -> Self:
        return self

    def __native_namespace__(self) -> ModuleType: ...
    @classmethod
    def from_native(cls, data: NativeSeriesT, /, *, context: _LimitedContext) -> Self: ...
    def to_narwhals(self) -> Series[NativeSeriesT]:
        return self._version.series(self, level="full")

    def _with_native(self, series: Any) -> Self: ...
    def _with_version(self, version: Version) -> Self: ...

    # NOTE: `polars`
    @property
    def dtype(self) -> DType: ...
    @property
    def name(self) -> str: ...
    def __array__(self, dtype: Any, *, copy: bool | None) -> _1DArray: ...
    def __contains__(self, other: Any) -> bool: ...
    def __getitem__(self, item: MultiIndexSelector[Self]) -> Any: ...
    def __iter__(self) -> Iterator[Any]: ...
    def __len__(self) -> int:
        return len(self.native)

    @classmethod
    def from_numpy(cls, data: Into1DArray, /, *, context: _LimitedContext) -> Self: ...
    @classmethod
    def from_iterable(
        cls,
        data: Iterable[Any],
        /,
        *,
        context: _LimitedContext,
        name: str = "",
        dtype: IntoDType | None = None,
    ) -> Self: ...
    def __radd__(self, other: Any) -> Self: ...
    def __rand__(self, other: Any) -> Self: ...
    def __rmul__(self, other: Any) -> Self: ...
    def __ror__(self, other: Any) -> Self: ...
    def all(self) -> bool: ...
    def any(self) -> bool: ...
    def arg_max(self) -> int: ...
    def arg_min(self) -> int: ...
    def arg_true(self) -> Self: ...
    def count(self) -> int: ...
    def filter(self, predicate: Any) -> Self: ...
    def first(self) -> PythonLiteral: ...
    def last(self) -> PythonLiteral: ...
    def gather_every(self, n: int, offset: int) -> Self: ...
    def head(self, n: int) -> Self: ...
    def is_empty(self) -> bool:
        return self.len() == 0

    def is_sorted(self, *, descending: bool) -> bool: ...
    def item(self, index: int | None) -> Any: ...
    def kurtosis(self) -> float | None: ...
    def len(self) -> int: ...
    def max(self) -> Any: ...
    def mean(self) -> float: ...
    def median(self) -> float: ...
    def min(self) -> Any: ...
    def n_unique(self) -> int: ...
    def null_count(self) -> int: ...
    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> float: ...
    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self: ...
    def scatter(self, indices: int | Sequence[int], values: Any) -> Self: ...
    def shift(self, n: int) -> Self: ...
    def skew(self) -> float | None: ...
    def sort(self, *, descending: bool, nulls_last: bool) -> Self: ...
    def std(self, *, ddof: int) -> float: ...
    def sum(self) -> float: ...
    def tail(self, n: int) -> Self: ...
    def to_arrow(self) -> pa.Array[Any]: ...
    def to_dummies(
        self, *, separator: str, drop_first: bool
    ) -> CompliantDataFrame[Self, Any, Any, Any]: ...
    def to_frame(self) -> CompliantDataFrame[Self, Any, Any, Any]: ...
    def to_list(self) -> list[Any]: ...
    def to_pandas(self) -> pd.Series[Any]: ...
    def to_polars(self) -> pl.Series: ...
    def unique(self, *, maintain_order: bool = False) -> Self: ...
    def value_counts(
        self, *, sort: bool, parallel: bool, name: str | None, normalize: bool
    ) -> CompliantDataFrame[Self, Any, Any, Any]: ...
    def var(self, *, ddof: int) -> float: ...
    def zip_with(self, mask: Any, other: Any) -> Self: ...

    # NOTE: *Technically* `polars`
    @unstable
    def hist_from_bins(
        self, bins: list[float], *, include_breakpoint: bool
    ) -> CompliantDataFrame[Self, Any, Any, Any]:
        """`Series.hist(bins=..., bin_count=None)`."""
        ...

    @unstable
    def hist_from_bin_count(
        self, bin_count: int, *, include_breakpoint: bool
    ) -> CompliantDataFrame[Self, Any, Any, Any]:
        """`Series.hist(bins=None, bin_count=...)`."""
        ...


class EagerSeries(CompliantSeries[NativeSeriesT], Protocol[NativeSeriesT]):
    _native_series: Any
    _implementation: Implementation
    _version: Version
    _broadcast: bool

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @classmethod
    def _align_full_broadcast(cls, *series: Self) -> Sequence[Self]:
        """Ensure all of `series` have the same length (and index if `pandas`).

        Scalars get broadcasted to the full length of the longest Series.

        This is useful when you need to construct a full Series anyway, such as:

            DataFrame.select(...)

        It should not be used in binary operations, such as:

            nw.col("a") - nw.col("a").mean()

        because then it's more efficient to extract the right-hand-side's single element as a scalar.
        """
        ...

    def _from_scalar(self, value: Any) -> Self:
        return self.from_iterable([value], name=self.name, context=self)

    def _with_native(
        self, series: NativeSeriesT, *, preserve_broadcast: bool = False
    ) -> Self:
        """Return a new `CompliantSeries`, wrapping the native `series`.

        In cases when operations are known to not affect whether a result should
        be broadcast, we can pass `preserve_broadcast=True`.
        Set this with care - it should only be set for unary expressions which don't
        change length or order, such as `.alias` or `.fill_null`. If in doubt, don't
        set it, you probably don't need it.
        """
        ...

    def __narwhals_namespace__(
        self,
    ) -> EagerNamespace[Any, Self, Any, Any, NativeSeriesT]: ...

    def _gather(self, rows: SizedMultiIndexSelector[NativeSeriesT]) -> Self: ...
    def _gather_slice(self, rows: _SliceIndex | range) -> Self: ...
    def __getitem__(self, item: MultiIndexSelector[Self]) -> Self:
        if isinstance(item, (slice, range)):
            return self._gather_slice(item)
        if is_compliant_series(item):
            return self._gather(item.native)
        elif is_sized_multi_index_selector(item):  # noqa: RET505
            return self._gather(item)
        assert_never(item)

    @property
    def str(self) -> EagerSeriesStringNamespace[Self, NativeSeriesT]: ...
    @property
    def dt(self) -> EagerSeriesDateTimeNamespace[Self, NativeSeriesT]: ...
    @property
    def cat(self) -> EagerSeriesCatNamespace[Self, NativeSeriesT]: ...
    @property
    def list(self) -> EagerSeriesListNamespace[Self, NativeSeriesT]: ...
    @property
    def struct(self) -> EagerSeriesStructNamespace[Self, NativeSeriesT]: ...


class _SeriesNamespace(  # type: ignore[misc]
    _StoresCompliant[CompliantSeriesT_co],
    _StoresNative[NativeSeriesT_co],
    Protocol[CompliantSeriesT_co, NativeSeriesT_co],
):
    _compliant_series: CompliantSeriesT_co

    @property
    def compliant(self) -> CompliantSeriesT_co:
        return self._compliant_series

    @property
    def implementation(self) -> Implementation:
        return self.compliant._implementation

    @property
    def backend_version(self) -> tuple[int, ...]:
        return self.implementation._backend_version()

    @property
    def version(self) -> Version:
        return self.compliant._version

    @property
    def native(self) -> NativeSeriesT_co:
        return self._compliant_series.native  # type: ignore[no-any-return]

    def with_native(self, series: Any, /) -> CompliantSeriesT_co:
        return self.compliant._with_native(series)


class EagerSeriesNamespace(
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    Generic[EagerSeriesT_co, NativeSeriesT_co],
):
    _compliant_series: EagerSeriesT_co

    def __init__(self, series: EagerSeriesT_co, /) -> None:
        self._compliant_series = series


class EagerSeriesCatNamespace(  # type: ignore[misc]
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    CatNamespace[EagerSeriesT_co],
    Protocol[EagerSeriesT_co, NativeSeriesT_co],
): ...


class EagerSeriesDateTimeNamespace(  # type: ignore[misc]
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    DateTimeNamespace[EagerSeriesT_co],
    Protocol[EagerSeriesT_co, NativeSeriesT_co],
): ...


class EagerSeriesListNamespace(  # type: ignore[misc]
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    ListNamespace[EagerSeriesT_co],
    Protocol[EagerSeriesT_co, NativeSeriesT_co],
): ...


class EagerSeriesStringNamespace(  # type: ignore[misc]
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    StringNamespace[EagerSeriesT_co],
    Protocol[EagerSeriesT_co, NativeSeriesT_co],
): ...


class EagerSeriesStructNamespace(  # type: ignore[misc]
    _SeriesNamespace[EagerSeriesT_co, NativeSeriesT_co],
    StructNamespace[EagerSeriesT_co],
    Protocol[EagerSeriesT_co, NativeSeriesT_co],
): ...


class EagerSeriesHist(Protocol[NativeSeriesT, _CountsT_co]):
    _series: EagerSeries[NativeSeriesT]
    _breakpoint: bool
    _data: HistData[NativeSeriesT, _CountsT_co]

    @property
    def native(self) -> NativeSeriesT:
        return self._series.native

    @classmethod
    def from_series(
        cls, series: EagerSeries[NativeSeriesT], *, include_breakpoint: bool
    ) -> Self:
        obj = cls.__new__(cls)
        obj._series = series
        obj._breakpoint = include_breakpoint
        return obj

    def to_frame(self) -> EagerDataFrameAny: ...
    def _linear_space(  # NOTE: Roughly `pl.linear_space`
        self,
        start: float,
        end: float,
        num_samples: int,
        *,
        closed: Literal["both", "none"] = "both",
    ) -> _1DArray: ...

    # NOTE: *Could* be handled at narwhals-level
    def is_empty_series(self) -> bool: ...

    # NOTE: **Should** be handled at narwhals-level
    def data_empty(self) -> HistData[NativeSeriesT, _CountsT_co]:
        return {"breakpoint": [], "count": []} if self._breakpoint else {"count": []}

    # NOTE: *Could* be handled at narwhals-level, **iff** we add `nw.repeat`, `nw.linear_space`
    # See https://github.com/narwhals-dev/narwhals/pull/2839#discussion_r2215630696
    def series_empty(
        self, arg: int | list[float], /
    ) -> HistData[NativeSeriesT, _CountsT_co]: ...

    def with_bins(self, bins: list[float], /) -> Self:
        if len(bins) <= 1:
            self._data = self.data_empty()
        elif self.is_empty_series():
            self._data = self.series_empty(bins)
        else:
            self._data = self._calculate_hist(bins)
        return self

    def with_bin_count(self, bin_count: int, /) -> Self:
        if bin_count == 0:
            self._data = self.data_empty()
        elif self.is_empty_series():
            self._data = self.series_empty(bin_count)
        else:
            self._data = self._calculate_hist(self._calculate_bins(bin_count))
        return self

    def _calculate_breakpoint(self, arg: int | list[float], /) -> list[float] | _1DArray:
        bins = self._linear_space(0, 1, arg + 1) if isinstance(arg, int) else arg
        return bins[1:]

    def _calculate_bins(self, bin_count: int) -> _1DArray: ...
    def _calculate_hist(
        self, bins: list[float] | _1DArray
    ) -> HistData[NativeSeriesT, _CountsT_co]: ...
