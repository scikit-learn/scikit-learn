from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING, Any, Callable, Final, Literal, cast, overload

import narwhals as nw
from narwhals import exceptions, functions as nw_f
from narwhals._exceptions import issue_warning
from narwhals._typing_compat import TypeVar, assert_never
from narwhals._utils import (
    Implementation,
    Version,
    deprecate_native_namespace,
    generate_temporary_column_name,
    inherit_doc,
    is_ordered_categorical,
    maybe_align_index,
    maybe_convert_dtypes,
    maybe_get_index,
    maybe_reset_index,
    maybe_set_index,
    validate_strict_and_pass_though,
)
from narwhals.dataframe import DataFrame as NwDataFrame, LazyFrame as NwLazyFrame
from narwhals.exceptions import InvalidIntoExprError
from narwhals.expr import Expr as NwExpr
from narwhals.functions import _new_series_impl, concat, show_versions
from narwhals.schema import Schema as NwSchema
from narwhals.series import Series as NwSeries
from narwhals.stable.v1 import dependencies, dtypes, selectors
from narwhals.stable.v1.dtypes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Field,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Object,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Unknown,
)
from narwhals.stable.v1.typing import IntoDataFrameT, IntoLazyFrameT
from narwhals.translate import _from_native_impl, get_native_namespace, to_py_scalar

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence
    from types import ModuleType

    from typing_extensions import ParamSpec, Self

    from narwhals._translate import IntoArrowTable
    from narwhals._typing import (
        Arrow,
        Backend,
        EagerAllowed,
        IntoBackend,
        LazyAllowed,
        Pandas,
        Polars,
    )
    from narwhals.dataframe import MultiColSelector, MultiIndexSelector
    from narwhals.dtypes import DType
    from narwhals.typing import (
        FileSource,
        IntoDType,
        IntoExpr,
        IntoFrame,
        IntoSchema,
        IntoSeries,
        NonNestedLiteral,
        SingleColSelector,
        SingleIndexSelector,
        _1DArray,
        _2DArray,
    )

    DataFrameT = TypeVar("DataFrameT", bound="DataFrame[Any]")
    LazyFrameT = TypeVar("LazyFrameT", bound="LazyFrame[Any]")
    SeriesT = TypeVar("SeriesT", bound="Series[Any]")
    T = TypeVar("T", default=Any)
    P = ParamSpec("P")
    R = TypeVar("R")

IntoSeriesT = TypeVar("IntoSeriesT", bound="IntoSeries", default=Any)


# NOTE legit
class DataFrame(NwDataFrame[IntoDataFrameT]):  # type: ignore[type-var]
    _version = Version.V1

    @inherit_doc(NwDataFrame)
    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        assert df._version is Version.V1  # noqa: S101
        super().__init__(df, level=level)

    # We need to override any method which don't return Self so that type
    # annotations are correct.

    @classmethod
    def from_arrow(
        cls, native_frame: IntoArrowTable, *, backend: IntoBackend[EagerAllowed]
    ) -> DataFrame[Any]:
        result = super().from_arrow(native_frame, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        schema: Mapping[str, DType] | Schema | None = None,
        *,
        backend: IntoBackend[EagerAllowed] | None = None,
    ) -> DataFrame[Any]:
        result = super().from_dict(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Any],
        schema: IntoSchema | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> DataFrame[Any]:
        result = super().from_dicts(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        schema: Mapping[str, DType] | Schema | Sequence[str] | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> DataFrame[Any]:
        result = super().from_numpy(data, schema, backend=backend)
        return cast("DataFrame[Any]", result)

    @property
    def _series(self) -> type[Series[Any]]:
        return cast("type[Series[Any]]", Series)

    @property
    def _lazyframe(self) -> type[LazyFrame[Any]]:
        return cast("type[LazyFrame[Any]]", LazyFrame)

    @overload
    def __getitem__(self, item: tuple[SingleIndexSelector, SingleColSelector]) -> Any: ...

    @overload
    def __getitem__(  # type: ignore[overload-overlap]
        self, item: str | tuple[MultiIndexSelector, SingleColSelector]
    ) -> Series[Any]: ...

    @overload
    def __getitem__(
        self,
        item: (
            SingleIndexSelector
            | MultiIndexSelector
            | MultiColSelector
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> Self: ...
    def __getitem__(
        self,
        item: (
            SingleIndexSelector
            | SingleColSelector
            | MultiColSelector
            | MultiIndexSelector
            | tuple[SingleIndexSelector, SingleColSelector]
            | tuple[SingleIndexSelector, MultiColSelector]
            | tuple[MultiIndexSelector, SingleColSelector]
            | tuple[MultiIndexSelector, MultiColSelector]
        ),
    ) -> Series[Any] | Self | Any:
        return super().__getitem__(item)

    def get_column(self, name: str) -> Series:
        # Type checkers complain that `nw.Series` is not assignable to `nw.v1.stable.Series`.
        # However the return type actually is `nw.v1.stable.Series`, check `tests/v1_test.py`.
        return super().get_column(name)  # type: ignore[return-value]

    def lazy(
        self,
        backend: IntoBackend[LazyAllowed] | None = None,
        *,
        session: Any | None = None,
    ) -> LazyFrame[Any]:
        return _stableify(super().lazy(backend=backend, session=session))

    @overload  # type: ignore[override]
    def to_dict(self, *, as_series: Literal[True] = ...) -> dict[str, Series[Any]]: ...
    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...
    @overload
    def to_dict(
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]: ...
    def to_dict(  # pyright: ignore[reportIncompatibleMethodOverride]
        self, *, as_series: bool = True
    ) -> dict[str, Series[Any]] | dict[str, list[Any]]:
        # Type checkers complain that `nw.Series` is not assignable to `nw.v1.stable.Series`.
        # However the return type actually is `nw.v1.stable.Series`, check `tests/v1_test.py::test_to_dict_as_series`.
        return super().to_dict(as_series=as_series)  # type: ignore[return-value]

    def is_duplicated(self) -> Series[Any]:
        return _stableify(super().is_duplicated())

    def is_unique(self) -> Series[Any]:
        return _stableify(super().is_unique())

    def _l1_norm(self) -> Self:
        # Private, just used to test the stable API.
        return self.select(all()._l1_norm())


class LazyFrame(NwLazyFrame[IntoLazyFrameT]):
    @inherit_doc(NwLazyFrame)
    def __init__(self, df: Any, *, level: Literal["full", "lazy", "interchange"]) -> None:
        assert df._version is Version.V1  # noqa: S101
        super().__init__(df, level=level)

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        return DataFrame

    def _extract_compliant(self, arg: Any) -> Any:
        # After v1, we raise when passing order-dependent, length-changing,
        # or filtration expressions to LazyFrame
        from narwhals.expr import Expr
        from narwhals.series import Series

        if isinstance(arg, Series):  # pragma: no cover
            msg = "Mixing Series with LazyFrame is not supported."
            raise TypeError(msg)
        if isinstance(arg, (Expr, str)):
            return self.__narwhals_namespace__().parse_into_expr(arg, str_as_lit=False)
        raise InvalidIntoExprError.from_invalid_type(type(arg))

    def collect(
        self, backend: IntoBackend[Polars | Pandas | Arrow] | None = None, **kwargs: Any
    ) -> DataFrame[Any]:
        return _stableify(super().collect(backend=backend, **kwargs))

    def _l1_norm(self) -> Self:
        # Private, just used to test the stable API.
        return self.select(all()._l1_norm())

    def tail(self, n: int = 5) -> Self:
        r"""Get the last `n` rows."""
        return super().tail(n)

    def gather_every(self, n: int, offset: int = 0) -> Self:
        r"""Take every nth row in the DataFrame and return as a new DataFrame.

        Arguments:
            n: Gather every *n*-th row.
            offset: Starting index.
        """
        return self._with_compliant(
            self._compliant_frame.gather_every(n=n, offset=offset)  # type: ignore[attr-defined]
        )

    def with_row_index(
        self, name: str = "index", *, order_by: str | Sequence[str] | None = None
    ) -> Self:
        order_by_ = [order_by] if isinstance(order_by, str) else order_by
        return self._with_compliant(
            self._compliant_frame.with_row_index(
                name=name,
                order_by=order_by_,  # type: ignore[arg-type]
            )
        )


class Series(NwSeries[IntoSeriesT]):
    _version = Version.V1

    @inherit_doc(NwSeries)
    def __init__(
        self, series: Any, *, level: Literal["full", "lazy", "interchange"]
    ) -> None:
        assert series._version is Version.V1  # noqa: S101
        super().__init__(series, level=level)

    # We need to override any method which don't return Self so that type
    # annotations are correct.

    @classmethod
    def from_numpy(
        cls,
        name: str,
        values: _1DArray,
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        result = super().from_numpy(name, values, dtype, backend=backend)
        return cast("Series[Any]", result)

    @classmethod
    def from_iterable(
        cls,
        name: str,
        values: Iterable[Any],
        dtype: IntoDType | None = None,
        *,
        backend: IntoBackend[EagerAllowed],
    ) -> Series[Any]:
        result = super().from_iterable(name, values, dtype, backend=backend)
        return cast("Series[Any]", result)

    @property
    def _dataframe(self) -> type[DataFrame[Any]]:
        return DataFrame

    def to_frame(self) -> DataFrame[Any]:
        return _stableify(super().to_frame())

    def value_counts(
        self,
        *,
        sort: bool = False,
        parallel: bool = False,
        name: str | None = None,
        normalize: bool = False,
    ) -> DataFrame[Any]:
        return _stableify(
            super().value_counts(
                sort=sort, parallel=parallel, name=name, normalize=normalize
            )
        )

    def hist(
        self,
        bins: list[float] | None = None,
        *,
        bin_count: int | None = None,
        include_breakpoint: bool = True,
    ) -> DataFrame[Any]:
        from narwhals.exceptions import NarwhalsUnstableWarning

        msg = (
            "`Series.hist` is being called from the stable API although considered "
            "an unstable feature."
        )
        issue_warning(msg, NarwhalsUnstableWarning)
        return _stableify(
            super().hist(
                bins=bins, bin_count=bin_count, include_breakpoint=include_breakpoint
            )
        )


class Expr(NwExpr):
    def _l1_norm(self) -> Self:
        return super()._taxicab_norm()

    def head(self, n: int = 10) -> Self:
        r"""Get the first `n` rows."""
        return self._with_orderable_filtration(
            lambda plx: self._to_compliant_expr(plx).head(n)  # type: ignore[attr-defined]
        )

    def tail(self, n: int = 10) -> Self:
        r"""Get the last `n` rows."""
        return self._with_orderable_filtration(
            lambda plx: self._to_compliant_expr(plx).tail(n)  # type: ignore[attr-defined]
        )

    def gather_every(self, n: int, offset: int = 0) -> Self:
        r"""Take every nth value in the Series and return as new Series.

        Arguments:
            n: Gather every *n*-th row.
            offset: Starting index.
        """
        return self._with_orderable_filtration(
            lambda plx: self._to_compliant_expr(plx).gather_every(n=n, offset=offset)  # type: ignore[attr-defined]
        )

    def unique(self, *, maintain_order: bool | None = None) -> Self:
        """Return unique values of this expression."""
        if maintain_order is not None:
            msg = (
                "`maintain_order` has no effect and is only kept around for backwards-compatibility. "
                "You can safely remove this argument."
            )
            issue_warning(msg, UserWarning)
        return self._with_filtration(lambda plx: self._to_compliant_expr(plx).unique())

    def sort(self, *, descending: bool = False, nulls_last: bool = False) -> Self:
        """Sort this column. Place null values first."""
        return self._with_window(
            lambda plx: self._to_compliant_expr(plx).sort(  # type: ignore[attr-defined]
                descending=descending, nulls_last=nulls_last
            )
        )

    def arg_max(self) -> Self:
        """Returns the index of the maximum value."""
        return self._with_orderable_aggregation(
            lambda plx: self._to_compliant_expr(plx).arg_max()  # type: ignore[attr-defined]
        )

    def arg_min(self) -> Self:
        """Returns the index of the minimum value."""
        return self._with_orderable_aggregation(
            lambda plx: self._to_compliant_expr(plx).arg_min()  # type: ignore[attr-defined]
        )

    def arg_true(self) -> Self:
        """Find elements where boolean expression is True."""
        return self._with_orderable_filtration(
            lambda plx: self._to_compliant_expr(plx).arg_true()  # type: ignore[attr-defined]
        )

    def sample(
        self,
        n: int | None = None,
        *,
        fraction: float | None = None,
        with_replacement: bool = False,
        seed: int | None = None,
    ) -> Self:
        """Sample randomly from this expression.

        Arguments:
            n: Number of items to return. Cannot be used with fraction.
            fraction: Fraction of items to return. Cannot be used with n.
            with_replacement: Allow values to be sampled more than once.
            seed: Seed for the random number generator. If set to None (default), a random
                seed is generated for each sample operation.
        """
        return self._with_filtration(
            lambda plx: self._to_compliant_expr(plx).sample(  # type: ignore[attr-defined]
                n, fraction=fraction, with_replacement=with_replacement, seed=seed
            )
        )


class Schema(NwSchema):
    _version = Version.V1

    @inherit_doc(NwSchema)
    def __init__(
        self, schema: Mapping[str, DType] | Iterable[tuple[str, DType]] | None = None
    ) -> None:
        super().__init__(schema)


@overload
def _stableify(obj: NwDataFrame[IntoDataFrameT]) -> DataFrame[IntoDataFrameT]: ...  # type: ignore[type-var]
@overload
def _stableify(obj: NwLazyFrame[IntoLazyFrameT]) -> LazyFrame[IntoLazyFrameT]: ...
@overload
def _stableify(obj: NwSeries[IntoSeriesT]) -> Series[IntoSeriesT]: ...
@overload
def _stableify(obj: NwExpr) -> Expr: ...


def _stableify(
    obj: NwDataFrame[IntoDataFrameT]  # type: ignore[type-var]
    | NwLazyFrame[IntoLazyFrameT]
    | NwSeries[IntoSeriesT]
    | NwExpr,
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT] | Series[IntoSeriesT] | Expr:
    if isinstance(obj, NwDataFrame):
        return DataFrame(obj._compliant_frame._with_version(Version.V1), level=obj._level)
    if isinstance(obj, NwLazyFrame):
        return LazyFrame(obj._compliant_frame._with_version(Version.V1), level=obj._level)
    if isinstance(obj, NwSeries):
        return Series(obj._compliant_series._with_version(Version.V1), level=obj._level)
    if isinstance(obj, NwExpr):
        return Expr(obj._to_compliant_expr, obj._metadata)
    assert_never(obj)


@overload
def from_native(native_object: SeriesT, **kwds: Any) -> SeriesT: ...


@overload
def from_native(native_object: DataFrameT, **kwds: Any) -> DataFrameT: ...


@overload
def from_native(native_object: LazyFrameT, **kwds: Any) -> LazyFrameT: ...


@overload
def from_native(
    native_object: DataFrameT | LazyFrameT, **kwds: Any
) -> DataFrameT | LazyFrameT: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeriesT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeriesT,
    *,
    strict: Literal[False],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[False],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    strict: Literal[False],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoLazyFrameT | IntoSeriesT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoLazyFrameT,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> LazyFrame[IntoLazyFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    strict: Literal[False],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoLazyFrameT,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> LazyFrame[IntoLazyFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoFrame | IntoSeries,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[Any] | LazyFrame[Any] | Series[Any]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    strict: Literal[True] | None = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeries,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    pass_through: Literal[True],
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoLazyFrameT | IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT] | Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT | IntoLazyFrameT,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT] | LazyFrame[IntoLazyFrameT]: ...


@overload
def from_native(
    native_object: T,
    *,
    pass_through: Literal[True],
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> T: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[True],
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[True],
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoFrame | IntoSeries,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: Literal[True],
) -> DataFrame[Any] | LazyFrame[Any] | Series[Any]: ...


@overload
def from_native(
    native_object: IntoSeriesT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[True],
    allow_series: None = ...,
) -> Series[IntoSeriesT]: ...


@overload
def from_native(
    native_object: IntoDataFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> DataFrame[IntoDataFrameT]: ...


@overload
def from_native(
    native_object: IntoLazyFrameT,
    *,
    pass_through: Literal[False] = ...,
    eager_only: Literal[False] = ...,
    eager_or_interchange_only: Literal[False] = ...,
    series_only: Literal[False] = ...,
    allow_series: None = ...,
) -> LazyFrame[IntoLazyFrameT]: ...


# All params passed in as variables
@overload
def from_native(
    native_object: Any,
    *,
    pass_through: bool,
    eager_only: bool,
    eager_or_interchange_only: bool = False,
    series_only: bool,
    allow_series: bool | None,
) -> Any: ...


def from_native(
    native_object: IntoDataFrameT
    | IntoLazyFrameT
    | IntoFrame
    | IntoSeriesT
    | IntoSeries
    | T,
    *,
    strict: bool | None = None,
    pass_through: bool | None = None,
    eager_only: bool = False,
    eager_or_interchange_only: bool = False,
    series_only: bool = False,
    allow_series: bool | None = None,
    **kwds: Any,
) -> LazyFrame[IntoLazyFrameT] | DataFrame[IntoDataFrameT] | Series[IntoSeriesT] | T:
    """Convert `native_object` to Narwhals Dataframe, Lazyframe, or Series.

    See `narwhals.from_native` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    # Early returns
    if isinstance(native_object, (DataFrame, LazyFrame)) and not series_only:
        return native_object
    if isinstance(native_object, Series) and (series_only or allow_series):
        return native_object

    pass_through = validate_strict_and_pass_though(
        strict, pass_through, pass_through_default=False
    )
    if kwds:
        msg = f"from_native() got an unexpected keyword argument {next(iter(kwds))!r}"
        raise TypeError(msg)

    return _from_native_impl(  # type: ignore[no-any-return]
        native_object,
        pass_through=pass_through,
        eager_only=eager_only,
        eager_or_interchange_only=eager_or_interchange_only,
        series_only=series_only,
        allow_series=allow_series,
        version=Version.V1,
    )


@overload
def to_native(
    narwhals_object: DataFrame[IntoDataFrameT], *, strict: Literal[True] = ...
) -> IntoDataFrameT: ...
@overload
def to_native(
    narwhals_object: LazyFrame[IntoLazyFrameT], *, strict: Literal[True] = ...
) -> IntoLazyFrameT: ...
@overload
def to_native(
    narwhals_object: Series[IntoSeriesT], *, strict: Literal[True] = ...
) -> IntoSeriesT: ...
@overload
def to_native(narwhals_object: Any, *, strict: bool) -> Any: ...
@overload
def to_native(
    narwhals_object: DataFrame[IntoDataFrameT], *, pass_through: Literal[False] = ...
) -> IntoDataFrameT: ...
@overload
def to_native(
    narwhals_object: LazyFrame[IntoLazyFrameT], *, pass_through: Literal[False] = ...
) -> IntoLazyFrameT: ...
@overload
def to_native(
    narwhals_object: Series[IntoSeriesT], *, pass_through: Literal[False] = ...
) -> IntoSeriesT: ...
@overload
def to_native(narwhals_object: Any, *, pass_through: bool) -> Any: ...


def to_native(
    narwhals_object: DataFrame[IntoDataFrameT]
    | LazyFrame[IntoLazyFrameT]
    | Series[IntoSeriesT],
    *,
    strict: bool | None = None,
    pass_through: bool | None = None,
) -> IntoLazyFrameT | IntoDataFrameT | IntoSeriesT | Any:
    """Convert Narwhals object to native one.

    See `narwhals.to_native` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    from narwhals._utils import validate_strict_and_pass_though

    pass_through = validate_strict_and_pass_though(
        strict, pass_through, pass_through_default=False
    )
    return nw.to_native(narwhals_object, pass_through=pass_through)


def narwhalify(
    func: Callable[..., Any] | None = None,
    *,
    strict: bool | None = None,
    pass_through: bool | None = None,
    eager_only: bool = False,
    eager_or_interchange_only: bool = False,
    series_only: bool = False,
    allow_series: bool | None = True,
) -> Callable[..., Any]:
    """Decorate function so it becomes dataframe-agnostic.

    See `narwhals.narwhalify` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    pass_through = validate_strict_and_pass_though(
        strict, pass_through, pass_through_default=True
    )

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            args = [
                from_native(
                    arg,
                    pass_through=pass_through,
                    eager_only=eager_only,
                    eager_or_interchange_only=eager_or_interchange_only,
                    series_only=series_only,
                    allow_series=allow_series,
                )
                for arg in args
            ]  # type: ignore[assignment]

            kwargs = {
                name: from_native(
                    value,
                    pass_through=pass_through,
                    eager_only=eager_only,
                    eager_or_interchange_only=eager_or_interchange_only,
                    series_only=series_only,
                    allow_series=allow_series,
                )
                for name, value in kwargs.items()
            }

            backends = {
                b()
                for v in (*args, *kwargs.values())
                if (b := getattr(v, "__native_namespace__", None))
            }

            if backends.__len__() > 1:
                msg = "Found multiple backends. Make sure that all dataframe/series inputs come from the same backend."
                raise ValueError(msg)

            result = func(*args, **kwargs)

            return to_native(result, pass_through=pass_through)

        return wrapper

    if func is None:
        return decorator
    # If func is not None, it means the decorator is used without arguments
    return decorator(func)


def all() -> Expr:
    return _stableify(nw.all())


def col(*names: str | Iterable[str]) -> Expr:
    return _stableify(nw.col(*names))


def exclude(*names: str | Iterable[str]) -> Expr:
    return _stableify(nw.exclude(*names))


def nth(*indices: int | Sequence[int]) -> Expr:
    return _stableify(nw.nth(*indices))


def len() -> Expr:
    return _stableify(nw.len())


def lit(value: NonNestedLiteral, dtype: IntoDType | None = None) -> Expr:
    return _stableify(nw.lit(value, dtype))


def min(*columns: str) -> Expr:
    return _stableify(nw.min(*columns))


def max(*columns: str) -> Expr:
    return _stableify(nw.max(*columns))


def mean(*columns: str) -> Expr:
    return _stableify(nw.mean(*columns))


def median(*columns: str) -> Expr:
    return _stableify(nw.median(*columns))


def sum(*columns: str) -> Expr:
    return _stableify(nw.sum(*columns))


def sum_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    return _stableify(nw.sum_horizontal(*exprs))


def all_horizontal(
    *exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool = False
) -> Expr:
    return _stableify(nw.all_horizontal(*exprs, ignore_nulls=ignore_nulls))


def any_horizontal(
    *exprs: IntoExpr | Iterable[IntoExpr], ignore_nulls: bool = False
) -> Expr:
    return _stableify(nw.any_horizontal(*exprs, ignore_nulls=ignore_nulls))


def mean_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    return _stableify(nw.mean_horizontal(*exprs))


def min_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    return _stableify(nw.min_horizontal(*exprs))


def max_horizontal(*exprs: IntoExpr | Iterable[IntoExpr]) -> Expr:
    return _stableify(nw.max_horizontal(*exprs))


def concat_str(
    exprs: IntoExpr | Iterable[IntoExpr],
    *more_exprs: IntoExpr,
    separator: str = "",
    ignore_nulls: bool = False,
) -> Expr:
    return _stableify(
        nw.concat_str(exprs, *more_exprs, separator=separator, ignore_nulls=ignore_nulls)
    )


def format(f_string: str, *args: IntoExpr) -> Expr:
    """Format expressions as a string."""
    return _stableify(nw.format(f_string, *args))


def coalesce(exprs: IntoExpr | Iterable[IntoExpr], *more_exprs: IntoExpr) -> Expr:
    return _stableify(nw.coalesce(exprs, *more_exprs))


def get_level(
    obj: DataFrame[Any] | LazyFrame[Any] | Series[IntoSeriesT],
) -> Literal["full", "lazy", "interchange"]:
    """Level of support Narwhals has for current object.

    Arguments:
        obj: Dataframe or Series.

    Returns:
        This can be one of

            - 'full': full Narwhals API support
            - 'lazy': only lazy operations are supported. This excludes anything
              which involves iterating over rows in Python.
            - 'interchange': only metadata operations are supported (`df.schema`)
    """
    return obj._level


class When(nw_f.When):
    @classmethod
    def from_when(cls, when: nw_f.When) -> When:
        return cls(when._predicate)

    def then(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Then:
        return Then.from_then(super().then(value))


class Then(nw_f.Then, Expr):
    @classmethod
    def from_then(cls, then: nw_f.Then) -> Then:
        return cls(then._to_compliant_expr, then._metadata)

    def otherwise(self, value: IntoExpr | NonNestedLiteral | _1DArray) -> Expr:
        return _stableify(super().otherwise(value))


def when(*predicates: IntoExpr | Iterable[IntoExpr]) -> When:
    return When.from_when(nw_f.when(*predicates))


@deprecate_native_namespace(required=True)
def new_series(
    name: str,
    values: Any,
    dtype: IntoDType | None = None,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
) -> Series[Any]:
    """Instantiate Narwhals Series from iterable (e.g. list or array).

    See `narwhals.new_series` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[EagerAllowed]", backend)
    return _stableify(_new_series_impl(name, values, dtype, backend=backend))


@deprecate_native_namespace(required=True)
def from_arrow(
    native_frame: IntoArrowTable,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
) -> DataFrame[Any]:
    """Construct a DataFrame from an object which supports the PyCapsule Interface.

    See `narwhals.from_arrow` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[EagerAllowed]", backend)
    return _stableify(nw_f.from_arrow(native_frame, backend=backend))


@deprecate_native_namespace()
def from_dict(
    data: Mapping[str, Any],
    schema: Mapping[str, DType] | Schema | None = None,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
) -> DataFrame[Any]:
    """Instantiate DataFrame from dictionary.

    See `narwhals.from_dict` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    return _stableify(nw_f.from_dict(data, schema, backend=backend))


from_dicts: Final = DataFrame.from_dicts


@deprecate_native_namespace(required=True)
def from_numpy(
    data: _2DArray,
    schema: Mapping[str, DType] | Schema | Sequence[str] | None = None,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
) -> DataFrame[Any]:
    """Construct a DataFrame from a NumPy ndarray.

    See `narwhals.from_numpy` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[EagerAllowed]", backend)
    return _stableify(nw_f.from_numpy(data, schema, backend=backend))


@deprecate_native_namespace(required=True)
def read_csv(
    source: FileSource,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
    **kwargs: Any,
) -> DataFrame[Any]:
    """Read a CSV file into a DataFrame.

    See `narwhals.read_csv` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[EagerAllowed]", backend)
    return _stableify(nw_f.read_csv(source, backend=backend, **kwargs))


@deprecate_native_namespace(required=True)
def scan_csv(
    source: FileSource,
    *,
    backend: IntoBackend[Backend] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
    **kwargs: Any,
) -> LazyFrame[Any]:
    """Lazily read from a CSV file.

    See `narwhals.scan_csv` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[Backend]", backend)
    return _stableify(nw_f.scan_csv(source, backend=backend, **kwargs))


@deprecate_native_namespace(required=True)
def read_parquet(
    source: FileSource,
    *,
    backend: IntoBackend[EagerAllowed] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
    **kwargs: Any,
) -> DataFrame[Any]:
    """Read into a DataFrame from a parquet file.

    See `narwhals.read_parquet` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[EagerAllowed]", backend)
    return _stableify(nw_f.read_parquet(source, backend=backend, **kwargs))


@deprecate_native_namespace(required=True)
def scan_parquet(
    source: FileSource,
    *,
    backend: IntoBackend[Backend] | None = None,
    native_namespace: ModuleType | None = None,  # noqa: ARG001
    **kwargs: Any,
) -> LazyFrame[Any]:
    """Lazily read from a parquet file.

    See `narwhals.scan_parquet` for full docstring. Note that `native_namespace` is
    an is the same as `backend` but only accepts module types - for new code, we
    recommend using `backend`, as that's available beyond just `narwhals.stable.v1`.
    """
    backend = cast("IntoBackend[Backend]", backend)
    return _stableify(nw_f.scan_parquet(source, backend=backend, **kwargs))


__all__ = [
    "Array",
    "Binary",
    "Boolean",
    "Categorical",
    "DataFrame",
    "Date",
    "Datetime",
    "Decimal",
    "Duration",
    "Enum",
    "Expr",
    "Field",
    "Float32",
    "Float64",
    "Implementation",
    "Int8",
    "Int16",
    "Int32",
    "Int64",
    "Int128",
    "LazyFrame",
    "List",
    "Object",
    "Schema",
    "Series",
    "String",
    "Struct",
    "Time",
    "UInt8",
    "UInt16",
    "UInt32",
    "UInt64",
    "UInt128",
    "Unknown",
    "all",
    "all_horizontal",
    "any_horizontal",
    "coalesce",
    "col",
    "concat",
    "concat_str",
    "dependencies",
    "dtypes",
    "exceptions",
    "exclude",
    "format",
    "from_arrow",
    "from_dict",
    "from_dicts",
    "from_native",
    "from_numpy",
    "generate_temporary_column_name",
    "get_level",
    "get_native_namespace",
    "is_ordered_categorical",
    "len",
    "lit",
    "max",
    "max_horizontal",
    "maybe_align_index",
    "maybe_convert_dtypes",
    "maybe_get_index",
    "maybe_reset_index",
    "maybe_set_index",
    "mean",
    "mean_horizontal",
    "median",
    "min",
    "min_horizontal",
    "narwhalify",
    "new_series",
    "nth",
    "read_csv",
    "read_parquet",
    "scan_csv",
    "scan_parquet",
    "selectors",
    "show_versions",
    "sum",
    "sum_horizontal",
    "to_native",
    "to_py_scalar",
    "when",
]
